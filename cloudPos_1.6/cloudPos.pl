#Perl -Sx "{0}"; Exit {Status}
#!perl

use strict;
use warnings;

use File::Glob ':glob';
use POSIX;
#use Math::Random;

use globals qw(ProgSetUp finish);
use kClouds qw(buildPClouds importClouds exportPClouds getCloudKmers revCompKmers);
use seqTools qw(readFasta makeSeqHash);

#globals
my %program =(
	     id => 1,
	     name => "cloudPos.pl",
         version => "1.6",
	     nickname => "cloudPos",
	     authors => "David D. Pollock, Corey Cox",
	     began => "11/29/2012",
	     modified => "02/04/2015", 	#updateable
	     uses => "Find locations of clouds in fasta file",
	     runrecord => "RunRecord_cloudPos.txt",
	     computer => "Craggenmore",	#updateable
	    );

#Settings for using control/default/factory system
my $globals = hash_ref();
my $factoryfile = "factory"; 	# factory file setting is meant to be an example, or fixed if *

#main memory
my $silent = 0; my $quiet = 1; my $murmur = 2; my $talk = 3; my $loud = 5; 
my $read = "<"; my $write = ">"; my $append = ">>"; my $read_append = "+>>"; my $no_record = 1;

### Your variables that will be used through the program (main memory) should go here.
my $clouds = {}; my $kmers = {}; my $seqs = {}; my $kseqs = {}; my $ptrs = {};
my $search_kmers = {}; my $locs = [];

#############################
##                         ##
##       Begin Main        ##
##                         ##
#############################

#system("purge");	# clear out inactive memory - on my system each purge takees about 8 seconds to run
ProgSetUp($factoryfile, $globals, \%program);	# read control files, import all values in $globals hash
my $flags = {
    kmerbits => 1,
    recordpos => 1,  ### currently positions are needed to build clouds!!!
    cloudbits => 1,
    cloudid => 1,    ### currently required for proper export and position annotation
    consensus => 1,
    old_himp_clouds => 0
};

### Get P-Clouds for further processing - import or build our own from fasta seqs.
if ($globals->{import_clouds}) {
    readx($globals->{import_clouds}, \&importClouds, "importing P-Clouds", $clouds, $kmers, $flags);
}
elsif ($globals->{cloud_fasta}) {
    readFasta($globals->{cloud_fasta}, \&makeSeqHash, "building P-Clouds", $seqs, $ptrs);
	buildPClouds($seqs, $kseqs, $clouds, $kmers, $globals->{core}, $globals->{outer}, $globals->{kmer_len}, $flags);
}
else { die "\n###\n\nNo clouds imported or built!!\nEnding run.\n\n###\n"; }

my $cloudkmers = getCloudKmers($clouds);

if ($globals->{export_clouds}) {
    printx($globals->{export_clouds}, \&exportPClouds, "exporting PClouds", $kmers);
}

### Annotate clouds onto sequences/genome
if ($globals->{revcomp_kmers}) { %$search_kmers = %$cloudkmers; revCompKmers($search_kmers); } else { $search_kmers = $cloudkmers; }
my $files = getFiles($globals->{seq_prefix}, $globals->{seq_postfix}, $globals->{seq_dir});

if (1 == $globals->{no_memorize}) {
    my $fp_out = openup($globals->{full_output}, $write, "output file");
    multireadx($files, \&getCloudPos_noMem, "reading sequence file", $search_kmers, $globals->{merge_window}, $globals->{print_ids}, $fp_out);
}
else {
    multireadx($files, \&getCloudPos, "reading sequence file", $search_kmers, $locs);
    my $print_locs = $locs;
    if ($globals->{merge_locations}) { $print_locs = mergeLocs($locs); }
    printx($globals->{full_output}, \&cloudLocsToBed, "exporting cloud locations to BED", $print_locs, $globals->{print_ids});
}

#system("purge"); # clear out inactive memory
finish($globals->{starttime}); # finish program and close out

#############################
##                         ##
##        End Main         ##
##                         ##
#############################

##					##
##   subroutines 	##
##					##
sub getCloudPos { my ($fp_in, $kmers, $locs) = @_;
    my $nextLoc = makeLocFinder($fp_in, $kmers);
    my $loc = $nextLoc->();
    while ($loc->{begin} > 0) {
        push (@$locs, $loc);
        $loc = $nextLoc->();
    }
}

sub mergeLocs { my $locs = shift; #my $merge = shift;
    my $newlocs = array_ref(); my $new_loc;
    foreach my $loc (@$locs) {
        if ($new_loc && $new_loc->{contig} eq $loc->{contig} && $new_loc->{end} > $loc->{begin}) {
            $new_loc->{end} = $loc->{end}; push (@{$new_loc->{kmers}}, $loc->{kinfo});
        }
        else { $new_loc = { %$loc }; #  clone current location to new hash with same info
            push (@$newlocs, $new_loc);
        }
    } return $newlocs;
}

sub getCloudPos_noMem { my ($fp_in, $kmers, $merge, $print_ids, $fp_out) = @_;
    if (! defined $merge) { $merge = 0; }
    my $nextLoc = makeLocFinder($fp_in, $kmers);
    my $prevLoc = $nextLoc->();
    while ($prevLoc->{begin} > 0) {
        my $loc = $nextLoc->();
        if ($prevLoc->{contig} eq $loc->{contig}) {
            if ($prevLoc->{end} + $merge >= $loc->{begin}) {
                $loc->{begin} = $prevLoc->{begin};
                push (@{$prevLoc->{kmers}}, $loc->{kinfo});
                $loc->{kmers} = $prevLoc->{kmers}; $loc->{kinfo} = $prevLoc->{kinfo};
            }
            else { locToBed($fp_out, $prevLoc, $print_ids); }
        }
        else { locToBed($fp_out, $prevLoc, $print_ids); }
        $prevLoc = $loc;
    }
}

sub makeLocFinder { my ($fp_in, $kmers) = @_;
    my $pos = 0; my $contig = ""; my $buffer = ""; my $i = 0; my $len = length((keys %{$kmers})[0]); my $ilen = 0; my $line = "";
    return sub {
        while (1) {
            for (; $i < $ilen; $i++) { my $kmer = substr ($buffer, $i, $len);
                if ($kmers->{$kmer}) {
                    my $loc = { begin => $pos, end => $pos + $len, kinfo => $kmers->{$kmer}, contig => $contig };
                    $pos++; $i++; push (@{$loc->{kmers}}, $kmers->{$kmer});
                    return $loc;
                } $pos++;
            }
            $i = 0; $buffer = substr ($buffer, $ilen);
            if ($line = <$fp_in>) { chomp $line; } else { return { begin => -1, end => -1, contig => "" }; }
            if ($line =~ s/^>//) { $contig = $line; $pos = 0; $buffer = "";
                talk("Found contig $contig; resetting position\n", $murmur);
                $line = <$fp_in>; chomp $line;
            }
            $buffer .= uc($line); $ilen = length($buffer) - $len;
        }
    }
}

sub cloudLocsToBed { my $fp_out = shift; my $locs = shift; my $print_ids = shift;
    foreach my $loc (@$locs) { locToBed($fp_out, $loc, $print_ids); }
}

sub locToBed { my $fp_out = shift; my $loc = shift; my $print_ids = shift;
    print $fp_out "$loc->{contig}\t$loc->{begin}\t$loc->{end}\tpcloud:$loc->{kinfo}->{cloud}->{id}\t0\t";
    my $num_plus = 0; my $num_minus = 0;
    foreach my $kmer (@{$loc->{kmers}}){ if ($kmer->{strand} eq "+") { $num_plus++; } else { $num_minus++; }; }
    if ($num_plus >= $num_minus ) { print $fp_out "+"; } else { print $fp_out "-"; }
    if ($print_ids && $loc->{kmers}) {
        foreach my $kmer (@{$loc->{kmers}}){ print $fp_out "\t$kmer->{cloud}->{id}"; }
    } print $fp_out "\n";
}

#### All of the subroutines below should probably be moved to globals.pm (or its successor) eventually.
#### However, I think it may be necessary to implement globals as an Object to do this without ugly hacks.

### sub-routines that call other subroutines should go here ###
sub readx { my ($file, $sub, $blurb, @args) = @_;
    if (my $fp = openup($file, $read, $blurb)) { $sub->($fp, @args); close $fp; }
    else { warn "Unable to read input file with readx!\n\n Aborting!!\n\n"; return -1; }
}

sub multireadx { my ($files, @args) = @_;
    # auto-detect single file, add it to the array
    if (!ref $files) { my $tempfile = $files; $files = array_ref(); push (@{$files}, $tempfile); }
    # re-write this to evaluate return values from readx and deal with them appropriately.
    foreach my $file (@{$files}){ readx($file, @args); }
}

sub readWritex { my ($flag, $sub, $infile, $outfile, $blurb, @args) = @_;
    my $fpin = openup($infile, $read, "input file for ".$blurb);
    my $fpout = \*STDOUT;
    if ($flag) { $fpout = openup($outfile, $write, "output file for ".$blurb); }
    if($fpin && $fpout) { $sub->($fpin, $fpout, @args); }
    else { print STDERR "Unable to read input/output file!\n\n Aborting!!\n\n"; die; }
    close $fpin; close $fpout; ### this could lead to closing STDOUT - which could have bad consequences
}

sub printx { my $target = shift; my $sub_ref = shift; my $blurb = shift; my @args = @_;
	my $file_handle; my $handle;
    if (!$target) { $handle = \*STDERR; }
    if (!$blurb) { $blurb = "unknown information"; }
    if ($target eq \*OUTFILE || $target eq \*STDERR || $target eq \*STDOUT) {
        talk("standard output ID: file pointer is $target\n", $murmur); $handle = $target;
    }
    else { $file_handle = openup($target, $write, "printx output file for ".$blurb); $handle = $file_handle; }
    if ($handle) { $sub_ref->($handle, @args); }
    if ($file_handle) { close($file_handle); }
}

#### read any number of files and write to same number or no output files -
sub m_readWritex { my ($infiles, $outfiles, $sub, $blurb, @args) = @_;
    if (!ref $infiles) { my $tempfile = $infiles; $infiles = array_ref(); push (@{$infiles}, $tempfile); }
    if (!ref $outfiles) { my $tempfile = $outfiles; $outfiles = array_ref(); push (@{$outfiles}, $tempfile); }
    my $return = hash_ref(); my $fpout;
    for (my $i = 0; $i < scalar @{$infiles}; $i++) { (my $name = $infiles->[$i]) =~ s/\.\/.*\///;
        my $fpin = openup($infiles->[$i], $read, "input file for ".$blurb);
        if (!$outfiles->[0]) { if ($fpin) { $return->{$name} = $sub->($fpin, @args); } else { $return->{$name} = -1; } next;}
        elsif (@$infiles == @$outfiles) { $fpout = openup($outfiles->[$i], $write, "output file for ".$blurb); }
        elsif (1 eq scalar @$outfiles) { if (!$fpout) { $fpout = openup($outfiles->[0], $write, "output file for ".$blurb); } }
        else { die "\n\nreadWritex called with multiple outfiles not equal to number of infiles\n\n"; }
        if ($fpin && $fpout) { $return->{$name} = $sub->($fpin, $fpout, @args); } else { $return->{$name} = -1; }
        close $fpin;
    } close $fpout;
}

##
##		Supporting functions: print
##

sub print_value { my $out = shift; my $value = shift; my $indent = shift; my $ptrs = shift;
    if (!$ptrs){ $ptrs = {}; } elsif ($ptrs->{$value}) { return; }
    my $val_ref = ref $value;
    if ($val_ref eq "") { print $out "$indent$value\n"; }
    elsif ($val_ref eq  "ARRAY") { $ptrs->{$value} = 1; print_array($out, $value, $indent." ", $ptrs); }
    elsif ($val_ref eq  "HASH") { $ptrs->{$value} = 1; print_hash($out, $value, $indent." ", $ptrs); }
}

sub print_hash { my $out = shift; my $hash = shift; my $indent = shift; my $ptrs = shift;
    while (my ($key, $value) = each %{$hash}) { print $out "$indent$key:\n"; print_value($out, $value, $indent." ", $ptrs); }
}

sub print_array { my $out = shift; my $array = shift; my $indent = shift; my $ptrs = shift;
    for (my $i = 0; $i <= $#{$array}; $i++) { print_value($out, $array->[$i], $indent, $ptrs); }
}

##
##		Basic subroutines needed for many apps
##

sub openup { my $file = shift; my $mode = shift; my $blurb = shift; my $no_record = shift; if (! $blurb) { $blurb = "unknown"; }
    open (my $file_handle, $mode, $file) or warn "Couldn't open file: $file. \n";
    if (!$no_record) { talk("Opened $file as $blurb.\n", $murmur); }
	return $file_handle;
}

sub getFiles { my $pre = shift; my $post = shift; my $folder = shift;
	talk("Searching for files in folder $folder starting with $pre and ending with $post.\n", $murmur);
	my @files = <./$folder/$pre*$post>;
	foreach my $file (@files){ talk("Source file $file detected.\n", $murmur); }
	return \@files;
}

sub scalar_ref { my $scalar; return \$scalar; }
sub array_ref { my @array; return \@array; }
sub hash_ref { my %hash; return \%hash; }
sub skip {}

##
##		Basic subroutines for program communication and recording
##

sub talk { my $phrase = shift; my $loudness = shift;
	if ($loudness <= $globals->{record_loudness}) { record($phrase, $globals->{record_run}); }
	if ($loudness <= $globals->{loudness}) {  print STDOUT $phrase; }
}

sub record { my $phrase = shift;
    my $recordfile = openup($globals->{record_run}, $append, "", $no_record);
    my $now = time(); my $delta = $now - $globals->{starttime};
	print $recordfile "$delta: $phrase"; close($recordfile);
}