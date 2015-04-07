# Replicate what the second half of PClouds does when it scans the entire
# sequence and tells whether the kmer at that position is in the set of kmers
# in a file.
# This looks for whole regions instead of individual sites

package Annotate;

use strict;
use warnings;

use Data::Dumper;

unless (caller) {
	
	AnnotateBigFasta("fasta", "", 50);
	exit;
	my $kmer_file = "kmers_above_threhold";
	my $sequence_file = "sequence";

	open my $fh, $sequence_file or die $!;
	my $sequence = <$fh>;
	chomp $sequence;

	my $kmer_assignments = ReadKmersFromFile($kmer_file);

	my $regions = AnnotateSequence($sequence, $kmer_assignments);

}

sub ReadKmersFromFile {
	my ($kmer_file) = @_;
	open my $fh, $kmer_file or die $!;
	<$fh>; # Burn header line

	my $kmers = {};
	while (<$fh>) {
		$kmers->{shift split $_} = 1;
	}
	return $kmers;
}

# Only immediately sequential kmers will be joined into regions
sub AnnotateSequence {
	my ($sequence, $kmer_assignments, $sequence_name, $initial_position) = @_;

	$initial_position = $initial_position // 0;

	my $regions = [];
	my $region_start = -1;
	my $previous_kmer_was_repetitive = 0;
	my $k = length((keys %$kmer_assignments)[0]);

	my $final_position = $initial_position + length($sequence) - $k;

	for my $position ($initial_position .. $final_position) {
		my $current_kmer_is_repetitive =
		  $kmer_assignments->{substr($sequence, $position, $k)};
		if ($current_kmer_is_repetitive) {
			print "Found kmer "
			  . substr($sequence, $position, $k)
			  . " at position $position";
			if(not $previous_kmer_was_repetitive) {
				$region_start = $position;
			}
		}
		else {

			# Current kmer is not repetitive
			if ($previous_kmer_was_repetitive) {

				# Ends are exclusive
				my $region_end = $position - 1 + $k;

				# Could print the region to BED here
				if ($sequence_name) {
					push @$regions,[$sequence_name, $region_start, $region_end];
				}
				else {
					push @$regions, [$region_start, $region_end];
				}
			}
			$previous_kmer_was_repetitive = 0;
		}
	}

	# If end of sequence is a repeat region
	if ($previous_kmer_was_repetitive) {

		# Ends are exclusive
		my $region_end = length($sequence);
		if ($sequence_name) {
			push @$regions,[$sequence_name, $region_start, $region_end];
		}
		else {
			push @$regions, [$region_start, $region_end];
		}
	}

	return $regions;
}

sub AnnotateFasta {
	my ($fasta_file, $kmer_assignments) = @_;
	
	my $regions = [];
	
	open my $fh, "<", $fasta_file or die $!;

	my $name = "";
	my $sequence = "";

	while (<$fh>) {
		chomp;
		if (s/>//) {
			if ($sequence) {
				push @$regions, @{AnnotateSequence($sequence, $kmer_assignments, $name)};
				$sequence = "";
			}
			$name = $_;
		}
		else { $sequence .= $_ }
	}

	# For the last sequence
	push @$regions, @{AnnotateSequence($sequence, $kmer_assignments, $name)};
	
	return $regions;
}

# Not complete. I'm thinking that we should just push this to C and not worry
# about it in perl for now...
sub AnnotateBigFasta {
	my ($big_fasta, $kmer_assignments, $chunk_size) = @_;
	
	open my $fh, "<", $big_fasta;
	
	my $size_read = read $fh, (my $chunk), $chunk_size;
	if ($size_read != $chunk_size) {
		warn "read less than chunk size";
		$chunk_size = $size_read;
	}
	print $chunk . "\n\n";
	
	my $previous_name_continuing = 0;
	my $previous_sequence_continuing = 0;
	my $previous_name = "";
	my $previous_sequence = "";
	
	my $sequence_name = "";
	
	if ($previous_name_continuing) {
		$chunk =~ /(.*)\n?/;
		$sequence_name = $previous_name . $1;
	}
	
	
	if ($chunk =~ /\n*>(.*)\n/) {
		$sequence_name = $1;
	}
	print "BEFORE: " . substr($chunk, 0, $-[0]) . "\n\n";
	print "MATCH: " . substr($chunk, $-[0], $+[0]-$-[0]) . "\n\n";
	
	print "AFTER: " . substr($chunk, $+[0]). "\n\n";
	print "\n\n\n";
	print $sequence_name;
		
	
		
}


sub RegionsToBed {
	my ($regions, $bed_file) = @_;

	open my $bed_fh, ">", $bed_file;
	if (@{$regions->[0]} == 2) {
		print $bed_fh "Start\tEnd\n";
	}
	elsif (@{$regions->[0]} == 3) {
		print $bed_fh "Chrom\tStart\tEnd\n";
	}
	else { die "Incorrect region format"; }
	for my $region (@$regions) {
		print $bed_fh join("\t", @$region) . "\n";
	}
}

1;
