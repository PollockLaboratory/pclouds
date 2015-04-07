# Make clouds from oligos and counts
# 
# The cloud structure in this program is not lab standard. It is minimal. 
# 
# This uses Jaime's method for forming clouds.
#
# It does not take into account reverse complement.  
#
# Stephen Pollard

use strict;
use warnings;

use Data::Dumper;

local $\ = "\n";

print "Making clouds and kmer assignments";


my $kmer_file = "kmer_counts";
my $cloud_out_file = "clouds";
my $kmer_assignment_out_file = "kmer_assignments";
my $core_threshold = $ARGV[0] // 5;
my $outer_threshold = $ARGV[1] // 4;
my $seed_file = $ARGV[2];

my $edges_file = "cloud_building_jumps";
open my $edges, '>', $edges_file;


die "Provide core threshold and outer threshold" if not $core_threshold or not $outer_threshold;




my $kmers  = {};
my $kmer_assignments = {};
my $clouds = {};
my $k;

ReadKmersAndCounts();
MakeClouds();
PrintClouds();
PrintKmerAssignments();


print "Made " . scalar(keys %$clouds) . " clouds";

sub ReadKmersAndCounts {
	open my $fh, $kmer_file or die $!;

	<$fh>; # Burn header line
	while (<$fh>) {
		my ($kmer, $count) = split;
		$kmers->{$kmer} = $count;
		$k = length($kmer) if not $k;
	}
}

sub MakeClouds {
	my $core_kmers = GetCoreKmers();
	
	my $seed_kmers = $core_kmers; 
	
	foreach my $seed_kmer (@$seed_kmers) {
		if (not $kmer_assignments->{$seed_kmer}) {
			my $cloud_id = $seed_kmer;

			$kmer_assignments->{$seed_kmer} = $cloud_id;
			$clouds->{$cloud_id} = {$seed_kmer => $kmers->{$seed_kmer}};
			ExpandCloudAroundKmer($cloud_id, $seed_kmer);
		}
	}
}

sub GetCoreKmers {
	my $core_kmers = [grep {$kmers->{$_} >= $core_threshold} keys %$kmers];
	return [ sort {$kmers->{$b} <=> $kmers->{$a}} @$core_kmers];	
}

sub ExpandCloudAroundKmer {
	my ($cloud_id, $core_kmer) = @_;
	print "Expanding cloud $cloud_id around $core_kmer";
	for my $position (0 .. length($core_kmer) - 1) {
		for my $nuc ("A", "C", "G", "T") {
			my $testmer = $core_kmer;
			substr($testmer, $position, 1) = $nuc;

			# Have already removed all kmers that aren't in outer or core
			if ($kmers->{$testmer} and not $kmer_assignments->{$testmer}) {
				# Keep track of edge from core kmer to testmer
				print $edges "$core_kmer\t$testmer";
				print "assigning $testmer to $cloud_id";
				
				$kmer_assignments->{$testmer} = $cloud_id;
				$clouds->{$cloud_id}->{$testmer} = $kmers->{$testmer};
				if ($kmers->{$testmer} >= $core_threshold) {
					ExpandCloudAroundKmer($cloud_id, $testmer);
				}
			}
		}
	}
}

sub PrintClouds {
	open my $fh, '>', $cloud_out_file or die $!;
	
	print "Printing clouds to \"$cloud_out_file\"";
	
	while (my ($cloud_id, $cloud) = each %$clouds) {
		# Print cloud id first
		# Can print level
#		print $fh "$cloud_id\tcore";
		print $fh "$cloud_id";
		# What is the general term for 'core' and 'outer'. I'm using 'level'.
		while (my ($kmer, $level) = each %$cloud) {
			next if $kmer eq $cloud_id;
			# Can print level
#			print $fh "$kmer\t$level";
			print $fh "$kmer";
		}
		# New line indicates the end of this cloud
		print $fh "";
	}
}

sub PrintKmerAssignments {
	open my $fh, '>', $kmer_assignment_out_file or die $!;
	
	print "Printing kmer assignments to \"$kmer_assignment_out_file\"";
	print $fh "Kmer\tCloud_id";
	while (my ($kmer, $cloud_id) = each %$kmer_assignments) {
		print $fh "$kmer\t$cloud_id";
	}
}
