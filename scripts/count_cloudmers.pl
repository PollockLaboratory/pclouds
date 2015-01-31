# Counts the number of times cloudmers in regions of a specific length occur
# in all regions. Prints the total number and the region length to a file.
# 
# 
# Usage
# perl count_cloudmers.pl <.region file> <kmer_length> <region_length> <output_file> 
# perl count_cloudmers.pl sequence.region 16 24 cloudmer_counts

use strict; 
use warnings;
use autodie;


use Data::Dumper;

local $\ = "\n";
local $, = " ";

my $regions = $ARGV[0] || "hg_tester.region";
my $kmer_length = $ARGV[1] || 14;
my $region_length = $ARGV[2] || 24;
my $counts_out = $ARGV[3] || "cloudmer_counts_by_length";

my $cloudmer_sequence_length = $region_length - $kmer_length;


my $time = time;

open my $fp_in, "<", $regions;
<$fp_in>;

my $cloudmers = {};

while(<$fp_in>) {
	next if not $_;

	my $fields = [split];
	my $start = shift @$fields;
	my $end = shift @$fields;

	if ($end - $start == $region_length) {
		my $clouds = $fields;
		my $cloudmer = join "_", @$clouds;
		$cloudmers->{$cloudmer} = 0;
	}
}

# Reset the file pointer
open $fp_in, "<", $regions;
<$fp_in>;

my $total_cloudmers = 0;
while(<$fp_in>) {
	my $fields = [split];
	# Remove the start and end numbers
	shift @$fields; shift @$fields;

	my $region_clouds = $fields;
	for my $position (0 .. (scalar(@$region_clouds) - $cloudmer_sequence_length)) {
		my $clouds = [@$region_clouds[$position..($position + $cloudmer_sequence_length - 1)]];
		my $cloudmer = join "_", @$clouds;
		if (exists $cloudmers->{$cloudmer}) {
			$total_cloudmers++;
		}
	} 
}

open my $fp_out, ">>", $counts_out;
print "$region_length\t$total_cloudmers";
print $fp_out "$region_length\t$total_cloudmers";
print time - $time;
