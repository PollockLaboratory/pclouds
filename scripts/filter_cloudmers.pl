# Counts the number of times cloudmers in regions of a specific length occur
# in all regions. Prints the total number and the region length to a file.

use strict; 
use warnings;
use autodie;


use Data::Dumper;

local $\ = "\n";
local $, = " ";


open my $fp_in, "<", "C:/Users/Stephen/Desktop/hg19_nomit.fasta.processed.C10.region";
<$fp_in>;

my $cloudmers = {};

open my $patterns, ">patterns";
open my $regions, ">regions"; 

while(<$fp_in>) {
	print $. if $. % 10000 == 0; 
	next if not $_;

	my $fields = [split];
	my $start = shift @$fields;
	my $end = shift @$fields;

	if (@$fields < 35) {
		print $patterns @$fields;
	}
	else {
		print $regions @$fields;
	}
}
