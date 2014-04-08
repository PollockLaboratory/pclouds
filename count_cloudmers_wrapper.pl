# Count the sequence of 10 clouds in a row


use strict; 
use warnings;
use autodie;


local $\ = "\n";
local $, = " ";

open my $fh, ">cloudmer_counts";
close $fh;

my $input = "hg_tester.region";
my $kmer_length = 15; 
my $file_out = "cloudmer_counts";

for my $region_length (23 .. 24) {
	system "perl count_cloudmers.pl $input $kmer_length $region_length $file_out &";	
}