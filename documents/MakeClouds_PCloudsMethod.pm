# Make clouds from oligos and counts the same way PClouds does.
#
# This makes the clouds using the same steps as PClouds so the clouds should be
# identical.
#
# Currently, the core expands exactly the same way as PClouds. I have not
# checked that the outer parts of the clouds are expanded the same way.
#
# Form clouds using the P-clouds method.
# I believe this is right but the original paper leaves some ambiguity. I am
# not sure that the way of connecting the core oligos is correct. I'm not sure
# if the paper said that any core oligos that have a distance of 3 or less
# are connected to any core oligo. The other option is that the distance
# cutoff is determined by the count of the core oligo. I think that the range
# of expansion is dependent on the count of the core oligo considered. But the
# paper does not completely support this option.
#
# This method uses the method that the range of expansion is dependent on the
# count of the core oligo being expanded around.
#
#
# So it turns out that the paper says that the count of the core oligo should
# influence the range of other core oligos it can reach. The source code and
# testing of PClouds does not support this. Expanding the cores is completely
# independent of the primary, secondary, and tertiary thresholds. There are
# two rounds of core expansion and then the outer kmers are added after core
# expansion is complete.
#
# PClouds does take into account reverse complement. This representation of it
# does also as of Oct 5, 2013.
#
# Pclouds does not combine the counts of the sense sequence and its reverse
# complement. Should it?
# Read "How should we handle reverse complements.odt".
#
#
# Stephen Pollard


package MakeClouds_PCloudsMethod;

use strict;
use warnings;

use Data::Dumper;

open my $edges, '>', "edges";
open my $testmer_fh, '>', "testmers";

unless (caller) {
	print "Making clouds and kmer assignments using the P-clouds method\n";

	#my $kmer_file = "small_temp.fasta.processed.out_lumped";
	my $kmer_file = "small.out_lumped";

	#my $kmer_file = "kmer_counts";
	my $cloud_file = "clouds";
	my $kmer_assignment_file = "kmer_assignments";
	my $expansion_network_edges_file = "expansion_network_edges";


	print $edges "Source_core\tDestination_core\n";

	my $c_parameters = [2, 5, 10, 100, 1000];

	my $kmer_counts = ReadKmersAndCounts($kmer_file);
	my ($clouds, $kmer_assignments) = MakeClouds($kmer_counts, $c_parameters);

	PrintClouds($clouds, $cloud_file);
	PrintKmerAssignments($kmer_assignments, $kmer_assignment_file);

	print "Made " . scalar(keys %$clouds) . " clouds\n";
}

sub ReadKmersAndCounts {
	my ($kmer_file) = @_;
	open my $fh, $kmer_file or die $!;

	my $kmer_counts = {};

	<$fh>; # Burn header line
	while (<$fh>) {
		my ($kmer, $count) = split;
		$kmer_counts->{$kmer} = $count;
	}
	return $kmer_counts;
}

sub MakeClouds {
	my ($kmer_counts, $c_parameters) = @_;
	
	my $core_kmers = DetermineCoreKmers($kmer_counts, $c_parameters->[1]);

	my $kmer_assignments = {};
	my $clouds = {};
	my $expanded_cores = {};

	for my $seed_kmer (@$core_kmers) {
		if (not $kmer_assignments->{$seed_kmer}) {
			MakeNewCloudFromSeedKmer($seed_kmer, $kmer_counts, $clouds,
				$kmer_assignments, $c_parameters->[1], $expanded_cores);
		}
	}

	# Now expand around all core kmers that were assigned to clouds. Not all
	# core kmers might be assigned if using seed clouds.
	
	for my $core_kmer (@$core_kmers) {
		ExpandOuterRegionAroundCoreKmer($kmer_assignments->{$core_kmer},
			$core_kmer, $kmer_counts, $c_parameters, $kmer_assignments, $clouds);
	}

	return $clouds, $kmer_assignments;
}

sub DetermineCoreKmers {
	my ($kmer_counts, $core_threshold) = @_;
	my $core_kmers =
	  [grep {$kmer_counts->{$_} >= $core_threshold} keys %$kmer_counts];

	# sort the core kmers because order might matter with the PClouds method
	# for making clouds.
	@$core_kmers =
	  sort {$kmer_counts->{$b} <=> $kmer_counts->{$a}} @$core_kmers;
	print "There are " . scalar(@$core_kmers) . " core kmers\n";
	return $core_kmers;
}

sub MakeNewCloudFromSeedKmer {
	my ($seed_kmer, $kmer_counts, $clouds, $kmer_assignments, $core_threshold, $expanded_cores) =
	  @_;

	print "Making new cloud from $seed_kmer\n";
	my $cloud_id = $seed_kmer;

	# Add seed to cloud
	$clouds->{$cloud_id}->{$seed_kmer} = $kmer_counts->{$seed_kmer};
	$kmer_assignments->{$seed_kmer} = $cloud_id;

	ExpandCoreAroundSeed($cloud_id, $seed_kmer, $kmer_counts, $kmer_assignments,
		$clouds,$core_threshold, $expanded_cores);
}

# Most efficient place to handle reverse complements is to call this function
# on the sense strand and the reverse complement.
sub ExpandCoreAroundSeed{
	my ($cloud_id, $seed_kmer, $kmer_counts, $kmer_assignments, $clouds,
		$core_threshold, $expanded_cores)
	  = @_;

	# print "Expanding seed of cloud $cloud_id from $seed_kmer\n";
	

	# Find all testmers 1 around, 2 around, and 3 around the seed in that order
	my $testmers = MakeTestmers($seed_kmer);

	# First round of expansion
	foreach my $testmer (@$testmers) {

		# Here we require that the count be higher than the core threshold
		if (    $kmer_counts->{$testmer}
			and $kmer_counts->{$testmer} >= $core_threshold
			and not $kmer_assignments->{$testmer})
		{

			#			print "Adding $testmer to cloud $cloud_id from $seed_kmer\n";
			$kmer_assignments->{$testmer} = $cloud_id;
			$clouds->{$cloud_id}->{$testmer} = $kmer_counts->{$testmer};

			# For expansion network
			print $edges "$seed_kmer\t$testmer\n";

		}

		# Second round of expansion
		# Notice that we expand even if the core has already been assigned
		# This is one advantage to keeping track of which cores have been
		# expanded already. Also an advantage to calculating all the possible
		# testmers in 3 substitutions first.
		if (    $kmer_counts->{$testmer}
			and $kmer_counts->{$testmer} >= $core_threshold
			and not $expanded_cores->{$testmer})
		{
			$expanded_cores->{$testmer} = 1;

			my $testmers_2 = MakeTestmers($testmer);
			foreach my $testmer_2 (@$testmers_2) {
				if (    $kmer_counts->{$testmer_2}
					and $kmer_counts->{$testmer_2} >= $core_threshold
					and not $kmer_assignments->{$testmer_2})
				{

			   #					print "Adding $testmer_2 to cloud $cloud_id from $testmer\n";
					$kmer_assignments->{$testmer_2} = $cloud_id;
					$clouds->{$cloud_id}->{$testmer_2} =
					  $kmer_counts->{$testmer_2};

					# For expansion network
					print $edges "$testmer\t$testmer_2\n";
				}
			}
		}
	}
}

sub MakeTestmers {
	my ($core_kmer) = @_;
	
	my $testmers  = [];
	# First add testmers with distance 1, then 2, then 3
	for my $position (0 .. (length($core_kmer) - 1)) {
		for my $nuc ("A", "C", "G", "T") {
			if ($nuc ne substr($core_kmer, $position, 1)){
				my $testmer = $core_kmer;
				substr($testmer, $position, 1) = $nuc;
				push @$testmers, $testmer;
				print $testmer_fh $testmer . "\n";

			}
		}
	}

	# Now add the reverse complement
	my $reverse_core_kmer = ReverseComplement($core_kmer);

	for my $position (0 .. length($reverse_core_kmer) - 1) {
		for my $nuc ("A", "C", "G", "T") {
			if ($nuc ne substr($reverse_core_kmer, $position, 1)){
				my $testmer = $reverse_core_kmer;
				substr($testmer, $position, 1) = $nuc;
				push @$testmers, $testmer;
				print $testmer_fh $testmer . "\n";
			}
		}
	}

	# Now add two away from core and reverse complement

	for my $position_1 (0 .. length($core_kmer) - 1 - 1) {
		for my $position_2 (($position_1 + 1) .. length($core_kmer) - 1){
			for my $nuc_1 ("A", "C", "G", "T") {
				for my $nuc_2 ("A", "C", "G", "T") {
					if (    $nuc_1 ne substr($core_kmer, $position_1, 1)
						and $nuc_2 ne substr($core_kmer, $position_2, 1))
					{
						my $testmer = $core_kmer;
						substr($testmer, $position_1, 1) = $nuc_1;
						substr($testmer, $position_2, 1) = $nuc_2;
						push @$testmers, $testmer;
						print $testmer_fh $testmer . "\n";
					}
				}
			}
		}
	}

	for my $position_1 (0 .. length($reverse_core_kmer) - 1 - 1) {
		for my $position_2 (($position_1 + 1) .. length($reverse_core_kmer) - 1)
		{
			for my $nuc_1 ("A", "C", "G", "T") {
				for my $nuc_2 ("A", "C", "G", "T") {
					if (    $nuc_1 ne substr($reverse_core_kmer, $position_1, 1)
						and $nuc_2 ne substr($reverse_core_kmer, $position_2, 1)
					  )
					{
						my $testmer = $reverse_core_kmer;
						substr($testmer, $position_1, 1) = $nuc_1;
						substr($testmer, $position_2, 1) = $nuc_2;
						push @$testmers, $testmer;
						print $testmer_fh $testmer . "\n";
					}
				}
			}
		}
	}

	# Add three away substitutions

	for my $position_1 (0 .. length($core_kmer) - 1 - 2) {
		for my $position_2 (($position_1 + 1) .. length($core_kmer) - 1 - 1){
			for my $position_3 (($position_2 + 1) .. length($core_kmer) - 1){
				for my $nuc_1 ("A", "C", "G", "T") {
					for my $nuc_2 ("A", "C", "G", "T") {
						for my $nuc_3 ("A", "C", "G", "T") {
							if (    $nuc_1 ne substr($core_kmer, $position_1, 1)
								and $nuc_2 ne substr($core_kmer, $position_2, 1)
								and $nuc_3 ne substr($core_kmer, $position_3, 1)
							  )
							{
								my $testmer = $core_kmer;
								substr($testmer, $position_1, 1) =$nuc_1;
								substr($testmer, $position_2, 1) =$nuc_2;
								substr($testmer, $position_3, 1) =$nuc_3;
								push @$testmers, $testmer;
								print $testmer_fh $testmer . "\n";
							}
						}
					}
				}
			}
		}
	}

	# And 3 away from reverse complement

	for my $position_1 (0 .. length($reverse_core_kmer) - 1 - 2) {
		for my $position_2 (
			($position_1 + 1) .. length($reverse_core_kmer) - 1 - 1)
		{
			for my $position_3 (
				($position_2 + 1) .. length($reverse_core_kmer) - 1)
			{
				for my $nuc_1 ("A", "C", "G", "T") {
					for my $nuc_2 ("A", "C", "G", "T") {
						for my $nuc_3 ("A", "C", "G", "T") {
							if ($nuc_1 ne
								substr($reverse_core_kmer, $position_1,1)
								and $nuc_2 ne
								substr($reverse_core_kmer, $position_2,1)
								and $nuc_3 ne
								substr($reverse_core_kmer, $position_3,1))
							{
								my $testmer = $reverse_core_kmer;
								substr($testmer, $position_1, 1) =$nuc_1;
								substr($testmer, $position_2, 1) =$nuc_2;
								substr($testmer, $position_3, 1) =$nuc_3;
								push @$testmers, $testmer;
								print $testmer_fh $testmer . "\n";
							}
						}
					}
				}
			}
		}
	}
	return $testmers;
}

# This function makes testmers in reverse order as the above function
# any takes into account the frequency count of the core
sub MakeTestmersForOuterExpansion {
	my ($core_kmer, $kmer_counts, $c_parameters) = @_;

	my $testmers = [];
	
	# Now add the reverse complement
	my $reverse_core_kmer = ReverseComplement($core_kmer);


	# First add testmers with distance 3, then 2, then 1 as high
	if ($kmer_counts->{$core_kmer} >= $c_parameters->[4]) {

		# Add three away substitutions
		for my $position_1 (0 .. length($core_kmer) - 1 - 2) {
			for my $position_2 (($position_1 + 1) .. length($core_kmer) - 1 - 1)
			{
				for my $position_3 (($position_2 + 1) .. length($core_kmer) - 1)
				{
					for my $nuc_1 ("A", "C", "G", "T") {
						for my $nuc_2 ("A", "C", "G", "T") {
							for my $nuc_3 ("A", "C", "G", "T") {
								if ($nuc_1 ne substr($core_kmer, $position_1, 1)
									and $nuc_2 ne
									substr($core_kmer, $position_2, 1)
									and $nuc_3 ne
									substr($core_kmer, $position_3, 1))
								{
									my $testmer = $core_kmer;
									substr($testmer, $position_1, 1) =$nuc_1;
									substr($testmer, $position_2, 1) =$nuc_2;
									substr($testmer, $position_3, 1) =$nuc_3;
									push @$testmers, $testmer;
									print $testmer_fh $testmer . "\n";
								}
							}
						}
					}
				}
			}
		}

		# And 3 away from reverse complement

		for my $position_1 (0 .. length($reverse_core_kmer) - 1 - 2) {
			for my $position_2 (
				($position_1 + 1) .. length($reverse_core_kmer) - 1 - 1)
			{
				for my $position_3 (
					($position_2 + 1) .. length($reverse_core_kmer) - 1)
				{
					for my $nuc_1 ("A", "C", "G", "T") {
						for my $nuc_2 ("A", "C", "G", "T") {
							for my $nuc_3 ("A", "C", "G", "T") {
								if ($nuc_1 ne
									substr($reverse_core_kmer, $position_1,1)
									and $nuc_2 ne
									substr($reverse_core_kmer, $position_2,1)
									and $nuc_3 ne
									substr($reverse_core_kmer, $position_3,1))
								{
									my $testmer = $reverse_core_kmer;
									substr($testmer, $position_1, 1) =$nuc_1;
									substr($testmer, $position_2, 1) =$nuc_2;
									substr($testmer, $position_3, 1) =$nuc_3;
									push @$testmers, $testmer;
									print $testmer_fh $testmer . "\n";
								}
							}
						}
					}
				}
			}
		}
	}

	if ($kmer_counts->{$core_kmer} >= $c_parameters->[3]) {

		# Now add two away from core and reverse complement

		for my $position_1 (0 .. length($core_kmer) - 1 - 1) {
			for my $position_2 (($position_1 + 1) .. length($core_kmer) - 1){
				for my $nuc_1 ("A", "C", "G", "T") {
					for my $nuc_2 ("A", "C", "G", "T") {
						if (    $nuc_1 ne substr($core_kmer, $position_1, 1)
							and $nuc_2 ne substr($core_kmer, $position_2, 1))
						{
							my $testmer = $core_kmer;
							substr($testmer, $position_1, 1) = $nuc_1;
							substr($testmer, $position_2, 1) = $nuc_2;
							push @$testmers, $testmer;
							print $testmer_fh $testmer . "\n";
						}
					}
				}
			}
		}

		for my $position_1 (0 .. length($reverse_core_kmer) - 1 - 1) {
			for my $position_2 (
				($position_1 + 1) .. length($reverse_core_kmer) - 1)
			{
				for my $nuc_1 ("A", "C", "G", "T") {
					for my $nuc_2 ("A", "C", "G", "T") {
						if ($nuc_1 ne substr($reverse_core_kmer, $position_1, 1)
							and $nuc_2 ne
							substr($reverse_core_kmer, $position_2, 1))
						{
							my $testmer = $reverse_core_kmer;
							substr($testmer, $position_1, 1) = $nuc_1;
							substr($testmer, $position_2, 1) = $nuc_2;
							push @$testmers, $testmer;
							print $testmer_fh $testmer . "\n";
						}
					}
				}
			}
		}
	}

	# And finally add 1 away testmers
	if ($kmer_counts->{$core_kmer} >= $c_parameters->[2]) {
		for my $position (0 .. (length($core_kmer) - 1)) {
			for my $nuc ("A", "C", "G", "T") {
				if ($nuc ne substr($core_kmer, $position, 1)){
					my $testmer = $core_kmer;
					substr($testmer, $position, 1) = $nuc;
					push @$testmers, $testmer;
					print $testmer_fh $testmer . "\n";

				}
			}
		}

		for my $position (0 .. length($reverse_core_kmer) - 1) {
			for my $nuc ("A", "C", "G", "T") {
				if ($nuc ne substr($reverse_core_kmer, $position, 1)){
					my $testmer = $reverse_core_kmer;
					substr($testmer, $position, 1) = $nuc;
					push @$testmers, $testmer;
					print $testmer_fh $testmer . "\n";
				}
			}
		}
	}
	return $testmers;
}

# Most efficient place to handle reverse complements is to call this function
# on the sense strand and the reverse complement.
sub ExpandOuterRegionAroundCoreKmer {
	my ($cloud_id, $core_kmer, $kmer_counts, $c_parameters, $kmer_assignments, $clouds) =
	  @_;

	my $testmers = MakeTestmersForOuterExpansion($core_kmer, $kmer_counts,
		$c_parameters);

	for my $testmer (@$testmers){

		# Here we require that the count be higher than the OUTER threshold
		if (    $kmer_counts->{$testmer}
			and $kmer_counts->{$testmer} >= $c_parameters->[0]
			and not $kmer_assignments->{$testmer})
		{

			#			print "Adding $testmer to cloud $cloud_id from $seed_kmer\n";
			$kmer_assignments->{$testmer} = $cloud_id;
			$clouds->{$cloud_id}->{$testmer} = $kmer_counts->{$testmer};

			# For expansion network
			print $edges "$core_kmer\t$testmer\n";
		}
	}
}

sub ReverseComplement {
	my ($kmer) = @_;
	$kmer = reverse $kmer;
	$kmer =~ tr/ACGTacgt/TGCAtgca/;
	return $kmer;
}

sub PrintClouds {
	my ($clouds, $cloud_file) = @_;
	open my $fh, '>', $cloud_file or die $!;

	print "Printing clouds to \"$cloud_file\"\n";

	print $fh "Clouds\n";

	while (my ($cloud_id, $cloud) = each %$clouds) {

		# Print cloud id first
		# Can print level
		#		print $fh "$cloud_id\tcore";
		print $fh "$cloud_id\n";

		# What is the general term for 'core' and 'outer'. I'm using 'level'.
		while (my ($kmer, $level) = each %$cloud) {
			next if $kmer eq $cloud_id;

			# Can print level
			#			print $fh "$kmer\t$level";
			print $fh "$kmer\n";
		}

		# New line indicates the end of this cloud
		print $fh "\n";
	}
}

sub PrintKmerAssignments {
	my ($kmer_assignments, $kmer_assignment_file) = @_;
	open my $fh, '>', $kmer_assignment_file or die $!;

	print"Printing kmer assignments to \"$kmer_assignment_file\"\n";
	print $fh "Kmer\tCloud_id\n";
	while (my ($kmer, $cloud_id) = each %$kmer_assignments) {
		print $fh "$kmer\t$cloud_id\n";
	}
}

sub MakeCloudFromSeedmer {
	my ($seed_kmer, $kmer_counts, $c_parameters) = @_;

	my $kmer_assignments = {};
	my $clouds = {};
	my $expanded_cores = {};

	MakeNewCloudFromSeedKmer($seed_kmer, $kmer_counts, $clouds,
		$kmer_assignments, $c_parameters->[1], $expanded_cores);

	ExpandOuterRegionAroundCoreKmer($kmer_assignments->{$seed_kmer},
		$seed_kmer, $kmer_counts, $c_parameters, $kmer_assignments, $clouds);

	return [keys $kmer_assignments];
}

1;