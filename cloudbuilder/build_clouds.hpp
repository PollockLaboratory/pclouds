#pragma once

#include <string>
#include <vector>

namespace pclouds {

struct Kmer {
	unsigned long number_pattern = 0;
	int cloud = 0;
};

struct CoreKmer {
	unsigned long number_pattern;
	int count = 0;
	int cloud = 0;
	bool has_been_extended = 0;
};

// transform the pattern sequence to the index of the array
// coding method: A=0, C=1, G=2, T=3,
void kmer_sequence_to_number_pattern(const char *kmer_sequence, unsigned long& number_pattern, const int& kmer_size);
void number_pattern_to_kmer_sequence(char *pchPattern, const unsigned long& index, const int& patternsize);
void get_one_substitutions(const char* core_kmer, unsigned long* array_of_testmer_number_patterns, int& number_of_testmers);
void get_two_substitutions(const char* core_kmer, unsigned long* array_of_testmer_number_patterns, int& number_of_testmers);
void get_three_substitutions(char* core_kmer, unsigned long* array_of_testmer_number_patterns, int& number_of_testmers);
void count_core_and_outer_kmers(std::string kmer_counts_file, int& number_of_core_kmers, int& number_of_outer_kmers, const int core_threshold, const int outer_threshold);
void read_core_kmers(std::string kmer_counts_file, CoreKmer* core_kmers, int& number_of_core_kmers, const int kmer_size, const int core_threshold);
void build_cloud_cores(CoreKmer* core_kmers, std::string cloud_summary_file, const int& number_of_core_kmers, const int& kmer_size, const int& core_threshold);
void output_cloud_cores(CoreKmer* core_kmers, std::string core_kmers_assign_file, const int& number_of_core_kmers, const int& size);
void read_cloud_cores(std::string core_assign_file, Kmer* tertiary_cores, int& number_of_tertiary_cores, Kmer* secondary_cores, int& number_of_secondary_cores, Kmer* primary_cores, int& number_of_primary_cores, const int& size, const int& primary_threshold, const int& secondary_threshold, const int& tertiary_threshold);
void build_cloud_outer(std::string kmer_counts_file, std::string outer_kmers_assign_file, Kmer* core_kmers_above_tertiary, const int number_of_cores_above_tertiary, Kmer* core_kmers_above_secondary, const int number_of_cores_above_secondary, Kmer* core_kmers_above_primary, const int number_of_cores_above_primary, const int kmer_size, const int core_threshold, const int outer_threshold);
void build_clouds(std::string controlfile);


}
