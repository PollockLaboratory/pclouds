#pragma once

#include <string>
#include <vector>

namespace pclouds {

struct Kmer {
	unsigned long number_pattern;
	int cloud;
};

struct CoreKmer : public Kmer {
	int count;
	bool has_been_extended;
};

// transform the pattern sequence to the index of the array
// coding method: A=0, C=1, G=2, T=3,
void kmer_sequence_to_number_pattern(const char *kmer_sequence, unsigned long& number_pattern, const int& kmer_size);

}
