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

bool greater_count (const CoreKmer& a, const CoreKmer& b) {
	return (a.count > b.count);
}

bool lesser_pattern (const CoreKmer& a, const CoreKmer& b) {
	return (a.number_pattern < b.number_pattern);
}

bool lesser_pattern (const Kmer& a, const Kmer& b) {
	return (a.number_pattern < b.number_pattern);
}


// transform the pattern sequence to the index of the array
// coding method: A=0, C=1, G=2, T=3,
void kmer_sequence_to_number_pattern(const char *kmer_sequence, unsigned long& number_pattern, const int& kmer_size);

}
