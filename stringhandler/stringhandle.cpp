/**
 * STP: This is clearly a cpp file. It used to be named *.c.
 *
 */

#include "../stringhandler/stringhandle.h"
#include "../include/macrodefine.h"
#include <cstring>

#include <string>

namespace pclouds {

void getreversecomplement(const char* pchSource, char* pchTarget) {
	int count;

	int size = strlen(pchSource);

	for (count = 0; count < static_cast<int>(strlen(pchSource)); count++) {
		if (*(pchSource + count) == 'A' || *(pchSource + count) == 'a')
			*(pchTarget + size - 1 - count) = 'T';
		else if (*(pchSource + count) == 'T' || *(pchSource + count) == 't')
			*(pchTarget + size - 1 - count) = 'A';
		else if (*(pchSource + count) == 'C' || *(pchSource + count) == 'c')
			*(pchTarget + size - 1 - count) = 'G';
		else if (*(pchSource + count) == 'G' || *(pchSource + count) == 'g')
			*(pchTarget + size - 1 - count) = 'C';
		else
			*(pchTarget + size - 1 - count) = 'N';
	}

	*(pchTarget + size) = '\0';
}

std::string getreversecomplement(std::string kmer) {
	std::string complement = "";
	for (int site = 0; site < static_cast<int>(kmer.length()); site++) {
		char nucleotide = kmer[site];
		if (nucleotide == 'A' or nucleotide == 'a')
			complement += 'T';
		else if (nucleotide == 'C' or nucleotide == 'c')
			complement += 'G';
		else if (nucleotide == 'G' or nucleotide == 'g')
			complement += 'C';
		else if (nucleotide == 'T' or nucleotide == 't')
			complement += 'A';
		else
			complement += 'N';
	}
	std::string reverse_complement = std::string(complement.rbegin(),
			complement.rend());
	return reverse_complement;
}

bool is_SSR(const char *kmer) {
	return is_one_SSR(kmer) or is_two_SSR(kmer) or is_three_SSR(kmer)
		or is_four_SSR(kmer);
}

bool is_one_SSR(const char *kmer) {
	int i;

	for (i = 0; i < static_cast<int>(strlen(kmer)) - 1; i++) {
		if (kmer[i + 1] != kmer[i])
			return (false);
	}

	return (true);
}

bool is_two_SSR(const char *kmer) {
	int i;

	for (i = 0; i < static_cast<int>(strlen(kmer)) - 2; i++) {
		if (kmer[i + 2] != kmer[i])
			return (false);
	}

	return (true);
}

bool is_three_SSR(const char *kmer) {
	int i;

	for (i = 0; i < static_cast<int>(strlen(kmer)) - 3; i++) {
		if (kmer[i + 3] != kmer[i])
			return (false);
	}

	return (true);
}

bool is_four_SSR(const char *kmer) {
	int i;

	for (i = 0; i < static_cast<int>(strlen(kmer)) - 4; i++) {
		if (kmer[i + 4] != kmer[i])
			return (false);
	}

	return (true);
}

bool issegmentvalid(char *pchPattern) {
	int count;

	//STP: Passing in size will be faster. Perhaps a global would work better.
	int size = strlen(pchPattern);
	if (size == 0) {
		return (0);
	}
	char temp;

	for (count = 0; count < size; count++) {
		temp = *(pchPattern + count);

		if (temp != 'A' && temp != 'a' && temp != 'C' && temp != 'c'
				&& temp != 'T' && temp != 't' && temp != 'G' && temp != 'g')
			return (0);
	}

	return (1);
}

bool is_segment_valid(char *sequence_start, int size) {
	bool is_valid = true;
	for (int site = 0; site < size; site++) {
		char character = sequence_start[site];

		if (character != 'A' && character != 'a' && character != 'C'
				&& character != 'c' && character != 'T' && character != 't'
				&& character != 'G' && character != 'g') {
			is_valid = false;
			break;
		}
	}

	return is_valid;
}

}
