/**
 * STP: This is clearly a cpp file. It used to be named *.c.
 *
 */

#include "../stringhandler/stringhandle.h"
#include "../include/macrodefine.h"
#include <cstring>

#include <string>

int strcompare(int i, char *source, char *target) {
	int count;
	for (count = 0; count < i; count++)
		if (*(source + count) != *(target + count))
			return 0;

	return 1;
}

void strright(int i, char *source, char *temp) {
	int count, j;

	j = i - 1;

	int size = strlen(source);

	for (count = size - 1; count > size - i - 1; count--) {
		*(temp + j) = *(source + count);
		j--;
	}

	*(temp + i) = '\0';
}

void strleft(int i, char* source, char* target) {
	int count, j;

	for (count = 0; count < i; count++)
		*(target + count) = *(source + count);

	*(target + i) = '\0';
}

int strreplaceright(char *source, char *target, char *replace) {
	int i;
	for (i = 0; i < (signed) strlen(source) - (signed) strlen(replace); i++)
		*(target + i) = *(source + i);
	for (i = (signed) strlen(source) - (signed) strlen(replace);
			i < (signed) strlen(source); i++)
		*(target + i) = *(replace + i - (signed) strlen(source)
				+ (signed) strlen(replace));
	return 1;
}

int straddright(char *source, char *target, char *add) {
	int i;
	for (i = 0; i < (signed) strlen(source); i++)
		*(target + i) = *(source + i);
	for (i = (signed) strlen(source);
			i < (signed) strlen(source) + (signed) strlen(add); i++)
		*(target + i) = *(add + i - (signed) strlen(source));
	return 1;
}

int stringinitial(char *pchSource, int iLength) {
	memset(pchSource, ' ', iLength);

	*(pchSource + iLength) = '\0';

	return (1);
}

void getsubstring(char* pchSource, char *pchTarget, int nStart, int nEnd) {
	int count;

	for (count = nStart; count <= nEnd; count++)
		*(pchTarget + count - nStart) = *(pchSource + count);
}

void getreversecomplement(const char* pchSource, char* pchTarget) {
	int count;

	int size = strlen(pchSource);

	for (count = 0; count < strlen(pchSource); count++) {
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
	for (int site = 0; site < kmer.length(); site++) {
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

int stringcompare(char *pchSource, char *pchTarget) {
	int count = 0;

	int nCompared = 1;

	for (count = 0; count < (signed) strlen(pchSource);) {
		if ((*(pchSource + count) == *(pchTarget + count))
				|| (*(pchSource + count) + 32 == *(pchTarget + count))
				|| (*(pchSource + count) - 32 == *(pchTarget + count)))
			count++;
		else
			return (1);
	}

	return (0);
}

int stringcompare(const char *pchSource, const char *pchTarget) {
	int count = 0;

	int nCompared = 1;

	for (count = 0; count < (signed) strlen(pchSource);) {
		if ((*(pchSource + count) == *(pchTarget + count))
				|| (*(pchSource + count) + 32 == *(pchTarget + count))
				|| (*(pchSource + count) - 32 == *(pchTarget + count)))
			count++;
		else
			return (1);
	}

	return (0);
}

int stringcopy(char *pchSource, char *pchTarget, int iLength) {
	int count;

	stringinitial(pchTarget, iLength);

	for (count = 0; count < iLength; count++)
		*(pchTarget + count) = *(pchSource + count);

	return (1);

}

int stringtonumber(const std::string& m_strTemp) {
	int temp = 0;

	int length = m_strTemp.length();

	for (int i = 0; i < length; i++)
		temp = temp * 10 + m_strTemp[i] - ZERO;

	return (temp);
}

unsigned int stringtolargenumber(const std::string& m_strTemp) {
	unsigned int temp = 0;

	int length = m_strTemp.length();

	for (int i = 0; i < length; i++)
		temp = temp * 10 + (unsigned) (m_strTemp[i] - ZERO);

	return (temp);
}

void stringtoarray(const std::string& src, char* target) {
	int i;

	int length = src.length();

	for (i = 0; i < length; i++)
		target[i] = src[i];

	target[i] = '\0';
}

bool is_SSR(const char *kmer) {
	return is_one_SSR(kmer) or is_two_SSR(kmer) or is_three_SSR(kmer)
			or is_four_SSR(kmer);
}

bool is_one_SSR(const char *kmer) {
	int i;

	for (i = 0; i < strlen(kmer) - 1; i++) {
		if (kmer[i + 1] != kmer[i])
			return (false);
	}

	return (true);
}

bool is_two_SSR(const char *kmer) {
	int i;

	for (i = 0; i < strlen(kmer) - 2; i++) {
		if (kmer[i + 2] != kmer[i])
			return (false);
	}

	return (true);
}

bool is_three_SSR(const char *kmer) {
	int i;

	for (i = 0; i < strlen(kmer) - 3; i++) {
		if (kmer[i + 3] != kmer[i])
			return (false);
	}

	return (true);
}

bool is_four_SSR(const char *kmer) {
	int i;

	for (i = 0; i < strlen(kmer) - 4; i++) {
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

