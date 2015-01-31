#ifndef stringhandle_h__
#define stringhandle_h__

#include <string>

namespace pclouds {

void getreversecomplement(const char* pchSource, char* pchTarget);
std::string getreversecomplement(std::string kmer);

bool is_SSR(const char *kmer);

bool is_one_SSR(const char *kmer);

bool is_two_SSR(const char *kmer);

bool is_three_SSR(const char *kmer);

bool is_four_SSR(const char *kmer);

bool issegmentvalid(char *pchPattern);

bool is_segment_valid(char *sequence_start, int size);

}

#endif // stringhandle_h__
