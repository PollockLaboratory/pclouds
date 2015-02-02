#include <string>
#include <vector>

#include "build_clouds.hpp"
#include "../include/macrodefine.h"

namespace pclouds {

template <typename t>
int greater_count (t a, t b) {
	return (a.count > b.count);
}

template <typename t>
bool lesser_pattern (t a, t b) {
	return (a.number_pattern < b.number_pattern);
}

template <typename t>
int search(int number_of_kmers, const std::vector <t> ts, const unsigned long key) {
	int m;
	int site = 0;

	while (number_of_kmers >= 1) {
		m = number_of_kmers / 2;
		if (ts[site + m].number_pattern == key)
			return (site + m);
		else if (ts[site + m].number_pattern > key)
			number_of_kmers = m;
		else {
			site = site + m + 1;
			number_of_kmers = number_of_kmers - m - 1;
		}
	}

	return (-1);
}

}
