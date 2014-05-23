#ifndef macrodefine_h__
#define macrodefine_h__


#include <vector>

#include <iostream>

typedef std::vector <bool> bitvector;


struct Kmer {
	unsigned long number_pattern;
	int cloud;
	Kmer() {
		number_pattern = 0;
		cloud = 0;
	}
};

struct CoreKmer {
	unsigned long number_pattern;
	int count;
	int cloud;
	bool has_been_extended;

	CoreKmer() {
		number_pattern = -1;
		count = 0;
		cloud = 0;
		has_been_extended = 0;
	}
};


/*definition of macros*/
#define MAX_LINE_LENGTH 6000		//maximum characters in each line of a file
#define MAX_FILENAME_LENGTH 5000		//maximum length of filename

#define MAXMAINMER 40000000		//maximum number of oligos in a cloud class
#define MAXACCMER 80000000		//maximum number of oligos in accessory cloud
#define MAXCLOUD 3000000		//maximum number of p cloud

#define ZERO 48				//ASCII code of character '0'
#define NINE 57				//ASCII code of character '9'

#endif //madrodefine_h__
