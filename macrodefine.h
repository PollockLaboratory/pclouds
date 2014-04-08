#ifndef macrodefine_h__
#define macrodefine_h__

/* define hash type and bit matrix type*/
#include <vector>
#include <ext/hash_map>

#include <iostream>

using namespace std;

using std::vector;
using __gnu_cxx::hash;
//using __gnu_cxx::hashtable;
using __gnu_cxx::hash_map;
//using __gnu_cxx::hash_multimap;

struct equalint
{
	bool operator()(unsigned long i1, unsigned long i2)
	{
		return i1 == i2;
	}
};

typedef hash_map<unsigned long, unsigned long> patternhash_type;
typedef vector<int> intvector;
typedef vector<intvector> intmatrix;
typedef vector <bool> bitvector;
typedef vector<bitvector> bit_matrix;

typedef struct _Pcloud
{
	unsigned long index;
	int number;
} cloud_type;

typedef struct _Pcloud1
{
	unsigned long index;
	int number;
	int cloud;
} cloud_type1;

typedef struct _Pcloud2
{
	unsigned long index;
	int cloud;
	_Pcloud2() {
		cloud = 0;
	}
} cloud_type2;

typedef struct _Pcloud3
{
	unsigned long number_pattern;
	int number;
	int cloud;
	bool has_been_extended;

	_Pcloud3() {
		cloud = 0;
		has_been_extended = 0;
	}
} cloud_type3;


/*definition of macros*/
#define MAXLINECHAR 6000		//maximum characters in each line of a file
#define MAXFILENAMELENGTH 5000		//maximum length of filename
#define BUFF_SIZE	1000000  	//size of the buffer when reading larger chromosome files

#define DIRECTLENGTH 13     	//maximum oligo length to use direct method
#define OVERLAPLENGTH 17		//minimum oligo length to use overlap method

#define HASHOLIGONUMBER 20000000	//maximum number of oligos to store with hash structure when using mixed method to calculate 16mers counts
#define LOWOLIGONUMBER 10000000		//maximum number of oligos to read with array structure when using overlapmethod to calcualter oligo counts

#define MAXMAINMER 40000000		//maximum number of oligos in a cloud class
#define MAXACCMER 80000000		//maximum number of oligos in accessory cloud
#define MAXCLOUD 3000000		//maximum number of p cloud

#define ZERO 48				//ASCII code of character '0'
#define NINE 57				//ASCII code of character '9'

#endif //madrodefine_h__
