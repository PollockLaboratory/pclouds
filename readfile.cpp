/**
 * STP: This is clearly a cpp file. It used to be named *.c.
 *
 */

#include <cstdio>
#include <cstdlib>
#include <stdlib.h>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

#include "readfile.h"
#include "stringhandle.h"
// debug

using namespace std;

// Read a segment of sequence to a buffer from sequence file
// input: pfSource, source file pointer; pchBuff, the buffer to store the read data;
//	 iOffset, the position to start to read data from the beginning of file;
//       iLength, the segment length to read in a time;
//       piReadCount, the number of data read in a time;	
// return value: 0, read error;
//              1, read correct, but read at the end of file;
//              2, read correct, still not at the end of file;

extern bool keep_SSRs;
extern bool print_clouds_in_regions;
extern bool expand_recursively;
extern bool dont_care_about_clouds;
extern bool print_testmers;

int readfromfile(FILE *pfSource, char *pchBuff, long long iOffset, int iLength,
		int *piReadCount) {
#ifdef _WIN32
	if (0 != fseek(pfSource, iOffset, SEEK_SET)) {
#else
		if (0 != fseeko(pfSource, iOffset, SEEK_SET)) {
#endif
		printf("can not go to position: %lld!\n", iOffset);
		return (0);
	}

	*piReadCount = fread(pchBuff, sizeof(char), iLength, pfSource);

	if (*piReadCount < (signed) (sizeof(char) * iLength)) {
		if (feof(pfSource))
			return (1);
		else
			return (0);
	}

	return (2);
}

void read_controlfile(string controlfile, int& kmer_size,
		int& outer_threshold, int& core_threshold_1, int& core_threshold_2,
		int& core_threshold_3, int& core_threshold_4, int& chunk_size,
		unsigned int& genome_size, int& window_size, int& percent,
		bool& build_clouds, bool& annotate_genome, string kmer_counts_file,
		string genome_file, string clouds_summary_file,
		string core_kmers_assign_file, string outer_kmers_assign_file,
		string annotation_file, string region_file) {

	ifstream ifControlfile(controlfile.c_str());

	if (not ifControlfile.good()) {
		cerr << "Cannot read control file" << endl;
		exit(-1);
	}

	bool in_comment = false;
	string option = "";
	string value = "";

	while (ifControlfile.good()) {
		ifControlfile >> option;

		if (option == "#") {
			in_comment = not in_comment;
		} else if (not in_comment) {
			ifControlfile >> value;
			cout << "Found option " << option << " : " << value << endl;

			if (option == "KmerSize")
				kmer_size = stringtonumber(value);
			else if (option == "OuterThreshold")
				outer_threshold = stringtonumber(value);
			else if (option == "CoreThreshold_1")
				core_threshold_1 = stringtonumber(value);
			else if (option == "CoreThreshold_2")
				core_threshold_2 = stringtonumber(value);
			else if (option == "CoreThreshold_3")
				core_threshold_3 = stringtonumber(value);
			else if (option == "CoreThreshold_4")
				core_threshold_4 = stringtonumber(value);
			else if (option == "GenomeChunkSize")
				chunk_size = stringtonumber(value);
			else if (option == "GenomeSize")
				genome_size = stringtolargenumber(value);
			else if (option == "WindowSize")
				window_size = stringtonumber(value);
			else if (option == "PercentCutoff")
				percent = stringtonumber(value);
			else if (option == "BuildClouds")
				build_clouds = stringtonumber(value);
			else if (option == "Annotate")
				annotate_genome = stringtonumber(value);
			else if (option == "KmerCounts")
				kmer_counts_file = value;
			else if (option == "Genome")
				genome_file = value;
			else if (option == "CloudSummaries")
				clouds_summary_file = value;
			else if (option == "CoreKmers")
				core_kmers_assign_file = value;
			else if (option == "OuterKmers")
				outer_kmers_assign_file = value;
			else if (option == "CloudAnnotation")
				annotation_file = value;
			else if (option == "RepeatRegion")
				region_file = value;
			else if (option == "KeepSSRs")
				keep_SSRs = stringtonumber(value);
			else if (option == "PrintCloudsInRegions")
				print_clouds_in_regions = stringtonumber(value);
			else if (option == "ExpandRecursively")
				expand_recursively = stringtonumber(value);
			else if (option == "DontCareAboutClouds")
				dont_care_about_clouds = stringtonumber(value);
			else if (option == "PrintTestmers")
				print_testmers = stringtonumber(value);
		}
	}
}
