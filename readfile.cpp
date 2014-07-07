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
extern bool genome_has_header;

/*
 * I do not understand the return status from this function.
 * return 2 means read was successful.
 * return 1 means ??
 * return 0 mean ?? - cannot go to position
 *
 */
int read_chunk_from_genome_file(FILE *genome_file, char *genome_chunk,
		long long genome_chunk_start, int genome_chunk_size,
		int& genome_chunk_actual_size) {
#ifdef _WIN32
	if (fseek(genome_file, genome_chunk_start, SEEK_SET)) {
#else
		if (fseeko(genome_file, genome_chunk_start, SEEK_SET) != 0) {
#endif
		printf("can not go to position: %lld!\n", genome_chunk_start);
		return (0);
	}

	genome_chunk_actual_size = fread(genome_chunk, sizeof(char),
			genome_chunk_size, genome_file);

	cerr << "chunk size requested " << genome_chunk_size << "\n";
	cerr << "chunk size received " << genome_chunk_actual_size << "\n";
//STP: If read was not successful
	if (genome_chunk_actual_size != genome_chunk_size) {
		genome_chunk[genome_chunk_actual_size] = '\0';
		cout << genome_chunk;
		if (feof(genome_file)) {
			cerr << "Reached the end of the genome file" << endl;
			return (1);
		} else {
			cerr << "Read error in genome file" << endl;
			return (0);
		}
	}
	return (2);
}

void read_chunk_from_genome_file2(FILE *genome_file, char *genome_chunk,
		long long genome_chunk_start, int& genome_chunk_size) {
#ifdef _WIN32
	if (fseek(genome_file, genome_chunk_start, SEEK_SET)) {
#else
		if (fseeko(genome_file, genome_chunk_start, SEEK_SET) != 0) {
#endif
		printf("Can not go to genome chunk start position: %lld!\n", genome_chunk_start);
		exit(-1);
	}
	int expected_genome_chunk_size = genome_chunk_size;

	genome_chunk_size = fread(genome_chunk, sizeof(char),
			genome_chunk_size, genome_file);

	//STP: If read did not return what you expected
	if (genome_chunk_size != expected_genome_chunk_size) {
		genome_chunk[genome_chunk_size] = '\0';
		cout << "This is the last chunk:\n" << genome_chunk << "\n";
		if (feof(genome_file)) {
			cout << "Reached the end of the genome file" << endl;
		} else {
			cerr << "Read error in genome file" << endl;
			exit(-1);
		}
	}
}

void read_controlfile(string controlfile, int& kmer_size, int& outer_threshold,
		int& core_threshold_1, int& core_threshold_2, int& core_threshold_3,
		int& core_threshold_4, int& chunk_size, unsigned int& genome_size,
		int& window_size, int& percent, bool& build_clouds,
		bool& annotate_genome, string& kmer_counts_file, string& genome_file,
		string& clouds_summary_file, string& core_kmers_assign_file,
		string& outer_kmers_assign_file, string& annotation_file,
		string& region_file) {

	ifstream ifControlfile(controlfile.c_str());

	if (not ifControlfile.good()) {
		cerr << "Cannot read control file: " << controlfile << "\n";
		exit(-1);
	}
	cout << "Reading control file: " << controlfile << "\n";

	bool in_comment = false;
	string option = "";
	string value = "";

	while (ifControlfile.good()) {
		ifControlfile >> option;

		if (option == "#") {
			in_comment = not in_comment;
		} else if (not in_comment) {
			ifControlfile >> value;
			// cout << "Found option " << option << " : " << value << endl;

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
			else if (option == "GenomeHasHeader")
				genome_has_header = stringtonumber(value);
		}
	}
}
