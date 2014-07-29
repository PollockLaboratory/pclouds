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

void read_chunk_from_genome_file(FILE *genome_file, char *genome_chunk,
		long long genome_chunk_start, int& genome_chunk_size) {
#ifdef _WIN32
	if (fseek(genome_file, genome_chunk_start, SEEK_SET)) {
#else
		if (fseeko(genome_file, genome_chunk_start, SEEK_SET) != 0) {
#endif
		printf("Can not go to genome chunk start position: %lld!\n",
				genome_chunk_start);
		exit(-1);
	}
	int expected_genome_chunk_size = genome_chunk_size;

	genome_chunk_size = fread(genome_chunk, sizeof(char), genome_chunk_size,
			genome_file);

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

void read_controlfile(string controlfile_name, int& kmer_size,
		int& outer_threshold, int& core_threshold_1, int& core_threshold_2,
		int& core_threshold_3, int& core_threshold_4, int& chunk_size,
		unsigned int& genome_size, int& window_size, int& percent,
		bool& build_clouds, bool& annotate_genome, string& kmer_counts_file,
		string& genome_file, string& clouds_summary_file,
		string& core_kmers_assign_file, string& outer_kmers_assign_file,
		string& annotation_file, string& region_file) {

	ifstream controlfile(controlfile_name.c_str());

	if (not controlfile.good()) {
		cerr << "Cannot read control file: " << controlfile_name << "\n";
		exit(-1);
	}
	cout << "Reading control file: " << controlfile_name << "\n";

	bool in_comment = false;
	string option = "";
	string value = "";

	while (controlfile.good()) {
		controlfile >> option;

		if (option == "#") {
			in_comment = not in_comment;
		} else if (not in_comment) {
			if (option == "KmerSize")
				controlfile >> kmer_size;
			else if (option == "OuterThreshold")
				controlfile >> outer_threshold;
			else if (option == "CoreThreshold_1")
				controlfile >> core_threshold_1;
			else if (option == "CoreThreshold_2")
				controlfile >> core_threshold_2;
			else if (option == "CoreThreshold_3")
				controlfile >> core_threshold_3;
			else if (option == "CoreThreshold_4")
				controlfile >> core_threshold_4;
			else if (option == "GenomeChunkSize")
				controlfile >> chunk_size;
			else if (option == "GenomeSize")
				controlfile >> genome_size;
			else if (option == "WindowSize")
				controlfile >> window_size;
			else if (option == "PercentCutoff")
				controlfile >> percent;
			else if (option == "BuildClouds")
				controlfile >> build_clouds;
			else if (option == "AnnotateGenome")
				controlfile >> annotate_genome;
			else if (option == "KmerCounts")
				controlfile >> kmer_counts_file;
			else if (option == "Genome")
				controlfile >> genome_file;
			else if (option == "CloudSummaries")
				controlfile >> clouds_summary_file;
			else if (option == "CoreKmers")
				controlfile >> core_kmers_assign_file;
			else if (option == "OuterKmers")
				controlfile >> outer_kmers_assign_file;
			else if (option == "CloudAnnotation")
				controlfile >> annotation_file;
			else if (option == "RepeatRegion")
				controlfile >> region_file;
			else if (option == "KeepSSRs")
				controlfile >> keep_SSRs;
			else if (option == "PrintCloudsInRegions")
				controlfile >> print_clouds_in_regions;
			else if (option == "ExpandRecursively")
				controlfile >> expand_recursively;
			else if (option == "DontCareAboutClouds")
				controlfile >> dont_care_about_clouds;
			else if (option == "PrintTestmers")
				controlfile >> print_testmers;
			else if (option == "GenomeHasHeader")
				controlfile >> genome_has_header;
			else if (option == "CutoffValues") {
				cout << "Found Cutoff Values" << endl;
				controlfile >> outer_threshold;
				controlfile >> core_threshold_1;
				controlfile >> core_threshold_2;
				controlfile >> core_threshold_3;
				controlfile >> core_threshold_4;
			}
		}
	}
	if (dont_care_about_clouds and print_clouds_in_regions) {
		cerr << "You have selected that you don't care about clouds\n"
				<< "but you asked to print clouds in regions.\n"
				<< "The cloud ids printed will not mean much.\n";
	}
}
