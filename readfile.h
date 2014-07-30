#ifndef readfile_h__
#define readfile_h__

#include <string>
using std::string;

// read a specific region of source file, pfSource - the file pointer of the source file;
// pchBuff - the target pointer to store the read sequence; 
// iOffset - the start point of the region from the start point of the beginning of the file;
// iLength - the length of the region to be read;
// piReadCount - the actual length of the region which have been read.
//               when it read the end of the file, piReadcount will be less than iLength;
//               or these two values should be the same.

void read_chunk_from_genome_file(FILE *genome_file, char *genome_chunk,
		long long genome_chunk_start, int& genome_chunk_size);


void read_controlfile(string controlfile, int& kmer_size,
		int& outer_threshold, int& core_threshold_1, int& core_threshold_2,
		int& core_threshold_3, int& core_threshold_4, int& chunk_size,
		unsigned int& genome_size, int& window_size, int& percent,
		bool& build_clouds, bool& annotate_genome, string& kmer_counts_file,
		string& genome_file, string& clouds_summary_file,
		string& core_kmers_assign_file, string& outer_kmers_assign_file,
		string& annotation_file, string& region_file);

#endif
