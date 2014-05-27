#ifndef readfile_h__
#define readfile_h__


// read a specific region of source file, pfSource - the file pointer of the source file;
// pchBuff - the target pointer to store the read sequence; 
// iOffset - the start point of the region from the start point of the beginning of the file;
// iLength - the length of the region to be read;
// piReadCount - the actual length of the region which have been read.
//               when it read the end of the file, piReadcount will be less than iLength;
//               or these two values should be the same.
//int readfromfile(FILE *pfSource, char *pchBuff, int iOffset, int iLength, int *piReadCount);
int readfromfile(FILE *pfSource, char *pchBuff, long long iOffset, int iLength, int *piReadCount);

bool ReadPcloudsControlfile(const char* controlfile, int& kmer_size,
		int& outer_threshold, int& core_1_threshold, int& core_2_threshold,
		int& core_3_threshold, int& core_4_threshold, int& chunk_size,
		unsigned int& genome_size, int& window_size, int& percent,
		bool& build_clouds, bool& annotate_genome, char* kmer_counts_file,
		char* genome_file, char* clouds_summary_file,
		char* core_kmers_assign_file, char* outer_kmers_assign_file,
		char* annotation_file, char* region_file);

bool ReadCountsControlfile(const char* pchControlfile, int& size, int&
m_nChunksize, unsigned int& m_nGenomesize, int& m_nCounts, char* pchGenome, char* pchCountfile, unsigned int& m_nMemory);

#endif // readfile_h__
