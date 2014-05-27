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

#include "readfile.h"
#include "stringhandle.h"
// debug

using std::ifstream;
using std::swap;
using std::ofstream;
using std::cout;
using std::cerr;
using std::endl;

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

//Read the controlfile to get the control settings for program
bool ReadPcloudsControlfile(const char* controlfile, int& kmer_size,
		int& outer_threshold, int& core_1_threshold, int& core_2_threshold,
		int& core_3_threshold, int& core_4_threshold, int& chunk_size,
		unsigned int& genome_size, int& window_size, int& percent,
		bool& build_clouds, bool& annotate_genome, char* kmer_counts_file,
		char* genome_file, char* clouds_summary_file,
		char* core_kmers_assign_file, char* outer_kmers_assign_file,
		char* annotation_file, char* region_file) {

	ifstream ifControlfile(controlfile);

	if (!ifControlfile)
		return (false);

	std::string m_strLine, option, value;

	int i;

	while (!ifControlfile.eof()) {
		getline(ifControlfile, m_strLine);

		if (m_strLine.find('#', 0) < m_strLine.length()) {
			i = 0;

			while (m_strLine[i] != '#')
				;
			i++;

			while (m_strLine[i] == ' ')
				i++;

			option = "";

			while (m_strLine[i] != ' ') {
				option.append(1, m_strLine[i]);
				i++;
			}

			while (m_strLine[i] != '>')
				i++;

			i++;

			while (m_strLine[i] == ' ')
				i++;

			value = "";
			int len = m_strLine.length();
			while (i < len) //&& ((m_strLine[i] != ' ') || (m_strLine[i] != '\n')))
			{
				if ((m_strLine[i] != ' ') || (m_strLine[i] != '\n')) {
					value.append(1, m_strLine[i]);
					i++;
				}
			}

			cout << "Found option " << option << " : " << value << endl;

			if (option == "OligoSize")
				kmer_size = stringtonumber(value);
			else if (option == "COPYTHRESHOLD")
				outer_threshold = stringtonumber(value);
			else if (option == "ENDTHRESHOLD")
				core_1_threshold = stringtonumber(value);
			else if (option == "STEP1THRESHOLD")
				core_2_threshold = stringtonumber(value);
			else if (option == "STEP2THRESHOLD")
				core_3_threshold = stringtonumber(value);
			else if (option == "STEP3THRESHOLD")
				core_4_threshold = stringtonumber(value);
			else if (option == "CALCHUNCKSIZE")
				chunk_size = stringtonumber(value);
			else if (option == "GENOMESIZE")
				genome_size = stringtolargenumber(value);
			else if (option == "WindowSize")
				window_size = stringtonumber(value);
			else if (option == "PercentCutoff")
				percent = stringtonumber(value);
			else if (option == "GETPCLOUDS")
				build_clouds = stringtonumber(value);
			else if (option == "DISSECTION")
				annotate_genome = stringtonumber(value);
			else if (option == "OligoSets")
				stringtoarray(value, kmer_counts_file);
			else if (option == "GenomeInput")
				stringtoarray(value, genome_file);
			else if (option == "MaincloudsInfo")
				stringtoarray(value, clouds_summary_file);
			else if (option == "MaincloudsAssign")
				stringtoarray(value, core_kmers_assign_file);
			else if (option == "AcccloudsAssign")
				stringtoarray(value, outer_kmers_assign_file);
			else if (option == "CloudAnnotation")
				stringtoarray(value, annotation_file);
			else if (option == "RepeatRegion")
				stringtoarray(value, region_file);
			else if (option == "KeepSSRs")
				keep_SSRs = stringtonumber(value);
			else if (option == "PrintCloudsInRegions")
				print_clouds_in_regions = stringtonumber(value);
			else if (option == "ExpandRecursively")
				expand_recursively = stringtonumber(value);
			else if (option == "DontCareAboutClouds")
				dont_care_about_clouds = stringtonumber(value);
		}
	}
	return (true);
}

bool ReadCountsControlfile(const char* pchControlfile, int& size,
		int& m_nChunksize, unsigned int& m_nGenomesize, int& m_nCounts,
		char* pchGenome, char* pchCountfile, unsigned int& m_nMemory) {
	ifstream ifControlfile(pchControlfile);

	if (!ifControlfile)
		return (false);

	std::string m_strLine, m_strTemp, m_strResult;

	int i;

	while (!ifControlfile.eof()) {
		getline(ifControlfile, m_strLine);

		if (m_strLine.find('#', 0) < m_strLine.length()) {
			i = 0;

			while (m_strLine[i] != '#')
				;
			i++;

			while (m_strLine[i] == ' ')
				i++;

			m_strTemp = "";

			while (m_strLine[i] != ' ') {
				m_strTemp.append(1, m_strLine[i]);
				i++;
			}

			while (m_strLine[i] != '>')
				i++;

			i++;

			while (m_strLine[i] == ' ')
				i++;

			m_strResult = "";

			int len = m_strLine.length();
			while (i < len) //&& ((m_strLine[i] != ' ') || (m_strLine[i] != '\n')))
			{
				if ((m_strLine[i] != ' ') || (m_strLine[i] != '\n')) {
					m_strResult.append(1, m_strLine[i]);
					i++;
				}
			}
			if (m_strTemp == "OligoSize")
				size = stringtonumber(m_strResult);
			else if (m_strTemp == "CALCHUNCKSIZE")
				m_nChunksize = stringtonumber(m_strResult);
			else if (m_strTemp == "GENOMESIZE")
				m_nGenomesize = stringtolargenumber(m_strResult);
			else if (m_strTemp == "CALCOUNTS")
				m_nCounts = stringtonumber(m_strResult);
			else if (m_strTemp == "CountGenome")
				stringtoarray(m_strResult, pchGenome);
			else if (m_strTemp == "OligoSets")
				stringtoarray(m_strResult, pchCountfile);
			else if (m_strTemp == "ExpectedAvailableMemory")
				m_nMemory = stringtolargenumber(m_strResult);
		}
	}
	return (true);
}
