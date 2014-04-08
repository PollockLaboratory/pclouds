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
extern bool print_legacy_regions;
extern bool expand_recursively;

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
bool ReadPcloudsControlfile(const char* pchControlfile, int& oligo_size,
		int& m_nCopy, int& m_nEndthreshold, int& m_nStep1, int& m_nStep2,
		int& m_nStep3, int& m_nChunksize, unsigned int& m_nGenomesize,
		int& windowsize, int& percent, int& m_nGetclouds, int& m_nDissection,
		char* pchrepeatfile, char* pchGenome, char* pchMainClouds,
		char* pchMainAssign, char* pchAccAssign, char* pchAnnotationfile,
		char* pchRegionfile) {
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

			cout << "Found option " << m_strTemp << " : " << m_strResult << endl;

			if (m_strTemp == "OligoSize")
				oligo_size = stringtonumber(m_strResult);
			else if (m_strTemp == "COPYTHRESHOLD")
				m_nCopy = stringtonumber(m_strResult);
			else if (m_strTemp == "ENDTHRESHOLD")
				m_nEndthreshold = stringtonumber(m_strResult);
			else if (m_strTemp == "STEP1THRESHOLD")
				m_nStep1 = stringtonumber(m_strResult);
			else if (m_strTemp == "STEP2THRESHOLD")
				m_nStep2 = stringtonumber(m_strResult);
			else if (m_strTemp == "STEP3THRESHOLD")
				m_nStep3 = stringtonumber(m_strResult);
			else if (m_strTemp == "CALCHUNCKSIZE")
				m_nChunksize = stringtonumber(m_strResult);
			else if (m_strTemp == "GENOMESIZE")
				m_nGenomesize = stringtolargenumber(m_strResult);
			else if (m_strTemp == "WindowSize")
				windowsize = stringtonumber(m_strResult);
			else if (m_strTemp == "PercentCutoff")
				percent = stringtonumber(m_strResult);
			else if (m_strTemp == "GETPCLOUDS")
				m_nGetclouds = stringtonumber(m_strResult);
			else if (m_strTemp == "DISSECTION")
				m_nDissection = stringtonumber(m_strResult);
			else if (m_strTemp == "OligoSets")
				stringtoarray(m_strResult, pchrepeatfile);
			else if (m_strTemp == "GenomeInput")
				stringtoarray(m_strResult, pchGenome);
			else if (m_strTemp == "MaincloudsInfo")
				stringtoarray(m_strResult, pchMainClouds);
			else if (m_strTemp == "MaincloudsAssign")
				stringtoarray(m_strResult, pchMainAssign);
			else if (m_strTemp == "AcccloudsAssign")
				stringtoarray(m_strResult, pchAccAssign);
			else if (m_strTemp == "CloudAnnotation")
				stringtoarray(m_strResult, pchAnnotationfile);
			else if (m_strTemp == "RepeatRegion")
				stringtoarray(m_strResult, pchRegionfile);
			else if (m_strTemp == "KeepSSRs")
				keep_SSRs = stringtonumber(m_strResult);
			else if (m_strTemp == "PrintCloudsInRegions")
				print_clouds_in_regions = stringtonumber(m_strResult);
			else if (m_strTemp == "ExpandRecursively")
				expand_recursively = stringtonumber(m_strResult);
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
