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

bool ReadPcloudsControlfile(const char* pchControlfile, int& size, int& m_nCopy, int& m_nEndthreshold, int& m_nStep1, int& m_nStep2, int& m_nStep3, int& m_nChunksize, unsigned int&
m_nGenomesize, int& windowsize, int& percent, int& m_nGetclouds, int& m_nDissection, char* pchrepeatfile, char* pchGenome, char* pchMainClouds, char* pchMainAssign, char*
pchAccAssign, char* pchAnnotationfile, char* pchRegionfile);

bool ReadCountsControlfile(const char* pchControlfile, int& size, int&
m_nChunksize, unsigned int& m_nGenomesize, int& m_nCounts, char* pchGenome, char* pchCountfile, unsigned int& m_nMemory);

#endif // readfile_h__
