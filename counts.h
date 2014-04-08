#ifndef counts_h__
#define counts_h__


int calculatecounts(const char* pchControlfile);
 bool getSmallCountsFile(char *pchCountfile, int size, int overlaplength, char *pchSmallCountfile);
 void getCutoffsForCountingMethods(int &overlaplength, int &directlength, int m_nMemory);



#endif // counts_h__
