#ifndef genomehandle_h__
#define genomehandle_h__



// *************************************************************************************
// handle the original sequence file to another file which only contains nucleotides
// pchSource - source file filename; pchTarget - target file filename
// *************************************************************************************
int handleformat(char *pchSource, char *pchTarget);

int handleuppercase(char *pchSource, char *pchTarget);

#endif // genomehandle_h__
