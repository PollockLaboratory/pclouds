#include <cstdio>
#include <cstdlib>
#include <stdlib.h>
#include <cstring>
#include <math.h>
#include <time.h>
#include <algorithm>
#include <vector>

#include "mixedmethod02.h"
#include "macrodefine.h"
#include "genomehandle.h"
#include "readfile.h"
#include "stringhandle.h"

static int repeatfinding(FILE *pfSource, FILE *pfResult, const char* pchCountfile, const int& m_nSize, const int& m_nChunck, const unsigned int& m_nGenome);

static int patterntoindex(char *pchPattern, unsigned long& index, int patternsize);

static int indextopattern(char *pchPattern, unsigned long index, int patternsize);

static int getallpatternstohash(char *pchSequence, int iLength, int patternsize, bit_matrix& PatternBitMatrix1, bit_matrix& PatternBitMatrix2, const int* matrixsize, long long *piTotalsegments);

static int outputpatternsequence(patternhash_type& Patterns, const int* matrixsize, const unsigned long* indexsize, FILE *pfResult, int patternsize);

static int getpatternmatrixindex(char* pchPattern, unsigned long& i, unsigned long& j);

static int getmatrixindextopattern(char* pchPattern, const unsigned long& i, const unsigned long& j, const int& patternsize);

static int outputpatterntemp(bit_matrix& PatternBitMatrix, const int* matrixsize, const unsigned long* indexsize, FILE *pfResult, int patternsize);

static int gettempoligos(FILE* pfSource, int patternsize, unsigned long* indexsize, int* matrixsize, const int& m_nChunck, const unsigned int& m_nGenome);

static int calculatenumbers(FILE* pfSource, const char* pchCountfile, unsigned long* indexsize, int* matrixsize, const int& m_nSize, const int& m_nChunck, const unsigned int& m_nGenome);

static int getnumberstohash(char *pchSequence, int iLength, int patternsize, patternhash_type& Patterns, unsigned long* indexsize, int* matrixsize, long long *piTotalsegments);

static int getmatrixindexsize(int patternsize, unsigned long* indexsize, int* matrixsize);

static int matrixindextopattern(unsigned long index, char* pchPattern, const unsigned long* indexsize, int patternsize);

static int sbsearch(int n, const unsigned long *indexsize, unsigned long &key);

// get the patterns from a source file *pchSource
// input: pchSource --  the source file name containing the source sequence, which is in fasta file format.
int mixedmethod02(char *pchSource, char* pchCountfile, const int& m_nSize, const int& m_nChunck, const unsigned int& m_nGenome)
{
	FILE * pfSource = NULL;
	FILE * pfResult = NULL;
	
	char *pchTarget, *pchResult;
	
	pchTarget = (char *)malloc(sizeof(char) * MAXFILENAMELENGTH);
	pchResult = (char *)malloc(sizeof(char) * MAXFILENAMELENGTH);
	strcpy(pchTarget, pchSource);
	
	strcpy(pchResult, pchTarget);
	
    char temp1[4] = ".up";
	strreplaceright(pchSource, pchTarget, temp1);
	
	char temp2[4] = ".rt";
	strreplaceright(pchSource, pchResult, temp2);
	
	pfSource = fopen(pchTarget, "rb");
       
	if ( NULL == pfSource )
	{
		printf("The file %s does not exist!\n", pchTarget);
		handleuppercase(pchSource, pchTarget);
		pfSource = fopen(pchTarget, "rb");
	}
	
	if (NULL == pfSource)
	{
		printf("The file %s could not be opened!\n", pchTarget);
		free(pchTarget);
		free(pchResult);
		return(0);
	}
	else
	{
		// *open the result file
		pfResult = fopen(pchResult, "wb");
		if ( NULL == pfResult ) 
		{
			printf("The file %s could not be opened!\n", pchResult);
			fclose(pfSource);
			free(pchTarget);
			free(pchResult);
			return(0);
		}
		if ( repeatfinding(pfSource, pfResult, pchCountfile, m_nSize, m_nChunck, m_nGenome) == 0)
		{  
			fclose(pfSource);
			fclose(pfResult);
			free(pchTarget);
			free(pchResult);
			return(0);
		}
		else
		{	
			fclose(pfSource);
			fclose(pfResult);
			free(pchTarget);
			free(pchResult);
		}
	}
	
	return (1);
}

//find the REPEATARR pattern from a source file pfSource, calculate the pattern in pHead, and retrieve the pHead
//to find the significant patterns and significant motifs, write to the pfResult file
int repeatfinding(FILE *pfSource, FILE *pfResult, const char* pchCountfile, const int& m_nSize, const int& m_nChunck, const unsigned int& m_nGenome)
{
	clock_t start, end;
	int seconds;

	start = clock();

	unsigned long* indexsize = (unsigned long*) malloc(sizeof(unsigned long) * (unsigned long) pow(4, m_nSize/2) ); 
        int*  matrixsize = (int*) malloc(sizeof(int) * (unsigned long) pow(4, m_nSize/2) ); 

	// get the oligos which are more observed and write to temp files;
	gettempoligos(pfSource, m_nSize, indexsize, matrixsize, m_nChunck, m_nGenome);
	calculatenumbers(pfSource, pchCountfile, indexsize, matrixsize, m_nSize, m_nChunck, m_nGenome);

	free(indexsize);
	free(matrixsize);

	end = clock();
	seconds = (int)((double)(end - start) / (double)CLOCKS_PER_SEC);
	fprintf(pfResult, "Retrieving time Taken Size %d: %d seconds\n", m_nSize, seconds);

	return(1);	
}

int getmatrixindexsize(int patternsize, unsigned long* indexsize, int* matrixsize)
{
	char* pchPattern = (char*) malloc(sizeof(char) * (patternsize/2+1));
        char* pchReverse = (char*) malloc(sizeof(char) * (patternsize/2+1));
        unsigned long index;

	unsigned long total  = 0;

	for (int i = 0; i < pow(4, patternsize/2); i++)
	{
		// index equal to i, get the last 8 characters of the pattern from index
		index = i;

		indextopattern(pchPattern, index, patternsize/2);
		
		// get the reverse complement of the last 8 characters, then get the index of the reverse complement LR
		//the size of the vector[i] is LR+1;

		 getreversecomplement(pchPattern, pchReverse);

		patterntoindex(pchReverse, index, patternsize/2);

		matrixsize[i] = index+1;

		indexsize[i] = total;

		total += index+1;
	}

	free(pchPattern);
        free(pchReverse);
        
	return(1);
}

int gettempoligos(FILE* pfSource, int patternsize, unsigned long* indexsize, int* matrixsize, const int& m_nChunck, const unsigned int& m_nGenome)
{
	// bit_matrix construction;
	bit_matrix PatternMatrix1((unsigned long)pow(4, patternsize/2));

	bit_matrix PatternMatrix2((unsigned long)pow(4, patternsize/2));

	unsigned long index;

	char* pchPattern = (char*) malloc(sizeof(char) * (patternsize/2+1));
        char* pchReverse = (char*) malloc(sizeof(char) * (patternsize/2+1));
        // each line correspond to each pattern with index of the last 8 chars.

	unsigned long total  = 0;
        for (int i = 0; i < pow(4, patternsize/2); i++)
	{
		// index equal to i, get the last 8 characters of the pattern from index
		index = i;

		indextopattern(pchPattern, index, patternsize/2);
		
		// get the reverse complement of the last 8 characters, then get the index of the reverse complement LR
		//the size of the vector[i] is LR+1;

		 getreversecomplement(pchPattern, pchReverse);

		patterntoindex(pchReverse, index, patternsize/2);

		PatternMatrix1[i].reserve(index+1); 

		PatternMatrix2[i].reserve(index+1); 

		matrixsize[i] = index+1;

		indexsize[i] = total;

		total += index+1;
	}

	free(pchPattern);
        free(pchReverse);
        
	// Read a sequence segment
		
	long long piTotalsegments;   //save the total pattern segments in the genome segments

	int iRead;
	long long iOffset, genomesize;
	int piReadCount;
	char * pchBuff = NULL;

	iOffset = 0;
	
	pchBuff = (char *)malloc(sizeof(char) * m_nChunck);
       
	genomesize = (long long) m_nGenome;
	do 
	{	
		

		//allocate memory to read the sequence, store the number of read nucleotide

		piReadCount = 0;

		if ( NULL == pchBuff ) 
		{
			printf("The system can not allocate the memory!\n");
			return(0);
		}

		if (iOffset + m_nChunck -1 <= genomesize)
			iRead = readfromfile(pfSource, pchBuff, iOffset, m_nChunck, &piReadCount);
		else 
			iRead = readfromfile(pfSource, pchBuff, iOffset, genomesize - iOffset, &piReadCount);

		if ( iRead == 0 ) 
		{
			free(pchBuff);
			return(0);
		}
		else
		{	
			if (!getallpatternstohash(pchBuff, piReadCount, patternsize, PatternMatrix1, PatternMatrix2, matrixsize, &piTotalsegments))
			{
				free(pchBuff);
				return(0);
			}
		}

		

		iOffset += m_nChunck- patternsize +1;
	}
	while ((iOffset < genomesize - patternsize +1) && (iRead == 2)); 
        free(pchBuff);

	FILE* pfTemp = fopen("temp.index", "wb");
	outputpatterntemp(PatternMatrix2, matrixsize, indexsize, pfTemp, patternsize);
	fclose(pfTemp);

	for (int i = 0; i < pow(4, patternsize/2); i++)
	{
		PatternMatrix1[i].clear(); 
		PatternMatrix2[i].clear(); 
	}

	PatternMatrix1.clear();
	PatternMatrix2.clear();
		

	return(1);
}

int calculatenumbers(FILE* pfSource, const char* pchCountfile, unsigned long* indexsize, int* matrixsize, const int& m_nSize, const int& m_nChunck, const unsigned int& m_nGenome)
{
	// read the temp files to get the oligos into an array or hash twble
	FILE* pfTemp = fopen("temp.index", "r");
	
	FILE *pfOut = fopen(pchCountfile, "wb");

	unsigned long index;

	unsigned long patternnumber = 0;

	//need revision to include all the oligos
	while (!feof(pfTemp))
	{
		patternhash_type Patterns;

		patternnumber = 0;

		while (!feof(pfTemp) &&(patternnumber <= HASHOLIGONUMBER))
		{
			if (patternnumber <= HASHOLIGONUMBER)
			{
				fscanf(pfTemp, "%ld", &index);
				Patterns.insert( patternhash_type::value_type(index, 0));
				patternnumber++;
			}
		}
		// search the whole genome and then calculate the numbers
		// Read a sequence segment
		
		long long piTotalsegments;   //save the total pattern segments in the genome segments

		int iRead, loop;
		long long iOffset, genomesize;
		//int *piReadCount = NULL;
		char * pchBuff = NULL;

		iOffset = 0;
	
		loop = 0;

		pchBuff = (char *)malloc(sizeof(char) * m_nChunck);
                genomesize = (long long) m_nGenome;

		do 
		{	
			loop++;
	
			int piReadCount = 0;
			if ( NULL == pchBuff ) 
			{
				printf("The system can not allocate the memory!\n");
				return(0);
			}

			if (iOffset + m_nChunck -1 <= genomesize)
				iRead = readfromfile(pfSource, pchBuff, iOffset, m_nChunck, &piReadCount);
			else 
				iRead = readfromfile(pfSource, pchBuff, iOffset, genomesize - iOffset, &piReadCount);

			if ( iRead == 0 ) 
			{
				free(pchBuff);
				return(0);
			}
			else
			{	
				if (!getnumberstohash(pchBuff, piReadCount, m_nSize, Patterns, indexsize, matrixsize, &piTotalsegments))
				{
					free(pchBuff);
					return(0);
				}
			}
	

			iOffset += m_nChunck- m_nSize +1;
		}
		while ((iOffset < genomesize - m_nSize +1) && (iRead == 2)); 

		free(pchBuff);	
		//output the patterns and numbers from hashtwble
		outputpatternsequence(Patterns, matrixsize, indexsize, pfOut, m_nSize);
		Patterns.clear();
	}

	fclose(pfOut);
	fclose(pfTemp);
	return(1);
}

int getnumberstohash(char *pchSequence, int iLength, int patternsize, patternhash_type& Patterns, unsigned long* indexsize, int *matrixsize, long long *piTotalsegments)
{
	int count;
	unsigned long index = 0;

	char* pchPattern = (char *) calloc(sizeof(char), patternsize+1);
        char* pchReverse = (char *) calloc(sizeof(char), patternsize+1);
        unsigned long left; 

	unsigned long right;

	patternhash_type::iterator pPatterns;
		
	for (count =0; count <= iLength - patternsize-1; count ++)
	{
		getsubstring(pchSequence, pchPattern, count, count + patternsize -1);

		pchPattern[patternsize] = '\0';
		
		if (issegmentvalid(pchPattern)==1)
		{
			(*piTotalsegments)++;

			getpatternmatrixindex(pchPattern, right, left);

			if (left >= matrixsize[right])
			{
				getreversecomplement(pchPattern, pchReverse);
				strcpy(pchPattern, pchReverse);
				getpatternmatrixindex(pchPattern, right, left);
			}
			
			index = indexsize[right] + left;

			pPatterns = Patterns.find(index);
			
			if (pPatterns != Patterns.end())
				((*pPatterns).second)++;
		}
	}

	free(pchPattern);
	free(pchReverse);

	return (1);
}

// get all patterns of a part of genome from pchSequence to a nodelist pPatternHead
// input: pchSequence -- source sequence;
//        iLength -- length of the source sequence;
//        patternsize -- length of the pattern
// output: pPatternHead -- Head pointer of the patterns nodelist
//         piTotalsegments -- total number of patterns in this sequence
 int getallpatternstohash(char *pchSequence, int iLength, int patternsize, bit_matrix& PatternBitMatrix1, bit_matrix& PatternBitMatrix2, const int* matrixsize, long long *piTotalsegments)
{
	int count;
	unsigned long index = 0;

	char* pchPattern = (char *) calloc(sizeof(char), patternsize+1);
        char* pchReverse = (char *) calloc(sizeof(char), patternsize+1);
        char* pchTemp = (char *) calloc(sizeof(char), patternsize);
        unsigned long left; 

	unsigned long right;
		
	for (count =0; count <= iLength - patternsize-1; count ++)
	{
		getsubstring(pchSequence, pchPattern, count, count + patternsize -1);

		pchPattern[patternsize] = '\0';
		
		if (issegmentvalid(pchPattern)==1)
		{
			(*piTotalsegments)++;

			getpatternmatrixindex(pchPattern, right, left);

			if (left >= matrixsize[right])
			{
				getreversecomplement(pchPattern, pchReverse);
				strcpy(pchPattern, pchReverse);
				getpatternmatrixindex(pchPattern, right, left);
			}

			if (!PatternBitMatrix1[right][left])
				PatternBitMatrix1[right][left] = true;
			else if (!PatternBitMatrix2[right][left])
				PatternBitMatrix2[right][left] = true;
		}
	}

	free(pchPattern);
	free(pchReverse);
	free(pchTemp);

	return (1);
}

int outputpatterntemp(bit_matrix& PatternBitMatrix, const int* matrixsize, const unsigned long* indexsize, FILE *pfResult, int patternsize)
{
	unsigned long temp;
	
	char *pchPattern = (char*) malloc(sizeof(char) * (patternsize+1));
        for (int i = 0; i < pow(4, patternsize/2); i++)
	{
		for (int j = 0; j < matrixsize[i]; j++)
		{
			if (PatternBitMatrix[i][j])
			{
				temp = indexsize[i] + j;
				fprintf(pfResult, "%ld ", temp);
			}
		}
	}

	free(pchPattern);
}

int outputpatternsequence(patternhash_type& Patterns, const int* matrixsize, const unsigned long* indexsize, FILE *pfResult, int patternsize)
{
	unsigned long temp;
	
	char *pchPattern = (char*) malloc(sizeof(char) * (patternsize+1));
        patternhash_type::iterator pPatterns;

	for (pPatterns = Patterns.begin(); pPatterns != Patterns.end(); pPatterns++)
	{
		matrixindextopattern((*pPatterns).first, pchPattern, indexsize, patternsize);
		fprintf(pfResult, "%s %ld\n", pchPattern, (*pPatterns).second);
	}


	free(pchPattern);
}

int matrixindextopattern(unsigned long index, char* pchPattern, const unsigned long* indexsize, int patternsize)
{
	int i, j;
	
	i = 0;

	i = sbsearch((int)pow(4, patternsize/2), indexsize, index);

	j = index - indexsize[i];

	char* pchTemp = (char*) malloc( sizeof(char) * (patternsize /2+1));

	// get the index of first 8 chars;
	indextopattern(pchTemp, j, patternsize/2);

	strcpy(pchPattern, pchTemp);
	
	// get the index of the last 8 chars;
	indextopattern(pchTemp, i, patternsize/2);

	strcat(pchPattern, pchTemp);
	free(pchTemp);
	return(1);
}

// transform the pattern sequence to the index of the array
// coding method: A=00, C=01, G= 10, T=11,
int patterntoindex(char *pchPattern, unsigned long& index, int patternsize)
{	
	int i;
	index = 0;

	for (i = 0; i< patternsize; i++)
	{
		index *= 4;
		if ( (*(pchPattern+i) == 'a') || (*(pchPattern+i) == 'A'))
			index += 0;
		else if ( (*(pchPattern+i) == 'c') || (*(pchPattern+i) == 'C'))
			index += 1;
		else if ( (*(pchPattern+i) == 'g') || (*(pchPattern+i) == 'G'))
			index += 2;
		else if ( (*(pchPattern+i) == 't') || (*(pchPattern+i) == 'T'))
			index += 3;
	}
	return (1);
}

int indextopattern(char *pchPattern, unsigned long index, int patternsize)
{
	int i, value;
	int temp = index;
	for (i = patternsize - 1; i>= 0 ; i--)
	{
		value = temp % 4;
		if (value == 0)
			*(pchPattern +i) = 'A';
		else if (value == 1)
			*(pchPattern +i) = 'C';
		else if (value == 2)
			*(pchPattern +i) = 'G';
		else if (value == 3)
			*(pchPattern +i) = 'T';
		temp = temp /4;
	}
	pchPattern[patternsize] = '\0';
	return (1);
}

// ------------------------------
// write a subfunction to get the (i, j) for a pattern
int getpatternmatrixindex(char* pchPattern, unsigned long& i, unsigned long& j)
{
	char* pchTemp = (char*) malloc( sizeof(char) * (strlen(pchPattern) /2+1));


	// get the index of last 8 chars;
	strright(strlen(pchPattern)/2, pchPattern, pchTemp);
	
	patterntoindex(pchTemp, i, strlen(pchPattern)/2);

	// get the index of the first 8 chars;
	strleft(strlen(pchPattern)/2, pchPattern, pchTemp);
	patterntoindex(pchTemp, j, strlen(pchPattern)/2);
	free(pchTemp);
	return(1);
}

int getmatrixindextopattern(char* pchPattern, const unsigned long& i, const unsigned long& j, const int& patternsize)
{
	char* pchTemp = (char*) malloc( sizeof(char) * (patternsize /2+1));
        // get the index of first 8 chars;
	indextopattern(pchTemp, j, patternsize/2);

	strcpy(pchPattern, pchTemp);
	
	// get the index of the last 8 chars;
	indextopattern(pchTemp, i, patternsize/2);

	strcat(pchPattern, pchTemp);
	free(pchTemp);
	return(1);
}

int sbsearch(int n, const unsigned long *indexsize, unsigned long &key)
{
	int m;
	int length = n;
	int site = 0;
	while (n >= 1) 
	{
		m = n/2;
		if ( site+m < length -1)
		{
	        	if ( ( key >= *(indexsize+site+m) ) && ( key < *(indexsize+site+m+1)))
				return(site+m);
			else if (key < *(indexsize+site+m))
		       		n = m;
			else 
			{
				site = site + m +1; 
				n = n - m - 1;
		        }
		}
		else
			return(length-1);	
	}

	return(-1);
}
