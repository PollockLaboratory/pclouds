#include <cstdio>
#include <cstdlib>
#include <malloc.h>
#include <cstring>
#include <math.h>
#include <time.h>
#include <algorithm>
#include <vector>

#include "mixedmethod01.h"
#include "macrodefine.h"
#include "genomehandle.h"
#include "readfile.h"
#include "stringhandle.h"

static int repeatfinding(FILE *pfSource, FILE *pfResult, const char* pchCountfile, const int& m_nSize, const int& m_nChunck, const unsigned int& m_nGenome);

static int gettempoligos(FILE* pfSource, int patternsize, const int& m_nChunck, const unsigned int& m_nGenome);

static int getallpatternstobitarray(char *pchSequence, int iLength, int patternsize, bitvector& pPatternA0, bitvector& pPatternC0, bitvector& pPatternG0, bitvector& pPatternT0,
bitvector& pPatternA1, bitvector& pPatternC1, bitvector& pPatternG1, bitvector& pPatternT1, long long *piTotalsegments);

static int outputpatterntemp(const bitvector& pPatternA1, const bitvector& pPatternC1, const bitvector& pPatternG1, const bitvector& pPatternT1, FILE *pfResult, int
patternsize);

static int calculatenumbers(FILE* pfSource, const char* pchCountfile, const int& m_nSize, const int& m_nChunck, const unsigned int& m_nGenome);

static int getnumberstohash(char *pchSequence, int iLength, int patternsize, patternhash_type& Patterns, long long *piTotalsegments);

static int outputpatternsequence(patternhash_type& Patterns, FILE *pfResult, int patternsize);

static int getreversedindexandchar(char* pchPattern, unsigned long& index, const int& patternsize, char& firstchar);

static int getreversedindex(char* pchPattern, unsigned long& index, const int& patternsize);

static int patterntoindex(char *pchPattern, unsigned long& index, int patternsize);

static int indextopattern(char *pchPattern, unsigned long index, int patternsize);

//get the patterns from a source file *pchSource
//input: pchSource --  the source file name containing the source sequence, which is in fasta file format.
int mixedmethod01(char *pchSource, char* pchCountfile, const int& m_nSize, const int& m_nChunck, const unsigned int& m_nGenome)
{
	FILE * pfSource = NULL;
	FILE * pfResult = NULL;
	
	char *pchTarget, *pchResult;
	
	pchTarget = (char *)malloc(sizeof(char) * MAXFILENAMELENGTH);
	pchResult = (char *)malloc(sizeof(char) * MAXFILENAMELENGTH);
	strcpy(pchTarget, pchSource);
	
	strcpy(pchResult, pchTarget);
	
	strreplaceright(pchSource, pchTarget, ".up");
	
	strreplaceright(pchSource, pchResult, ".rt");
	
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

	// get the oligos which are more observed and write to temp files;
	gettempoligos(pfSource, m_nSize, m_nChunck, m_nGenome);

	calculatenumbers(pfSource, pchCountfile, m_nSize, m_nChunck, m_nGenome);

	end = clock();
	seconds = (int)((double)(end - start) / (double)CLOCKS_PER_SEC);
	fprintf(pfResult, "Retrieving time Taken Size %d: %d seconds\n", m_nSize, seconds);

	return(1);	
}

int gettempoligos(FILE* pfSource, int patternsize, const int& m_nChunck, const unsigned int& m_nGenome)
{
	// bit_matrix construction;
	bitvector pPatternA0((unsigned long)pow(4,patternsize-1), false);

	bitvector pPatternC0((unsigned long)pow(4,patternsize-1), false);

	bitvector pPatternT0((unsigned long)pow(4,patternsize-1), false);

	bitvector pPatternG0((unsigned long)pow(4,patternsize-1), false);

	bitvector pPatternA1((unsigned long)pow(4,patternsize-1), false);

	bitvector pPatternC1((unsigned long)pow(4,patternsize-1), false);

	bitvector pPatternT1((unsigned long)pow(4,patternsize-1), false);

	bitvector pPatternG1((unsigned long)pow(4,patternsize-1), false);

	long long piTotalsegments=0;   //save the total pattern segments in the genome segments

	int iRead, loop;
	long long iOffset, genomesize;
	int *piReadCount = NULL;
	char * pchBuff = NULL;

	iOffset = 0;
	
	loop = 0;

	pchBuff = (char *)malloc(sizeof(char) * m_nChunck);
       
	genomesize = (long long) m_nGenome;
	


	do 
	{	
		loop++;

		//allocate memory to read the sequence, store the number of read nucleotide

		piReadCount = (int*) malloc(sizeof(int));
                if ( NULL == pchBuff ) 
		{
			printf("The system can not allocate the memory!\n");
			return(0);
		}

		if (iOffset + m_nChunck -1 <= genomesize)
			iRead = readfromfile(pfSource, pchBuff, iOffset, m_nChunck, piReadCount);
		else 
			iRead = readfromfile(pfSource, pchBuff, iOffset, genomesize - iOffset, piReadCount);

		if ( iRead == 0 ) 
		{
			free(pchBuff);
			free(piReadCount);
			return(0);
		}
		else
		{	
			if (!getallpatternstobitarray(pchBuff, *piReadCount, patternsize, pPatternA0, pPatternC0, pPatternG0, pPatternT0,pPatternA1, pPatternC1, pPatternG1,  
			pPatternT1, &piTotalsegments))
			{
				free(pchBuff);
				free(piReadCount);
				return(0);
			}
		}

		free(piReadCount);

		iOffset += m_nChunck- patternsize +1;
	}
	while ((iOffset < genomesize - patternsize +1) && (iRead == 2)); 

	FILE* pfTemp = fopen("temp.index", "wb");
	outputpatterntemp(pPatternA1, pPatternC1, pPatternG1, pPatternT1, pfTemp, patternsize);
	fclose(pfTemp);

	free(pchBuff);

	return(1);
}

// get all patterns of a part of genome from pchSequence to a nodelist pPatternHead
// input: pchSequence -- source sequence;
//        iLength -- length of the source sequence;
//        patternsize -- length of the pattern
// output: pPatternHead -- Head pointer of the patterns nodelist
//         piTotalsegments -- total number of patterns in this sequence
int getallpatternstobitarray(char *pchSequence, int iLength, int patternsize, bitvector& pPatternA0, bitvector& pPatternC0, bitvector& pPatternG0, bitvector& pPatternT0,
bitvector& pPatternA1, bitvector& pPatternC1, bitvector& pPatternG1, bitvector& pPatternT1, long long *piTotalsegments)
{
	int count;
	unsigned long index = 0;
	unsigned long temp;
	char *pchPattern;
	char pchTemp;

	pchPattern = (char *) calloc(sizeof(char), patternsize+1);
	for (count =0; count <= iLength - patternsize; count++)
	{
		getsubstring(pchSequence, pchPattern, count, count + patternsize-1);
		
		if ((issegmentvalid(pchPattern)==1) && (pchSequence[count] != 'N') && (pchSequence[count] != 'n'))
		{
			(*piTotalsegments)++;

			getreversedindexandchar(pchPattern, index, patternsize, pchTemp);

			switch (pchTemp)
			{
				case 'A' :
					if (!pPatternA0[index])		//oligos end with 'A'
						pPatternA0[index] = true;
					else if (!pPatternA1[index])
						pPatternA1[index] = true;
					break;
				case 'C' :
					if (!pPatternC0[index])		//oligos end with 'C'
						pPatternC0[index] = true;
					else if (!pPatternC1[index])
						pPatternC1[index] = true;
					break;
				case 'G' :
					if (!pPatternG0[index])		//oligos end with 'G'
						pPatternG0[index] = true;
					else if (!pPatternG1[index])
						pPatternG1[index] = true;
					break;
				case 'T' :
					if (!pPatternT0[index])		//oligos end with 'T'
						pPatternT0[index] = true;
					else if (!pPatternT1[index])
						pPatternT1[index] = true;
					break;
				default:
					break;
			}
		}
	}
	
	free(pchPattern);

	return (1);
}

int outputpatterntemp(const bitvector& pPatternA1, const bitvector& pPatternC1, const bitvector& pPatternG1, const bitvector& pPatternT1, FILE *pfResult, int patternsize)
{
	unsigned long temp;
	
	char *pchPattern = (char*) malloc(sizeof(char) * (patternsize));
        for (int i = 0; i < pow(4, patternsize-1); i++)
	{
		if (pPatternA1[i])
		{
			indextopattern(pchPattern, i, patternsize-1);
			fprintf(pfResult, "%sA\n", pchPattern);
		}
		if (pPatternC1[i])
		{
			indextopattern(pchPattern, i, patternsize-1);
			fprintf(pfResult, "%sC\n", pchPattern);
		}
		if (pPatternG1[i])
		{
			indextopattern(pchPattern, i, patternsize-1);
			fprintf(pfResult, "%sG\n", pchPattern);
		}
		if (pPatternT1[i])
		{
			indextopattern(pchPattern, i, patternsize-1);
			fprintf(pfResult, "%sT\n", pchPattern);
		}
	}

	free(pchPattern);
	
	return(1);
}

int calculatenumbers(FILE* pfSource, const char* pchCountfile, const int& m_nSize, const int& m_nChunck, const unsigned int& m_nGenome)
{
	// read the temp files to get the oligos into an array or hash twble
	FILE* pfTemp = fopen("temp.index", "r");

	FILE *pfOut = fopen(pchCountfile, "wb");

	unsigned long index;

	unsigned long patternnumber = 0;
	
	char* pchPattern = (char*) malloc(sizeof(char) * (m_nSize+1));
        while (!feof(pfTemp))
	{
		patternhash_type Patterns;

		patternnumber = 0;

		while (!feof(pfTemp) &&(patternnumber <= HASHOLIGONUMBER))
		{
			if (patternnumber <= HASHOLIGONUMBER)
			{
				fscanf(pfTemp, "%s", pchPattern);
				getreversedindex(pchPattern, index, m_nSize);
				Patterns.insert( patternhash_type::value_type(index, 0));
				patternnumber++;
			}
		}

		// search the whole genome and then calculate the numbers
		// Read a sequence segment
		
		long long piTotalsegments=0;   //save the total pattern segments in the genome segments

		int iRead, loop;
		long long iOffset, genomesize;
		int *piReadCount = NULL;
		char * pchBuff = NULL;

		iOffset = 0;
	
		loop = 0;

		pchBuff = (char *)malloc(sizeof(char) * m_nChunck);
               genomesize = (long long) m_nGenome;

		do 
		{	
			loop++;
	
			//allocate memory to read the sequence, store the number of read nucleotide

			piReadCount = (int*) malloc(sizeof(int));
                        if ( NULL == pchBuff ) 
			{
				printf("The system can not allocate the memory!\n");
				return(0);
			}

			if (iOffset + m_nChunck -1 <= genomesize)
				iRead = readfromfile(pfSource, pchBuff, iOffset, m_nChunck, piReadCount);
			else 
				iRead = readfromfile(pfSource, pchBuff, iOffset, genomesize - iOffset, piReadCount);

			if ( iRead == 0 ) 
			{
				free(pchBuff);
				free(piReadCount);
				return(0);
			}
			else
			{	
				if (!getnumberstohash(pchBuff, *piReadCount, m_nSize, Patterns, &piTotalsegments))
				{
					free(pchBuff);
					free(piReadCount);
					return(0);
				}
			}
	
			free(piReadCount);

			iOffset += m_nChunck- m_nSize +1;
		}
		while ((iOffset < genomesize - m_nSize +1) && (iRead == 2)); 
                free(pchBuff);
		//output the patterns and numbers from hashtwble
		outputpatternsequence(Patterns, pfOut, m_nSize);
	}

	fclose(pfOut);
	fclose(pfTemp);

	return(1);
}

int getnumberstohash(char *pchSequence, int iLength, int patternsize, patternhash_type& Patterns, long long *piTotalsegments)
{
	int count;
	unsigned long index = 0;

	char* pchPattern = (char *) malloc(sizeof(char)*(patternsize+1));
        patternhash_type::iterator pPatterns;
		
	for (count =0; count <= iLength - patternsize-1; count ++)
	{
		getsubstring(pchSequence, pchPattern, count, count + patternsize -1);

		pchPattern[patternsize] = '\0';
		
		if (issegmentvalid(pchPattern)==1)
		{
			(*piTotalsegments)++;

			getreversedindex(pchPattern, index, patternsize);

			pPatterns = Patterns.find(index);
	
			if (pPatterns != Patterns.end())
				((*pPatterns).second)++;
		}
	}

	free(pchPattern);
	return (1);
}

int outputpatternsequence(patternhash_type& Patterns, FILE *pfResult, int patternsize)
{
	
	char *pchPattern = (char*) malloc(sizeof(char) * (patternsize+1));
        patternhash_type::iterator pPatterns;

 	for (pPatterns = Patterns.begin(); pPatterns != Patterns.end(); pPatterns++)
	{
		indextopattern(pchPattern, (*pPatterns).first, patternsize);
		fprintf(pfResult, "%s %u\n", pchPattern, (*pPatterns).second);
	}

	free(pchPattern);
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

int getreversedindexandchar(char* pchPattern, unsigned long& index, const int& patternsize, char& firstchar)
{
	unsigned long index1, index2;
	
	patterntoindex(pchPattern, index1, patternsize);
	
	char* pchReversed = (char*) malloc(sizeof(char) * (patternsize+1));
	getreversecomplement(pchPattern, pchReversed);
	
	patterntoindex(pchReversed, index2, patternsize);
	
	if (index1 < index2)
	{
		patterntoindex(pchPattern, index, patternsize-1);
		firstchar = pchPattern[patternsize-1];
	}
	else
	{
		patterntoindex(pchReversed, index, patternsize-1);
		firstchar = pchReversed[patternsize-1];
	}
	
	free(pchReversed);
	return(1);
}

int getreversedindex(char* pchPattern, unsigned long& index, const int& patternsize)
{
	unsigned long index1, index2;
	
	patterntoindex(pchPattern, index1, patternsize);
	
	char* pchReversed = (char*) malloc(sizeof(char) * (patternsize+1));
	getreversecomplement(pchPattern, pchReversed);
	
	patterntoindex(pchReversed, index2, patternsize);
	
	if (index1 < index2)
		index = index1;
	else
		index = index2;
			
	free(pchReversed);
	return(1);
}

int indextopattern(char *pchPattern, unsigned long index, int patternsize)
{
	int i, value;
	unsigned long temp = index;
	
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
	
	*(pchPattern + patternsize) = '\0';
	return (1);
}
