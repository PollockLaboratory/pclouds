#include <cstdio>
#include <cstdlib>
#include <stdlib.h>
#include <cstring>
#include <math.h>
#include <time.h>
#include <algorithm>

#include "directmethod.h"
#include "macrodefine.h"
#include "genomehandle.h"
#include "readfile.h"
#include "stringhandle.h"

static int repeatfinding(FILE *pfSource, FILE *pfResult, const char* pchCountfile, const int& m_nSize, const int& m_nChunck, const unsigned int& m_nGenome);

static int patterntoindex(char *pchPattern, int *index, int patternsize);

static int indextopattern(char *pchPattern, int index, int patternsize);

static int getallpatternstoarray(char *pchSequence, int iLength, int patternsize, unsigned long *pPattern, long long *piTotalsegments);

static int outputpattern( unsigned long *pPattern, FILE *pfOutput);

static int outputpatternsequence(unsigned long *pPattern, FILE *pfOutput, int patternsize);

// get the patterns from a source file *pchSource
// input: pchSource --  the source file name containing the source sequence, which is in fasta file format.
int directmethod(char *pchSource, char* pchCountfile, const int& m_nSize, const int& m_nChunck, const unsigned int& m_nGenome)
{
	FILE * pfSource = NULL;
	FILE * pfResult = NULL;
	
	char *pchTarget, *pchResult;
	
	pchTarget = (char *)malloc(sizeof(char) * MAXFILENAMELENGTH);
	pchResult = (char *)malloc(sizeof(char) * MAXFILENAMELENGTH);
	strcpy(pchTarget, pchSource);
	
	strcpy(pchResult, pchTarget);

	char temp1[4] = ".nu";
	strreplaceright(pchSource, pchTarget, temp1);
	
	char temp2[4] = ".rt";
	strreplaceright(pchSource, pchResult, temp2);
	
	pfSource = fopen(pchTarget, "rb");

	if ( NULL == pfSource )
	{
		printf("The file %s does not exist!\n", pchTarget);
		handleformat(pchSource, pchTarget);
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
	int iRead, loop;
	long long iOffset, genomesize;
	int *piReadCount = NULL;
	char * pchBuff = NULL;
	int *piPatterns;
	FILE *pfOut = NULL;
	char temp;

	clock_t start, end;
	int seconds;

	unsigned long *pPattern = NULL;

	pPattern = (unsigned long *)malloc(sizeof(unsigned long)*(unsigned int)pow(4,m_nSize));
        for (loop = 0; loop < pow(4, m_nSize); loop ++)
		pPattern[loop] = 0;
		
	long long *piTotalsegments;   //save the total pattern segments in the genome segments

	start = clock();

	// Read a sequence segment
		
	iOffset = 0;
	
	loop = 0;

	pchBuff = (char *)malloc(sizeof(char) * m_nChunck);
        piTotalsegments = (long long *)malloc(sizeof(long long));
        genomesize = (long long) m_nGenome;

	do 
	{	
		loop++;

		//allocate memory to read the sequence, store the number of read nucleotide

		piReadCount = (int*) malloc(sizeof(int));
                

		if ( NULL == pchBuff ) 
		{
			printf("The system can not allocate the memory!\n");
                        free(piTotalsegments);	
                        free(pchBuff);
                        free(pPattern);
                        free(piReadCount);
                        return(0);
		}

		if (iOffset + m_nChunck -1 <= genomesize)
			iRead = readfromfile(pfSource, pchBuff, iOffset, m_nChunck, piReadCount);
		else 
			iRead = readfromfile(pfSource, pchBuff, iOffset, genomesize - iOffset, piReadCount);

		if ( iRead == 0 ) 
		{
			free(piTotalsegments);	
                        free(pchBuff);
                        free(pPattern);
                        free(piReadCount);
			return(0);
		}
		else
		{	
			if (getallpatternstoarray(pchBuff, *piReadCount, m_nSize, pPattern, piTotalsegments))
			{
				printf("Read Patterns completed.........\n");
				fflush( stdout);
			}
			else 
			{
				free(piTotalsegments);	
                                free(pchBuff);
                                free(pPattern);
                                free(piReadCount);
				return(0);
			}
		}

		free(piReadCount);

		iOffset += m_nChunck- m_nSize +1;
	}
	while ((iOffset < genomesize - m_nSize +1) && (iRead == 2)); 

	pfOut = fopen(pchCountfile, "wb");

	outputpatternsequence(pPattern, pfOut, m_nSize);

	end = clock();

	seconds = (int)((double)(end - start) / (double)CLOCKS_PER_SEC);

	fprintf(pfResult, "Retrieving time Taken Size %d: %d seconds\n", m_nSize, seconds);

	free(piTotalsegments);	
	free(pchBuff);
	free(pPattern);
        fclose(pfOut);
	return(1);	
}

// get all patterns of a part of genome from pchSequence to a nodelist pPatternHead
// input: pchSequence -- source sequence;
//        iLength -- length of the source sequence;
//        patternsize -- length of the pattern
// output: pPatternHead -- Head pointer of the patterns nodelist
//         piTotalsegments -- total number of patterns in this sequence
int getallpatternstoarray(char *pchSequence, int iLength, int patternsize, unsigned long *pPattern, long long *piTotalsegments)
{
	int count, index = 0;
	char *pchPattern;
	
	pchPattern = (char *) calloc(sizeof(char), patternsize+1);
	for (count =0; count <= iLength - patternsize; count ++)
	{
		getsubstring(pchSequence, pchPattern, count, count + patternsize-1);
		
		if (issegmentvalid(pchPattern)==1)
		{
			(*piTotalsegments)++;
		
			patterntoindex(pchPattern, &index, patternsize);

			pPattern[index] += 1;
		}
	}
	
	free(pchPattern);
        return (1);
}

int outputpattern( unsigned long *pPattern, FILE *pfOutput, int patternsize)
{
	unsigned long temp;

	for (int count = 0; count < pow(4, patternsize); count ++)
	{
		temp = pPattern[count];
		if (temp>1)
			fprintf(pfOutput, "%d %ld\n", count, temp);
	}
}

int outputpatternsequence(unsigned long *pPattern, FILE *pfOutput, int patternsize)
{
	unsigned long temp;
	
	char *pchPattern = (char*) malloc(sizeof(char) * (patternsize+1));
        for (int count = 0; count < pow(4, patternsize); count ++)
	{
		temp = pPattern[count];
		if (temp>1)
		{
			indextopattern(pchPattern, count, patternsize);

			fprintf(pfOutput, "%s %ld\n", pchPattern, temp);
		}
	}

	free(pchPattern);
        
}

// transform the pattern sequence to the index of the array
// coding method: A=00, C=01, G= 10, T=11,
int patterntoindex(char *pchPattern, int *index, int patternsize)
{	
	int i;
	*index = 0;
	for (i = 0; i< patternsize; i++)
	{
		*index *= 4;
		if ( (*(pchPattern+i) == 'a') || (*(pchPattern+i) == 'A'))
			*index += 0;
		else if ( (*(pchPattern+i) == 'c') || (*(pchPattern+i) == 'C'))
			*index += 1;
		else if ( (*(pchPattern+i) == 'g') || (*(pchPattern+i) == 'G'))
			*index += 2;
		else if ( (*(pchPattern+i) == 't') || (*(pchPattern+i) == 'T'))
			*index += 3;
	}
	return (1);
}

int indextopattern(char *pchPattern, int index, int patternsize)
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
	
	*(pchPattern+patternsize) = '\0';
	
	return (1);
}
