#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <time.h>

#include "macrodefine.h"
#include "genomehandle.h"
#include "readfile.h"
#include "overlapmethod.h"
#include "stringhandle.h"

static int compare( char **arg1, char **arg2 );
static int sort(int argc, char **argv);
static int search( int argc, char **argv, char *key);
static int readlowrepeats(char *pchLowfile, char** indextwble, const int& offset, const int& chunck, int& number, const int& size);
static int getpatterns(char *pchSequence,  char** indextwble, const int& number, unsigned long *pPattern, const int& size, const int& genomesize);
static int parsepatternsandoutput(unsigned long *pPattern, char** indextwble, char *pchOut, const int& number, const int& size);
static int repeatfinding(FILE *pfSource, FILE *pfResult, const char* pchCountfile, const int& m_nSize, const int& m_nChunck, const unsigned int& m_nGenome);

int overlapmethod(char *pchSource, char* pchCountfile, const int& m_nSize, const int& m_nChunck, const unsigned int& m_nGenome)
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
	int iRead, iPatternoffset, iReadPattern;
	long long iOffset, genomesize;
	int *piReadCount = NULL;
	char * pchBuff = NULL;

	char temp;
	int number;

	unsigned long *pPattern;

	char* pchInput = (char*) malloc(sizeof(char) * 100);
        char* pchOutput = (char*) malloc(sizeof(char) * 100);
        clock_t start, end;
	int seconds;

	// Read a sequence segment
		
	pchBuff = (char *)malloc(sizeof(char) * m_nChunck);
        genomesize = (long long) m_nGenome;

	for (int size = OVERLAPLENGTH; size <= m_nSize; size++)
	{
		start = clock();

		strcpy(pchInput, pchCountfile);
		temp = 48 + (size-1) / 10;
		strncat(pchInput, &temp, 1);
		temp = 48 + (size-1) %10;
		strncat(pchInput, &temp, 1);
		strcat(pchInput, ".txt");
	
		strcpy(pchOutput, pchCountfile);
		temp = 48 + size / 10;
		strncat(pchOutput, &temp, 1);
		temp = 48 + size %10;
		strncat(pchOutput, &temp, 1);
		strcat(pchOutput, ".txt");

		iPatternoffset = 0;

		//begin loop of different patterns in pattern file;
		do 
		{
			char** indextwble  = (char**) malloc(sizeof(char) * LOWOLIGONUMBER* (size));
                        iReadPattern = readlowrepeats(pchInput, indextwble, iPatternoffset, LOWOLIGONUMBER, number, size-1);

			if (iReadPattern != 0)

			{
				pPattern = (unsigned long*) malloc(sizeof(unsigned long) * number * 4);
                                iOffset = 0;

				do 
				{	
					//allocate memory to read the sequence, store the number of read nucleotide
					piReadCount = (int*) malloc(sizeof(int));
                                        if ( NULL == pchBuff ) 
					{
						printf("The system can not allocate the memory!\n");
						free(pchInput);
                                                free(pchOutput);
						free(pchBuff);
						free(piReadCount);
                                                for (int count = 0; count < number; count++)
                                                {
                                                        free(*(indextwble+count));
                                                }
                                                free(indextwble);
                                                        
                                                free(pPattern);
                                                return(0);
					}
		
					if (iOffset + m_nChunck -1 <= genomesize)
						iRead = readfromfile(pfSource, pchBuff, iOffset, m_nChunck, piReadCount);
					else 
						iRead = readfromfile(pfSource, pchBuff, iOffset, genomesize - iOffset, piReadCount);

					if ( iRead == 0 ) 
					{
                                                free(pchInput);
                                                free(pchOutput);
						free(pchBuff);
						free(piReadCount);
                                                for (int count = 0; count < number; count++)
                                                {
                                                        free(*(indextwble+count));
                                                       
                                                }
                                                free(indextwble);
                                                        
                                                free(pPattern);
						return(0);
					}
					else
						getpatterns(pchBuff, indextwble, number, pPattern, size-1, *piReadCount);
 
					free(piReadCount);
		
					iOffset += m_nChunck - m_nSize +1;
				}
				while ((iOffset < genomesize - m_nSize +1) && (iRead == 2)); 

				parsepatternsandoutput(pPattern, indextwble, pchOutput, number, size-1);

				for (int count = 0; count < number; count++)
                                {
                                    free(*(indextwble+count));
                                }

				free(indextwble);
				free(pPattern);
				iPatternoffset += LOWOLIGONUMBER;
			}	
		}
		while(iReadPattern == 2);

		end = clock();
		seconds = (int)((double)(end - start) / (double)CLOCKS_PER_SEC);
		fprintf(pfResult, "Time Taken for size %d: %d seconds\n\n", size, seconds);
	}

	free(pchInput);
	free(pchOutput);
	free(pchBuff);

	return(1);	
}

int compare( char **arg1, char **arg2 )
{
   /* Compare all of both strings: */
   return strcmp( *arg1, *arg2 );
}

int sort(int argc, char **argv)
{
   int i;

	/* Sort using Quicksort algorithm: */
   qsort( (void *)argv, (size_t)argc, sizeof( char * ), (int (*)(const void*, const void*))compare );

   return(1);
}

int search( int argc, char **argv, char *key)
{
	char **result;

	/* Find the word in key using a binary search algorithm: */
	result = (char **)bsearch( (char *) &key, (char *)argv, argc, sizeof( char * ), (int (*)(const void*, const void*))compare );

	if( result)
		return(result-argv);
	else
		return(-1);
}

int readlowrepeats(char *pchLowfile, char** indextwble, const int& offset, const int& chunck, int& number, const int& size)
{
	FILE *pfLowfile;
	char *pchLine;

	number = 0;
	
	int count = 0;
	if( ( pfLowfile = fopen( pchLowfile, "r" )) == NULL )
	{
		printf("Can not find the lower repeat file: %s.\n", pchLowfile);
		return(0);
	}
	else 
	{
		pchLine = (char *)malloc(sizeof(char) * MAXLINECHAR);
                while ( (!feof(pfLowfile)) && (count< offset) )
		{
			if (count < offset)
			{
				fgets(pchLine, MAXLINECHAR, pfLowfile);
				count++;
			}
		}

		while ((!feof(pfLowfile)) && (number < chunck))
		{
			if( ( fgets( pchLine, MAXLINECHAR, pfLowfile ) != NULL) && (number < chunck))
			{
				// according to each line, build the sequence index array
				*(indextwble + number) = (char*) malloc(sizeof(char)*(size+1));
				strncpy( *(indextwble + number), pchLine, size); 
				*(*(indextwble + number) + size) = '\0';
				number += 1;
			}
		}

		free(pchLine);
		fclose( pfLowfile );
	}

	if (number == 0) 
		return(0);	//no patterns read
	else 
		sort(number, indextwble);

	if (number < chunck) 
		return(1);	//read the patterns, and till to the end of file
	else 
		return(2);	//read the patterns and not to the end of file
}

int getpatterns(char *pchSequence,  char** indextwble, const int& number, unsigned long *pPattern, const int& size, const int& genomesize)
{
	char *pchPattern;

	char *pchTemp;
	
	int *piReadCount = (int*) malloc(sizeof(int));;
        pchPattern = (char *) malloc(sizeof(char) * (size+2));
        pchTemp = (char*) malloc(sizeof(char)* (size+1));
	for (int count =0; count <= genomesize - size -1 ; count ++)
	{
		getsubstring(pchSequence, pchPattern, count, count + size);
		
		pchPattern[size+1] = '\0';

		if (issegmentvalid(pchPattern)==1)
		{
			// get the left substring of the pattern (Length = size )
			strncpy(pchTemp, pchPattern, size);
			pchTemp[size] = '\0';

			// search the indextwble to find whether the left substring is a significant pattern
			int result = search(number, indextwble, pchTemp);

			// if it is a significant pattern of length (size)
			if (result >= 0 )
			{
				// according to the index, get the index number of the pattern of length size in the array
				char temp = pchPattern[strlen(pchPattern) - 1];
				
				int i;

				if ((temp == 'A') || (temp == 'a'))
					i = 0;
				else if ((temp == 'C') || (temp == 'c'))
					i = 1;
				else if ((temp == 'G') || (temp == 'g'))
					i = 2;
				else if ((temp == 'T') || (temp == 't'))
					i = 3;

				int index = 4* result + i;
	
				// increment the number of the pattern
				pPattern[index] += 1;
			}
		}
	}

	free(pchPattern);
	free(pchTemp);
	free(piReadCount);
}

int parsepatternsandoutput(unsigned long *pPattern, char** indextwble, char *pchOut, const int& number, const int& size)
{
	FILE *pfOut = fopen(pchOut, "wb");

	char *pchPattern = (char*) malloc(sizeof(char) * size+2);	
        for (int count = 0; count < number*4; count++)
	{
		if (pPattern[count] > 1)
		{
			int index = count / 4;
			int offset = count %4;
			strncpy(pchPattern, *(indextwble+index), size);
			if (offset == 0)
				pchPattern[size] = 'A';
			else if (offset == 1)
				pchPattern[size] = 'C';
			else if (offset == 2)
				pchPattern[size] = 'G';
			else if (offset == 3)
				pchPattern[size] = 'T';
			pchPattern[size+1] = '\0';
			fprintf(pfOut, "%s %lu\n", pchPattern, pPattern[count]);
		}
	}

	free(pchPattern);
	fclose(pfOut);
}


