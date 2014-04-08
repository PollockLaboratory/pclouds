/**
 * STP: This is clearly a cpp file. It used to be named *.c.
 *
 */



#include <cstdio>
#include <cstdlib>
#include <stdlib.h>
#include <cstring>

#include "genomehandle.h"
#include "macrodefine.h"
#include "readfile.h"

// handle the read original genome sequence into a specific format which contains only the nucleotides,
// in which the first line of the FASTA format and '\n' character has been removed.
// pchSource - the original sequence in the original FASTA format file;
// piSource - the length of the original sequence; 
// pchTarget - the handled sequence in the format of pure nucleotides;
// piTarget - the length of the handled sequence;
static int sequenceprehandle(char *pchSource, int *piSource, char *pchTarget, int *piTarget);

// write the handled sequence in pure nucleotides format in a target file
// pfTarget - the file pointer of the target file to store the handled sequence;
// pchTarget -  the sequence pointer which contains the sequence in pure nucleotide format to be stored;
// piTargetCount - the length of the target nucleotide sequence to be stored
// piWriteTotal - the total length of the target sequence which have been stored in the target file
static int writesequencefile(FILE *pfTarget, char *pchTarget, int *piTargetCount, long long *piWriteTotal);

// Transform the original genome file which is in FASTA format to the target genome file in pure nucleotides format
// the step is first split the large original genome file into several segments in a specific length (defined by BUFF_SIZE)
// secondly, read each segment of the original file into pchBuffer; prehandle each original sequence; write each handled 
// sequence into the target file; repeat the wbove steps till the original sequence have all been handled, and the handled
// sequence are stored in the target file
// pfSource - file pointer of the origina genome file in FASTA format
// pfTarget - file pointer of the target genome file in pure nucleotide format
static int filetransform(FILE *pfSource, FILE *pfTarget);

static int fileuppercase(FILE *pfSource, FILE *pfTarget);

static int sequenceuppercase(char *pchSource, int *piSource, char *pchTarget, int *piTarget);

//prehandle the sequence format to pure nucleotide sequence, omit the first line of each file
//and the "\n" character
//Input: pchSource, source sequence; piSource, the number of characters in source sequence; 
//       pchTarget, target sequence of handled; piTarget, the number of characters in target seq;
//Output: 1, correctly handled
//        0, error in handling

int sequenceprehandle(char *pchSource, int *piSource, char *pchTarget, int *piTarget)
{
	int countSource = 0;
	int countTarget = 0;
	int flag = 0;

	while ( countSource < *piSource)
	{

		//delete the first line of each file;
		if (*(pchSource + countSource) == '>')
		{
			flag = 1;	
		}

		else if ( (*(pchSource + countSource) == '\n') && (flag == 1))
		{
			flag = 0;
		}

		//delete the "\n" symbol in sequence data;

		else if ( (flag == 0) && (*(pchSource + countSource) != '\n') )
		{
			*(pchTarget + countTarget) = *(pchSource + countSource);
			countTarget++;
		}

		countSource++;
	}
	
	*piTarget = countTarget;
	
	return(1);
}

int sequenceuppercase(char *pchSource, int *piSource, char *pchTarget, int *piTarget)
{

	int countSource = 0;
	int countTarget = 0;
	int flag = 0;

	while ( countSource < *piSource)
	{

		//delete the first line of each file;
		if (*(pchSource + countSource) == '>')
		{
			flag = 1;	
		}

		else if ( (*(pchSource + countSource) == '\n') && (flag == 1))
		{
			flag = 0;
		}

		//delete the "\n" symbol in sequence data;

		else if ( (flag == 0) && (*(pchSource + countSource) != '\n') )
		{
			if ( *(pchSource + countSource) >96)
				*(pchSource+countSource) -= 32;
			*(pchTarget + countTarget) = *(pchSource + countSource);
			countTarget++;
		}

		countSource++;
	}
	
	*piTarget = countTarget;
	
	return(1);
}

//write the target sequence to file pfTarget;
//input: pfTarget: file pointer of target file;  pchTarget: char pointer of target sequence;
//       piTargetCount: number of Target sequence;
//       piWriteTotal: total number of writed sequence in the target file;
//Output:0 if write error; 1 if write accurate;

int writesequencefile(FILE *pfTarget, char *pchTarget, int *piTargetCount, long long *piWriteTotal)
{

	int iWriteCount;

	iWriteCount = fwrite(pchTarget, sizeof(char), *piTargetCount, pfTarget);
			
	*piWriteTotal += iWriteCount;

	if ( iWriteCount < (signed) (sizeof(char) * (*piTargetCount)) ) 
	{
		printf("Only write %d data because of some error!\n", iWriteCount);
		return(0);
	}
	else 
	{
		printf("write %d data.\n", iWriteCount);
	}

	return(1);
}


//Tranform the original sequence file to another file only contains the nucleotides sequences.
//Input: pfSource: file handler of the source sequence file;
//       pfTarget: file handler of the target handled sequence file;
//Output:0 if error; 1 if accurate;

int filetransform(FILE *pfSource, FILE *pfTarget)
{
	int iRead;
	long long iOffset;
	long long piWriteTotal = 0;
	int *piReadCount = NULL;
	int *piTargetCount = NULL;
	int loop = 0;
	char * pchBuff = NULL;
	char * pchTarget = NULL;

	// Read a sequence segment
		
	iOffset = 0;


	printf("\n");

	printf("***************************\n");

	printf("\n");
	
	do 
	{	
		loop++;

		printf("Loop %d:\n", loop);

		//allocate memory to read the sequence, store the number of read nucleotide

		pchBuff = (char *)malloc(sizeof(char) * BUFF_SIZE);
                piReadCount = (int*) malloc(sizeof(int));
                if ( NULL == pchBuff ) 
		{
			printf("The system can not allocate the memory!\n");
			fclose(pfSource);
			return(0);
		}

		iRead = readfromfile(pfSource, pchBuff, iOffset, BUFF_SIZE, piReadCount);

		if ( iRead == 0 ) 
		{
			fclose(pfSource);
			free(pchBuff);
                        free(piReadCount);
                        return(0);
		}
		else
		{
			if ( *piReadCount < (sizeof(char) * BUFF_SIZE) ) 
			{
				if ( 0 != feof(pfSource) )
				{
					printf("only read %d data because the end of the file!\n", *piReadCount);
				}
			}
			else 
			{
				printf("read %d data.\n", *piReadCount);
			}

			pchTarget = (char *)malloc(sizeof(char) * BUFF_SIZE);
                        piTargetCount = (int*) malloc(sizeof(int));
                        // handle the sequence data
			
			sequenceprehandle(pchBuff, piReadCount, pchTarget, piTargetCount);

			if (writesequencefile(pfTarget, pchTarget, piTargetCount, &piWriteTotal))
			{
				free(pchTarget);
                                free(piTargetCount);
			}
			else 
			{
				fclose(pfTarget);
				free(pchTarget);
				free(piTargetCount);
				free(pchBuff);
				free(piReadCount);
				return(0);
			}
		}

		free(pchBuff);
		free(piReadCount);

		iOffset += BUFF_SIZE;
		
		printf("\n");

	}
	while (iRead == 2);

	printf("\n");

	printf("***************************\n");

	printf("\n");

	printf("Total %d Loops. \n", loop);

	printf("Total write %lld data. \n", piWriteTotal);

	printf("\n");

	printf("\n");


	return(1);
}

// *************************************************************************************
// handle the original sequence file to another file which only contains nucleotides
// *************************************************************************************
int handleformat(char *pchSource, char *pchTarget)
{

	FILE * pfSource = NULL;
	FILE * pfTarget = NULL;
	
	//open the source file
	pfSource = fopen(pchSource, "rb");
	if ( NULL == pfSource )
	{
		printf("%s could not be opened!\n", pchSource);
		return(0);
	}
	
	//open the target file
	pfTarget = fopen(pchTarget, "wb");
	if ( NULL == pfTarget ) 
	{
		printf("%s could not be opened!\n",pchTarget);
		fclose(pfSource);
		return(0);
	}

	//handle the original sequence file to another file which only contains nucleotides
	if ( filetransform(pfSource, pfTarget) == 0)
	{
		fclose(pfSource);
		fclose(pfTarget);
		return(0);
	}
	else
	{	
		fclose(pfSource);
		fclose(pfTarget);
	}
	return(1);
}


int fileuppercase(FILE *pfSource, FILE *pfTarget)
{
	int iRead;
	long long iOffset;
	long long piWriteTotal = 0;
	int *piReadCount = NULL;
	int *piTargetCount = NULL;
	int loop = 0;
	char * pchBuff = NULL;
	char * pchTarget = NULL;

	// Read a sequence segment
		
	iOffset = 0;

	

	printf("\n");

	printf("***************************\n");

	printf("\n");
	
	do 
	{	
		loop++;

		printf("Loop %d:\n", loop);

		//allocate memory to read the sequence, store the number of read nucleotide

		pchBuff = (char *)malloc(sizeof(char) * BUFF_SIZE);
                piReadCount = (int*) malloc(sizeof(int));
                if ( NULL == pchBuff ) 
		{
			printf("The system can not allocate the memory!\n");
			fclose(pfSource);
			return(0);
		}

		iRead = readfromfile(pfSource, pchBuff, iOffset, BUFF_SIZE, piReadCount);

		if ( iRead == 0 ) 
		{
			fclose(pfSource);
			free(pchBuff);
                        free(piReadCount);
                        return(0);
		}
		else
		{
			if ( *piReadCount < (sizeof(char) * BUFF_SIZE) ) 
			{
				if ( 0 != feof(pfSource) )
				{
					printf("only read %d data because the end of the file!\n", *piReadCount);
				}
			}
			else 
			{
				printf("read %d data.\n", *piReadCount);
			}

			pchTarget = (char *)malloc(sizeof(char) * BUFF_SIZE);
                        piTargetCount = (int*) malloc(sizeof(int));
                        // handle the sequence data
			
			sequenceuppercase(pchBuff, piReadCount, pchTarget, piTargetCount);

			if (writesequencefile(pfTarget, pchTarget, piTargetCount, &piWriteTotal))
			{
				free(pchTarget);
                                free(piTargetCount);
			}
			else 
			{
				fclose(pfTarget);
				free(pchTarget);
                                free(piTargetCount);
                                free(pchBuff);
                                free(piReadCount);
                                return(0);
			}
		}

		free(pchBuff);
                free(piReadCount);
                
		iOffset += BUFF_SIZE;
		
		printf("\n");

	}
	while (iRead == 2);

	printf("\n");

	printf("***************************\n");

	printf("\n");

	printf("Total %d Loops. \n", loop);

	printf("Total write %lld data. \n", piWriteTotal);

	printf("\n");

	printf("\n");

	return(1);
}

// *************************************************************************************
// handle the original sequence file to another file which only contains nucleotides
// *************************************************************************************
int handleuppercase(char *pchSource, char *pchTarget)
{

	FILE * pfSource = NULL;
	FILE * pfTarget = NULL;
	
	//open the source file
	pfSource = fopen(pchSource, "rb");
	if ( NULL == pfSource )
	{
		printf("%s could not be opened!\n", pchSource);
		return(0);
	}
	
	//open the target file
	pfTarget = fopen(pchTarget, "wb");
	if ( NULL == pfTarget ) 
	{
		printf("%s could not be opened!\n", pchTarget);
		fclose(pfSource);
		return(0);
	}

	//handle the original sequence file to another file which only contains nucleotides
	if ( fileuppercase(pfSource, pfTarget) == 0)
	{
		fclose(pfSource);
		fclose(pfTarget);
		return(0);
	}
	else
	{	
		fclose(pfSource);
		fclose(pfTarget);
	}
	return(1);
}
