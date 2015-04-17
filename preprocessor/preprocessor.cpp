#include <cstdlib>
#include <stdio.h>


#define MAXANNOLINE_SIZE 600000

bool isannotationline(char* pchLine)
{
	return(pchLine[0] == '>');
}

int combine454fna(char* pch454fa, char* pchOutput)
{
	FILE *pf454fa, *pfOutput;
	
	pfOutput = fopen(pchOutput, "ab");
	
	fprintf(pfOutput, ">Processed reads\n");

	bool include = false;
	
	long int reads = 0;
	
	// read the fa file for each read 
	if( ( pf454fa = fopen( pch454fa, "r" )) == NULL )
	{
		printf("Can not find the Sequence file.\n");
		return(0);
	}
	else 
	{
		char* pchLine = (char *)malloc(sizeof(char) * MAXANNOLINE_SIZE);

		while (!feof(pf454fa))
		{
			if (fgets( pchLine, MAXANNOLINE_SIZE, pf454fa ) != NULL)
			{
				if (isannotationline(pchLine)) //check if we want to include this read
				{
					if (include)
						fprintf(pfOutput, "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN");

					reads++;
					
					include = true;		// we set this to be true all the time for now
				}
				else if (include)	//put all the sequence in upppercase
				{
					int i = 0;
					
					while (pchLine[i] != '\n' && pchLine[i] != '\r')
					{
						if ((pchLine[i] <= 122 ) && (pchLine[i] >= 97))
							pchLine[i] = pchLine[i] - 32;
							
						i++;
					}
					
					fprintf(pfOutput, "%s", pchLine);
				}
			}
		}
		
		free(pchLine);
		
		fclose(pf454fa);
	
		fclose(pfOutput);	
		
		printf("total number of reads: %ld. \n", reads);
	}

	return(1);
}

int main(int argc, char *argv[])
{

	if(argc != 3)
	{
		printf("Usage: preprocessor <fastafile> <output filename>\n");
		return(1);
	}

	// File name of original fasta
	char* pch454fa = *++argv;
	// Pre-processed filename for pClouds input
	char* pchOutput = *++argv;
	printf("Converting %s to %s.\n",pch454fa,pchOutput);	
	combine454fna(pch454fa, pchOutput);
	
	return(1);
}

