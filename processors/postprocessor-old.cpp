/* 
 * File:   npostproc.cpp
 * Author: Kathryn Hall, Wanjun Gu
 *
 * Created on February 27, 2012, 1:50 PM
 */

#include <cstdlib>
#include <cstring>
#include <stdio.h>

using namespace std;
#define MAXANNOLINE_SIZE 600000
#define SPACER_SIZE 30

/*
 * 
 */

unsigned int findseqlength(char* pchLine)
{
	
	unsigned int i = 0;

	while ((pchLine[i] >= 'A') && (pchLine[i] <= 'Z') || (pchLine[i] >= 'a') && (pchLine[i] <= 'z'))
        {
                i++;
        }

	return(i);
}

bool isannotationline(char* pchLine)
{
	return(pchLine[0] == '>');
}

void parseid(char* pchLine, char* readname)
{
	unsigned int i = 0;
	
	int k = 0;	
	while ((pchLine[i] != ' ') && (pchLine[i] != '\t') && (pchLine[i] != '\n') && (pchLine[i] != '\0'))
	{
		if (pchLine[i] != '>')
		{
			memcpy((char*)(readname+k),(char*)(pchLine+i),1);   
			k++;
		}

		i++;
	}
	
	*(readname+k)='\0';
}

bool parseannotation(char* pchLine, long long& start, long long& end)
{
	unsigned int i = 0;
	
	start = 0;
	
	end = 0;
	
	while ((pchLine[i] != ' ') && (pchLine[i] != '\t') && (pchLine[i] != '\n') && (pchLine[i] != '\0'))
	{
		start = start *10 + pchLine[i] - 48;
		i++;
	}

	i++;
	
	while ((pchLine[i] != ' ') && (pchLine[i] != '\t') && (pchLine[i] != '\n') && (pchLine[i] != '\0'))
	{
		end = end *10 + pchLine[i] - 48;
		i++;
	}
        // make sure the region makes sense
        // the last line of a region file is usually garbage
        return (start < end);
}

int convertRegionToAnnotation(char *pchFasta,char *pchRegion, char *pchOutput)
{
	FILE *pfFasta, *pfRegion, *pfOutput;
	
	pfOutput = fopen(pchOutput, "w");
	
	// read the fasta file for each read 
	if( ( pfFasta = fopen( pchFasta, "r" )) == NULL )
	{
		printf("Cannot find the sequence file.\n");
                fclose(pfOutput);
		return(0);
	}
        
        if( ( pfRegion = fopen(pchRegion, "r" )) == NULL )
        {
                printf("Can not find the region file.\n");
                fclose(pfOutput);
                fclose(pfFasta);
                return(0);
        }
        
        long long readstart=0;
        long long readend=0;
        long long regionstart=0;
        long long regionend=0;
        long long total=0;
        char *readname = (char *)malloc(sizeof(char) * MAXANNOLINE_SIZE);
        char* pchLine = (char *)malloc(sizeof(char) * MAXANNOLINE_SIZE);
        char* pchLine2 = (char *)malloc(sizeof(char) * MAXANNOLINE_SIZE);
        bool header = false;
        
        // The goal here is to keep as little in memory as possible
        while(!feof(pfRegion) && (regionstart <= regionend))
        {
            if (fgets( pchLine, MAXANNOLINE_SIZE, pfRegion ) != NULL)
            {
                // get region info
                if(!parseannotation(pchLine, regionstart, regionend)) 
                    break;
            }
            // Get the next read if needed
            while (!((regionstart >= readstart-SPACER_SIZE)&&(regionend <= readend+SPACER_SIZE)))
            {
                // Get fasta header
                if (!header) // Take care of the initial starting read
                { 
                   fgets( pchLine2, MAXANNOLINE_SIZE, pfFasta ); 
                }
                if (isannotationline(pchLine2)) //it better be true
                {
                    parseid(pchLine2,readname);
                    readstart = total+1;
                    readend = total; // will be incremented
                    header = false;
                }
                // Find the sequence start and end.
                while (!feof(pfFasta) && !header)
                {
                    fgets( pchLine2, MAXANNOLINE_SIZE, pfFasta );
                    if (isannotationline(pchLine2)) //header for the next read, stop
                    {
                        header = true;
                        total += SPACER_SIZE; // spacing between reads
                    }
                    else // increment readend by the length of the sequence
                    {
                        unsigned int len = findseqlength(pchLine2);
                        readend += len;
                        total += len;
                    }
                }
            }
            int start = regionstart - readstart + 1;
            int end = regionend - readstart + 1;
            // Fix ends - Not sure why pClouds is like this.
            // The pClouds regions don't entirely overlap the reads,
            // likely why the combined file has to have N spacers.
            if (start < 1) 
                start = 1;
            if (end > readend) 
                end = readend;
				
            // Print to output file
            fprintf(pfOutput, "%s\t%d\t%d\n", readname, start, end);

            fflush(pfOutput);
				
        }

        fclose(pfRegion);
        fclose(pfOutput);
        fclose(pfFasta);
	return(1);
}

int main(int argc, char *argv[]) {
        if(argc != 4)
	{
		printf("Usage: postprocessor <fastafile> <region file> <annotation file (output)>\n");
		return(1);
	}

	// File name of original fasta
	char* pchFasta = *++argv;
	// region file
	char* pchRegion =  *++argv;
	// Annotation output file
	char* pchOutput = *++argv;
	printf("Converting %s using %s to %s\n",pchFasta,pchRegion, pchOutput);	
	
	convertRegionToAnnotation(pchFasta, pchRegion, pchOutput);
    return 0;
}

