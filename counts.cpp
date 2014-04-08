/**
 * STP: This is clearly a cpp file. It used to be named *.c.
 *
 */

#include "macrodefine.h"
#include "counts.h"
#include "directmethod.h"
#include "overlapmethod.h"
#include "mixedmethod01.h"
#include "mixedmethod02.h"
#include "readfile.h"
#include <iostream>
#include <cstring>
#include <cstdio>

using std::cout;
using std::cerr;
using std::endl;
int calculatecounts(const char* pchControlfile)
{
	cout << "Calculating counts" << endl;
	int size, m_nChunksize, m_nCounts;
	unsigned int m_nGenomesize, m_nMemory;

	// defaults
	int directlength = DIRECTLENGTH;
	int overlaplength = OVERLAPLENGTH;

	char pchGenome[MAXFILENAMELENGTH], pchCountfile[MAXFILENAMELENGTH];
				
	// read the control file the get the parameters of the program
	if (!ReadCountsControlfile(pchControlfile, size, m_nChunksize, m_nGenomesize,
			m_nCounts, pchGenome, pchCountfile, m_nMemory))
	{
		cerr << "Open Control file error." <<endl;
		return(-1);
	}

	if (m_nCounts)
	{
		getCutoffsForCountingMethods(overlaplength,directlength,m_nMemory);
// debug
directlength = 16;
overlaplength = 18;
// debug
		if (size <= directlength)
			directmethod(pchGenome, pchCountfile, size, m_nChunksize, m_nGenomesize);	
		else if (size >= overlaplength) // Generate small oligo count files for input into the overlap method.
		{
			int precountsize = overlaplength-1;

			cerr    << "Caution: due to the limited amount of memory available on this system, "
				<< "obtaining oligo counts for oligos of size " << size << " will take a "
				<< "VERY long time.  You may wish to use alternative software for this "
				<< "phase of the analysis." << endl;

			char *pchSmallCountfile = (char *)malloc(sizeof(char) * MAXFILENAMELENGTH);
			if ( !getSmallCountsFile(pchCountfile,size,overlaplength,pchSmallCountfile) )
			{
				// Generate the small oligo count file since it doesn't exist
				if ((precountsize % 2) == 0)
					mixedmethod02(pchGenome, pchSmallCountfile, precountsize, m_nChunksize, m_nGenomesize);
				else
					mixedmethod01(pchGenome, pchSmallCountfile, precountsize, m_nChunksize, m_nGenomesize);
			}

			// Perform the overlap method
			overlapmethod(pchGenome, pchCountfile, size, m_nChunksize, m_nGenomesize);
		}
		else if ((size % 2) == 0)
			mixedmethod02(pchGenome, pchCountfile, size, m_nChunksize, m_nGenomesize);
		else
			mixedmethod01(pchGenome, pchCountfile, size, m_nChunksize, m_nGenomesize);
	}
	
	return(1);
}		

bool getSmallCountsFile(char *pchCountfile, int size, int overlaplength, char *pchSmallCountfile)
// Returns true if the file already exists, otherwise returns false.  The name of the small count file is stored in pchSmallCountfile either way.
{
	// Check for the existence of the proper small oligo count file
	FILE * pfSource = NULL;
	strcpy(pchSmallCountfile, pchCountfile);
	char temp = 48 + (overlaplength-1) / 10;
	strncat(pchSmallCountfile, &temp, 1);
	temp = 48 + (overlaplength-1) %10;
	strncat(pchSmallCountfile, &temp, 1);
	strcat(pchSmallCountfile, ".txt");
	pfSource = fopen(pchSmallCountfile, "rb");

	if(NULL == pfSource) 
		return false; 
	else {
		fclose(pfSource);
		return true;
	}
}

void getCutoffsForCountingMethods(int &overlaplength, int &directlength, int m_nMemory)
{
	// If expected available memory is > 1Gb, then we might switch from the default parameters
	if(m_nMemory > 1000000000)
	{
		// 2*4^size
		int i=0;
		unsigned long memreq = 2;
		// Find maximum direct length oligo size for the available memory
		while(memreq < m_nMemory && memreq != 0)
		{
			memreq=memreq * 4;
			i++;

		}
		directlength = i-1;
		cerr << "Directmethod cutoff =" << directlength << endl;

		// Find the smallest oligo size that necessitates using the overlap method.
		memreq =  74000000; // empirical observation for mixedmethod01, oligosize 14.  Mixedmethod01 uses more memory than mixedmethod02.
		i = 14;
		while(memreq < m_nMemory)
		{
			memreq = memreq * 4.2; // The observed factor is actually closer to 4, but this is a more conservative estimate
			i++;
		}
		overlaplength = i;
		
		cerr << "Overlapmethod cutoff = " << overlaplength << endl;
		// cerr << "Mixedmethod will use variant " << ((size%2)==0?2:1)<< endl;
	}
}
