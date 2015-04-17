#include "macrodefine.h"

#include "directmethod.h"
#include "overlapmethod.h"
#include "mixedmethod01.h"
#include "mixedmethod02.h"
#include "readfile.h"
#include <iostream>
using std::cerr;
using std::endl;
int calculatecounts(const char* pchControlfile)
{
	int size, m_nChunksize, m_nCounts;
	unsigned int m_nGenomesize;

	char pchGenome[MAXFILENAMELENGTH], pchCountfile[MAXFILENAMELENGTH];
				
	// read the control file the get the parameters of the program
	if (!ReadCountsControlfile(pchControlfile, size, m_nChunksize, m_nGenomesize, m_nCounts, pchGenome, pchCountfile))
	{
		cerr << "Open Control file error." <<endl;
		return(-1);
	}

	if (m_nCounts)
	{
		if (size <= DIRECTLENGTH)
			directmethod(pchGenome, pchCountfile, size, m_nChunksize, m_nGenomesize);	
		else if (size >= OVERLAPLENGTH)
			overlapmethod(pchGenome, pchCountfile, size, m_nChunksize, m_nGenomesize);
		else if (size == 16)
			mixedmethod02(pchGenome, pchCountfile, size, m_nChunksize, m_nGenomesize);
		else
			mixedmethod01(pchGenome, pchCountfile, size, m_nChunksize, m_nGenomesize);
	}
	
	return(1);
}		
