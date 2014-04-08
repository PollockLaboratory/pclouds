#include <cstdlib>
#include <stdio.h>
#include <time.h>
#include <algorithm>
#include <ext/algorithm>
#include <math.h>
#include <iostream>
#include <cstring>
#include "macrodefine.h"
#include "readfile.h"
#include "stringhandle.h"

using std::swap;
using std::cerr;
using std::endl;
using std::sort;

static int sbsearch(int n, const cloud_type *argv, const unsigned long& key)
{
	int  m;
	int site = 0;

	while (n >= 1) 
	{
		m = n/2;
        	if (argv[site+m].index == key)
			return(site+m);
		else if (argv[site+m].index > key)
	       		n = m;
		else 
		{
			site = site + m +1; 
			n = n - m - 1;
	        }	
	}

	return(-1);
}

static int sbsearch1(int n, const cloud_type1 *argv, const unsigned long& key)
{
	int  m;
	int site = 0;

	while (n >= 1) 
	{
		m = n/2;
        	if (argv[site+m].index == key)
			return(site+m);
		else if (argv[site+m].index > key)
	       		n = m;
		else 
		{
			site = site + m +1; 
			n = n - m - 1;
	        }	
	}

	return(-1);
}

static int sbsearch2(int n, const cloud_type2 *argv, const unsigned long& key)
{
	int  m;
	int site = 0;

	while (n >= 1) 
	{
		m = n/2;
        	if (argv[site+m].index == key)
			return(site+m);
		else if (argv[site+m].index > key)
	       		n = m;
		else 
		{
			site = site + m +1; 
			n = n - m - 1;
	        }	
	}

	return(-1);
}

static int sbsearch3(int n, const cloud_type3 *argv, const unsigned long& key)
{
	int  m;
	int site = 0;

	while (n >= 1) 
	{
		m = n/2;
        	if (argv[site+m].index == key)
			return(site+m);
		else if (argv[site+m].index > key)
	       		n = m;
		else 
		{
			site = site + m +1; 
			n = n - m - 1;
	        }	
	}

	return(-1);
}

static bool highnumber(cloud_type a, cloud_type b)
{
	if ( a.number > b.number)
		return(1);
	else return(0);
}

static bool highnumber1(cloud_type1 a, cloud_type1 b)
{
	if ( a.number > b.number)
		return(1);
	else return(0);
}

static bool highnumber3(cloud_type3 a, cloud_type3 b)
{
	if ( a.number > b.number)
		return(1);
	else return(0);
}

static bool lowsequence1(cloud_type1 a, cloud_type1 b)
{
	if (a.index < b.index)
		return(1);
	else return(0);
}

static bool lowsequence2(cloud_type2 a, cloud_type2 b)
{
	if (a.index < b.index)
		return(1);
	else return(0);
}

static bool lowsequence3(cloud_type3 a, cloud_type3 b)
{
	if (a.index < b.index)
		return(1);
	else return(0);
}

static int getnumber(char* pchLine)
{
	int i = 0;
	int temp = 0;

	while (pchLine[i] != ' ')
		i++;

	i++;

	while ((pchLine[i] >= ZERO) && (pchLine[i] <= NINE))
	{
		temp = temp * 10 + pchLine[i] - ZERO;
		i++;
	}

	return(temp);
}

static void getcloudandnumber(char* pchLine, int& cloud, int& number)
{
	int i = 0;

	cloud = 0;

	while (pchLine[i] != ' ')
		i++;

	i++;

	while ((pchLine[i] >= ZERO) && (pchLine[i] <= NINE))
	{
		cloud = cloud * 10 + pchLine[i] - ZERO;
		i++;
	}

	while (pchLine[i] != ' ')
		i++;

	i++;

	number = 0;

	while ((pchLine[i] >= ZERO) && (pchLine[i] <= NINE))
	{
		number = number * 10 + pchLine[i] - ZERO;
		i++;
	}
}

static void swapcloud(cloud_type3& a, cloud_type3& b)
{
	swap(a.index, b.index);
	swap(a.number, b.number);
	swap(a.cloud, b.cloud);
	swap(a.extension, b.extension);
}

static void q_sort(cloud_type3* repeats, int left, int right)
{
	int pivot, l, r;
	l = left;
	r = right;
	pivot = repeats[(left+right)/2].number;
	while (l < r)
	{
		while (repeats[l].number > pivot)
			++l;

		while (repeats[r].number < pivot)
			--r;

		if (l >= r)
			break;

		swapcloud(repeats[l], repeats[r]);
		++l;
		--r;		
	}
	
	if (l == r) l++;

	if (left < r)
		q_sort(repeats, left, l-1);
	if (l < right)
		q_sort(repeats, r+1, right);
}

// transform the pattern sequence to the index of the array
// coding method: A=00, C=01, G= 10, T=11,
static void patterntoindex(const char *pchPattern, unsigned long& index, const int& patternsize)
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
}

static void indextopattern(char *pchPattern, const unsigned long& index, const int& patternsize)
{
	int  i;
	unsigned long value;
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
    
    // This needs to be null terminated
    *(pchPattern + patternsize) = '\0';
}

static void getonesubstitutions(const char* pchCore, unsigned long* piSubstitution, int& number)
{
	char alphwbet[5] = "ACGT";

	int size = strlen(pchCore);
	
	char* pchTemp = (char*) malloc(sizeof(char)*(size+1));

	for (int i = 0; i < size; i++)
	{	
		for (int k = 0; k<4; k++)
		{
			if (alphwbet[k] != pchCore[i])
			{
				strcpy(pchTemp,pchCore);
				pchTemp[i] = alphwbet[k];
				patterntoindex(pchTemp, piSubstitution[number], size);
				number++;
			}
		}
	}
	
	char* pchReverse = (char*) malloc(sizeof(char) * (size+1));

	getreversecomplement(pchCore, pchReverse);
	
	if (strcmp(pchReverse, pchCore) != 0)
	{
		for (int i = 0; i < strlen(pchReverse); i++)
		{	
			for (int k = 0; k<4; k++)
			{
				if (alphwbet[k] != pchReverse[i])
				{
					strcpy(pchTemp, pchReverse);
					pchTemp[i] = alphwbet[k];
					patterntoindex(pchTemp, piSubstitution[number], size);
					number++;
				}
			}
		}
	}

	free(pchReverse);
	free(pchTemp);
}

static void gettwosubstitutions(const char* pchCore, unsigned long* piSubstitution, int& number)
{
	char alphwbet[5] = "ACGT";

	int size = strlen(pchCore);
	
	char* pchTemp = (char*) malloc(sizeof(char)*(size+1));

	for (int i = 0; i < size-1; i++)
	{	
		for (int j = i+1;  j< size;  j++)
		{
			for (int k = 0; k<4; k++)
			{
				for (int l = 0; l<4; l++)
				{
					if ( (alphwbet[k] != pchCore[i]) && (alphwbet[l] != pchCore[j]) )
					{
						strcpy(pchTemp, pchCore);
						pchTemp[i] = alphwbet[k];
						pchTemp[j] = alphwbet[l];
						patterntoindex(pchTemp, piSubstitution[number], size);
						number++;
					}	
				}
			}
		}
	}

	char* pchReverse = (char*) malloc(sizeof(char) * (size+1));

	getreversecomplement(pchCore, pchReverse);
	
	for (int i = 0; i < strlen(pchReverse)-1; i++)
	{	
		for (int j = i+1;  j< strlen(pchReverse);  j++)
		{
			for (int k = 0; k<4; k++)
			{
				for (int l = 0; l<4; l++)
				{
					if ( (alphwbet[k] != pchReverse[i]) && (alphwbet[l] != pchReverse[j]) )
					{
						strcpy(pchTemp, pchReverse);
						pchTemp[i] = alphwbet[k];
						pchTemp[j] = alphwbet[l];
						patterntoindex(pchTemp, piSubstitution[number], size);
						number++;
					}	
				}
			}
		}
	}

	free(pchReverse);
	free(pchTemp);
}

static void getthreesubstitutions(char* pchCore, unsigned long* piSubstitution, int& number)
{
	char alphwbet[5] = "ACGT";

	int size = strlen(pchCore);
	
	char* pchTemp = (char*) malloc(sizeof(char)*(size+1));

	for (int i = 0; i < size-2; i++)
	{	
		for (int j = i+1;  j< size-1;  j++)
		{
			for (int m = j+1; m < size; m++)
			{
				for (int k = 0; k<4; k++)
				{
					for (int l = 0; l<4; l++)
					{
						for (int n = 0; n<4; n++)
						{
							if ( (alphwbet[k] != pchCore[i]) && (alphwbet[l] != pchCore[j]) &&((alphwbet[n] != pchCore[m])))
							{
								strcpy(pchTemp, pchCore);
								pchTemp[i] = alphwbet[k];
								pchTemp[j] = alphwbet[l];
								pchTemp[m] = alphwbet[n];
								patterntoindex(pchTemp, piSubstitution[number], size);
								number++;
							}
						}	
					}
				}
			}
		}
	}
	
	char* pchReverse = (char*) malloc(sizeof(char) * (size+1));

	getreversecomplement(pchCore, pchReverse);
	
	for (int i = 0; i < size-2; i++)
	{	
		for (int j = i+1;  j< strlen(pchReverse)-1;  j++)
		{
			for (int m = j+1; m < strlen(pchReverse); m++)
			{
				for (int k = 0; k<4; k++)
				{
					for (int l = 0; l<4; l++)
					{
						for (int n = 0; n<4; n++)
						{
							if ( (alphwbet[k] != pchReverse[i]) && (alphwbet[l] != pchReverse[j]) &&((alphwbet[n] != pchReverse[m])))
							{
								strcpy(pchTemp, pchReverse);
								pchTemp[i] = alphwbet[k];
								pchTemp[j] = alphwbet[l];
								pchTemp[m] = alphwbet[n];
								patterntoindex(pchTemp, piSubstitution[number], size);
								number++;
							}
						}	
					}
				}
			}
		}
	}

	free(pchReverse);
	free(pchTemp);
}

// specific to use 16mers, use some different kind of coding method for each 16mer
static void buildbitmatrix(bit_matrix& PatternMatrix, int*  matrixsize, const int& patternsize)
{
	unsigned long index;

	char* pchPattern = (char*) malloc(sizeof(char) * (patternsize/2+1));

	char* pchReverse = (char*) malloc(sizeof(char) * (patternsize/2+1));

	// each line correspond to each pattern with index of the last 8 chars.
	for (int i = 0; i < pow(4, patternsize/2); i++)
	{
		// index equal to i, get the last 8 characters of the pattern from index
		index = i;

		indextopattern(pchPattern, index, patternsize/2);
		
		// get the reverse complement of the last 8 characters, then get the index of the reverse complement LR
		//the size of the vector[i] is LR+1;

		getreversecomplement(pchPattern, pchReverse);

		patterntoindex(pchReverse, index, patternsize/2);

		PatternMatrix[i].reserve(index+1); 

		matrixsize[i] = index+1;
	}

	free(pchPattern);
	free(pchReverse);

}

// specific to use 16mers, use some different kind of coding method for each 16mer
static void getpatternmatrixindex(char* pchPattern, unsigned long& i, unsigned long& j)
{
        int size = strlen(pchPattern);

	char* pchTemp = (char*) malloc( sizeof(char) * (size/2 + 1));

	// get the index of last 8 chars;
	strright(size/2, pchPattern, pchTemp);
	
	patterntoindex(pchTemp, i, size/2);
        
	// get the index of the first 8 chars;
	strleft(size/2, pchPattern, pchTemp);

	patterntoindex(pchTemp, j, size/2);

	free(pchTemp);
}

static bool isrepeatregion(vector<int>& occurrence, const int& percent)
{
	vector<int>::iterator pOccurrence;

	int count = 0;

	for (pOccurrence = occurrence.begin(); pOccurrence < occurrence.end(); pOccurrence++)
	{
		if (*(pOccurrence) > 0)
			count++;
	}

	if ( (int) ((float)count / (float)(occurrence.size()) *100) >= percent)
		return(1);
	else
		return(0);
}

// for 16mers calculation, oligos and its reverse complement are set to one
// read this kind of oligo count sets, get number1 and number2
static int readrepeatnumber(char* pchrepeatfile, int& number1, int& number2, const int& size, const int& m_nStep1,
const int& m_nEndthreshold, const int& m_nCopy)
{
	FILE *pfRepeatfile;
	char *pchLine;
	char *pchTemp;

	number1 = 0;
	number2 = 0;
	
	int threshold;

	if (m_nStep1 < m_nEndthreshold)
		threshold = m_nStep1;
	else threshold = m_nEndthreshold;

	if( ( pfRepeatfile = fopen( pchrepeatfile, "r" )) == NULL )
	{
		printf("Can not find the repeat file: %s.\n", pchrepeatfile);
		return(0);
	}
	else 
	{
		pchLine = (char *)malloc(sizeof(char) * MAXLINECHAR);

		pchTemp = (char*) malloc(sizeof(char) * (size+1));

		while (!feof(pfRepeatfile))
		{
			if( fgets( pchLine, MAXLINECHAR, pfRepeatfile ) != NULL)
			{
				strncpy( pchTemp, pchLine, size); 
				*(pchTemp + size) = '\0';

				if ( (isonessr(pchTemp) == 0) && (istwossr(pchTemp) == 0) && (isthreessr(pchTemp) == 0) && (isfourssr(pchTemp) == 0) )
				{
					int occur = getnumber(pchLine);
					
					if (occur >= threshold)
						number1++;
					else if (occur >= m_nCopy)
						number2++;
				}
			}
		}

		free(pchLine);
		free(pchTemp);
		fclose( pfRepeatfile);
	}

	return(1);
}

// read the oligos which is in mainclouds into repeats1 array
static int readmainoligos(char *pchrepeatfile, cloud_type3* repeats1, int& number1, const int& size, const int& m_nStep1,
const int& m_nEndthreshold)
{
	FILE *pfRepeatfile;
	char *pchLine;
	char *pchTemp, *pchReverse;
	unsigned long index;

	number1 = 0;
	
	int threshold;

	if (m_nStep1 < m_nEndthreshold)
		threshold = m_nStep1;
	else threshold = m_nEndthreshold;

	if( ( pfRepeatfile = fopen( pchrepeatfile, "r" )) == NULL )
	{
		printf("Can not find the repeat file: %s.\n", pchrepeatfile);
		return(0);
	}
	else 
	{
		pchLine = (char *)malloc(sizeof(char) * MAXLINECHAR);
		pchTemp = (char*) malloc(sizeof(char) * (size+1));


		while (!feof(pfRepeatfile))
		{
			if( fgets( pchLine, MAXLINECHAR, pfRepeatfile ) != NULL)
			{
				strncpy( pchTemp, pchLine, size); 
				*(pchTemp + size) = '\0';

				if ( (isonessr(pchTemp) == 0) && (istwossr(pchTemp) == 0) && (isthreessr(pchTemp) == 0) && (isfourssr(pchTemp) == 0) )
				{
					int occur = getnumber(pchLine);
					
					if (occur >= threshold)
					{
						patterntoindex(pchTemp, repeats1[number1].index, size);
						repeats1[number1].number = occur;
						number1++;
					}
				}
			}
		}

		free(pchLine);
		free(pchTemp);
		fclose( pfRepeatfile);
	}

	return(1);
}

// build main clouds based on the algorithm
static void buildmainpcloud(cloud_type3* repeats1, char* pchOutput, const int& number1, const int& size, const int& m_nEndthreshold)
{
	//set the pointer to the most often observed repeat;
        int count = 0;
	char* pchCore = (char*) malloc(sizeof(char) * (size+1));
        char* pchRepeat = (char*) malloc(sizeof(char) * (size+1));
        int *totalnumber = (int *) malloc(sizeof(int) * MAXCLOUD);
        char** core = (char**) malloc(sizeof(char) * MAXCLOUD * (size+1));
        int *totalsize = (int*) malloc(sizeof(int) * MAXCLOUD);
        int pcloudnumber = 0;
	int distance;	

	unsigned long *piRepeatThreesubstitution =  (unsigned long*) malloc(2*sizeof(unsigned long)* (( 3 * size) + (9*size*(size-1)/2)+ (27*size*(size-1)*(size-2)/2)));
        unsigned long *piCoreThreesubstitution =  (unsigned long*) malloc(2*sizeof(unsigned long)* (( 3 * size) + (9*size*(size-1)/2)+ (27*size*(size-1)*(size-2)/2)));
        int iRepeatThreesubstitution=0;
	int iCoreThreesubstitution=0;
	int result=0;
        

	unsigned long* repeats3 = (unsigned long*) malloc(sizeof(unsigned long) * number1);
        sort(repeats1, repeats1+number1, highnumber3);

	int k;

	for (k = 0; k< number1; k++)
		repeats3[k] = repeats1[k].index;

	sort(repeats1, repeats1+number1, lowsequence3);

	int coreindex = sbsearch3(number1, repeats1, repeats3[count]);

	while ((count< number1) && (repeats1[coreindex].number >= m_nEndthreshold))
	{
		pcloudnumber++;

		// get the core sequence of the pcloud
		core[pcloudnumber -1] = (char*) malloc(sizeof(char) * (size+1));
                
		indextopattern(core[pcloudnumber -1], repeats3[count], size);
		strcpy(pchCore, core[pcloudnumber -1]);

		repeats1[coreindex].cloud = pcloudnumber;
		repeats1[coreindex].extension = 1;
		totalsize[pcloudnumber -1] = 1;
		totalnumber[pcloudnumber -1] = repeats1[coreindex].number;

		iCoreThreesubstitution = 0;
		getonesubstitutions(pchCore, piCoreThreesubstitution, iCoreThreesubstitution);
		gettwosubstitutions(pchCore, piCoreThreesubstitution, iCoreThreesubstitution);
		getthreesubstitutions(pchCore, piCoreThreesubstitution, iCoreThreesubstitution);
                // get all the repeats in this pclouds within 3 subs
		for (int i = 0; i < iCoreThreesubstitution; i++)
		{
			result = sbsearch3(number1, repeats1, piCoreThreesubstitution[i]);

			if (result >=0)
			{
				if (repeats1[result].cloud == 0)
				{
					repeats1[result].cloud = pcloudnumber;
					indextopattern(pchRepeat, repeats1[result].index, size);
					totalsize[pcloudnumber -1]++;
					totalnumber[pcloudnumber -1] += repeats1[result].number;
				}

				// get another 1-sub, 2-sub and 3-sub with repeat wbove some threshold
				if (repeats1[result].extension == 0)
				{
					repeats1[result].extension = 1;

					iRepeatThreesubstitution = 0;
					indextopattern(pchRepeat, repeats1[result].index, size);
					getonesubstitutions(pchRepeat, piRepeatThreesubstitution, iRepeatThreesubstitution);
                                        
                                        gettwosubstitutions(pchRepeat, piRepeatThreesubstitution, iRepeatThreesubstitution);
					                                        
                                        getthreesubstitutions(pchRepeat, piRepeatThreesubstitution, iRepeatThreesubstitution);

					for (int j = 0; j<iRepeatThreesubstitution; j++)
					{
						result = sbsearch3(number1, repeats1, piRepeatThreesubstitution[j]);

						if (result >= 0) 
						{
							if (repeats1[result].cloud == 0)
							{
								repeats1[result].cloud = pcloudnumber;
								indextopattern(pchRepeat, repeats1[result].index, size);
								totalsize[pcloudnumber -1]++;
								totalnumber[pcloudnumber -1] += repeats1[result].number;
							}
						}
					}
				}
			}
		}

		// move to next core repeats
		while((repeats1[coreindex].cloud !=0) && (count < number1))
		{
			count++;
			
			if (count <= number1-1)
				coreindex =  sbsearch3(number1, repeats1, repeats3[count]);
		}
	}

	free(piRepeatThreesubstitution);
	free(piCoreThreesubstitution);
	//output the main clouds information
	FILE *pfOutput = fopen(pchOutput, "wb");

	for (k = 0; k < pcloudnumber; k++)
		fprintf(pfOutput, "%d\t%d\t%d\t%s\n", k+1, totalsize[k], totalnumber[k], core[k]);

	free(repeats3);
	free(pchRepeat);
	free(pchCore);
	fclose(pfOutput);

	for (k = 0; k < pcloudnumber; k++)
        {
		free(core[k]);
        }

	free(core);
        free(totalnumber);
        free(totalsize);
        
}

//output the mainclouds assignments of each oligo in main clouds
static void outputmainclouds(cloud_type3* Repeat1, const char* pchResult, const int& number1, const int& size)
{
	FILE *pfResult = fopen(pchResult, "wb");
	
	char *pchPattern = (char*) malloc(sizeof(char) * (size+1));
        
	sort(Repeat1, Repeat1 + number1, highnumber3);

	for (int i=0; i< number1; i++)
	{
		indextopattern(pchPattern, Repeat1[i].index, size);
        	fprintf(pfResult, "%s %d %d\n", pchPattern, Repeat1[i].cloud, Repeat1[i].number);
	}

	fclose(pfResult);
	free(pchPattern);
}

// read main clouds assignments information into three different sets
static int readmainclouds(const char* pchMainCloudassign, cloud_type2* pCloudsA, int& numberA, cloud_type2* pCloudsB, int&
numberB, cloud_type2* pCloudsC, int& numberC, const int& size, const int& m_nStep1, const int& m_nStep2, const 
int& m_nStep3)
{
	FILE *pfCloudassign;
	char *pchLine;
	char *pchTemp;

	numberA = 0;
	numberB = 0;
	numberC = 0;

	int count;
	int cloud;
	unsigned long index;
	
	if( ( pfCloudassign = fopen( pchMainCloudassign, "r" )) == NULL )
	{
		printf("Can not find the cloud assign file: %s.\n", pchMainCloudassign);
		return(0);
	}
	else 
	{
		pchLine = (char *)malloc(sizeof(char) * MAXLINECHAR);
                pchTemp = (char*) malloc(sizeof(char) * (size+1));
                while (!feof(pfCloudassign))
		{
			if( fgets( pchLine, MAXLINECHAR, pfCloudassign ) != NULL)
			{
				getcloudandnumber(pchLine, cloud, count);

				if (count >= m_nStep3)
				{
					pCloudsA[numberA].cloud = cloud;
					strncpy( pchTemp, pchLine, size); 
					*(pchTemp + size) = '\0';
					patterntoindex(pchTemp, index, size); 
					pCloudsA[numberA].index = index;
					numberA++;
				}
				else if (count >= m_nStep2)
				{
					pCloudsB[numberB].cloud = cloud;
					strncpy( pchTemp, pchLine, size); 
					*(pchTemp + size) = '\0';
					patterntoindex(pchTemp, index, size); 
					pCloudsB[numberB].index = index;
					numberB++;
				}
				else if (count >= m_nStep1)
				{
					pCloudsC[numberC].cloud = cloud;
					strncpy( pchTemp, pchLine, size); 
					*(pchTemp + size) = '\0';
					patterntoindex(pchTemp, index, size); 
					pCloudsC[numberC].index = index;
					numberC++;
				}
			}
		}
		free(pchLine);
		free(pchTemp);
		fclose( pfCloudassign);
	}

	sort(pCloudsA, pCloudsA+numberA, lowsequence2);
	sort(pCloudsB, pCloudsB+numberB, lowsequence2);
	sort(pCloudsC, pCloudsC+numberC, lowsequence2);

	return(1);
}

// assign the oligos in accessory regions into P clouds constructed in former step
static int buildaccessarypcloud(char *pchrepeatfile, char* pchOutput, cloud_type2* pCloudsA, const int& numberA, cloud_type2*
pCloudsB, const int& numberB, cloud_type2* pCloudsC, const int& numberC, int& number2, const int& size, const int& m_nStep1,
const int& m_nEndthreshold, const int& m_nCopy)
{
	FILE *pfRepeatfile;
	char *pchLine;
	char *pchTemp, *pchReverse;

	unsigned long *piRepeatOnesubstitution =  (unsigned long*) malloc(2*sizeof(unsigned long) * ( 3 * size));
        unsigned long *piRepeatTwosubstitution =  (unsigned long*) malloc(2*sizeof(unsigned long)* ( ( 3 * size) + (9*size*(size-1)/2)));
        unsigned long *piRepeatThreesubstitution =  (unsigned long*) malloc(2*sizeof(unsigned long)* (( 3 * size) + (9*size*(size-1)/2)+ (27*size*(size-1)*(size-2)/(3*2))));
        int iRepeatOnesubstitution;
	int iRepeatTwosubstitution;
	int iRepeatThreesubstitution;
	int number;

	number2 = 0;

	int threshold;

	if (m_nStep1 < m_nEndthreshold)
		threshold = m_nStep1;
	else threshold = m_nEndthreshold;

	int copynumber;
	int result;
	int cloudindex;

	FILE* pfOutput= fopen(pchOutput, "wb");

	if( ( pfRepeatfile = fopen( pchrepeatfile, "r" )) == NULL )
	{
		printf("Can not find the repeat file: %s.\n", pchrepeatfile);
		return(0);
	}
	else 
	{
		pchLine = (char *)malloc(sizeof(char) * MAXLINECHAR);
                pchTemp = (char*) malloc(sizeof(char) * (size+1));
                while (!feof(pfRepeatfile))
		{
			number = 0;

			cloud_type2* repeatacc = (cloud_type2*) malloc(sizeof(cloud_type2) * MAXACCMER);
                        // read in a chunk of accessory repeats
			while ((!feof(pfRepeatfile)) && (number < MAXACCMER))
			{
				if( fgets( pchLine, MAXLINECHAR, pfRepeatfile ) != NULL)
				{
					strncpy( pchTemp, pchLine, size); 
					*(pchTemp + size) = '\0';

					if ( (isonessr(pchTemp) == 0) && (istwossr(pchTemp) == 0) && (isthreessr(pchTemp) == 0) && (isfourssr(pchTemp) == 0) )
					{
						int occur = getnumber(pchLine);
					
						if ( (occur < threshold) && (occur >= m_nCopy))
						{
							// insert the acc into repeats[number]
							patterntoindex(pchTemp, repeatacc[number].index, size); 
							number++;
						}
					}
				}
			}

			sort(repeatacc, repeatacc + number, lowsequence2);

			// build the accessory p clouds for these repeats chunck
			for (int i =0; i < numberA; i++)
			{
				indextopattern(pchTemp, pCloudsA[i].index, size);

				// get the 3-mutation sets of pchTemp;
				iRepeatThreesubstitution = 0;

				getthreesubstitutions(pchTemp, piRepeatThreesubstitution, iRepeatThreesubstitution);
				gettwosubstitutions(pchTemp, piRepeatThreesubstitution, iRepeatThreesubstitution);
				getonesubstitutions(pchTemp, piRepeatThreesubstitution, iRepeatThreesubstitution);

				for (int count = 0; count < iRepeatThreesubstitution; count++)
				{
					result = sbsearch2(number, repeatacc, piRepeatThreesubstitution[count]);
		
					if (result >= 0) 
					{
						if (repeatacc[result].cloud == 0)
							repeatacc[result].cloud = pCloudsA[i].cloud;
					}
				}
			}

			for (int i =0; i < numberB; i++)
			{
				indextopattern(pchTemp, pCloudsB[i].index, size);

				// get the 2-mutation sets of pchTemp;
				iRepeatTwosubstitution = 0;

				gettwosubstitutions(pchTemp, piRepeatTwosubstitution, iRepeatTwosubstitution);
				getonesubstitutions(pchTemp, piRepeatTwosubstitution, iRepeatTwosubstitution);

				for (int count = 0; count < iRepeatTwosubstitution; count++)
				{
					result = sbsearch2(number, repeatacc, piRepeatTwosubstitution[count]);
			
					if (result >= 0) 
					{
						if (repeatacc[result].cloud == 0)
							repeatacc[result].cloud = pCloudsB[i].cloud;
					}
				}
			}

			for (int i = 0; i < numberC; i++)
			{
				indextopattern(pchTemp, pCloudsC[i].index, size);

				// get the 1-mutation sets of pchTemp;
				iRepeatOnesubstitution = 0;

				getonesubstitutions(pchTemp, piRepeatOnesubstitution, iRepeatOnesubstitution);

				for (int count = 0; count < iRepeatOnesubstitution; count++)
				{
					result = sbsearch2(number, repeatacc, piRepeatOnesubstitution[count]);
			
					if (result >= 0) 
					{
						if (repeatacc[result].cloud == 0)
							repeatacc[result].cloud = pCloudsC[i].cloud;
					}
				}
			}

			//output the accessory p clouds information for these repeats
			for (int i = 0; i < number; i++)
			{
				if (repeatacc[i].cloud != 0)
				{
					indextopattern(pchTemp, repeatacc[i].index, size);
					fprintf(pfOutput, "%s %d\n", pchTemp, repeatacc[i].cloud);
				}
			}

			free(repeatacc);
		}

		free(pchLine);
		free(pchTemp);
		fclose(pfRepeatfile);
		fclose(pfOutput);
	}
        free(piRepeatOnesubstitution);
        free(piRepeatTwosubstitution);
        free(piRepeatThreesubstitution);

	return(1);
}

//for 16mers, build the bool patternmatrix to ascertain which oligo is in the P clouds
static int readclouds(const char* pchMainCloudassign, const char* pchAccCloudassign, bit_matrix& PatternMatrix, int* matrixsize, const int& size)
{
	FILE *pfMainCloudassign;
	FILE *pfAccCloudassign;
	char *pchLine;
	char *pchTemp;
	char* pchReverse;
	unsigned long left, right; 
	
	if( ( pfMainCloudassign = fopen( pchMainCloudassign, "r" )) == NULL )
	{
		printf("Can not find the cloud assign file: %s.\n", pchMainCloudassign);
		return(0);
	}
	else 
	{
		pchLine = (char *)malloc(sizeof(char) * MAXLINECHAR);
                pchTemp = (char*) malloc(sizeof(char) * (size+1));
                pchReverse = (char*) malloc(sizeof(char) * (size +1));
                
		while (fgets( pchLine, MAXLINECHAR, pfMainCloudassign ) != NULL)
		{
			if (strlen(pchLine) < size){ printf("String not long enough %s\n",pchLine);}
			strncpy( pchTemp, pchLine, size); 
			*(pchTemp + size) = '\0';
                       getpatternmatrixindex(pchTemp, right, left);

			if (left >= matrixsize[right])
			{
                                getreversecomplement(pchTemp, pchReverse);
				strcpy(pchTemp, pchReverse);
				getpatternmatrixindex(pchTemp, right, left);
                                
			}
                        PatternMatrix[right][left] = true;
			
		}
		free(pchLine);
                free(pchTemp);
                free(pchReverse);
                fclose(pfMainCloudassign);
	}

	if( ( pfAccCloudassign = fopen( pchAccCloudassign, "r" )) == NULL )
	{
		printf("Can not find the cloud assign file: %s.\n", pchAccCloudassign);
		return(0);
	}
	else 
	{
		pchLine = (char *)malloc(sizeof(char) * MAXLINECHAR);
                pchTemp = (char*) malloc(sizeof(char) * (size+1));
                pchReverse = (char*) malloc(sizeof(char) * (size +1));
                
		while (fgets( pchLine, MAXLINECHAR, pfAccCloudassign ) != NULL)
		{
				strncpy( pchTemp, pchLine, size); 
				*(pchTemp + size) = '\0';

				getpatternmatrixindex(pchTemp, right, left);

				if (left >= matrixsize[right])
				{ 
					getreversecomplement(pchTemp, pchReverse);
					strcpy(pchTemp, pchReverse);
					getpatternmatrixindex(pchTemp, right, left); 
				}

				PatternMatrix[right][left] = true;
			
		}
		free(pchLine);
                free(pchTemp);
                free(pchReverse);
                fclose(pfAccCloudassign);
	}
	
	return(1);
}

//for 16mers, annotate the genome based upon whether it is in the p clouds, and get the repeat regions
//but no p cloud pattern included in that annotation.
static int GenomeScanAndIdentify(const char* pchGenome, const char* pchOutfile, const char* pchRegionfile, const bit_matrix&
PatternMatrix, const int* matrixsize, const int& size, const int& windowsize, const int& percent, const int& m_nChunksize,
const unsigned int& m_nGenomesize)
{
	FILE *pfGenome;

	char *pchPattern;
	int *piReadCount = (int*) malloc(sizeof(int));
        char *pchSequence = (char*) malloc(sizeof(char) * m_nChunksize);
        pchPattern = (char *) malloc(sizeof(char) * (size+1));
        char *pchReverse = (char *) malloc(sizeof(char) * (size+1));
        long long iOffset = 0;
	int iRead;

	FILE *pfOut;
	
	pfOut = fopen(pchOutfile, "wb");

	int patternnumber;
	unsigned long index;
	unsigned long left; 
	unsigned long right;

	FILE *pfRegion;

	vector<int> occurrence;

	long long icount = 0;
	long long iStart = 0, iEnd = 0;
	long long iFormerStart = 0, iFormerEnd = 0;
	long long totalsize = 0;

	pfRegion = fopen(pchRegionfile, "wb");

	if( ( pfGenome = fopen( pchGenome, "r" )) == NULL )
	{
		printf("Can not find the genome file: %s.\n", pchGenome);
		return(0);
	}
	else 
	{ 
		do 
		{	
			if ( NULL == pchSequence ) 
			{
				printf("The system can not allocate the memory!\n");
				free(pchPattern);
                                free(pchReverse);
                                free(piReadCount);
                                return(0);
			}
		
			if (iOffset + m_nChunksize -1 <= m_nGenomesize)
				iRead = readfromfile(pfGenome, pchSequence, iOffset, m_nChunksize, piReadCount);
			else 
				iRead = readfromfile(pfGenome, pchSequence, iOffset, m_nGenomesize - iOffset, piReadCount);
                        
			if ( iRead == 0 ) 
			{
                                
                                free(pchPattern);
				free(pchSequence);
				free(pchReverse);
				free(piReadCount);
				return(0);
			}
			else
			{
				for (int count =0; count <= *piReadCount - size; count ++)
				{
					if ( (count % 50 == 0) && (count>0))
						fprintf(pfOut, "\n");

					getsubstring(pchSequence, pchPattern, count, count + size - 1);
					pchPattern[size] = '\0';
					patternnumber = 0;

					if (issegmentvalid(pchPattern)==1)
					{
						getpatternmatrixindex(pchPattern, right, left);

						if (left >= matrixsize[right])
						{
							getreversecomplement(pchPattern, pchReverse);
							strcpy(pchPattern, pchReverse);
							getpatternmatrixindex(pchPattern, right, left);
						}

						if (PatternMatrix[right][left])
							patternnumber = 1;
					}

					//output the p cloud annotation file;
					fprintf(pfOut, "%d ", patternnumber);

					//get the annotated repeat region
					icount++;

					if (icount <= windowsize)
					{
						occurrence.push_back(patternnumber);

						if (icount == windowsize)
						{
							if (isrepeatregion(occurrence, percent))
							{
								iStart = icount - windowsize +1;
								iFormerStart = iStart;
								iEnd = 0;
							}
						}
					}
					else
					{
						occurrence.erase(occurrence.begin());
						occurrence.push_back(patternnumber);

						if (isrepeatregion(occurrence, percent))
						{
							if (iStart == 0)
							{
								iStart = icount - windowsize+1;

								if (iStart > iFormerEnd + 1)
								{
									if ((iFormerEnd != 0) && (iFormerStart !=0))
									{
										fprintf(pfRegion, "%lld %lld\n", iFormerStart, iFormerEnd);
										totalsize += iFormerEnd - iFormerStart + 1;
									}

									iFormerStart = iStart;
								}

								iEnd = 0;
							}
						}
						else
						{
							if ((iStart != 0) && (iEnd == 0))
							{
								int i = windowsize -1;

								while (occurrence[i] == 0)
									i--;
		
								iEnd = icount + size - 2 - (windowsize-1 - i);
								iFormerEnd = iEnd;
								iStart = 0;
							}
						}
					}
				}
			}
			iOffset += m_nChunksize - size +1;
		}
		while ((iOffset < m_nGenomesize - size +1) && (iRead == 2)); 

		fclose( pfGenome);
	}

	fprintf(pfRegion, "%lld %lld\n", iFormerStart, iFormerEnd);

	totalsize += iFormerEnd - iFormerStart + 1;

	fprintf(pfRegion, "%lld", totalsize);
	fflush(pfRegion);
	fclose(pfRegion);
	fclose(pfOut);
        free(pchPattern);
	free(pchSequence);
	free(piReadCount);
	free(pchReverse);

	return(1);
}

//for 16mers, build the bool patternmatrix to ascertain which oligo is in the P clouds
static int readclouds1(const char* pchMainCloudassign, const char* pchAccCloudassign, bitvector& PatternVector, const int& size)
{
	FILE *pfMainCloudassign;
	FILE *pfAccCloudassign;
	char *pchLine;
	char *pchTemp;
	char* pchReverse;
	unsigned long index; 
	
	if( ( pfMainCloudassign = fopen( pchMainCloudassign, "r" )) == NULL )
	{
		printf("Can not find the cloud assign file: %s.\n", pchMainCloudassign);
		return(0);
	}
	else 
	{
		pchLine = (char *)malloc(sizeof(char) * MAXLINECHAR);
                pchTemp = (char*) malloc(sizeof(char) * (size+1));
                pchReverse = (char*) malloc(sizeof(char) * (size +1));
                
	
		while (!feof(pfMainCloudassign))
		{
			if( fgets( pchLine, MAXLINECHAR, pfMainCloudassign ) != NULL)
			{
				strncpy( pchTemp, pchLine, size); 
				*(pchTemp + size) = '\0';

				patterntoindex(pchTemp, index, size);
				PatternVector[index] = true;

				getreversecomplement(pchTemp, pchReverse);
				patterntoindex(pchReverse, index, size);

				PatternVector[index] = true;
			}
		}

		free(pchLine);
		free(pchTemp);
		free(pchReverse);
		fclose( pfMainCloudassign);
	}

	if( ( pfAccCloudassign = fopen( pchAccCloudassign, "r" )) == NULL )
	{
		printf("Can not find the cloud assign file: %s.\n", pchAccCloudassign);
		return(0);
	}
	else 
	{
		pchLine = (char *)malloc(sizeof(char) * MAXLINECHAR);
                pchTemp = (char*) malloc(sizeof(char) * (size+1));
                pchReverse = (char*) malloc(sizeof(char) * (size +1));
                
		while (!feof(pfAccCloudassign))
		{
			if( fgets( pchLine, MAXLINECHAR, pfAccCloudassign ) != NULL)
			{
				strncpy( pchTemp, pchLine, size); 
				*(pchTemp + size) = '\0';

				patterntoindex(pchTemp, index, size);
				PatternVector[index] = true;

				getreversecomplement(pchTemp, pchReverse);
				patterntoindex(pchReverse, index, size);

				PatternVector[index] = true;
			}
		}

		free(pchLine);
		free(pchTemp);
		free(pchReverse);
		fclose(pfAccCloudassign);
	}
	
	return(1);
}

//for 16mers, annotate the genome based upon whether it is in the p clouds, and get the repeat regions
//but no p cloud pattern included in that annotation.
static int GenomeScanAndIdentify1(const char* pchGenome, const char* pchOutfile, const char* pchRegionfile, const bitvector& PatternVector, const int& size, const int& windowsize, const int& percent, const int& m_nChunksize, const unsigned int& m_nGenomesize)
{
	FILE *pfGenome;

	char *pchPattern;
	int *piReadCount = (int*) malloc(sizeof(int));
        char *pchSequence = (char*) malloc(sizeof(char) * m_nChunksize);
        pchPattern = (char *) malloc(sizeof(char) * (size+1));
        char *pchReverse = (char *) malloc(sizeof(char) * (size+1));
        
        long long iOffset = 0;
	int iRead;

	FILE *pfOut;
	
	pfOut = fopen(pchOutfile, "wb");

	int patternnumber;
	unsigned long index;

	FILE *pfRegion;

	vector<int> occurrence;

	long long icount = 0;
	long long iStart = 0, iEnd = 0;
	long long iFormerStart = 0, iFormerEnd = 0;
	long long totalsize = 0;

	pfRegion = fopen(pchRegionfile, "wb");

	if( ( pfGenome = fopen( pchGenome, "r" )) == NULL )
	{
		printf("Can not find the genome file: %s.\n", pchGenome);
		return(0);
	}
	else 
	{
		do 
		{	
			if ( NULL == pchSequence ) 
			{
				printf("The system can not allocate the memory!\n");
				free(pchPattern);
				free(pchReverse);
				free(piReadCount);
                                return(0);
			}
		
			if (iOffset + m_nChunksize -1 <= m_nGenomesize)
				iRead = readfromfile(pfGenome, pchSequence, iOffset, m_nChunksize, piReadCount);
			else 
				iRead = readfromfile(pfGenome, pchSequence, iOffset, m_nGenomesize - iOffset, piReadCount);

			if ( iRead == 0 ) 
			{
				free(pchPattern);
				free(pchSequence);
				free(pchReverse);
				free(piReadCount);
				return(0);
			}
			else
			{
				for (int count =0; count <= *piReadCount - size; count ++)
				{
					if ( (count % 50 == 0) && (count>0))
						fprintf(pfOut, "\n");

					getsubstring(pchSequence, pchPattern, count, count + size - 1);
					pchPattern[size] = '\0';
					patternnumber = 0;

					if (issegmentvalid(pchPattern)==1)
					{
						patterntoindex(pchPattern, index, size);

						if (PatternVector[index])
							patternnumber = 1;
					}

					//output the p cloud annotation file;
					fprintf(pfOut, "%d ", patternnumber);

					//get the annotated repeat region
					icount++;

					if (icount <= windowsize)
					{
						occurrence.push_back(patternnumber);

						if (icount == windowsize)
						{
							if (isrepeatregion(occurrence, percent))
							{
								iStart = icount - windowsize +1;
								iFormerStart = iStart;
								iEnd = 0;
							}
						}
					}
					else
					{
						occurrence.erase(occurrence.begin());
						occurrence.push_back(patternnumber);

						if (isrepeatregion(occurrence, percent))
						{
							if (iStart == 0)
							{
								iStart = icount - windowsize+1;

								if (iStart > iFormerEnd + 1)
								{
									if ((iFormerEnd != 0) && (iFormerStart !=0))
									{
										fprintf(pfRegion, "%lld %lld\n", iFormerStart, iFormerEnd);
										totalsize += iFormerEnd - iFormerStart + 1;
									}

									iFormerStart = iStart;
								}

								iEnd = 0;
							}
						}
						else
						{
							if ((iStart != 0) && (iEnd == 0))
							{
								int i = windowsize -1;

								while (occurrence[i] == 0)
									i--;
		
								iEnd = icount + size - 2 - (windowsize-1 - i);
								iFormerEnd = iEnd;
								iStart = 0;
							}
						}
					}
				}
			}
			iOffset += m_nChunksize - size +1;
		}
		while ((iOffset < m_nGenomesize - size +1) && (iRead == 2)); 

		fclose( pfGenome);
	}

	fprintf(pfRegion, "%lld %lld\n", iFormerStart, iFormerEnd);

	totalsize += iFormerEnd - iFormerStart + 1;

	fprintf(pfRegion, "%lld", totalsize);

	fflush(pfRegion);
	fclose(pfRegion);
	fclose(pfOut);
	free(pchPattern);
	free(pchSequence);
	free(piReadCount);
        free(pchReverse);

	return(1);
}

// for oligo count sets which oligos and its reverse complement are not calculated as the same
// firstly read the oligo sets and put those into one temporary file "handleformats.txt", and get number1 and number2
static int readallrepeat(char *pchrepeatfile, int& number1, int& number2, const int& size, const int& m_nStep1,
const int& m_nEndthreshold, const int& m_nCopy)
{
	FILE *pfRepeatfile;
	char *pchLine;
	char *pchTemp, *pchReverse;
	unsigned long index;

	int threshold;

	if (m_nStep1 < m_nEndthreshold)
		threshold = m_nStep1;
	else threshold = m_nEndthreshold;

	int number = 0;
	
	number1 = 0;
	
	if( ( pfRepeatfile = fopen( pchrepeatfile, "r" )) == NULL )
	{
		printf("Can not find the repeat file: %s.\n", pchrepeatfile);
		return(0);
	}
	else 
	{
		pchLine = (char *)malloc(sizeof(char) * MAXLINECHAR);
                pchTemp = (char*) malloc(sizeof(char) * (size+1));
                pchReverse = (char*) malloc(sizeof(char) * (size+1));
                
		while (!feof(pfRepeatfile))
		{
			if( fgets( pchLine, MAXLINECHAR, pfRepeatfile ) != NULL)
			{
				strncpy( pchTemp, pchLine, size); 
				*(pchTemp + size) = '\0';

				if ( (isonessr(pchTemp) == 0) && (istwossr(pchTemp) == 0) && (isthreessr(pchTemp) == 0) && (isfourssr(pchTemp) == 0) )
				{
					int occur = getnumber(pchLine);

					if (occur >= m_nCopy)
					{
						if (occur >= threshold)
							number1++;
						
						number++;
					}
				}
			}
		}

		free(pchLine);
		free(pchTemp);
                free(pchReverse);
		fclose( pfRepeatfile);
	}

	number2 = number - number1;

	return(1);
}

static int outputrepeat(const cloud_type* Repeat, const char *pchResult, const int& number, const int& size)
{
	FILE *pfResult = fopen(pchResult, "wb");
	
	char *pchPattern = (char*) malloc(sizeof(char) * (size+1));
        for (int i=0; i< number; i++)
	{
		indextopattern(pchPattern, Repeat[i].index, size);
		fprintf(pfResult, "%s %d\n", pchPattern, Repeat[i].number);
	}

	fclose(pfResult);

	return(1);
}

static int readalloligos(char *pchrepeatfile, cloud_type3* repeats1, cloud_type1* repeats2, int& number1, int& number2, const
int& size, const int& m_nStep1, const int& m_nEndthreshold, const int& m_nCopy)
{
	FILE *pfRepeatfile;
	char *pchLine;
	char *pchTemp, *pchReverse;
	int index;

	number1 = 0;
	number2 = 0;
	
	int threshold1 = m_nStep1;

	int threshold2 = m_nEndthreshold;

	int threshold;

	if (threshold1 < threshold2)
		threshold = threshold1;
	else threshold = threshold2;

	if( ( pfRepeatfile = fopen( pchrepeatfile, "r" )) == NULL )
	{
		printf("Can not find the repeat file: %s.\n", pchrepeatfile);
		return(0);
	}
	else 
	{
		pchLine = (char *)malloc(sizeof(char) * MAXLINECHAR);
                pchTemp = (char*) malloc(sizeof(char) * (size+1));
                while (!feof(pfRepeatfile))
		{
			if( fgets( pchLine, MAXLINECHAR, pfRepeatfile ) != NULL)
			{
				strncpy( pchTemp, pchLine, size); 
				*(pchTemp + size) = '\0';

				if ( (isonessr(pchTemp) == 0) && (istwossr(pchTemp) == 0) && (isthreessr(pchTemp) == 0) && (isfourssr(pchTemp) == 0) )
				{
					int occur = getnumber(pchLine);
				
					if (occur >= threshold)
					{
						patterntoindex(pchTemp, repeats1[number1].index, size);
						repeats1[number1].number = occur;
						number1++;
					}
					else if (occur >= m_nCopy)
					{
						patterntoindex(pchTemp, repeats2[number2].index, size);
						repeats2[number2].number = occur;
						number2++;
					}
				}
			}
		}

		free(pchLine);
		free(pchTemp);
		fclose( pfRepeatfile);
	}

	return(1);
}

static int buildpcloud(cloud_type3* repeats1, cloud_type1* repeats2, char* pchOutput, const int& number1, const int& number2,
const int& size, const int& m_nStep3, const int& m_nStep2, const int& m_nStep1, const int& m_nEndthreshold)
{
	//set the pointer to the most often observed repeat;
	int count = 0;
	
	char* pchCore = (char*) malloc(sizeof(char) * (size+1));
        char* pchRepeat = (char*) malloc(sizeof(char) * (size+1));
        int *totalnumber = (int *) malloc(sizeof(int) * MAXCLOUD);
        char** core = (char**) malloc(sizeof(char) * MAXCLOUD * (size+1));
        int *totalsize = (int*) malloc(sizeof(int) * MAXCLOUD);
        int pcloudnumber = 0;

	int distance;	

	unsigned long *piRepeatOnesubstitution =  (unsigned long*) malloc(2*sizeof(unsigned long) * ( 3 * size));
        unsigned long *piRepeatTwosubstitution =  (unsigned long*) malloc(2*sizeof(unsigned long)* ( ( 3 * size) + (9*size*(size-1)/2)));
        unsigned long *piRepeatThreesubstitution =  (unsigned long*) malloc(2*sizeof(unsigned long)* (( 3 * size) + (9*size*(size-1)/2)+ (27*size*(size-1)*(size-2)/(3*2))));
        unsigned long *piCoreThreesubstitution =  (unsigned long*) malloc(2*sizeof(unsigned long)* (( 3 * size) + (9*size*(size-1)/2)+ (27*size*(size-1)*(size-2)/(3*2))));
        
	int k;

	int iRepeatOnesubstitution;

	int iRepeatTwosubstitution;
	
	int iRepeatThreesubstitution;

	int iCoreThreesubstitution;

	int result;

	int* repeats3 = (int*) malloc(sizeof(int) * number1);
	
	sort(repeats1, repeats1+number1, highnumber3);

	for (k = 0; k< number1; k++)
	{
		repeats3[k] = repeats1[k].index;
	}

	sort(repeats1, repeats1+number1, lowsequence3);

	sort(repeats2, repeats2+number2, lowsequence1);

	int coreindex = sbsearch3(number1, repeats1, repeats3[count]);

	while ((count< number1) && (repeats1[coreindex].number >= m_nEndthreshold))
	{
		pcloudnumber++;

		// get the core sequence of the pcloud
		core[pcloudnumber -1] = (char*) malloc(sizeof(char) * (size+1));
                indextopattern(core[pcloudnumber -1], repeats3[count], size);
		strcpy(pchCore, core[pcloudnumber -1]);

		if (repeats1[coreindex].number >= m_nStep3)
		{
			iCoreThreesubstitution = 0;
			getonesubstitutions(pchCore, piCoreThreesubstitution, iCoreThreesubstitution);
			gettwosubstitutions(pchCore, piCoreThreesubstitution, iCoreThreesubstitution);
			getthreesubstitutions(pchCore, piCoreThreesubstitution, iCoreThreesubstitution);
		}
		else if (repeats1[coreindex].number >= m_nStep2)
		{
			iCoreThreesubstitution = 0;
			getonesubstitutions(pchCore, piCoreThreesubstitution, iCoreThreesubstitution);
			gettwosubstitutions(pchCore, piCoreThreesubstitution, iCoreThreesubstitution);
		}
		else
		{
			iCoreThreesubstitution = 0;
			getonesubstitutions(pchCore, piCoreThreesubstitution, iCoreThreesubstitution);
		}

/*		iCoreThreesubstitution = 0;

		getonesubstitutions(pchCore, piCoreThreesubstitution, iCoreThreesubstitution);

		gettwosubstitutions(pchCore, piCoreThreesubstitution, iCoreThreesubstitution);
	
		getthreesubstitutions(pchCore, piCoreThreesubstitution, iCoreThreesubstitution);
		*/
		
		repeats1[coreindex].cloud = pcloudnumber;
		repeats1[coreindex].extension = 1;
		totalsize[pcloudnumber -1] = 1;
		totalnumber[pcloudnumber -1] = repeats1[coreindex].number;

		// get all the repeats in this pclouds within 3 subs
		for (int i = 0; i < iCoreThreesubstitution; i++)
		{
			result = sbsearch3(number1, repeats1, piCoreThreesubstitution[i]);

			if (result >=0)
			{
				if (repeats1[result].cloud == 0)
				{
					repeats1[result].cloud = pcloudnumber;
					totalsize[pcloudnumber -1]++;
					totalnumber[pcloudnumber -1] += repeats1[result].number;
				}

				// get another 1-sub, 2-sub and 3-sub with repeat wbove some threshold
				if (repeats1[result].extension == 0)
				{
					repeats1[result].extension = 1;

					if (repeats1[result].number >=  m_nStep3)
					{
						iRepeatThreesubstitution = 0;

						indextopattern(pchRepeat, repeats1[result].index, size);

						getonesubstitutions(pchRepeat, piRepeatThreesubstitution, iRepeatThreesubstitution);

						gettwosubstitutions(pchRepeat, piRepeatThreesubstitution, iRepeatThreesubstitution);
	
						getthreesubstitutions(pchRepeat, piRepeatThreesubstitution, iRepeatThreesubstitution);

						for (int j = 0; j<iRepeatThreesubstitution; j++)
						{
							result = sbsearch3(number1, repeats1, piRepeatThreesubstitution[j]);

							if (result >= 0) 
							{
								if (repeats1[result].cloud == 0)
								{
									repeats1[result].cloud = pcloudnumber;
									totalsize[pcloudnumber -1]++;
									totalnumber[pcloudnumber -1] += repeats1[result].number;
								}
							}
							else
							{
								result = sbsearch1(number2, repeats2, piRepeatThreesubstitution[j]);
								if ((result >= 0) && (repeats2[result].cloud == 0))
								{
									repeats2[result].cloud = pcloudnumber;
									totalsize[pcloudnumber -1]++;
									totalnumber[pcloudnumber -1] += repeats2[result].number;
								}
							}
						}
					}
					else if (repeats1[result].number >= m_nStep2)
					{
						iRepeatTwosubstitution = 0;
	
						indextopattern(pchRepeat, repeats1[result].index, size);

						getonesubstitutions(pchRepeat, piRepeatTwosubstitution, iRepeatTwosubstitution);
	
						gettwosubstitutions(pchRepeat, piRepeatTwosubstitution, iRepeatTwosubstitution);
		
						for (int j = 0; j<iRepeatTwosubstitution; j++)
						{
							result = sbsearch3(number1, repeats1, piRepeatTwosubstitution[j]);
	
							if (result >= 0) 
							{
								if (repeats1[result].cloud == 0)
								{
									repeats1[result].cloud = pcloudnumber;
									totalsize[pcloudnumber -1]++;
									totalnumber[pcloudnumber -1] += repeats1[result].number;
								}
							}
							else
							{
								result = sbsearch1(number2, repeats2, piRepeatTwosubstitution[j]);
								if ((result >= 0) && (repeats2[result].cloud == 0))
								{
									repeats2[result].cloud = pcloudnumber;
									totalsize[pcloudnumber -1]++;
									totalnumber[pcloudnumber -1] += repeats2[result].number;
								}
							}
						}
					}
					else if (repeats1[result].number >= m_nStep1)
					{
						iRepeatOnesubstitution = 0;
	
						indextopattern(pchRepeat, repeats1[result].index, size);
	
						getonesubstitutions(pchRepeat, piRepeatOnesubstitution, iRepeatOnesubstitution);
	
						for (int j = 0; j<iRepeatOnesubstitution; j++)
						{
							result = sbsearch3(number1, repeats1, piRepeatOnesubstitution[j]);
	
							if (result >= 0) 
							{
								if (repeats1[result].cloud == 0)
								{
									repeats1[result].cloud = pcloudnumber;
									totalsize[pcloudnumber -1]++;
									totalnumber[pcloudnumber -1] += repeats1[result].number;
								}
							}
							else
							{
								result = sbsearch1(number2, repeats2, piRepeatOnesubstitution[j]);
								if ((result >= 0) && (repeats2[result].cloud == 0))
								{
									repeats2[result].cloud = pcloudnumber;
									totalsize[pcloudnumber -1]++;
									totalnumber[pcloudnumber -1] += repeats2[result].number;
								}
							}
						}
					}
				}
			}
			else
			{
				result = sbsearch1(number2, repeats2, piCoreThreesubstitution[i]);

				if ((result >= 0) && (repeats2[result].cloud == 0))
				{
					repeats2[result].cloud = pcloudnumber;
					totalsize[pcloudnumber -1]++;
					totalnumber[pcloudnumber -1] += repeats2[result].number;
				}
			}
		}

		// move to next core repeats
		while((repeats1[coreindex].cloud !=0) && (count < number1))
		{
			count++;
			
			if (count <= number1-1)
				coreindex =  sbsearch3(number1, repeats1, repeats3[count]);

			// extension of repeats

			if ((repeats1[coreindex].extension == 0) && (repeats1[coreindex].cloud != 0))
			{
				repeats1[coreindex].extension = 1;

				if (repeats1[coreindex].number >=  m_nStep3)
				{
					iRepeatThreesubstitution = 0;

					indextopattern(pchRepeat, repeats1[coreindex].index, size);

					getonesubstitutions(pchRepeat, piRepeatThreesubstitution, iRepeatThreesubstitution);

					gettwosubstitutions(pchRepeat, piRepeatThreesubstitution, iRepeatThreesubstitution);
	
					getthreesubstitutions(pchRepeat, piRepeatThreesubstitution, iRepeatThreesubstitution);

					for (int j = 0; j<iRepeatThreesubstitution; j++)
					{
						result = sbsearch3(number1, repeats1, piRepeatThreesubstitution[j]);

						if (result >= 0) 
						{
							if (repeats1[result].cloud == 0)
							{
								repeats1[result].cloud = repeats1[coreindex].cloud;
								totalsize[repeats1[coreindex].cloud -1]++;
								totalnumber[repeats1[coreindex].cloud -1] += repeats1[result].number;
							}
						}
						else
						{
							result = sbsearch1(number2, repeats2, piRepeatThreesubstitution[j]);
							if ((result >= 0) && (repeats2[result].cloud == 0))
							{
								repeats2[result].cloud = repeats1[coreindex].cloud;
								totalsize[repeats1[coreindex].cloud -1]++;
								totalnumber[repeats1[coreindex].cloud -1] += repeats2[result].number;
							}
						}
					}
				}
				else if (repeats1[coreindex].number >= m_nStep2)
				{
					iRepeatTwosubstitution = 0;
					indextopattern(pchRepeat, repeats1[coreindex].index, size);
					getonesubstitutions(pchRepeat, piRepeatTwosubstitution, iRepeatTwosubstitution);
					gettwosubstitutions(pchRepeat, piRepeatTwosubstitution, iRepeatTwosubstitution);
		
					for (int j = 0; j<iRepeatTwosubstitution; j++)
					{
						result = sbsearch3(number1, repeats1, piRepeatTwosubstitution[j]);
						if (result >= 0) 
						{
							if (repeats1[result].cloud == 0)
							{
								repeats1[result].cloud =  repeats1[coreindex].cloud;
								totalsize[ repeats1[coreindex].cloud -1]++;
								totalnumber[ repeats1[coreindex].cloud -1] += repeats1[result].number;
							}
						}
						else
						{
							result = sbsearch1(number2, repeats2, piRepeatTwosubstitution[j]);
							if ((result >= 0) && (repeats2[result].cloud == 0))
							{
								repeats2[result].cloud = repeats1[coreindex].cloud;
								totalsize[repeats1[coreindex].cloud -1]++;
								totalnumber[repeats1[coreindex].cloud -1] += repeats2[result].number;
							}
						}
					}
				}
				else if (repeats1[coreindex].number >= m_nStep1)
				{
					iRepeatOnesubstitution = 0;
					indextopattern(pchRepeat, repeats1[coreindex].index, size);
	
					getonesubstitutions(pchRepeat, piRepeatOnesubstitution, iRepeatOnesubstitution);
	
					for (int j = 0; j<iRepeatOnesubstitution; j++)
					{
						result = sbsearch3(number1, repeats1, piRepeatOnesubstitution[j]);
	
						if (result >= 0) 
						{
							if (repeats1[result].cloud == 0)
							{
								repeats1[result].cloud =  repeats1[coreindex].cloud;
								totalsize[ repeats1[coreindex].cloud -1]++;
								totalnumber[ repeats1[coreindex].cloud -1] += repeats1[result].number;
							}
						}
						else
						{
							result = sbsearch1(number2, repeats2, piRepeatOnesubstitution[j]);
							if ((result >= 0) && (repeats2[result].cloud == 0))
							{
								repeats2[result].cloud =  repeats1[coreindex].cloud;
								totalsize[ repeats1[coreindex].cloud -1]++;
								totalnumber[ repeats1[coreindex].cloud -1] += repeats2[result].number;
							}
						}
					}
				}
			}
		}
	}

	free(piRepeatOnesubstitution);
	free(piRepeatTwosubstitution);
	free(piRepeatThreesubstitution);
	free(piCoreThreesubstitution);
	FILE *pfOutput = fopen(pchOutput, "wb");

	for (k = 0; k < pcloudnumber; k++)
		fprintf(pfOutput, "%d\t%d\t%d\t%s\n", k+1, totalsize[k], totalnumber[k], core[k]);

        

	free(repeats3);
	free(pchRepeat);
	free(pchCore);
	fclose(pfOutput);

	for (k = 0; k < pcloudnumber; k++)
        {    
		free(core[k]);
        }

	free(core);
	free(totalnumber);
	free(totalsize);
	return(1);
}

static int outputclouds(cloud_type3* Repeat1, cloud_type1* Repeat2, const char *pchResult, const int& number1, const int&
number2, const int& size, const int& m_nChunksize, const unsigned int& m_nGenomesize)
{
	FILE *pfResult = fopen(pchResult, "wb");
	
	char *pchPattern = (char*) malloc(sizeof(char) * (size+1));
        sort(Repeat1, Repeat1 + number1, highnumber3);

	sort(Repeat2, Repeat2+number2, highnumber1);

	for (int i=0; i< number1; i++)
	{
		if (Repeat1[i].cloud != 0)
		{
			indextopattern(pchPattern, Repeat1[i].index, size);
			fprintf(pfResult, "%s %d\n", pchPattern, Repeat1[i].cloud);
		}
	}

	for (int i=0; i< number2; i++)
	{
		if (Repeat2[i].cloud != 0)
		{
			indextopattern(pchPattern, Repeat2[i].index, size);
			fprintf(pfResult, "%s %d\n", pchPattern, Repeat2[i].cloud);
		}
	}

	fclose(pfResult);
	
	free(pchPattern);
	return(1);
}

static int readallclouds(const char* pchCloudassign, cloud_type2* pClouds, int& number, const int& size)
{
	FILE *pfCloudassign;
	char *pchLine;
	char *pchTemp;

	number = 0;
	
	if( ( pfCloudassign = fopen( pchCloudassign, "r" )) == NULL )
	{
		printf("Can not find the cloud assign file: %s.\n", pchCloudassign);
		return(0);
	}
	else 
	{
		pchLine = (char *)malloc(sizeof(char) * MAXLINECHAR);
                pchTemp = (char*) malloc(sizeof(char) * (size+1));
                while (!feof(pfCloudassign))
		{
			if( fgets( pchLine, MAXLINECHAR, pfCloudassign ) != NULL)
			{
				strncpy( pchTemp, pchLine, size); 
				*(pchTemp + size) = '\0';
				int cloud = getnumber(pchLine);
				patterntoindex(pchTemp, pClouds[number].index, size); 
				pClouds[number].cloud = cloud;
				number++;
			}
		}

		free(pchLine);
		free(pchTemp);
		fclose(pfCloudassign);
	}

	sort(pClouds, pClouds+number, lowsequence2);

	return(1);
}


// annotate the genome with the P cloud assignment
static int genomecloudannotation(const char* pchGenome, char* pchOutfile, const cloud_type2* pClouds, const int& number,
const int& size, const int& m_nChunksize, const unsigned int& m_nGenomesize)
{
	FILE *pfGenome;

	char *pchPattern;
	int *piReadCount = (int*) malloc(sizeof(int));
        char *pchSequence = (char*) malloc(sizeof(char) * m_nChunksize);
        pchPattern = (char *) malloc(sizeof(char) * (size+1));
        char *pchReverse = (char *) malloc(sizeof(char) * (size+1));
        long long iOffset = 0;
	int iRead;

	FILE *pfOut;
	
	pfOut = fopen(pchOutfile, "wb");

	int patternnumber;
	unsigned long index;

	if( ( pfGenome = fopen( pchGenome, "r" )) == NULL )
	{
		printf("Can not find the genome file: %s.\n", pchGenome);
		return(0);
	}
	else 
	{
		do 
		{	
			if ( !pchSequence ) 
			{
				printf("The system can not allocate the memory!\n");
				free(pchPattern);
				free(pchReverse);
				free(piReadCount);
                                return(0);
			}
		
			if (iOffset + m_nChunksize -1 <= m_nGenomesize)
				iRead = readfromfile(pfGenome, pchSequence, iOffset, m_nChunksize, piReadCount);
			else 
				iRead = readfromfile(pfGenome, pchSequence, iOffset, m_nGenomesize - iOffset, piReadCount);

			if ( iRead == 0 ) 
			{
				free(pchPattern);
				free(pchSequence);
				free(pchReverse);
				free(piReadCount);
				return(0);
			}
			else
			{
				for (int count =0; count <= *piReadCount - size; count ++)
				{
					if ( (count % 50 == 0) && (count>0))
						fprintf(pfOut, "\n");

					getsubstring(pchSequence, pchPattern, count, count + size - 1);
	
					pchPattern[size] = '\0';
					patternnumber = 0;

					if (issegmentvalid(pchPattern)==1)
					{
						// search the indextwble to find whether the left substring is a significant pattern
						patterntoindex(pchPattern, index, size);

						int result = sbsearch2(number, pClouds, index);
						
						if (result <0)
						{
							getreversecomplement(pchPattern, pchReverse);
							patterntoindex(pchReverse, index, size);
							result = sbsearch2(number, pClouds, index);
						}

						// if it is a significant pattern of length (size)
						if (result >= 0 )
						{
							patternnumber = pClouds[result].cloud;
						}
					}
					fprintf(pfOut, "%d ", patternnumber);
				}
			}
			iOffset += m_nChunksize - size +1;
		}
		while ((iOffset < m_nGenomesize - size +1) && (iRead == 2)); 

		fclose( pfGenome);
	}

	fclose(pfOut);
	free(pchPattern);
        free(pchSequence);
        free(pchReverse);
        free(piReadCount);

	return(1);
}

// get the repeat regions based on the p cloud annotation
static int getrepeatregion(char* pchOccurencefile, char* pchRegionfile, const int& windowsize, const int& percent, const int& size)
{
	FILE *pfOccurence;

	FILE *pfRegion;

	vector<int> occurrence;

	int number;
	long long count = 0;
	long long iStart = 0, iEnd = 0;
	long long iFormerStart = 0, iFormerEnd = 0;
	long long totalsize = 0;

	if( ( pfOccurence = fopen( pchOccurencefile, "r" )) == NULL )
	{
		printf("Can not find the occurrence file: %s.\n", pchOccurencefile);
		return(0);
	}
	else 
	{
		pfRegion = fopen(pchRegionfile, "wb");

		while (!feof(pfOccurence))
		{
			fscanf(pfOccurence, "%d", &number);

			count++;

			if (count <= windowsize)
			{
				occurrence.push_back(number);

				if (count == windowsize)
				{
					if (isrepeatregion(occurrence, percent))
					{
						iStart = count - windowsize +1;
						iFormerStart = iStart;
						iEnd = 0;
					}
				}
			}
			else
			{
				occurrence.erase(occurrence.begin());
				occurrence.push_back(number);

				if (isrepeatregion(occurrence, percent))
				{
					if (iStart == 0)
					{
						iStart = count - windowsize+1;

						if (iStart > iFormerEnd + 1)
						{
							if ((iFormerEnd != 0) && (iFormerStart !=0))
							{
								fprintf(pfRegion, "%lld %lld\n", iFormerStart, iFormerEnd);
								totalsize += iFormerEnd - iFormerStart + 1;
							}

							iFormerStart = iStart;
						}

						iEnd = 0;
					}
				}
				else
				{
					if ((iStart != 0) && (iEnd == 0))
					{
						int i = windowsize -1;

						while (occurrence[i] == 0)
							i--;

						iEnd = count + size - 2 - (windowsize-1 - i);

						iFormerEnd = iEnd;

						iStart = 0;
					}
				}
			}
		}

		fprintf(pfRegion, "%lld %lld\n", iFormerStart, iFormerEnd);
	
		totalsize += iFormerEnd - iFormerStart + 1;

		fprintf(pfRegion, "%lld", totalsize);
		fclose(pfOccurence);
		fflush(pfRegion);
		fclose(pfRegion);
	}

	return(1);
}

int pcloudsdissection(const char* pchControlfile)
{
	int size, windowsize, percent, m_nGetclouds, m_nDissection;
	
	int m_nStep1, m_nStep2, m_nStep3, m_nEndthreshold, m_nCopy, m_nChunksize;
	unsigned int m_nGenomesize;

	char pchrepeatfile[MAXFILENAMELENGTH], pchMainClouds[MAXFILENAMELENGTH], pchMainAssign[MAXFILENAMELENGTH], pchAccAssign[MAXFILENAMELENGTH];
	char pchGenome[MAXFILENAMELENGTH], pchAnnotationfile[MAXFILENAMELENGTH], pchRegionfile[MAXFILENAMELENGTH];

	int number1, number2;
				
	// read the control file the get the parameters of the program
	if (!ReadPcloudsControlfile(pchControlfile, size, m_nCopy, m_nEndthreshold, m_nStep1, m_nStep2, m_nStep3, m_nChunksize, m_nGenomesize,      
	windowsize, percent, m_nGetclouds, m_nDissection, pchrepeatfile, pchGenome, pchMainClouds, pchMainAssign, pchAccAssign, pchAnnotationfile, pchRegionfile))
	{
		cerr << "Open Control file error." <<endl;
		return(-1);
	}
	
	if (size == 16)
	{
		if (m_nGetclouds)
		{
			if (readrepeatnumber(pchrepeatfile, number1, number2, size, m_nStep1, m_nEndthreshold, m_nCopy))
			{
				cloud_type3* repeats1 = (cloud_type3*) malloc(sizeof(cloud_type3) * number1);
                                readmainoligos(pchrepeatfile, repeats1, number1, size, m_nStep1, m_nEndthreshold);

				buildmainpcloud(repeats1, pchMainClouds, number1, size, m_nEndthreshold);
				outputmainclouds(repeats1, pchMainAssign, number1, size);

				free(repeats1);
                                cloud_type2* pMainCloudsA = (cloud_type2*) malloc(sizeof(cloud_type2) * MAXMAINMER);
				cloud_type2* pMainCloudsB = (cloud_type2*) malloc(sizeof(cloud_type2) * MAXMAINMER);
				cloud_type2* pMainCloudsC = (cloud_type2*) malloc(sizeof(cloud_type2) * MAXMAINMER);
                                int numberA, numberB, numberC;
	
				readmainclouds(pchMainAssign, pMainCloudsA, numberA, pMainCloudsB, numberB, pMainCloudsC,
				numberC, size, m_nStep1, m_nStep2, m_nStep3);
				buildaccessarypcloud(pchrepeatfile, pchAccAssign, pMainCloudsA, numberA, pMainCloudsB,
				numberB, pMainCloudsC, numberC, number2, size, m_nStep1, m_nEndthreshold, m_nCopy);
	
				free(pMainCloudsA);
                                free(pMainCloudsB);
                                free(pMainCloudsC);
			}       
			else
			{
				cerr << "Open word count file error." <<endl;
				return(-1);
			}
		}


		if (m_nDissection)
		{
			bit_matrix PatternMatrix((unsigned long)pow(4, size/2));

			int*  matrixsize = (int*) malloc(sizeof(int) * (unsigned long) pow(4, size/2) ); 
                        buildbitmatrix(PatternMatrix, matrixsize, size);
	
			readclouds(pchMainAssign, pchAccAssign, PatternMatrix, matrixsize, size);

			if (!GenomeScanAndIdentify(pchGenome, pchAnnotationfile, pchRegionfile, PatternMatrix, matrixsize,
			size, windowsize, percent, m_nChunksize, m_nGenomesize))
			{
				cerr << "Genome Annotation and Repeat Region Identification Error." <<endl;
				free(matrixsize);
                                return(-1);
			}

			free(matrixsize);
                        
		}
	}
	else
	{
		if (m_nGetclouds)
		{
			if (readrepeatnumber(pchrepeatfile, number1, number2, size, m_nStep1, m_nEndthreshold, m_nCopy))
			{
				cloud_type3* repeats1 = (cloud_type3*) malloc(sizeof(cloud_type3) * number1);
                                readmainoligos(pchrepeatfile, repeats1, number1, size, m_nStep1, m_nEndthreshold);

				buildmainpcloud(repeats1, pchMainClouds, number1, size, m_nEndthreshold);
				outputmainclouds(repeats1, pchMainAssign, number1, size);

				free(repeats1);
                                        
				cloud_type2* pMainCloudsA = (cloud_type2*) malloc(sizeof(cloud_type2) * MAXMAINMER);
				cloud_type2* pMainCloudsB = (cloud_type2*) malloc(sizeof(cloud_type2) * MAXMAINMER);
				cloud_type2* pMainCloudsC = (cloud_type2*) malloc(sizeof(cloud_type2) * MAXMAINMER);
                                int numberA, numberB, numberC;
	
				readmainclouds(pchMainAssign, pMainCloudsA, numberA, pMainCloudsB, numberB, pMainCloudsC,
				numberC, size, m_nStep1, m_nStep2, m_nStep3);
				buildaccessarypcloud(pchrepeatfile, pchAccAssign, pMainCloudsA, numberA, pMainCloudsB,
				numberB, pMainCloudsC, numberC, number2, size, m_nStep1, m_nEndthreshold, m_nCopy);
	
				free(pMainCloudsA);
                                free(pMainCloudsB);
                                free(pMainCloudsC);
                                
			}
			else
			{
				cerr << "Open word count file error." <<endl;
				return(-1);
			}
		}


		if (m_nDissection)
		{
			bitvector PatternVector((unsigned long)pow(4, size), false);

			readclouds1(pchMainAssign, pchAccAssign, PatternVector, size);

			if (!GenomeScanAndIdentify1(pchGenome, pchAnnotationfile, pchRegionfile, PatternVector, 
			size, windowsize, percent, m_nChunksize, m_nGenomesize))
			{
				cerr << "Genome Annotation and Repeat Region Identification Error." <<endl;
				PatternVector.clear();
				return(-1);
			}

			PatternVector.clear();
		}

	}
	
	return(1);
}
