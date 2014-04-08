
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