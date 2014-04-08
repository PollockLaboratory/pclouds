

// build main clouds based on the algorithm
//STP: Which algorithm??
static void buildmainpcloud(cloud_type3* core_kmers, char* pchOutput,
		const int& number_of_core_kmers, const int& kmer_size,
		const int& core_threshold) {

	// Sort to descending order
	sort(core_kmers, core_kmers + number_of_core_kmers, highnumber3);

	unsigned long* repeats3 = (unsigned long*) malloc(
			sizeof(unsigned long) * number_of_core_kmers);

	for (int k = 0; k < number_of_core_kmers; k++)
		repeats3[k] = core_kmers[k].number_pattern;

	//STP: How are we sorting these kmers?
	sort(core_kmers, core_kmers + number_of_core_kmers, lowsequence3);

	//STP: What does this do?
	int count = 0;
	// This finds the index for the core oligo with the index main_oligos_index[count]
	// This maps from index in main_oligos_index to index in main_oligos
	int core_kmer_index = sbsearch3(number_of_core_kmers, core_kmers,
			repeats3[count]);

	char* seed_sequence = (char*) malloc(sizeof(char) * (kmer_size + 1));
	char* kmer_sequence = (char*) malloc(sizeof(char) * (kmer_size + 1));
	int *total_count_of_members_for_each_cloud = (int *) malloc(
			sizeof(int) * MAXCLOUD);
	char** core_kmer_sequences = (char**) malloc(
			sizeof(char) * MAXCLOUD * (kmer_size + 1));
	int *number_of_members_for_each_cloud = (int*) malloc(
			sizeof(int) * MAXCLOUD);
	int cloud_number_id = 0;

	int iRepeatThreesubstitution = 0;
	int iCoreThreesubstitution = 0;
	int result = 0;

	unsigned long *piRepeatThreesubstitution = (unsigned long*) malloc(
			2 * sizeof(unsigned long)
					* ((3 * kmer_size) + (9 * kmer_size * (kmer_size - 1) / 2)
							+ (27 * kmer_size * (kmer_size - 1)
									* (kmer_size - 2) / 2)));
	unsigned long *piCoreThreesubstitution = (unsigned long*) malloc(
			2 * sizeof(unsigned long)
					* ((3 * kmer_size) + (9 * kmer_size * (kmer_size - 1) / 2)
							+ (27 * kmer_size * (kmer_size - 1)
									* (kmer_size - 2) / 2)));

	//STP: For expansion network
	ofstream edges_out("edges_out");
	edges_out << "Source_core\tDestination_core" << endl;

	while ((count < number_of_core_kmers)
			&& (core_kmers[core_kmer_index].number >= core_threshold)) {
		// Cloud assignments must begin at 1 NOT 0
		cloud_number_id++;
		cout << "Cloud " << cloud_number_id << endl;
		number_pattern_to_kmer_sequence(kmer_sequence,
				core_kmers[core_kmer_index].number_pattern, kmer_size);

		cout << "Seed is " << kmer_sequence << endl;

		// get the core sequence of the pcloud
		core_kmer_sequences[cloud_number_id - 1] = new char[kmer_size + 1];

		number_pattern_to_kmer_sequence(
				core_kmer_sequences[cloud_number_id - 1],
				repeats3[count], kmer_size);

		cout << "Expanding around seed "
				<< core_kmer_sequences[cloud_number_id - 1] << endl;

		strcpy(seed_sequence, core_kmer_sequences[cloud_number_id - 1]);

		core_kmers[core_kmer_index].cloud = cloud_number_id;
		number_of_members_for_each_cloud[cloud_number_id - 1] = 1;
		total_count_of_members_for_each_cloud[cloud_number_id - 1] =
				core_kmers[core_kmer_index].number;

		core_kmers[core_kmer_index].extension = 1;

		//STP: iCoreThreesubstitution is the number of patterns that are 3 subs
		// away from the seed kmer.
		iCoreThreesubstitution = 0;
		// Determine all the patterns within 3 substitutions
		// Put them in piCoreThreesubstitution
		getonesubstitutions(seed_sequence, piCoreThreesubstitution,
				iCoreThreesubstitution);
		gettwosubstitutions(seed_sequence, piCoreThreesubstitution,
				iCoreThreesubstitution);

		getthreesubstitutions(seed_sequence, piCoreThreesubstitution,
				iCoreThreesubstitution);
//		exit(1);
		//STP: For outputting all the edges between cores for a network
		// representation of the building method.

		// For all the repeats in this cloud within 3 subs
		for (int i = 0; i < iCoreThreesubstitution; i++) {
			//STP: result is now the index in main_oligos for the pattern given
			// by the 'index' piCoreThreesubstitution[i] if found or -1 if not
			// found
			result = sbsearch3(number_of_core_kmers, core_kmers,
					piCoreThreesubstitution[i]);

			//If you find the pattern above in main_oligos
			if (result >= 0) {

				//If the main_oligo is not assigned, assign it
				// The unassigned value of .cloud must be 0
				// Cloud assignments must begin at 1 NOT 0
				if (core_kmers[result].cloud == 0) {
					core_kmers[result].cloud = cloud_number_id;
					number_pattern_to_kmer_sequence(kmer_sequence,
							core_kmers[result].number_pattern, kmer_size);
//					cout << "Adding " << kmer_sequence << endl;

					//STP: Added for expansion network
					edges_out << seed_sequence << "\t" << kmer_sequence << endl;

					number_of_members_for_each_cloud[cloud_number_id - 1]++;
					total_count_of_members_for_each_cloud[cloud_number_id - 1] +=
							core_kmers[result].number;
				}

				//NOTICE: how assignment and extension (expansion) are
				// completely independent
				// get another 1-sub, 2-sub and 3-sub with repeat above core
				// threshold
				// STP: I think 'extension' means has already been extended.
				if (core_kmers[result].extension == 0) {
					// STP: These lines are the same as above!
					core_kmers[result].extension = 1;

					// Find all the oligos that are 3 or less distance
					iRepeatThreesubstitution = 0;
					number_pattern_to_kmer_sequence(kmer_sequence,
							core_kmers[result].number_pattern, kmer_size);

					cout << "Expanding around core " << kmer_sequence << endl;

					getonesubstitutions(kmer_sequence,
							piRepeatThreesubstitution,
							iRepeatThreesubstitution);

					gettwosubstitutions(kmer_sequence,
							piRepeatThreesubstitution,
							iRepeatThreesubstitution);

					getthreesubstitutions(kmer_sequence,
							piRepeatThreesubstitution,
							iRepeatThreesubstitution);

					//STP: For expansion network
					string core_sequence(kmer_sequence);

					// Do the second round of expansion
					for (int j = 0; j < iRepeatThreesubstitution; j++) {
						result = sbsearch3(number_of_core_kmers, core_kmers,
								piRepeatThreesubstitution[j]);

						if (result >= 0) {
							if (core_kmers[result].cloud == 0) {
								core_kmers[result].cloud = cloud_number_id;
								number_pattern_to_kmer_sequence(kmer_sequence,
										core_kmers[result].number_pattern,
										kmer_size);
//								cout << "Adding " << kmer_sequence << endl;

								//STP: For making the expansion network
								edges_out << core_sequence << "\t"
										<< kmer_sequence << endl;

								//edges_out << distanceBetween(core_sequence, kmer_sequence);

								number_of_members_for_each_cloud[cloud_number_id
										- 1]++;
								total_count_of_members_for_each_cloud[cloud_number_id
										- 1] += core_kmers[result].number;
							}
						}
					}
					// This does not work correctly because pchRepeat is
					// overwritten.
					cout << "Done expanding around core " << core_sequence
							<< endl;
				}
			}
		}

		cout << "Done expanding around seed "
				<< core_kmer_sequences[cloud_number_id - 1] << endl;

		// move to next core repeats
		// Keep counting up until you find a main oligo that has not been
		// assigned
		// If this line is changed to checking if the main oligo has not been
		// /extended/ then the method would continue to expand the cloud after
		// the second expansion. Notice how we have not expanded around the
		// core kmers that
		while ((core_kmers[core_kmer_index].cloud != 0
		//STP: Added to make every cloud undergo more than 2 rounds of
		// expansion
				and core_kmers[core_kmer_index].extension != 0)
		//STP:  ^^Now the clouds expand more than twice
				&& (count < number_of_core_kmers)) {
			count++;

			// This is doing extra work. It should be outside this while loop.
			if (count <= number_of_core_kmers - 1)
				core_kmer_index = sbsearch3(number_of_core_kmers, core_kmers,
						repeats3[count]);
		}

	}
//STP: This number is now wrong and so is the assignment above
	cout << "total clouds formed is " << cloud_number_id << endl;
	free(piRepeatThreesubstitution);
	free(piCoreThreesubstitution);
	//output the main clouds information
	FILE *pfOutput = fopen(pchOutput, "wb");

	for (int k = 0; k < cloud_number_id; k++)
		fprintf(pfOutput, "%d\t%d\t%d\t%s\n", k + 1,
				number_of_members_for_each_cloud[k],
				total_count_of_members_for_each_cloud[k],
				core_kmer_sequences[k]);

	free(repeats3);
	free(kmer_sequence);
	free(seed_sequence);
	fclose(pfOutput);

	for (int k = 0; k < cloud_number_id; k++) {
		delete[] (core_kmer_sequences[k]);
	}

	free(core_kmer_sequences);
	free(total_count_of_members_for_each_cloud);
	free(number_of_members_for_each_cloud);

}
