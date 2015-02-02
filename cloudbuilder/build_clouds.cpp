#include <iostream>
#include <fstream> // for printing edges in expansion network
#include <cstring>
#include <algorithm>
#include <vector>

#include "./build_clouds.hpp"
#include "../macrodefine.h"
#include "../filereader/readfile.h"
#include "../stringhandler/stringhandle.h"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::cerr;
using std::sort;
using std::ofstream;
using std::ifstream;

namespace pclouds {

bool keep_SSRs = false;
bool expand_recursively = false;
bool dont_care_about_clouds = false;
bool print_testmers = false;

ofstream testmers("testmers");
ofstream edges_out("edges_out");

template <typename t>
int search(int number_of_kmers, const std::vector <t> ts, const unsigned long key) {
	int m;
	int site = 0;

	while (number_of_kmers >= 1) {
		m = number_of_kmers / 2;
		if (ts[site + m].number_pattern == key)
			return (site + m);
		else if (ts[site + m].number_pattern > key)
			number_of_kmers = m;
		else {
			site = site + m + 1;
			number_of_kmers = number_of_kmers - m - 1;
		}
	}

	return (-1);
}

//STP: Is this a binary search algorithm?
void number_pattern_to_kmer_sequence(char *pchPattern, const unsigned long& index, const int& patternsize) {
	int i;
	unsigned long value;
	unsigned long temp = index;

	for (i = patternsize - 1; i >= 0; i--) {
		value = temp % 4;
		if (value == 0)
			*(pchPattern + i) = 'A';
		else if (value == 1)
			*(pchPattern + i) = 'C';
		else if (value == 2)
			*(pchPattern + i) = 'G';
		else if (value == 3)
			*(pchPattern + i) = 'T';
		temp = temp / 4;
	}

	// This needs to be null terminated
	*(pchPattern + patternsize) = '\0';
}

// Generates all the possible kmers (and their reverse complements) that are
// one substitution away from core_kmer
// This fails to get the reverse complement
void get_one_substitutions(const char* core_kmer, unsigned long* array_of_testmer_number_patterns, int& number_of_testmers) {
	char alphabet[5] = "ACGT";

	int size = strlen(core_kmer);

	char* testmer = new char[size + 1];

	for (int i = 0; i < size; i++) {
		for (int k = 0; k < 4; k++) {
			if (alphabet[k] != core_kmer[i]) {
				strcpy(testmer, core_kmer);
				testmer[i] = alphabet[k];

				//STP: Added for reverse prototype
				if (print_testmers) {
					testmers << testmer << '\n';
				}

				kmer_sequence_to_number_pattern(testmer,
						array_of_testmer_number_patterns[number_of_testmers],
						size);
				number_of_testmers++;
			}
		}
	}

	char* reverse_core_kmer = new char[size + 1];

	getreversecomplement(core_kmer, reverse_core_kmer);

	if (strcmp(reverse_core_kmer, core_kmer) != 0) {
		for (int i = 0; i < static_cast<int>(strlen(reverse_core_kmer)); i++) {
			for (int k = 0; k < 4; k++) {
				if (alphabet[k] != reverse_core_kmer[i]) {
					strcpy(testmer, reverse_core_kmer);
					testmer[i] = alphabet[k];

					//STP: Added for reverse prototype
					if (print_testmers) {
						testmers << testmer << '\n';
					}

					kmer_sequence_to_number_pattern(testmer,
							array_of_testmer_number_patterns[number_of_testmers],
							size);
					number_of_testmers++;
				}
			}
		}
	}

	delete[] (reverse_core_kmer);
	delete[] (testmer);
}

void get_two_substitutions(const char* core_kmer, unsigned long* array_of_testmer_number_patterns, int& number_of_testmers) {
	char alphabet[5] = "ACGT";

	int size = strlen(core_kmer);

	char* testmer = (char*) malloc(sizeof(char) * (size + 1));

	for (int i = 0; i < size - 1; i++) {
		for (int j = i + 1; j < size; j++) {
			for (int k = 0; k < 4; k++) {
				for (int l = 0; l < 4; l++) {
					if ((alphabet[k] != core_kmer[i])
							&& (alphabet[l] != core_kmer[j])) {
						strcpy(testmer, core_kmer);
						testmer[i] = alphabet[k];
						testmer[j] = alphabet[l];

						//STP: Added for reverse prototype
						if (print_testmers) {
							testmers << testmer << '\n';
						}
						kmer_sequence_to_number_pattern(testmer,
								array_of_testmer_number_patterns[number_of_testmers],
								size);
						number_of_testmers++;
					}
				}
			}
		}
	}

	char* reverse_core_kmer = (char*) malloc(sizeof(char) * (size + 1));

	getreversecomplement(core_kmer, reverse_core_kmer);

	for (int i = 0; i < static_cast<int>(strlen(reverse_core_kmer)) - 1; i++) {
		for (int j = i + 1; j < static_cast<int>(strlen(reverse_core_kmer)); j++) {
			for (int k = 0; k < 4; k++) {
				for (int l = 0; l < 4; l++) {
					if ((alphabet[k] != reverse_core_kmer[i])
							&& (alphabet[l] != reverse_core_kmer[j])) {
						strcpy(testmer, reverse_core_kmer);
						testmer[i] = alphabet[k];
						testmer[j] = alphabet[l];

						//STP: Added for reverse prototype
						if (print_testmers) {
							testmers << testmer << '\n';
						}

						kmer_sequence_to_number_pattern(testmer,
								array_of_testmer_number_patterns[number_of_testmers],
								size);
						number_of_testmers++;
					}
				}
			}
		}
	}

	free(reverse_core_kmer);
	free(testmer);
}

void get_three_substitutions(char* core_kmer, unsigned long* array_of_testmer_number_patterns, int& number_of_testmers) {
	char alphabet[5] = "ACGT";

	int size = strlen(core_kmer);

	char* testmer = (char*) malloc(sizeof(char) * (size + 1));

	for (int i = 0; i < size - 2; i++) {
		for (int j = i + 1; j < size - 1; j++) {
			for (int m = j + 1; m < size; m++) {
				for (int k = 0; k < 4; k++) {
					for (int l = 0; l < 4; l++) {
						for (int n = 0; n < 4; n++) {
							if ((alphabet[k] != core_kmer[i])
									&& (alphabet[l] != core_kmer[j])
									&& ((alphabet[n] != core_kmer[m]))) {
								strcpy(testmer, core_kmer);
								testmer[i] = alphabet[k];
								testmer[j] = alphabet[l];
								testmer[m] = alphabet[n];

								//STP: Added for reverse prototype
								if (print_testmers) {
									testmers << testmer << '\n';
								}

								kmer_sequence_to_number_pattern(testmer,
										array_of_testmer_number_patterns[number_of_testmers],
										size);
								number_of_testmers++;
							}
						}
					}
				}
			}
		}
	}

	char* reverse_core_kmer = (char*) malloc(sizeof(char) * (size + 1));

	getreversecomplement(core_kmer, reverse_core_kmer);

	for (int i = 0; i < size - 2; i++) {
		for (int j = i + 1; j < static_cast<int>(strlen(reverse_core_kmer)) - 1; j++) {
			for (int m = j + 1; m < static_cast<int>(strlen(reverse_core_kmer)); m++) {
				for (int k = 0; k < 4; k++) {
					for (int l = 0; l < 4; l++) {
						for (int n = 0; n < 4; n++) {
							if ((alphabet[k] != reverse_core_kmer[i])
									&& (alphabet[l] != reverse_core_kmer[j])
									&& ((alphabet[n] != reverse_core_kmer[m]))) {
								strcpy(testmer, reverse_core_kmer);
								testmer[i] = alphabet[k];
								testmer[j] = alphabet[l];
								testmer[m] = alphabet[n];

								//STP: Added for reverse prototype
								if (print_testmers) {
									testmers << testmer << '\n';
								}

								kmer_sequence_to_number_pattern(testmer,
										array_of_testmer_number_patterns[number_of_testmers],
										size);
								number_of_testmers++;
							}
						}
					}
				}
			}
		}
	}

	free(reverse_core_kmer);
	free(testmer);
}

void count_core_and_outer_kmers(string kmer_counts_file, int& number_of_core_kmers, int& number_of_outer_kmers, const int core_threshold, const int outer_threshold) {
	number_of_core_kmers = 0;
	number_of_outer_kmers = 0;

	ifstream kmer_counts(kmer_counts_file.c_str());
	if (not kmer_counts.good()) {
		cerr << "Can not find the kmer counts file: " << __LINE__
			<< kmer_counts_file << "\n";
		//STP: This is a fatal error
		exit(-1);
	}

	string kmer = "";
	int kmer_count = 0;
	while (kmer_counts.good()) {
		kmer_counts >> kmer >> kmer_count;
		if (not is_SSR(kmer.c_str()) or keep_SSRs) {
			if (kmer_count >= core_threshold)
				number_of_core_kmers++;
			else if (kmer_count >= outer_threshold)
				number_of_outer_kmers++;
		}
	}
	cout << "Read " << number_of_core_kmers << " cores and "
		<< number_of_outer_kmers << " outer kmers\n";
	if (number_of_core_kmers == 0) {
		cerr << "Did not find any kmers above the core threshold.\n";
		cerr << "Cannot continue. Set core threshold lower.\n";
		exit(-1);
	}
}

// read the oligos which is in mainclouds into repeats1 array
//STP: This does not take into account reverse complement of kmers.
void read_core_kmers(string kmer_counts_file, vector <CoreKmer> core_kmers, int& number_of_core_kmers, const int kmer_size, const int core_threshold) {
	number_of_core_kmers = 0;
	ifstream kmer_counts(kmer_counts_file.c_str());

	if (not kmer_counts.good()) {
		cerr << "Can not find the repeat file: " << __LINE__ << kmer_counts_file
			<< "\n";
		//STP: This is a fatal error
		exit(-1);
	}

	while (kmer_counts.good()) {
		string kmer; // These are in here to limit scope.
		int kmer_count;
		kmer_counts >> kmer >> kmer_count;
		if (not is_SSR(kmer.c_str()) or keep_SSRs) {
			if (kmer_count >= core_threshold) {
				kmer_sequence_to_number_pattern(kmer.c_str(),
						core_kmers[number_of_core_kmers].number_pattern,
						kmer_size);
				core_kmers[number_of_core_kmers].count = kmer_count;
				number_of_core_kmers++;
			}
		}
	}
}

// build main clouds based on the algorithm
//STP: Which algorithm??
void build_cloud_cores(vector <CoreKmer> core_kmers, string cloud_summary_file, const int number_of_core_kmers, const int& kmer_size, const int& core_threshold) {
	// Sort to descending order of kmer count
	sort(core_kmers.begin(), core_kmers.begin() + number_of_core_kmers, greater_core_count);

	unsigned long* core_kmer_number_patterns = (unsigned long*) malloc(
			sizeof(unsigned long) * number_of_core_kmers);

	for (int k = 0; k < number_of_core_kmers; k++)
		core_kmer_number_patterns[k] = core_kmers[k].number_pattern;

	// Sort by number pattern (used to be called 'index')
	sort(core_kmers.begin(), core_kmers.begin() + number_of_core_kmers, lesser_core_pattern);

	int core_kmer_index_from_top = 0;
	// This finds the index for the core oligo with the index main_oligos_index[count]
	// This maps from index in main_oligos_index to index in main_oligos
	int core_kmer_index = search(number_of_core_kmers, core_kmers,
			core_kmer_number_patterns[core_kmer_index_from_top]);

	char* seed_sequence = (char*) malloc(sizeof(char) * (kmer_size + 1));
	char* kmer_sequence = (char*) malloc(sizeof(char) * (kmer_size + 1));

	int *total_count_of_members_for_each_cloud = (int *) malloc(
			sizeof(int) * MAX_CLOUDS);
	char** seed_sequences = (char**) malloc(
			sizeof(char) * MAX_CLOUDS * (kmer_size + 1));
	int *number_of_members_for_each_cloud = (int*) malloc(
			sizeof(int) * MAX_CLOUDS);
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
	edges_out << "Source_core\tDestination_core\n";

	int total_number_of_clouds = 0;
	if (dont_care_about_clouds) {
		cout << "basically ignoring clouds\n";
		total_number_of_clouds = number_of_core_kmers;
		for (int i = 0; i < number_of_core_kmers; i++) {
			core_kmers[i].cloud = i;
		}
	} else {
		while ((core_kmer_index_from_top < number_of_core_kmers)
				&& (core_kmers[core_kmer_index].count >= core_threshold)) {
			// Cloud assignments must begin at 1 NOT 0
			// If the core kmer has not been assigned
			if (core_kmers[core_kmer_index].cloud == 0) {
				total_number_of_clouds++;
				core_kmers[core_kmer_index].cloud = total_number_of_clouds;
			}
			cloud_number_id = core_kmers[core_kmer_index].cloud;

			//			cout << "Cloud " << cloud_number_id << "\n";
			number_pattern_to_kmer_sequence(kmer_sequence,
					core_kmers[core_kmer_index].number_pattern, kmer_size);

			//			cout << "Seed is " << kmer_sequence << "\n";

			// get the core sequence of the pcloud
			seed_sequences[cloud_number_id - 1] = new char[kmer_size + 1];

			number_pattern_to_kmer_sequence(seed_sequences[cloud_number_id - 1],
					core_kmer_number_patterns[core_kmer_index_from_top],
					kmer_size);

			//cout << "Expanding around seed "
			//	<< seed_sequences[cloud_number_id - 1] << "\n";

			strcpy(seed_sequence, seed_sequences[cloud_number_id - 1]);

			number_of_members_for_each_cloud[cloud_number_id - 1] = 1;
			total_count_of_members_for_each_cloud[cloud_number_id - 1] =
				core_kmers[core_kmer_index].count;

			core_kmers[core_kmer_index].has_been_extended = true;

			//STP: iCoreThreesubstitution is the number of patterns that are 3 subs
			// away from the seed kmer.
			iCoreThreesubstitution = 0;
			// Determine all the patterns within 3 substitutions
			// Put them in piCoreThreesubstitution
			//			cout <<"stuff"<< "\n";
			testmers << "EXPANDING AROUND " << seed_sequence << "\n";
			char* reverse_core_kmer = new char[kmer_size + 1];

			getreversecomplement(seed_sequence, reverse_core_kmer);
			testmers << "REV COMP is  " << reverse_core_kmer << "\n";

			get_one_substitutions(seed_sequence, piCoreThreesubstitution,
					iCoreThreesubstitution);
			//			exit(1);
			get_two_substitutions(seed_sequence, piCoreThreesubstitution,
					iCoreThreesubstitution);

			get_three_substitutions(seed_sequence, piCoreThreesubstitution,
					iCoreThreesubstitution);
			//					exit(1);
			//STP: For outputting all the edges between cores for a network
			// representation of the building method.

			// For all the testmers within 3 subs of the seedmer
			for (int i = 0; i < iCoreThreesubstitution; i++) {
				//STP: result is now the index in main_oligos for the pattern given
				// by the 'index' piCoreThreesubstitution[i] if found or -1 if not
				// found
				result = search(number_of_core_kmers, core_kmers,
						piCoreThreesubstitution[i]);

				//If you find the pattern above in main_oligos
				if (result >= 0) {

					//If the main_oligo is not assigned, assign it.
					// The unassigned value of .cloud must be 0
					// Cloud assignments must begin at 1 NOT 0
					if (core_kmers[result].cloud == 0) {
						core_kmers[result].cloud = cloud_number_id;
						number_pattern_to_kmer_sequence(kmer_sequence,
								core_kmers[result].number_pattern, kmer_size);
						//cout << "Adding " << kmer_sequence << "\n";

						//STP: Added for expansion network
						edges_out << seed_sequence << "\t" << kmer_sequence
							<< "\n";

						number_of_members_for_each_cloud[cloud_number_id - 1]++;
						total_count_of_members_for_each_cloud[cloud_number_id
							- 1] += core_kmers[result].count;
					}

					//NOTICE: how assignment and extension (expansion) are
					// completely independent

					// get another 1-sub, 2-sub and 3-sub with repeat above core
					// threshold
					// STP: I think 'extension' means has already been extended.
					if (not core_kmers[result].has_been_extended) {
						// STP: These lines are the same as above!
						core_kmers[result].has_been_extended = true;

						// Find all the oligos that have a distance of 3 or less
						iRepeatThreesubstitution = 0;
						number_pattern_to_kmer_sequence(kmer_sequence,
								core_kmers[result].number_pattern, kmer_size);

						//cout << "Expanding around core " << kmer_sequence << "\n";

						get_one_substitutions(kmer_sequence,
								piRepeatThreesubstitution,
								iRepeatThreesubstitution);

						get_two_substitutions(kmer_sequence,
								piRepeatThreesubstitution,
								iRepeatThreesubstitution);

						get_three_substitutions(kmer_sequence,
								piRepeatThreesubstitution,
								iRepeatThreesubstitution);

						//STP: For expansion network
						string core_sequence(kmer_sequence);

						// Do the second round of expansion
						for (int j = 0; j < iRepeatThreesubstitution; j++) {
							result = search(number_of_core_kmers, core_kmers,
									piRepeatThreesubstitution[j]);

							if (result >= 0) {
								if (core_kmers[result].cloud == 0) {
									core_kmers[result].cloud = cloud_number_id;
									number_pattern_to_kmer_sequence(
											kmer_sequence,
											core_kmers[result].number_pattern,
											kmer_size);
									//cout << "Adding " << kmer_sequence << "\n";

									//STP: For making the expansion network
									edges_out << core_sequence << "\t"
										<< kmer_sequence << "\n";

									//edges_out << distanceBetween(core_sequence, kmer_sequence);

									number_of_members_for_each_cloud[cloud_number_id
										- 1]++;
									total_count_of_members_for_each_cloud[cloud_number_id
										- 1] += core_kmers[result].count;
								}
							}
						}
						//cout << "Done expanding around core " << core_sequence
						//	<< "\n";
					}
				}
			}

			//cout << "Done expanding around seed "
			//	<< seed_sequences[cloud_number_id - 1] << "\n";

			// move to next core repeats
			// Keep counting up until you find a main oligo that has not been
			// assigned
			// If this line is changed to checking if the main oligo has not been
			// /extended/ then the method would continue to expand the cloud after
			// the second expansion. Notice how we have not expanded around the
			// core kmers that are 2 steps from a seedmer.
			// Keep counting if
			while (
					// We haven't seen the last core. Stop if we have.
					(core_kmer_index_from_top < number_of_core_kmers) and (
						// EITHER
						// We are not expanding recursively and the core we are
						// looking at has been assigned already. Stop if it has not.
						(not expand_recursively
						 and core_kmers[core_kmer_index].cloud != 0)
						// OR
						or
						// We are expanding recursively and the core we are looking at
						// has been assigned AND extended. Stop if it has not.
						(expand_recursively
						 and core_kmers[core_kmer_index].cloud != 0
						 and core_kmers[core_kmer_index].has_been_extended))) {
				core_kmer_index_from_top++;

				// This is doing extra work. It should be outside this while loop.
				if (core_kmer_index_from_top <= number_of_core_kmers - 1)
					core_kmer_index =
						search(number_of_core_kmers, core_kmers,
								core_kmer_number_patterns[core_kmer_index_from_top]);

			}
		}
		FILE *pfOutput = fopen(cloud_summary_file.c_str(), "wb");

		for (int k = 0; k < total_number_of_clouds; k++)
			fprintf(pfOutput, "%d\t%d\t%d\t%s\n", k + 1,
					number_of_members_for_each_cloud[k],
					total_count_of_members_for_each_cloud[k],
					seed_sequences[k]);
		fclose(pfOutput);

	}
	cout << total_number_of_clouds << " clouds formed\n";

	free(piRepeatThreesubstitution);
	free(piCoreThreesubstitution);
	for (int k = 0; k < total_number_of_clouds; k++) {
		delete[] (seed_sequences[k]);
	}
	free(core_kmer_number_patterns);
	free(kmer_sequence);
	free(seed_sequence);

	free(seed_sequences);
	free(total_count_of_members_for_each_cloud);
	free(number_of_members_for_each_cloud);

}

//output the mainclouds assignments of each oligo in main clouds
void output_cloud_cores(vector <CoreKmer> core_kmers, string core_kmers_assign_file, const int& number_of_core_kmers, const int& size) {
	FILE *core_kmers_assign = fopen(core_kmers_assign_file.c_str(), "wb");

	char *pchPattern = (char*) malloc(sizeof(char) * (size + 1));

	sort(core_kmers.begin(), core_kmers.begin() + number_of_core_kmers, greater_core_count);

	for (int i = 0; i < number_of_core_kmers; i++) {
		number_pattern_to_kmer_sequence(pchPattern,
				core_kmers[i].number_pattern, size);
		if (core_kmers[i].cloud > 9000)
			cerr << "Cloud id is over 9000. This is probably an error.\n"
				"Unless you don't care about clouds." << endl;
		fprintf(core_kmers_assign, "%s %d %d\n", pchPattern, core_kmers[i].cloud,
				core_kmers[i].count);
	}

	fclose(core_kmers_assign);
	free(pchPattern);
}

// read main clouds assignments information into three different sets
void read_cloud_cores(string core_assign_file, vector <CoreKmer> tertiary_cores, int& number_of_tertiary_cores, vector <CoreKmer> secondary_cores, int& number_of_secondary_cores, vector <CoreKmer> primary_cores, int& number_of_primary_cores, const int& size, const int& primary_threshold, const int& secondary_threshold, const int& tertiary_threshold) {
	number_of_tertiary_cores = 0;
	number_of_secondary_cores = 0;
	number_of_primary_cores = 0;

	ifstream cloud_assign(core_assign_file.c_str());

	if (not cloud_assign.good()) {
		cerr << "Can not find the cloud assign file: " << core_assign_file
			<< "\n";
		exit(-1);
	}

	while (cloud_assign.good()) {
		string kmer; // These declarations are here to limit the scope
		int cloud;
		int count;
		unsigned long number_pattern;

		cloud_assign >> kmer >> cloud >> count;

		if (count >= tertiary_threshold) {
			tertiary_cores[number_of_tertiary_cores].cloud = cloud;
			kmer_sequence_to_number_pattern(kmer.c_str(), number_pattern, size);
			tertiary_cores[number_of_tertiary_cores].number_pattern =
				number_pattern;
			number_of_tertiary_cores++;
		} else if (count >= secondary_threshold) {
			secondary_cores[number_of_secondary_cores].cloud = cloud;
			kmer_sequence_to_number_pattern(kmer.c_str(), number_pattern, size);
			secondary_cores[number_of_secondary_cores].number_pattern =
				number_pattern;
			number_of_secondary_cores++;
		} else if (count >= primary_threshold) {
			primary_cores[number_of_primary_cores].cloud = cloud;
			kmer_sequence_to_number_pattern(kmer.c_str(), number_pattern, size);
			primary_cores[number_of_primary_cores].number_pattern =
				number_pattern;
			number_of_primary_cores++;
		}
	}

	sort(tertiary_cores.begin(), tertiary_cores.begin() + number_of_tertiary_cores, lesser_core_pattern);
	sort(secondary_cores.begin(), secondary_cores.begin()+ number_of_secondary_cores, lesser_core_pattern);
	sort(primary_cores.begin(), primary_cores.begin() + number_of_primary_cores, lesser_core_pattern);
}

// assign the oligos in accessory regions into P clouds constructed in former step
void build_cloud_outer(string kmer_counts_file, string outer_kmers_assign_file, vector <CoreKmer> core_kmers_above_tertiary, const int number_of_cores_above_tertiary, vector <CoreKmer> core_kmers_above_secondary, const int number_of_cores_above_secondary, vector <CoreKmer> core_kmers_above_primary, const int number_of_cores_above_primary, const int kmer_size, const int core_threshold, const int outer_threshold) {

	ofstream outer_assign(outer_kmers_assign_file.c_str());

	//STP: For expansion network


	ifstream kmer_counts(kmer_counts_file.c_str());
	if (not kmer_counts.good()) {
		cerr << "Can not find the kmer counts file: " << __LINE__
			<< kmer_counts_file << "\n";
		//STP: This is a fatal error
		exit(-1);
	}

	vector <Kmer> outer_kmers;

	string kmer = "";
	int kmer_count = 0;
	int number_of_outer_kmers = 0;
	while (kmer_counts.good()) {
		kmer_counts >> kmer >> kmer_count;
		if (not is_SSR(kmer.c_str()) or keep_SSRs) {
			if (kmer_count < core_threshold and kmer_count >= outer_threshold) {
				kmer_sequence_to_number_pattern(kmer.c_str(),
						outer_kmers[number_of_outer_kmers].number_pattern,
						kmer_size);

				number_of_outer_kmers++;
			}
		}
	}
	cout << "Found " << number_of_outer_kmers << " outer kmers\n";
	sort(outer_kmers.begin(), outer_kmers.begin() + number_of_outer_kmers, lesser_pattern);

	//STP: It would be nice to sort the core kmers before expanding
	//STP: This is not possible because the counts are not included in
	// this type of cloud.
	//			sort(core_kmers_above_tertiary, core_kmers_above_tertiary + number_of_cores_above_tertiary, highnumber2);
	//			sort(core_kmers_above_secondary, core_kmers_above_secondary+ number_of_cores_above_secondary, greater_count);
	//			sort(core_kmers_above_primary, core_kmers_above_primary + number_of_cores_above_primary, greater_count);

	char *core_kmer = (char*) malloc(sizeof(char) * (kmer_size + 1));
	char *testmer = (char*) malloc(sizeof(char) * (kmer_size + 1));

	//STP: The sizes of these are the number of possible kmers with 1, 2, and 3,
	// substitutions
	// sum ( 3^n * size * (3(size - n))! / n!) for n = {1, 2, 3}

	unsigned long *distance_3_tesmers = (unsigned long*) malloc(
			2 * sizeof(unsigned long)
			* ((3 * kmer_size) + (9 * kmer_size * (kmer_size - 1) / 2)
				+ (27 * kmer_size * (kmer_size - 1)
					* (kmer_size - 2) / (3 * 2))));

	// build the accessory p clouds for these repeats chunk
	for (int i = 0; i < number_of_cores_above_tertiary; i++) {
		number_pattern_to_kmer_sequence(core_kmer,
				core_kmers_above_tertiary[i].number_pattern, kmer_size);

		int number_of_distance_3_tesmers = 0;

		get_three_substitutions(core_kmer, distance_3_tesmers,
				number_of_distance_3_tesmers);
		get_two_substitutions(core_kmer, distance_3_tesmers,
				number_of_distance_3_tesmers);
		get_one_substitutions(core_kmer, distance_3_tesmers,
				number_of_distance_3_tesmers);

		for (int count = 0; count < number_of_distance_3_tesmers; count++) {
			int index = search(number_of_outer_kmers, outer_kmers,
					distance_3_tesmers[count]);

			if (index >= 0) {
				if (outer_kmers[index].cloud == 0) {
					outer_kmers[index].cloud =
						core_kmers_above_tertiary[i].cloud;

					//STP: For making the expansion network
					number_pattern_to_kmer_sequence(testmer,
							distance_3_tesmers[count], kmer_size);
					edges_out << core_kmer << "\t" << testmer << "\n";
				}
			}
		}
	}

	// NOTICE: the distance_2_tesmers are not cleared between
	// rounds. This is OK because we keep track of how many repeat
	// two substitutions we have in number_of_distance_2_tesmers.

	unsigned long *distance_2_tesmers =
		(unsigned long*) malloc(
				2 * sizeof(unsigned long)
				* ((3 * kmer_size)
					+ (9 * kmer_size * (kmer_size - 1) / 2)));

	for (int i = 0; i < number_of_cores_above_secondary; i++) {
		number_pattern_to_kmer_sequence(core_kmer,
				core_kmers_above_secondary[i].number_pattern, kmer_size);

		// get the 2-mutation sets of core_kmer;
		int number_of_distance_2_tesmers = 0;

		get_two_substitutions(core_kmer, distance_2_tesmers,
				number_of_distance_2_tesmers);
		get_one_substitutions(core_kmer, distance_2_tesmers,
				number_of_distance_2_tesmers);

		for (int count = 0; count < number_of_distance_2_tesmers; count++) {
			int index = search(number_of_outer_kmers, outer_kmers,
					distance_2_tesmers[count]);

			if (index >= 0) {
				if (outer_kmers[index].cloud == 0) {
					outer_kmers[index].cloud =
						core_kmers_above_secondary[i].cloud;

					//STP: For making the expansion network
					number_pattern_to_kmer_sequence(testmer,
							distance_2_tesmers[count], kmer_size);
					edges_out << core_kmer << "\t" << testmer << "\n";
				}
			}
		}
	}

	unsigned long *distance_1_testmers = (unsigned long*) malloc(
			2 * sizeof(unsigned long) * (3 * kmer_size));

	for (int i = 0; i < number_of_cores_above_primary; i++) {
		number_pattern_to_kmer_sequence(core_kmer,
				core_kmers_above_primary[i].number_pattern, kmer_size);

		// get the 1-mutation sets of core_kmer;
		int number_of_distance_1_tesmers = 0;

		get_one_substitutions(core_kmer, distance_1_testmers,
				number_of_distance_1_tesmers);

		for (int count = 0; count < number_of_distance_1_tesmers; count++) {
			int index = search(number_of_outer_kmers, outer_kmers,
					distance_1_testmers[count]);

			if (index >= 0) {
				if (outer_kmers[index].cloud == 0) {
					outer_kmers[index].cloud =
						core_kmers_above_primary[i].cloud;

					//STP: For making the expansion network
					number_pattern_to_kmer_sequence(testmer,
							distance_1_testmers[count], kmer_size);
					edges_out << core_kmer << "\t" << testmer << "\n";
				}
			}
		}
	}


	free(distance_1_testmers);
	free(distance_2_tesmers);
	free(distance_3_tesmers);

	free(testmer);
	free(core_kmer);


	char *outer_kmer = (char*) malloc(sizeof(char) * (kmer_size + 1));
	for (int i = 0; i < number_of_outer_kmers; i++) {
		if (outer_kmers[i].cloud != 0) {
			number_pattern_to_kmer_sequence(outer_kmer,
					outer_kmers[i].number_pattern, kmer_size);
			outer_assign << outer_kmer << " " << outer_kmers[i].cloud << "\n";
		}
	}
	free(outer_kmer);
}

void build_clouds(string controlfile) {
	int kmer_size, window_size, percent;
	bool build_clouds, annotate_genome;

	int outer_threshold, core_threshold_1, core_threshold_2, core_threshold_3,
	    core_threshold_4;
	int chunk_size;

	string kmer_counts_file, clouds_summary_file, core_kmers_assign_file,
	       outer_kmers_assign_file;
	string genome_file, annotation_file, region_file;

	int number_of_core_kmers, number_of_outer_kmers;

	// Could use an "Options" object to store all these options rather than
	// giving a variable for each. This function definition is ugly. It has
	// way too many arguments.
	read_controlfile(controlfile, kmer_size, outer_threshold, core_threshold_1,
			core_threshold_2, core_threshold_3, core_threshold_4, chunk_size,
			window_size, percent, build_clouds, annotate_genome,
			kmer_counts_file, genome_file, clouds_summary_file,
			core_kmers_assign_file, outer_kmers_assign_file,
			region_file);

	if (build_clouds) {
		cout << "Building clouds" << endl;
		count_core_and_outer_kmers(kmer_counts_file, number_of_core_kmers,
				number_of_outer_kmers, core_threshold_1,
				outer_threshold);

		vector <CoreKmer> core_kmers;

		read_core_kmers(kmer_counts_file, core_kmers, number_of_core_kmers,
				kmer_size, core_threshold_1);

		build_cloud_cores(core_kmers, clouds_summary_file, number_of_core_kmers,
				kmer_size, core_threshold_1);

		output_cloud_cores(core_kmers, core_kmers_assign_file,
				number_of_core_kmers, kmer_size);

		vector <CoreKmer> tertiary_cores;
		/*
		 * STP: We could also use vectors. I'm not sure they will be slower but
		 * it would be much more readable.
		 *
		 * vector<Kmer> tertiary_coresv;
		 * tertiary_coresv.reserve(MAX_CORE_KMERS);
		 */

		vector <CoreKmer> secondary_cores;
		vector <CoreKmer> primary_cores;
		int number_of_tertiary_cores, number_of_secondary_cores,
		    number_of_primary_cores;

		read_cloud_cores(core_kmers_assign_file, tertiary_cores,
				number_of_tertiary_cores, secondary_cores,
				number_of_secondary_cores, primary_cores,
				number_of_primary_cores, kmer_size, core_threshold_2,
				core_threshold_3, core_threshold_4);

		build_cloud_outer(kmer_counts_file, outer_kmers_assign_file,
				tertiary_cores, number_of_tertiary_cores, secondary_cores,
				number_of_secondary_cores, primary_cores,
				number_of_primary_cores, kmer_size, core_threshold_1,
				outer_threshold);

	}
}

}
