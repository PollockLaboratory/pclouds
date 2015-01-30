#include <iostream>
#include <fstream> // for printing edges in expansion network
#include <cstring>
#include <algorithm>
#include <cmath>
#include <vector>

#include "../include/macrodefine.h"
#include "../filereader/readfile.h"
#include "../stringhandler/stringhandle.h"

using namespace std;

// Added by STP
#include <deque> // For clouds in regions
bool print_clouds_in_regions = false;
bool genome_has_header = false;

// transform the pattern sequence to the index of the array
// coding method: A=00, C=01, G= 10, T=11,
//STP: This is the only function shared between build_clouds.cpp and
// annotate_genome.cpp
void kmer_sequence_to_number_pattern(const char *kmer_sequence,
		unsigned long& number_pattern, const int& kmer_size) {
	number_pattern = 0;
	for (int i = 0; i < kmer_size; i++) {
		number_pattern *= 4;
		if (kmer_sequence[i] == 'a' or kmer_sequence[i] == 'A')
			number_pattern += 0;
		else if (kmer_sequence[i] == 'c' or kmer_sequence[i] == 'C')
			number_pattern += 1;
		else if (kmer_sequence[i] == 'g' or kmer_sequence[i] == 'G')
			number_pattern += 2;
		else if (kmer_sequence[i] == 't' or kmer_sequence[i] == 'T')
			number_pattern += 3;
	}
}

void update_index(unsigned long& index, char front, char back, int kmer_size) {
	// subtract 4^(k-1) * value of front character
	if (front == 'a' or front == 'A')
		index -= 0 * pow(4, kmer_size - 1);
	else if (front == 'c' or front == 'C')
		index -= 1 * pow(4, kmer_size - 1);
	else if (front == 'g' or front == 'G')
		index -= 2 * pow(4, kmer_size - 1);
	else if (front == 't' or front == 'T')
		index -= 3 * pow(4, kmer_size - 1);

	// Multiply by 4 (this is like a left shift)
	index *= 4;

	// Subtract out 4^0 * value of back char
	if (back == 'a' or back == 'A')
		index += 0;
	else if (back == 'c' or back == 'C')
		index += 1;
	else if (back == 'g' or back == 'G')
		index += 2;
	else if (back == 't' or back == 'T')
		index += 3;
}

bool is_repeat_region(vector<int>& occurrence, const int& percent) {
	vector<int>::iterator pOccurrence;

	int count = 0;

	for (pOccurrence = occurrence.begin(); pOccurrence < occurrence.end();
			pOccurrence++) {
		if (*(pOccurrence) > 0)
			count++;
	}

	if ((int) ((float) count / (float) (occurrence.size()) * 100) >= percent)
		return (1);
	else
		return (0);
}

void UpdatePatternAndCloudIdVectors(string kmer_assign_file,
		vector<bool>& PatternVector, vector<int>& CloudIdVector,
		const int& size) {

	ifstream kmers(kmer_assign_file.c_str());
	if (not kmers.good()) {
		cerr << "Can not read the kmer assign file: " << kmer_assign_file
				<< "\n";
		exit(-1);
	}

	unsigned long index = 0;
	while (kmers.good()) {
		string kmer = "";
		int assigned_cloud = 0;

		kmers >> kmer >> assigned_cloud;
		kmer_sequence_to_number_pattern(kmer.c_str(), index, size);

		PatternVector[index] = true;
		CloudIdVector[index] = assigned_cloud;

		string reverse_complement = getreversecomplement(kmer);
		kmer_sequence_to_number_pattern(reverse_complement.c_str(), index,
				size);

		PatternVector[index] = true;
		CloudIdVector[index] = assigned_cloud;
	}
}

/*
 * STP: I was talking with Aaron and he mentioned using a rolling hash for
 * the index rather than recalculating the entire index when only the end and
 * beginning actually change. I want to write a function update_index(index,
 * genome_chunk[site], genome_chunk[site+kmer_size(-1?)]) that removes the value
 * of the character at genome_chunk[site] from the index, multiplies the index
 * by 4, and then adds the value of the character at
 * genome_chunk[site+kmer_size(-1?)]. I don't want to think about if I need the
 * -1 there or not.
 *
 *
 * This will reduce the annotation time by 1/K because right now calculating
 * the index scales linearly with K while updating the index will be constant
 * in K.
 *
 * I was wrong. It will not reduce the annotation time by 1/K. Annotating
 * includes more than just calculating the index. Therefore the speed-up will
 * be much slower than 1/K.
 *
 * Also, we don't need to copy the kmer from the genome_chunk, instead we are
 * simply looking at locations along the genome_chunk.
 */

void find_repeat_regions(string genome_file, string region_file,
		const vector<bool>& PatternVector, const vector<int>& CloudIdVector,
		const int kmer_size, const int windowsize, const int percent,
		int& genome_chunk_size) {

	char *kmer_sequence = (char *) malloc(sizeof(char) * (kmer_size + 1));
	kmer_sequence[kmer_size] = '\0';
	char *reverse_complement_sequence = (char *) malloc(
			sizeof(char) * (kmer_size + 1));

	char *genome_chunk = (char*) malloc(sizeof(char) * genome_chunk_size);
	if (genome_chunk == NULL) {
		cerr
				<< "The system can not allocate the memory for the genome chunk!\n";
		exit(-1);
	}

	/*
	 * STP: The preprocessed fasta file has ">Processed reads" in the first line
	 * This removes that header so the regions actually correspond to
	 * the regions in the genome if the genome_has_header.
	 */
	long long genome_chunk_start = genome_has_header ? 18 : 0;

	bool patternnumber = 0;
	int cloud_id = 0;
	unsigned long index = 0;

	ofstream regions(region_file.c_str());

	vector<int> occurrence(windowsize, 0);

	//STP: Added for annotating which clouds are in which repeat regions
	// Since queue has no 'erase' function, I'll use a deque which is the
	// implementation anyway
	deque<int> cloud_ids_in_region;

	long long genome_position = 0;
	long long start = 0;
	long long former_start = 0;
	long long former_end = -1; // This must be initialized to < 0 in case the

	bool previous_window_was_repetitive = false;
	bool have_annotated_first_region = false;

	if (print_clouds_in_regions) {
		regions << "Start\tEnd\tClouds\n";
	} else {
		regions << "Start\tEnd\n";
	}
	FILE *genome = fopen(genome_file.c_str(), "r");
	if (genome == NULL) {
		cerr << "Can not find the genome file: " << genome_file << "\n";
		exit(-1);
	}

	while (not feof(genome)) {

		/* Does this handle the stretch of sequence at the junction
		 between reads? Those need to be handled specially.
		 They are handled specially below when the new iOffset is
		 calculated. */

		read_chunk_from_genome_file(genome, genome_chunk, genome_chunk_start,
				genome_chunk_size);

		for (int site = 0; site <= genome_chunk_size - kmer_size;
				site++, genome_position++) {

			//STP: We don't need to copy the sequence here.
			// We could simply have a pointer to the chunk
			strncpy(kmer_sequence, genome_chunk + site, kmer_size);
			/*
			 * STP: This is where we use the rolling function rather than
			 * recalculating the index every time. Simply update it.
			 */

			if (site == 0)
				kmer_sequence_to_number_pattern(kmer_sequence, index,
						kmer_size);
			else
				update_index(index, genome_chunk[site - 1],
						genome_chunk[site + kmer_size - 1], kmer_size);

			//if (issegmentvalid(kmer_sequence)) {
			if (is_segment_valid(genome_chunk + site, kmer_size)
					and PatternVector[index]) {
				patternnumber = true;
				cloud_id = CloudIdVector[index];
			} else {
				patternnumber = false;
				cloud_id = 0;
			}

			occurrence.erase(occurrence.begin());
			occurrence.push_back(patternnumber);

			cloud_ids_in_region.push_back(cloud_id);

			if (genome_position < (windowsize - 1)) {
				// for example, when genome_position is 9, we have
				// looked at 10 kmers and filled a window of 10
				continue;
				// Don't even consider if it is repetitive or not
			}

			if (is_repeat_region(occurrence, percent)) {
				if (not previous_window_was_repetitive) {

					start = genome_position - (windowsize - 1);

					if (start > former_end) {
						// If the next repeat starts past the end
						// of the former repeat
						// If not, merge.
						// I'm not sure this if statement works with
						// the legacy regions.

						if (have_annotated_first_region) {
							// Since we print the previous region,
							// we need to have seen the first region
							// before we print
							regions << former_start << "\t" << former_end;

							if (print_clouds_in_regions) {
								int number_of_clouds_to_print = former_end
										- kmer_size - former_start + 1;
								regions << "\t";
								for (int cloud = 0;
										cloud < number_of_clouds_to_print;
										cloud++) {
									if (cloud == 0) {
										regions << cloud_ids_in_region.front();
									}
									else {
										regions << " " << cloud_ids_in_region.front();
									}

									cloud_ids_in_region.pop_front();
								}
							}
							regions << "\n";

							// Erase all clouds up to a window size before the end
							cloud_ids_in_region.erase(
									cloud_ids_in_region.begin(),
									cloud_ids_in_region.end() - windowsize);
						} else {
							have_annotated_first_region = true;
						}

						former_start = start;
					}
				}
				previous_window_was_repetitive = true;
			} else {
				// Current window is not repetitive
				if (previous_window_was_repetitive) {
					/*
					 * The new way to calculate the former end is
					 * more true to the method published in the
					 * pclouds paper by wanjun.
					 */
					former_end = genome_position - 1 + kmer_size; // ends are exclusive
				} else {
					if (not have_annotated_first_region) {
						cloud_ids_in_region.pop_front();
					}
				}
				previous_window_was_repetitive = false;
			}
		}

		// Subtract kmer size in order to handle the break between chunks
		// properly
		genome_chunk_start += genome_chunk_size - kmer_size + 1;
	}

	fclose(genome);

	// If the last window was repetitive, update the former end
	// This essentially makes the former end to be the genome_file size
	if (previous_window_was_repetitive) {
		/*
		 * The new way to calculate the former end is
		 * more true to the method published in the
		 * pclouds paper by wanjun.
		 */
		former_end = genome_position - 1 + kmer_size; // ends are exclusive
	}

	if (have_annotated_first_region) {
		regions << former_start << "\t" << former_end;
		if (print_clouds_in_regions) {
			int number_of_clouds_to_print = former_end - kmer_size
					- former_start + 1;
			regions << "\t";

			for (int cloud = 0; cloud < number_of_clouds_to_print; cloud++) {
				if (cloud == 0) {
					regions << cloud_ids_in_region.front();
				}
				else {
					regions << " " << cloud_ids_in_region.front();
				}
				cloud_ids_in_region.pop_front();
			}
		}
		regions << "\n";
	}


	free(kmer_sequence);
	free(reverse_complement_sequence);
	free(genome_chunk);
}

void annotate_genome(string controlfile) {
	int kmer_size, window_size, percent;
	bool build_clouds, annotate_genome;

	int outer_threshold, core_threshold_1, core_threshold_2, core_threshold_3,
			core_threshold_4;
	int chunk_size;

	string kmer_counts_file, clouds_summary_file, core_kmers_assign_file,
			outer_kmers_assign_file;
	string genome_file, annotation_file, region_file;

	read_controlfile(controlfile, kmer_size, outer_threshold, core_threshold_1,
			core_threshold_2, core_threshold_3, core_threshold_4, chunk_size,
			window_size, percent, build_clouds, annotate_genome,
			kmer_counts_file, genome_file, clouds_summary_file,
			core_kmers_assign_file, outer_kmers_assign_file, region_file);

	if (annotate_genome) {
		cout << "Annotating the genome" << endl;
		vector<bool> pattern_vector(pow(4, kmer_size), false);
		vector<int> cloud_id_vector(pow(4, kmer_size), 0);

		UpdatePatternAndCloudIdVectors(core_kmers_assign_file, pattern_vector,
				cloud_id_vector, kmer_size);
		UpdatePatternAndCloudIdVectors(outer_kmers_assign_file, pattern_vector,
				cloud_id_vector, kmer_size);

		find_repeat_regions(genome_file, region_file, pattern_vector,
				cloud_id_vector, kmer_size, window_size, percent, chunk_size);
	}
}
