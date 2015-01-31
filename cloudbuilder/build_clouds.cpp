#include <iostream>
#include <fstream> // for printing edges in expansion network
#include <cstring>
#include <algorithm>
#include <math.h>
#include <vector>
#include <map>

#include "../include/macrodefine.h"
#include "../filereader/readfile.h"
#include "../stringhandler/stringhandle.h"

using namespace std;

// Added by STP
bool keep_SSRs = false;
bool expand_recursively = false;
bool dont_care_about_clouds = false;
bool print_testmers = false;

ofstream testmers("testmers");
ofstream edges_out("edges_out");

map<string, int> kmers;
map<string, string> kmer_assignments;

map<string, int> read_kmers(string kmer_counts_file,
		int outer_threshold) {

	map<string, int> kmers;
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
		if ((not is_SSR(kmer.c_str()) or keep_SSRs)
				and kmer_count >= outer_threshold) {
			kmers[kmer] = kmer_count;
		}
	}
	return kmers;
}

bool compare_kmer_counts(pair<string, int> kmer_1, pair<string, int> kmer_2) {
	return (kmer_1.second < kmer_2.second);
}

vector<string> get_core_kmers(map<string, int> kmers, int core_threshold) {

	// Collect the core kmers and their counts
	vector<pair<string, int> > core_kmers_pairs;
	for (map<string, int>::iterator kmer = kmers.begin(); kmer != kmers.end();
			kmer++) {
		if (kmer->second >= core_threshold) {
			core_kmers_pairs.push_back(make_pair(kmer->first, kmer->second));
		}
	}

	// sort core kmers
	sort(core_kmers_pairs.begin(), core_kmers_pairs.end(), compare_kmer_counts);

	// Extract just the kmers and ignore the counts
	vector<string> core_kmers;
	for (vector<pair<string, int> >::iterator kmer = core_kmers_pairs.begin();
			kmer < core_kmers_pairs.end(); kmer++) {
		core_kmers.push_back(kmer->first);
	}

	return core_kmers;
}

void ExpandCloudAroundKmer(string cloud_id, string core_kmer,
		int core_threshold) {
	char nucs[] = "ACGT";
	for (int position = 0; position < core_kmer.length(); position++) {
		for (int nuc = 0; nuc < 4; nuc++) {
			string testmer = core_kmer;
			testmer.at(position) = nucs[nuc];
			testmers << testmer << "\n";
			// If the testmer was seen in the genome and it has not been
			// assigned yet
			if (kmers.find(testmer) != kmers.end()
					and kmer_assignments.find(testmer)
					== kmer_assignments.end()) {
				// Keep track of edge from core kmer to testmer
				edges_out << core_kmer << "\t" << testmer << "\n";

				kmer_assignments[testmer] = cloud_id;
				if (kmers[testmer] >= core_threshold) {
					ExpandCloudAroundKmer(cloud_id, testmer, core_threshold);
				}
			}
		}
	}
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

	read_controlfile(controlfile, kmer_size, outer_threshold, core_threshold_1,
			core_threshold_2, core_threshold_3, core_threshold_4, chunk_size,
			window_size, percent, build_clouds, annotate_genome,
			kmer_counts_file, genome_file, clouds_summary_file,
			core_kmers_assign_file, outer_kmers_assign_file, region_file);

	if (build_clouds) {
		cout << "Building clouds" << endl;
		kmers = read_kmers(kmer_counts_file, outer_threshold);

		vector<string> core_kmers = get_core_kmers(kmers, core_threshold_1);

		for (vector<string>::iterator core_kmer = core_kmers.begin();
				core_kmer < core_kmers.end(); core_kmer++) {
			// If the kmer has not been assigned.
			if (kmer_assignments.find(*core_kmer) == kmer_assignments.end()) {
				string cloud_id = *core_kmer;

				kmer_assignments[*core_kmer] = cloud_id;
				ExpandCloudAroundKmer(cloud_id, *core_kmer, core_threshold_1);
			}
		}
	}
}
