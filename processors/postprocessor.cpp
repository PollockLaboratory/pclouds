//
//  Post-processor for P-Clouds 
//
//  Created by A.P. Jason de Koning on 7/5/12.
//  Copyright (c) 2012 University of Colorado Denver, School of Medicine. All rights reserved.

// This program takes an input FASTA formatted file and calculates the start position in a concatenated version
// using the same spacer length as in the pre-processor.  The output file insures annotations don't go off the
// end of a sequence (i.e., into the spacer sequence).

// This version supersedes previous versions.

// This program is meant for use with P-Clouds (Gu et al., 2008; de Koning et al., 2011)

// jason.de.koning@gmail.com
// http://jasondk.org

// Does not correct for the off by 18 error in the regions file
// caused by the preproccessor inserting ">Processed reads" before the genome.

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <map>
#include <utility>

#include <algorithm> 
#include <functional> 
#include <locale>

#include <time.h>
#include <string.h>

// Delete any non-canonical nucleotide characters?
const int totalDeleteNonCanonical = 0;
const int spacerLength = 30;

using namespace std;

// primary functions
void getSequenceNamesAndLengths(string fname, vector<string>& headers,
		vector<long long>& sizes);
void ReadRegionsFile(string regions_file,
		vector<pair<long long int, long long int> >& regions);

// utility functions
static inline std::string &chomp(std::string &s);
void findAndReplace(std::string& source, const char* find, const char* replace);
void cleanseNonCanonical(string &in);
bool _predicate(char c);

int main(int argc, const char * argv[]) {
	vector<pair<long long int, long long int> > regions;

	// WARNING: This version reads the contigs into memory one at a time
	// Check for expected command line args
	if (argc != 4) {
		cerr
				<< "Usage: postprocessor <fastafile> <regionfile> <output annotation filename>"
				<< endl;
		exit(-1);
	}

	string genome_file = argv[1];
	vector<string> sequence_names;
	vector<long long int> sequence_lengths;

	// Read the sequence_names and sequence_lengths of the contigs
	getSequenceNamesAndLengths(genome_file, sequence_names, sequence_lengths);

	// Assuming 0-based counting for the region annotations
	long long int positionPointer = 0;

	vector<long long int> starting_positions(sequence_names.size());
	for (int i = 0; i < sequence_names.size(); i++) {
		starting_positions.at(i) = positionPointer;

		positionPointer += sequence_lengths.at(i) + spacerLength;
	}

	string genome_regions_file = argv[2];
	// Now read in the regions file
	ReadRegionsFile(genome_regions_file, regions);

	string contig_regions_file = argv[3];

	// Open an output file
	ofstream contig_regions_out(contig_regions_file.c_str());

	if (!contig_regions_out.good()) {
		cerr << "Error.  Could not initialize output file." << endl;
		exit(-1);
	}

	// Determine the contig for each region
	for (int i = 0; i < regions.size(); i++) {
		//STP: Default to the last contig
		int correct_contig = starting_positions.size() - 1;
		for (int contig = 1; contig < starting_positions.size(); contig++) {
			if (regions.at(i).first < starting_positions.at(contig)) {
				correct_contig = contig - 1;
				break;
			}
		}
		contig_regions_out << sequence_names.at(correct_contig) << "\t"
				<< (regions.at(i).first - starting_positions.at(correct_contig))
				<< "\t"
				<< (regions.at(i).second - starting_positions.at(correct_contig))
				<< "\n";
	}

	return 0;
}

void ReadRegionsFile(string regions_file,
		vector<pair<long long int, long long int> >& regions) {

	ifstream regions_in(regions_file.c_str());

	if (not regions_in.good()) {
		cerr << "Could not read regions file" << endl;
		exit(-1);
	}

	while (regions_in.good()) {
		pair<long long int, long long int> region;
		regions_in >> region.first;
		regions_in >> region.second;

		//STP: Skip blank lines (the last line is often blank)
		if (region.first == 0 and region.second == 0) {
			continue;
		}

		if (region.first >= region.second) {
			cerr << "Found a region with beginning equal to or after the end. "
					<< "Start: " << region.first << " End: " << region.second
					<< endl;
			exit(-1);
		}
		regions.push_back(region);
	}
}

void getSequenceNamesAndLengths(string fname, vector<string>& headers,
		vector<long long>& sizes) {
	ifstream in;
	string delim = ">";

	in.open(fname.c_str(), ios::in);

	if (!in.is_open()) {
		cerr << "Could not open file" << endl;
		exit(-1);
	}

	// Clear the buffer
	headers.clear();
	sizes.clear();

	// Temp strings
	string header;
	string seq;

	int openSeq = 0;

	// Dump the file to buffer, line by line
	string line;
	while (in.good()) {
		getline(in, line);
		chomp(line);

		// New sequence?
		if (std::string::npos != line.find(delim)) {
			if (openSeq == 1) {
				if (totalDeleteNonCanonical)
					cleanseNonCanonical(seq);
				headers.push_back(header);
				sizes.push_back(seq.size());

				header.clear();
				seq.clear();
				openSeq = 0;
			}

			// Prep the next one
			findAndReplace(line, ">", "");
			header = line;

			openSeq = 1;
		} else {
			seq += line;
		}
	} // while

	// Last one
	if (openSeq == 1) {
		if (totalDeleteNonCanonical)
			cleanseNonCanonical(seq);
		headers.push_back(header);
		sizes.push_back(seq.size());

		header.clear();
		seq.clear();
		openSeq = 0;
	}

	in.close();
}

// Utility functions

// trim whitespace from end
static inline std::string &chomp(std::string &s) {
	s.erase(
			std::find_if(s.rbegin(), s.rend(),
					std::not1(std::ptr_fun<int, int>(std::isspace))).base(),
			s.end());
	return s;
}

// simple find and replace alowing literals
void findAndReplace(std::string& source, const char* find,
		const char* replace) {
	size_t findLen = strlen(find);
	size_t replaceLen = strlen(replace);
	size_t pos = 0;

	//search for the next occurrence of find within source
	while ((pos = source.find(find, pos)) != std::string::npos) {
		source.replace(pos, findLen, replace);
		pos += replaceLen;
	}
}

// return true if character should be removed, false otherwise
bool _predicate(char c) {
	char valid[9] = "ATCGatcg";

	for (int i = 0; i < 9; i++) {
		if (c == valid[i])
			return false;
	}

	return true;
}

void cleanseNonCanonical(string &in) {
	// Remove all non-canonical nucleotide characters
	in.erase(std::remove_if(in.begin(), in.end(), _predicate), in.end());
}
