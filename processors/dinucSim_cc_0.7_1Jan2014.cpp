//
//  main.cpp
//  DinucSimulator
//
//  Created by A.P. Jason de Koning on 5/15/12.
//  Modified by Corey Cox on 6-Dec-2013
//  Copyright (c) 2012 University of Colorado Denver, School of Medicine. All rights reserved.

// This program takes an input FASTA formatted file and window size, and simulates new sequences
// having the same dinucleotide frequencies as the template.  Sampling is with replacement.
// All non-canonical nucleotides are removed prior to simulation.  Windows at the ends of sequences
// will be less than requested window size.

// IF USING MANY SHORT CONTIGS, one could either run this program on the sequences directly, or
// could concatenate these in random order and then run.

// This program is meant for use with P-Clouds (Gu et al., 2008; de Koning et al., 2011)

// jason.de.koning@gmail.com
// http://jasondk.org

#include <iostream>
#include <fstream>
#include <sstream>
#include <array>

#include <algorithm>
#include <functional>
#include <locale>

#include <time.h>
#include <string.h>

// Turn this off to resample the sequences including ambiguities
const int totalDeleteNonCanonical = 1;
using namespace std;

// set up arrays to hold our samples
char nucs[8][100];  // array to hold samples
int s_begin[8]; // counter first posiiton of samples in sample array
int s_end[8]; // counter for last position of samples in sample array

// Utility functions
// Note: this is ASCII dependant, Unicode will break this!!! FASTA is defined as ASCII so shoud be OK.
int getNucIndex (char & nuc_char) { return (nuc_char - 'A') % 16; } // A:0, C:2, G:3, T:6

// trim whitespace from end
static inline std::string &chomp(std::string &s) {
	s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
	return s;
}

bool _all(char c) { char valid_nucs[9] = "ATCGatcg";
	for (int i=0; i < 9; i++) { if ( c == valid_nucs[i] ) return false; }
	return true;
}

bool _masked_only (char c) { char valid_nucs[5] = "atcg";
	for (int i=0; i < 5; i++) { if ( c == valid_nucs[i] ) return false; }
	return true;
}

bool _unmasked_only (char c) { char valid_nucs[5] = "ATCG";
	for (int i=0; i < 5; i++) { if ( c == valid_nucs[i] ) return false; }
	return true;
}

bool (* _predicate)(char) = _all;

void cleanseNonCanonical( string &in ) {
	// Remove all nucleotide characters not in the appropriate set
	in.erase(std::remove_if(in.begin(), in.end(), _predicate), in.end());
}

// primary functions
void setup (int argc, const char * argv[], ifstream& fp_in, long int& winLength, ofstream& fp_out ) {
	// get filename from argv and open input stream
	string fname;
	unsigned int seed = 0;
	string masking;

	if (argc < 2 || argc > 5) {
		cerr << "Usage:\n\tdinucSim filename\nOR\n\t dinucSim filename [windowsize] [masked/unmasked/all] [seed]" << endl; exit(-1);
	}
	// check if window length, seed and masking are in argv and set.
	if (argc > 2) { winLength = atol(argv[2]); }
	if (argc > 3) { masking.assign(argv[3]); }
	if (argc > 4) { seed = (unsigned)atoi(argv[4]); }

	if (winLength < 1000) {
		cout << "Window length < 1000; set to default 1000000" << endl;
		winLength = 1000000;
	}
	if (seed < 1) {
		seed = time(NULL);
		cout << "Seed < 1 setting random seed" << endl;
	}
	if (masking == "masked") { _predicate = _masked_only;
		cout << "Performing simulation on masked portion of the genome only" << endl;
	}
	else if (masking == "unmasked") { _predicate = _unmasked_only;
		cout << "Performing simulation on unmasked portion of the genome only" << endl;
	}
	else { cout << "Performing simulation without accounting for masking" << endl; }

	// setup input and output files
	fname = argv[1];
	fp_in.open( fname.c_str(), ios::in );
	if (!fp_in.is_open()) { cerr << "Could not open file: fname" << endl; exit(-1); }

	string outname = fname + ".rand";
	fp_out.open( outname.c_str(), ios::out );
	if (!fp_out.is_open()) { cerr << "Error.  Could not open output file." << endl; exit(-1); }
	fp_out << ">Sequence simulated from " << fname;

	srand( seed ); // Initialize random number gen
	cout << "Random seed is " << seed << endl;
}

void getGenome( ifstream& fp_in, string& genome, long int& winLength ) {
	string delim = ">";

	// Clear the buffer
	genome.clear();
	string seq;

	// Dump the file to buffer, line by line
	string line;
	while (fp_in.good()) {
		getline( fp_in, line );
		chomp( line );

		// Check for start of new sequence
		if (std::string::npos != line.find( delim )) {
			if ( totalDeleteNonCanonical ) cleanseNonCanonical( seq );
			genome += seq; seq.clear();
		}
		else { seq += line; }
	}
	if ( totalDeleteNonCanonical ) cleanseNonCanonical( seq );
	genome += seq;  // add last sequence to genome
	fp_in.close();

	if ( genome.size() < winLength*2 ) {
		cout << "   Window length " << winLength << "is > 1/2 genome size " << genome.size() << endl;
		cout << "   Aborting " << endl;
		exit(-1);
	}
	cout << "   Genome length is " << genome.size() << "; window size is " << winLength << endl;

	// circularize genome
	genome += genome.substr( 0, winLength/2 );
}

void simulateSequence( string& newSequence, string& theWindow, int& prevnuc ) {
	// Sample new sequence across moving window 1/2 window size
	long int windowSize = theWindow.size()/2;
	newSequence.clear();
	newSequence.resize( windowSize );

	for (long int i=0; i < windowSize; i++) { int count = 0;
		while (s_begin[prevnuc] == s_end[prevnuc]) {
			long int u = rand() % windowSize + i;
			int nuc = getNucIndex(theWindow[u]);
			nucs[nuc][s_end[nuc]] = theWindow[u+1];
			if (s_end[nuc] < 99) { s_end[nuc]++; } else { s_end[nuc] = 0; }
			// handle the unusual case of too few samples in the window - could happen due to isochores
			if ( count > 100 ) { s_begin[prevnuc] = rand() % 99; }
		}
		newSequence.at(i) = nucs[prevnuc][s_begin[prevnuc]];
		if (s_begin[prevnuc] < 99) { s_begin[prevnuc]++; } else { s_begin[prevnuc] = 0; }
		prevnuc = getNucIndex(newSequence.at(i));
	}
}

int main (int argc, const char * argv[]) {
	ifstream fp_in;
	ofstream fp_out;

	string genome;
	string theWindow;
	long int winLength = 0;

	setup (argc, argv, fp_in, winLength, fp_out );

	// Read the FASTA file, concatenate, circularize - keep in memory
	getGenome (fp_in, genome, winLength);

	string newSeq;
	int prevnuc = getNucIndex(genome[0]);  // seed with first nucleotide in the genome

	// Get input sequence in chunks and randomized dinucs
	for (long int j=0; j < genome.size(); j+=winLength) {
		fp_out << newSeq << endl;
		theWindow = genome.substr( j, winLength*2 );
		simulateSequence( newSeq, theWindow, prevnuc );
	}
	// add last nucleotide if genome.size() is odd
	if (genome.size() % 2 > 0) { newSeq += nucs[prevnuc][rand() % 99]; }
	fp_out << newSeq << endl;
	fp_out.close();
	return 0;
}
