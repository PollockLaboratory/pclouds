//
//  main.cpp
//  DinucSimulator
//
//  Created by A.P. Jason de Koning on 5/15/12.
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
#include <vector>
#include <sstream>

#include <algorithm> 
#include <functional> 
#include <locale>

#include <time.h>
#include <string.h>

// Turn this off to resample the sequences including ambiguities
const int totalDeleteNonCanonical = 1;

using namespace std;

// primary functions
void getFileContents( string fname, vector<string>& headers, vector<string>& seqs );
void simulateSequence( string& newSequence, string& theWindow );

// utility functions
static inline std::string &chomp(std::string &s);
void findAndReplace( std::string& source, const char* find, const char* replace );
void cleanseNonCanonical( string &in );
bool _predicate(char c);

int main (int argc, const char * argv[]) {
    vector<string> headers;
    vector<string> seqs;
    
    string theWindow;
    long int windowSize;

    ifstream in;

    // Check for expected command line args
    if (argc < 3) {
        cerr << "Error.  Expecting input FASTA filename and window size." << endl;
        return -1;
    }
    
    string fname = argv[1];
    long int winLength = atol( argv[2] );    

    // Check that window length is divisible by 2
    if ( winLength % 2 != 0 ) { cerr << "Error. Window length must be divisible by two." << endl; exit(-1); }
    
    // Open an output file
    ofstream out;
    string outname = fname + ".rand";
    
    out.open( outname.c_str(), ios::out );
    if (!out.is_open()) {
        cerr << "Error.  Could not open output file." << endl;
        exit(-1);
    }

    // Initialize random number gen
    srand( time(NULL) );
    
    // Read the FASTA file
    getFileContents( fname, headers, seqs ); 

    // For each sequence...
    for (int i=0; i<headers.size(); i++) {
        cout << "Processing sequence " << (i+1) << "..." << endl;
        windowSize = winLength;
        
        // Handle case where window is longer than sequence
        if ( seqs.at(i).size() < winLength ) {
            windowSize = seqs.at(i).size();
            cout << "   Using window of length " << windowSize << ", since sequence length is < " << winLength << endl;
        }
        
        cout << "   Sequence length is " << seqs.at(i).size() << "; win size is " << windowSize << endl;
        
        int done=0;
        string newSeq;
        
        // Get input sequence in chunks and randomized dinucs
        for (long int j=0; j < seqs.at(i).size(); j+=windowSize) {

            if ( (j + windowSize) > seqs.at(i).size() ) {
                windowSize = seqs.at(i).size() - j;
                done = 1;
                
                cout << "   Using window of length " << windowSize << ", since we'd otherwise fall off the end..." << endl;
            }
            
            // Get the sequence in our current window
            theWindow = seqs.at(i).substr( j, windowSize );
            cout << "   Window starting at " << j << " of length " << windowSize << endl;

            // Sample new sequence with same dinucleotide frequencies
            simulateSequence( newSeq, theWindow );
        
            // Output
            stringstream tempHead;
            
            tempHead << ">" << headers.at(i) << " dinuc win from " << j << ", len " << windowSize;
            out << tempHead.str() << endl;
            out << newSeq << endl;
            
            if (done == 1) continue;
        }
        
    }
    
    out.close();
    return 0;
}

void simulateSequence( string& newSequence, string& theWindow ) {
    
    newSequence.clear();
    newSequence.resize( theWindow.size() );
    
    // Sample new sequence with replacement having same dinucleotide freqs
    for (long int i=0; i < theWindow.size()-1; i+=2) {
        // Sample dinucs from either phase, at random (with replacement)
        long int u = rand() % (theWindow.size()-1);
        newSequence.at(i  ) = theWindow.at( u   );
        newSequence.at(i+1) = theWindow.at( u+1 );
    }

}

void getFileContents( string fname, vector<string>& headers, vector<string>& seqs ) {
    ifstream in;
    string delim = ">";
    
    in.open( fname.c_str(), ios::in );
    
    if (!in.is_open()) {
        cerr << "Could not open file" << endl;
        exit(-1);
    }

    // Clear the buffer
    headers.clear();
    seqs.clear();
    
    // Temp strings
    string header;
    string seq;
    
    int openSeq = 0;
    
    // Dump the file to buffer, line by line
    string line;
    while (in.good()) {
        getline( in, line );
        chomp( line );
        
        // New sequence?
        if (std::string::npos != line.find( delim )) {
            if ( openSeq == 1 ) {
                if ( totalDeleteNonCanonical ) cleanseNonCanonical( seq );
                headers.push_back( header );
                seqs.push_back( seq );
                
                header.clear(); seq.clear();
                openSeq = 0;
            }
            
            // Prep the next one
            findAndReplace( line, ">", "" );
            header = line;

            openSeq = 1;
        } else {
            seq += line;
        }
    } // while

    // Last one
    if ( openSeq == 1 ) {
        if ( totalDeleteNonCanonical ) cleanseNonCanonical( seq );
        headers.push_back( header );
        seqs.push_back( seq );
        
        header.clear(); seq.clear();
        openSeq = 0;
    }

    in.close();
}


// Utility functions

// trim whitespace from end
static inline std::string &chomp(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}

// simple find and replace alowing literals
void findAndReplace( std::string& source, const char* find, const char* replace ) {
    size_t findLen = strlen(find);
    size_t replaceLen = strlen(replace);
    size_t pos = 0;
    
    //search for the next occurrence of find within source
    while ((pos = source.find( find, pos)) != std::string::npos) {
        source.replace( pos, findLen, replace );
        pos += replaceLen; 
    }
}

// return true if character should be removed, false otherwise
bool _predicate(char c) {
    char valid[9] = "ATCGatcg";
    
    for (int i=0; i < 9; i++) {
        if ( c == valid[i] ) return false;
    }
    
    return true;
}

void cleanseNonCanonical( string &in ) {
    // Remove all non-canonical nucleotide characters
    in.erase(std::remove_if(in.begin(), in.end(), _predicate), in.end());
}
