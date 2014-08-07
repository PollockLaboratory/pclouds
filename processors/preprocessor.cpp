//
//  Pre-processor for P-Clouds 
//
//  Created by A.P. Jason de Koning on 6/11/12.
//  Copyright (c) 2012 University of Colorado Denver, School of Medicine. All rights reserved.

// This program takes an input FASTA formatted file and concatenates multiple sequences together
// with an 'NNNN..' spacer.  This technique was introduced in the original P-Clouds distribution.

// This program is meant for use with P-Clouds (Gu et al., 2008; de Koning et al., 2011)

// jason.de.koning@gmail.com
// http://jasondk.org

// The off-by-18 error introduced by prepending the genome with ">Processed reads\n"
// has been fixed
// STP 2014-8-5

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#include <algorithm> 
#include <functional> 
#include <locale>

#include <time.h>
#include <string.h>

// Delete any non-canonical nucleotide characters?
const int totalDeleteNonCanonical = 0;

using namespace std;

// primary functions
void getFileContents( string fname, vector<string>& headers, vector<string>& seqs );

// utility functions
static inline std::string &chomp(std::string &s);
void findAndReplace( std::string& source, const char* find, const char* replace );
void cleanseNonCanonical( string &in );
bool _predicate(char c);

int main (int argc, const char * argv[]) {
    vector<string> headers;
    vector<string> seqs;
    
    ifstream in;

    // Check for expected command line args
    if (argc < 3) {
        cerr << "Usage: preprocessor <fastafile> <output filename>" << endl;
        return -1;
    }
    
    string fname = argv[1];
    string outname = argv[2];

    // Open an output file
    ofstream out;
    
    out.open( outname.c_str(), ios::out );
    if (!out.is_open()) {
        cerr << "Error.  Could not open output file." << endl;
        exit(-1);
    }

    // Read the FASTA file
    getFileContents( fname, headers, seqs ); 

    stringstream newSeq;
    // For each sequence...
    for (int i=0; i<headers.size(); i++) {
        
	    if (i>0) {
		newSeq << "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
	    }
	    
 	    newSeq << seqs.at(i);
    }

    cout << "Processed " << headers.size() << " sequences" << endl << endl;

    // Output
    /*STP:
     * tempHead is what causes the off by 18 problem in old p-clouds
     */
    stringstream tempHead;

    tempHead << ">Processed reads";
//STP: tempHead is no longer printed.
//    out << tempHead.str() << endl;
    out << newSeq.str() << endl;
    
    out.close();
    return 0;
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
