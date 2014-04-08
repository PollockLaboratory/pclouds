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
void getFileContents( string fname, vector<string>& headers, vector<string>& seqs );
void getRegionsFile( string fname, vector< pair< long long int, long long int > >& regions );

// utility functions
static inline std::string &chomp(std::string &s);
void findAndReplace( std::string& source, const char* find, const char* replace );
void cleanseNonCanonical( string &in );
bool _predicate(char c);
void tokenize( string myString, char separator, vector<string>& tokens, int theSize );

int main (int argc, const char * argv[]) {
    vector< pair< long long int, long long int > > regions;
    
    // WARNING: This version reads the entire FASTA file into memory
    
    vector<string> headers;
    vector<string> seqs;
    
    // Create a hash for the start of each sequence in the concatenated 'processed reads'
    map<string, long long int> sequenceHash;

    // Assuming 0-based counting for the region annotations
    long long int positionPointer = 0;
    
    ifstream in;

    // Check for expected command line args
    if (argc < 4) {
        cerr << "Usage: postprocessor <fastafile> <regionfile> <output annotation filename>" << endl;
        return -1;
    }
    
    string fname = argv[1];
    string regionFile = argv[2];
    string outname = argv[3];

    // Open an output file
    ofstream out;
    
    out.open( outname.c_str(), ios::out );
    if (!out.is_open()) {
        cerr << "Error.  Could not open output file." << endl;
        exit(-1);
    }

    // Read the FASTA file
    getFileContents( fname, headers, seqs ); 

    // For each sequence, build starting positions hash...
    
    for (int i=0; i<headers.size(); i++) {
        sequenceHash[ headers.at(i) ] = positionPointer;
        
        positionPointer += seqs.at(i).size() + spacerLength;
    }

    // Now read in the regions file
    getRegionsFile( regionFile, regions );
    
    // Assume that hits are in order; don't have to keep re-searching (i think)
    int lastHit = 0;
    
    // Now output the correct coords
    for (int i=0; i < regions.size(); i++) {
        
        int closestSeq = 0;
        
        // Find the corresponding sequence
        for (int j=lastHit; j < headers.size(); j++) {
            if ( sequenceHash[ headers.at(j) ] < regions.at(i).first ) {
                closestSeq = j;
            } else {
                lastHit = closestSeq;
                break;
            }
        }

        // Don't allow the end annotation to go into the spacer (this happens due to the annotation algo)
        long long int end = regions.at(i).second;
        if ( closestSeq < (headers.size()-1) ) {
            if ( ( end + spacerLength ) > sequenceHash[ headers.at(closestSeq+1) ] ) {
                end = ( sequenceHash[ headers.at(closestSeq+1) ] - spacerLength) ;
                // [diagnostic] cout << "Annotation falls off the end.  End of sequence should be " << end << ", but is listed as " << regions.at(i).second << endl;
            }
        }
        
        out << headers.at(closestSeq) << "\t" << ( regions.at(i).first - sequenceHash[ headers.at(closestSeq) ] ) << "\t" << ( end - sequenceHash[ headers.at(closestSeq) ] ) << endl; //
    }
    
    out.close();
    return 0;
}

void getRegionsFile( string fname, vector< pair< long long int, long long int > >& regions ) {
    vector<string> tok;
    pair< long long int, long long int > tempRegion;

    ifstream in;
    
    in.open( fname.c_str(), ios::in );
    
    if (!in.is_open()) {
        cerr << "Could not open file" << endl;
        exit(-1);
    }

    // Temp strings
    string header;
    string seq;

    // Fix the number of columns for speed
    tok.resize( 2 );
    
    int openSeq = 0;
    
    stringstream convert;
    long long int first, second;
    
    // Dump the file to buffer, line by line
    string line;
    while (in.good()) {
    
        getline( in, line );
        chomp( line );

        tokenize( line, '\t', tok, 2 );

        convert.clear();
        convert << tok[0];
        convert >> first;

        second = -1;
        convert.clear();
        convert << tok[1];
        convert >> second;

        if ( ( second > first ) && ( second != -1 ) ) {
            tempRegion = make_pair( first, second );
            regions.push_back( tempRegion );
        }        

    } // while
    
    in.close();
   
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

// tokenize a string
void tokenize( string myString, char separator, vector<string>& tokens, int theSize ) {
    // specifying size of tokens vector up front to avoid dynamic re-allocation
    static stringstream sstring;
    string token;

    sstring.clear();
    sstring << myString;
    
    // Clear the tokens vector passed in
    if ( tokens.size() != theSize ) {
        tokens.resize( theSize );
    }
    
    int count = 0;
    while ( getline( sstring, token, separator ) ) {
        tokens.at( count++ ) = token;
        if (count > theSize) break;
    }

}

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
