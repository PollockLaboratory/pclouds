//  main.cpp
//  cloudPos
//
//  Created by Corey Cox on 5-Feb-2015
//  Copyright (c) 2015 University of Colorado Denver, School of Medicine. All rights reserved.
//  Last updated on 17-Mar-2015

// This program takes an input FASTA formatted file and a set of P-Clouds and finds the locations of the
// P-Cloud kmers in the fasta sequence. Output is in bed file format.

#include <string>
#include <unordered_map>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
using std::unordered_map;
using std::string;
using std::cerr;
using std::invalid_argument;
using std::out_of_range;
using std::ifstream;
using std::ios;
using std::cout;
using std::endl;

namespace grizzly {
namespace pclouds {

struct kmer_map {
    unordered_map<string,string> all;
    unordered_map<string,string> rev;
    int length;
};

// get command line options
void getOptions (int argc, const char * argv[], unordered_map<string,string>& options) {
    string option; string value;
    for (int i = 1; i < argc; i++) {
        option = argv[i];
        if (option[0] == '-') { i++;
            if (i < argc ) { value = argv[i]; } else { cerr << option << " requires value\n"; exit(0); }
            if (option[1] == 'c') { options.emplace("cloud_file", value); };
            if (option[1] == 'f') { options.emplace("fasta_file", value); };
            if (option[1] == 'm') { options.emplace("merge", value); };
            if (option[1] == 'p') { options.emplace("print", value); };
            if (option[1] == 'r') { options.emplace("rev_comp", value); };
        }
    }
}

int getInt (unordered_map<string,string>& options, string key) {  string value = options.at(key);
    try { return stoi(value); }
    catch (invalid_argument&) { cerr<< "Aborting -- " << key <<" value must be a number.\n"; exit(1); }
    catch (out_of_range&) { cerr<< "Aborting -- " << key <<" value out of integer range.\n"; exit(1); }
}

void openfile (ifstream& fp, string& fname) {
    fp.open( fname.c_str(), ios::in );
    if (!fp.is_open()) { cerr << "Could not open file: " << fname << endl; exit(-1); }
}

char compliment (char nuc) {
    if (nuc == 'A') { return 'T'; }
    if (nuc == 'C') { return 'G'; }
    if (nuc == 'G') { return 'C'; }
    if (nuc == 'T') { return 'A'; }
    return nuc;
}

string reverseCompliment (const string& seq) {
    string revseq = ""; const int len = seq.length();
    for (int i = 1; i <= len; i++) { revseq.push_back(compliment(seq[len-i])); }
    return revseq;
}

// get search kmers from .assign files
void getKmers (ifstream& fp_in, kmer_map& kmers) {
//    int test = 10;
    while (fp_in.good()) {
        string kmer; string cloud; string count;
        fp_in >> kmer >> cloud >> count;
        kmers.all.emplace (kmer, cloud);
    }
    kmers.all.erase("");
    kmers.length = kmers.all.begin()->first.length();
}

// get reverse compliment of kmers, place in search map and separate map.
void inline getRevKmers (kmer_map& kmers) {
    for (unordered_map<string,string>::iterator it = kmers.all.begin(); it != kmers.all.end(); it++) {
        kmers.rev.emplace(reverseCompliment(it->first), it->second);
    }
    kmers.all.insert(kmers.rev.begin(), kmers.rev.end());
}

// find and merge locations matching kmers with strand orientation, then print.
void annotateContig (string& contig, string& seq, kmer_map& kmers, int& merge, bool& cloud_print) {
    const int len = seq.length() - kmers.length;
    int begin = 0; int end = 0; int num_plus = 0; int num_minus = 0;
    bool print = false; string start_cloud; string clouds = "";
    for (int pos = 0; pos < len; pos++) {
        string kmer = seq.substr(pos, kmers.length);
        if (1 == kmers.all.count(kmer)) { print = true;
            if (end + merge >= pos) {
                clouds += "\t" + kmers.all.at(kmer);
                end = pos + kmers.length;
                if (kmers.rev.count(kmer)) { num_minus++; } else { num_plus++; }
            }
            else {
                begin = pos; end = pos + kmers.length;
                start_cloud = kmers.all.at(kmer);
                clouds = kmers.all.at(kmer);
                if (kmers.rev.count(kmer)) { num_minus++; } else { num_plus++; }
            }
        }
        else {
            if (print && end + merge <= pos) {
                cout << contig <<"\t"<< begin <<"\t"<< end << "\tpcloud:" << start_cloud << "\t0\t";
                if (num_plus >= num_minus) { cout << "+"; } else { cout << "-"; }
                if (cloud_print) { cout <<"\t"<< clouds; }
                cout << endl;
                print = false; num_plus = 0; num_minus = 0;
            }
        }
    }
    if (print) {
        cout << contig <<"\t"<< begin <<"\t"<< end << "\tpcloud:" << start_cloud << "\t0\t";
        if (num_plus >= num_minus) { cout << "+"; } else { cout << "-"; }
        if (cloud_print) { cout <<"\t"<< clouds; }
        cout << endl;
    }
}

// traverse fasta file by contig and get annotation with annotateContig.
void annotateFasta (ifstream& fp_in, kmer_map& kmers, int merge, bool cloud_print) {
    string contig = ""; int pos = 0; string line = ""; string seq = "";
    while (fp_in.good()) {
        fp_in >> line;
        if ( '>' == line[0]) {
            annotateContig(contig, seq, kmers, merge, cloud_print);
            contig = line.substr(1);
            pos = 0; seq.clear();
        }
        else { transform(line.begin(), line.end(), line.begin(), ::toupper);
            seq += line;
        }
    }
    annotateContig(contig, seq, kmers, merge, cloud_print);
}

int main (int argc, const char * argv[]) {
    unordered_map<string,string> options;
    kmer_map kmers;
    ifstream fp_in;
    
    getOptions(argc, argv, options);
    
    openfile(fp_in, options.at("cloud_file"));
    getKmers(fp_in, kmers);
    if (getInt(options, "rev_comp")) { getRevKmers(kmers); }
    fp_in.close();
    
    openfile(fp_in, options.at("fasta_file"));
    annotateFasta(fp_in, kmers, getInt(options, "merge") , getInt(options, "print"));
    fp_in.close();
    
    return 0;
}

}
}
