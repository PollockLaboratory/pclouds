//  main.cpp
//  cloudPos
//
//  Created by Corey Cox on 5-Feb-2015
//  Copyright (c) 2015 University of Colorado Denver, School of Medicine. All rights reserved.
//  Last updated on 17-Mar-2015

// This program takes an input FASTA formatted file and a set of P-Clouds and finds the locations of the
// P-Cloud kmers in the fasta sequence. Output is in bed file format.

#include <iostream>
#include <fstream>
#include <sstream>
#include <array>
#include <vector>

#include <algorithm>
#include <functional>
#include <locale>

//#include <time.h>
#include <string.h>
#include <unordered_map>
#include <map>

using namespace std;

// trim whitespace from end
//string &chomp(string &s) {
//    s.erase(find_if(s.rbegin(), s.rend(), not1(ptr_fun<int, int>(isspace))).base(), s.end());
//    return s;
//}

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
//        if (test > 0) { test--;
//            cerr << "Kmer: " << kmer << " cloud:" << cloud << "\n";
//        }
        kmers.all.emplace (kmer, cloud);
    }
    kmers.all.erase("");
    kmers.length = kmers.all.begin()->first.length();
//    cerr << "Found " << kmers.all.size() << " kmers.\n";
//    cerr << "kmer length: " << kmers.length << "\nFirst kmer: " << kmers.all.begin()->first << endl;
}

// get reverse compliment of kmers, place in search map and separate map.
void inline getRevKmers (kmer_map& kmers) {
//    if (kmers.all.count(kmers.all.begin()->first)) { return; }
    for (unordered_map<string,string>::iterator it = kmers.all.begin(); it != kmers.all.end(); it++) {
        kmers.rev.emplace(reverseCompliment(it->first), it->second);
    }
    kmers.all.insert(kmers.rev.begin(), kmers.rev.end());
//    cerr << "Found " << kmers.rev.size() << " reverse kmers.\n";
}

// find and merge locations matching kmers with strand orientation, then print.
void annotateContig (string& contig, string& seq, kmer_map& kmers, int& merge, bool& cloud_print) {
//    cerr << "Annotating contig: " << contig << " Sequence: " << seq << endl;
    const int len = seq.length() - kmers.length;
    int begin = 0; int end = 0; int num_plus = 0; int num_minus = 0;
    bool print = false; string start_cloud; string clouds = "";
//    int test = 10;
    for (int pos = 0; pos < len; pos++) {
        string kmer = seq.substr(pos, kmers.length);
//        if (test > 0) { test--;
//            cerr << "Kmer: " << kmer << "\n";
//        }
//
        if (1 == kmers.all.count(kmer)) { print = true;
//            cerr << "Kmer: " << kmer << " merge: " << merge << "\n"; exit(1);
            if (end + merge >= pos) {
//            if (pos < end + merge) {
//                cerr << "Internal kmer found: " << kmer << "Print value: " << print << endl;
                clouds += "\t" + kmers.all.at(kmer);
                end = pos + kmers.length;
                if (kmers.rev.count(kmer)) { num_minus++; } else { num_plus++; }
            }
            else {
//                cerr << "Starting kmer found: " << kmer << "Print value: " << print << endl;
                begin = pos; end = pos + kmers.length;
                start_cloud = kmers.all.at(kmer);
                clouds = kmers.all.at(kmer);
                if (kmers.rev.count(kmer)) { num_minus++; } else { num_plus++; }
            }
        }
        else {
//            if (print && pos > end + merge) {
            if (print && end + merge <= pos) {
//                cerr << "Entered printing section\n";
                cout << contig <<"\t"<< begin <<"\t"<< end << "\tpcloud:" << start_cloud << "\t0\t";
                if (num_plus >= num_minus) { cout << "+"; } else { cout << "-"; }
                if (cloud_print) { cout <<"\t"<< clouds; }
                cout << endl;
                print = false; num_plus = 0; num_minus = 0;
            }
        }
    }
    if (print) {
//        cerr << "Entered printing section\n";
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
//        getline( fp_in, line );
//        chomp( line );
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