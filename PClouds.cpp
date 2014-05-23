

#include <iostream>
#include <string>

// In clouds.cpp
int build_clouds_and_annotate_genome(const char* pchControlfile);


int main() {
	std::cout << "STEPHEN's hacked version of P-clouds" << std::endl;

	std::string controlfile = "Controlfile";

	//STP: Count kmers
	//calculatecounts(controlfile.c_str());
	// Use Jellyfish or Perl for kmer counting

	//STP: By 'dissect' he means 'possible make clouds and possible annotate
	// the genome'. Most of the time 'dissect' means 'annotate'.
	build_clouds_and_annotate_genome(controlfile.c_str());

	std::cout << std::endl << "P-clouds finished!" << std::endl;

	return (0);
}
