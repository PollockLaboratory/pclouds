/**
 * STP: This is clearly a cpp file. It used to be named *.c.
 *
 * This used to be named "genomepro.c". I have renamed it to "PClouds.cpp".
 *
 *
 */


#include <iostream>
#include <string>


#include "counts.h"
#include "clouds.h"

int main() {
	std::cout << "STEPHEN's hacked version of P-clouds" << std::endl;

	std::string controlfile = "Controlfile";

	//STP: Count kmers
	calculatecounts(controlfile.c_str());

	//STP: By 'dissect' he means 'possible make clouds and possible annotate
	// the genome'. Most of the time 'dissect' means 'annotate'.
	pcloudsdissection(controlfile.c_str());

	std::cout << std::endl << "P-clouds finished!" << std::endl;

	return (0);
}
