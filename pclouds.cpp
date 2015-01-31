#include <iostream>
#include <string>

using std::string;

namespace pclouds {
void build_clouds(string controlfile); // In build_clouds.cpp
void annotate_genome(string controlfile); // In annotate_genome.cpp
}

int main() {
	string controlfile = "controlfile";

	pclouds::build_clouds(controlfile);

	pclouds::annotate_genome(controlfile);

	return (0);
}
