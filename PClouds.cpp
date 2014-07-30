

#include <iostream>
#include <string>

using namespace std;


// I didn't want to have a header file for two functions
void build_clouds(string controlfile); // In build_clouds.cpp
void annotate_genome(string controlfile); // In annotate_genome.cpp


int main() {
	cout << "STEPHEN's hacked version of P-clouds\n" << endl;

	string controlfile = "Controlfile";

	build_clouds(controlfile);

	annotate_genome(controlfile);

	cout << "\nSTEPHEN's hacked version of P-clouds finished!\n";

	return (0);
}
