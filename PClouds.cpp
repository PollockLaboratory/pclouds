

#include <iostream>
#include <string>

using namespace std;


// In clouds.cpp
// I didn't want to have a header file for two functions
void build_clouds(string pchControlfile);
void annotate_genome(string pchControlfile);


int main() {
	cout << "STEPHEN's hacked version of P-clouds" << endl;

	string controlfile = "Controlfile";

	build_clouds(controlfile);

	annotate_genome(controlfile);

	std::cout << std::endl <<
			"STEPHEN's hacked version of P-clouds finished!" << std::endl;

	return (0);
}
