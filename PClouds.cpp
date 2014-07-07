

#include <iostream>
#include <string>

using namespace std;


// In clouds.cpp
// I didn't want to have a header file for two functions
void build_clouds(string controlfile);
void annotate_genome(string controlfile);


int main() {
	cout << "STEPHEN's hacked version of P-clouds\n" << endl;

	string controlfile = "Controlfile";

	build_clouds(controlfile);

	annotate_genome(controlfile);

	cout << "\nSTEPHEN's hacked version of P-clouds finished!\n";

	return (0);
}
