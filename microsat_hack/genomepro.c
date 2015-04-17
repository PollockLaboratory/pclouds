#include "counts.h"
#include "clouds.h"
// debug
#include <cstdio>
int main()
{
	char* pchControlfile = "Controlfile";

	calculatecounts(pchControlfile);
        
        printf("Calculating...\n");
	pcloudsdissection(pchControlfile);
        // debug
        printf("DONE!\n");

	return (0);
}
