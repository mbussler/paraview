#include <stdio.h>
#include <float.h>
#include <stdlib.h>
#include <math.h>

#include "../datRaw_half.h"

#define CHECK(F) printf("%f -> %f (%f)\n", (F), halfToFloat(floatToHalf(F)), fabs(F - halfToFloat(floatToHalf(F))));

int main(int argc, char **argv)
{
	int i, j;
	
	CHECK(0.0);
	CHECK(-0.0);
	CHECK(1.0);
	CHECK(-1.0);
	CHECK(1000.0);
	CHECK(-1000.0);
	CHECK((float)0xFFFF);
	CHECK(-(float)0xFFFF);

	for (j = 0; j < 6; j++) {
		for (i = 0; i < 20; i++) {
			float f = (float)pow(10, j)*(float)rand()/(float)RAND_MAX;
			CHECK(f)
		}
	}
	return 0;
}
