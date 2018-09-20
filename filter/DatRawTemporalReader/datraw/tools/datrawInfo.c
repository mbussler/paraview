#include <stdio.h>
#include <string.h>
#include <datRaw.h>
#include <malloc.h>
#include <sys/time.h>
#include <stdlib.h>
#include <unistd.h>
#include <float.h>
#include <limits.h>
#include <ctype.h>

#include "helpers.h"

void usage(const char*pname)
{
	printf("Usage: %s datfile\n", pname);
}

int main(int argc, char *argv[])
{
    DatRawFileInfo info;
	float *input, *fp;
	float *min, *max, *avrg;
	unsigned long s, size;
	int t, i;

    if (argc != 2) {
		usage(argv[0]);
		exit(1);
    }

	if (!datRaw_readHeader(argv[1], &info, NULL)) {
		fprintf(stderr, "Loading header file ""%s"" failed\n", argv[1]);
		exit(1);
	}

	datRaw_printInfo(&info);

	if (!(input = malloc(datRaw_getBufferSize(&info, DR_FORMAT_FLOAT)))) {
		fprintf(stderr, "Failed to allocate memory for input data\n");
		exit(1);
	}
	if (!(min = malloc(info.numComponents*sizeof(float)))) {
		fprintf(stderr, "Failed to allocate memory for min data\n");
		exit(1);
	}
	if (!(max = malloc(info.numComponents*sizeof(float)))) {
		fprintf(stderr, "Failed to allocate memory for max data\n");
		exit(1);
	}
	if (!(avrg = malloc(info.numComponents*sizeof(float)))) {
		fprintf(stderr, "Failed to allocate memory for average data\n");
		exit(1);
	}
	
	size = 1;
	for (i = 0; i < info.dimensions; i++) {
		size *= info.resolution[i];
	}
	
	for (i = 0; i < info.numComponents; i++) {
		min[i] = FLT_MAX;
		max[i] = -FLT_MAX;
		avrg[i] = 0.0f;
	}
	
	for (t = 0; t < info.timeSteps; t++) {
		if (datRaw_getNext(&info, (void*)&input, DR_FORMAT_FLOAT) <= 0) {
			fprintf(stderr, "Error reading timestep %d\n", t);
			exit(1);
		}
		fp = input;
		for (s = 0; s < size; s++) {
			int j;
			for (j = 0; j < info.numComponents; j++) {
				if (*fp > max[j]) {
					max[j]  = *fp;
				}
				if (*fp < min[j]) {
					min[j] = *fp;
				}
				avrg[j] += *fp;
				fp++;
			}
		}

		printf("Timestep %i:\n", t);
		for (i = 0; i < info.numComponents; i++) {
			avrg[i] /= size;
			printf("\tmin[%i] = %f  max[%i] = %f average[%i] = %f\n", i, min[i],
					i, max[i], i, avrg[i]);
		}
	}

	
	datRaw_close(&info);
	datRaw_freeInfo(&info);

	free(input);
	free(min);
	free(max);
	free(avrg);

	return 0;
}

