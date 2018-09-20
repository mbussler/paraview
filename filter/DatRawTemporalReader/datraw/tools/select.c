#include <stdio.h>
#include <string.h>
#include <datRaw.h>
#include <malloc.h>
#include <sys/time.h>
#include <stdlib.h>
#include <ctype.h>

#include "helpers.h"

static void usage(const char *pname)
{
	fprintf(stderr, "Usage: %s <start0-end0:start1-end1:...> <inFile.dat>"
	                " <outFile>\n"
					"\tstart-end pairs denotes the first and last slice to extract\n"
					"\tfor each dimension.\n",
			pname);
	exit(1);
}

int main(int argc, char *argv[])
{
    DatRawFileInfo info, info2;
    void *ibuffer = NULL, *obuffer = NULL;
    int n, *selection, i, *indices;
	unsigned long size, numElements, j, elementSize;
	char *s, *selectionString, *dp, *df;
	
    if (argc != 4) {
		usage(argv[0]);
    }

	if (strchr(argv[3], '.')) {
		fprintf(stderr, "Output file name <%s> already has an extension!\n"
						"Continue? (y/n)", argv[3]);
		if (toupper(getchar()) != 'Y') {
			exit(0);
		}
	}

	if (!datRaw_readHeader(argv[2], &info, NULL)) {
		fprintf(stderr, "Loading header file ""%s"" failed\n", argv[2]);
		exit(1);
	}

	datRaw_printInfo(&info);

	if (info.gridType != DR_GRID_CARTESIAN) {
		fprintf(stderr, "No cartesian grid!\n"); 
		exit(1);
	}

	if (!(selection = malloc(info.dimensions * 2 * sizeof(int)))) {
		fprintf(stderr, "Failed to allocate selection info\n");
		exit(1);
	}

	s = selectionString = strdup(argv[1]);
	for (i = info.dimensions - 1; i >= 0; i--) {
		if (!s) {
			fprintf(stderr, "Selection string to short: %d ranges given but %d dimensional data set\n",
					info.dimensions - i - 1, info.dimensions);
			exit(1);
		} else if ((s = strrchr(selectionString, ':'))) {
			if (sscanf(s + 1, "%d-%d", &selection[2*i], &selection[2*i + 1]) != 2) {
				fprintf(stderr, "Error parsing selection string: %s\n", argv[1]);
				exit(1);
			}
			*s = '\0';
		} else {
			if (sscanf(selectionString, "%d-%d", &selection[2*i], &selection[2*i + 1]) != 2) {
				fprintf(stderr, "Error parsing selection string: %s\n", argv[1]);
				exit(1);
			}
		}
	}
	if (s != NULL) {
		fprintf(stderr, "Selection string to long: more than %d ranges given for %d dimensional data\n",
				info.dimensions, info.dimensions);
	}
	free(selectionString);
	
	datRaw_copyInfo(&info2, &info);
	
	for (i = 0; i < info.dimensions; i++) {
		if (selection[2*i] < 0 ||
			selection[2*i + 1] >= info.resolution[i] ||
			selection[2*i] > selection[2*i + 1]) {
			fprintf(stderr, "invalid range %d-%d given for dimension %d\n", 
					selection[2*i], selection[2*i + 1], i);
			exit(1);
		}
		info2.resolution[i] = selection[2*i + 1] - selection[2*i] + 1;
	}

	replaceFileName(&info2, argv[3]);

	datRaw_printInfo(&info2);

	if (!datRaw_writeHeaderFile(&info2, NULL)) {
		fprintf(stderr, "Writing header file failed\n");	
		exit(1);
	}
	
	numElements = datRaw_getElementCount(&info);
	elementSize = datRaw_getRecordSize(&info, DR_FORMAT_RAW);

	size = datRaw_getBufferSize(&info, DR_FORMAT_RAW);
	if (!(ibuffer = malloc(size))) {
		fprintf(stderr, "Failed to allocate buffer for input data\n");
		exit(1);
	}
	size = datRaw_getBufferSize(&info2, DR_FORMAT_RAW);
	if (!(obuffer = malloc(size))) {
		fprintf(stderr, "Failed to allocate buffer for selected data\n");
		exit(1);
	}

	if (!(indices = malloc(info.dimensions*sizeof(int)))) {
		fprintf(stderr, "Failed to allocate indices\n");
		exit(1);
	}

	for (n = 0; n < info.timeSteps; n++) {

		if (datRaw_getNext(&info, &ibuffer, DR_FORMAT_RAW) <= 0) {
			fprintf(stderr, "Failed to load timestep %d\n", n);
			exit(1);
		}
		
		dp = (char*)ibuffer;
		df = (char*)obuffer;
		memset(indices, 0, info.dimensions*sizeof(int));
	
		for (j = 0; j < numElements; j++) {
			int inside = 1;
			for (i = 0; i < info.dimensions; i++) {
				if (indices[i] < selection[2*i] || indices[i] > selection[2*i + 1]) {
					inside = 0;
					break;
				}
			}
			if (inside) {
				memcpy(df, dp, elementSize);
				df += elementSize;
			}

			indices[0]++;
			for (i = 0; i < info.dimensions - 1; i++) {
				if (indices[i] == info.resolution[i]) {
					indices[i + 1]++;
					indices[i] = 0;
				}
			}
			dp += elementSize;
		}
		
		if (!datRaw_writeTimestep(&info2, obuffer, DR_FORMAT_RAW, 1, n)) {
			fprintf(stderr, "Writing of timestep %d failed\n", n);
			exit(1);
		}
	}

	datRaw_close(&info);	
	datRaw_freeInfo(&info);
	datRaw_freeInfo(&info2);

	free(ibuffer);
	free(obuffer);
	free(indices);
	free(selection);
	return 0;
}
