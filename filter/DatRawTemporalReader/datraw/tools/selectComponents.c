#include <stdio.h>
#include <string.h>
#include <datRaw.h>
#include <malloc.h>
#include <sys/time.h>
#include <stdlib.h>
#include <ctype.h>

#include "helpers.h"

void usage(const char *pname)
{
	fprintf(stderr, "Usage: %s (c | cstart-cend){,(c | cstart-cend)}*> <inFile.dat> <outFile>\n"
					"\tstart-end pairs denotes the first and last component to"
					" include in the output\n",
			pname);
	exit(1);
}

int main(int argc, char *argv[])
{
    DatRawFileInfo info, info2;
    void *ibuffer = NULL, *obuffer = NULL;
    int n, *selection, i, numSelections, rangeSelect;
	unsigned long size, numElements, j, elementSize, valueSize;
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


	selection = NULL;
	numSelections = 0;
	rangeSelect = 0;

	s = selectionString = strdup(argv[1]);
	while (s && *s) {
		int select, range;

		if (sscanf(s, "%d", &select) != 1) {
			fprintf(stderr, "Error in selection string: ...%s\n", s);
			exit(1);
		}
		if (select < 0) {
			fprintf(stderr, "Invalid component: %d\n", select);
			exit(1);
		}
		if (rangeSelect) {
			if (select <= selection[numSelections - 1]) {
				fprintf(stderr, "Error invalid range: start > end\n");
				exit(1);
			}
			range = select - selection[numSelections - 1];
			fprintf(stderr, "RANGE: %d - %d\n", selection[numSelections - 1],
			select);
		} else {
			range = 1;
		}
		if (!(selection = realloc(selection, (numSelections + range)*sizeof(int)))) {
			fprintf(stderr, "Failed to allocate selection info\n");
			exit(1);
		}
		if (rangeSelect) {
			for (i = 0; i < range; i++) {
				selection[numSelections + i] = selection[numSelections - 1] + i + 1;
			}
			numSelections += range;
		} else {
			selection[numSelections++] = select;
		}

		while (*s && isdigit(*s)) {
			s++;
		}
		if (!*s) {
			break;
		} else if (*s == ',') {
			s++;
			rangeSelect = 0;
		} else if ((*s == '-')) {
			if (rangeSelect) {
				fprintf(stderr, "Invalid range in range\n");
				exit(1);
			}
			s++;
			rangeSelect = 1;
		} else {
			fprintf(stderr, "Error in selection string\n");
			exit(1);
		}
		
	}
	
	for (i = 0; i < numSelections; i++) {
		if (selection[i] >= info.numComponents) {
			fprintf(stderr, "Selecting componnet %d not possible, only %d"
			                " components in dataset\n",
							selection[i],info.numComponents); 
			exit(1);
		}
	}

	fprintf(stderr, "SELECTION:");
	for (i = 0; i < numSelections; i++) {
		fprintf(stderr, "%d,", selection[i]);
	}
	fprintf(stderr, "\n");
	
	free(selectionString);
	
	datRaw_copyInfo(&info2, &info);
	
	replaceFileName(&info2, argv[3]);

	info2.numComponents = numSelections;

	datRaw_printInfo(&info2);

	if (!datRaw_writeHeaderFile(&info2, NULL)) {
		fprintf(stderr, "Writing header file failed\n");	
		exit(1);
	}
	
	numElements = datRaw_getElementCount(&info);
	elementSize = datRaw_getRecordSize(&info, DR_FORMAT_RAW);
	valueSize = elementSize/info.numComponents;

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

	for (n = 0; n < info.timeSteps; n++) {

		if (datRaw_getNext(&info, &ibuffer, DR_FORMAT_RAW) <= 0) {
			fprintf(stderr, "Failed to load timestep %d\n", n);
			exit(1);
		}
		
		dp = (char*)ibuffer;
		df = (char*)obuffer;
	
		for (j = 0; j < numElements; j++) {
			for (i = 0; i < numSelections; i++) {
				memcpy(df, dp + selection[i]*valueSize, valueSize);
				df += valueSize;
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
	free(selection);
	return 0;
}
