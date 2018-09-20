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
	fprintf(stderr, "Usage: %s <sliceNum> <inFile.dat> <outFile>\n",
	        pname);
	exit(1);
}

int main(int argc, char *argv[])
{
    DatRawFileInfo info, info2;
    void *ibuffer = NULL;
    int i, slice;
	unsigned long size, offset;

	
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

	fprintf(stderr, "Input ");
	datRaw_printInfo(&info);

	if (info.gridType != DR_GRID_CARTESIAN) {
		fprintf(stderr, "No cartesian grid!\n"); 
		exit(1);
	}

	if (sscanf(argv[1], "%d", &slice) != 1) {
		fprintf(stderr, "Error parsing sliceNum argument: %s\n", argv[1]);
		exit(1);
	}

	if (slice < 0 || slice >= info.resolution[info.dimensions - 1]) {
		fprintf(stderr, "Selected slice %d out of bound\n", slice);
		exit(1);
	}

	datRaw_copyInfo(&info2, &info);

	replaceFileName(&info2, argv[3]);
	info2.dimensions--;  

	fprintf(stderr, "Output ");
	datRaw_printInfo(&info2);
	
	if (!datRaw_writeHeaderFile(&info2, NULL)) {
		fprintf(stderr, "Writing header file failed\n");	
		exit(1);
	}

	size = datRaw_getBufferSize(&info, DR_FORMAT_RAW);
	if (!(ibuffer = malloc(size))) {
		fprintf(stderr, "Failed to allocate buffer for input data\n");
		exit(1);
	}

	offset = 1;
	for (i = 0; i < info.dimensions - 1; i++) {
		offset *= info.resolution[i];
	}
	offset *= slice * datRaw_getRecordSize(&info, DR_FORMAT_RAW);

	for (i = 0; i < info2.timeSteps; i++) {
		
		if (datRaw_getNext(&info, &ibuffer, DR_FORMAT_RAW) <= 0) {
			fprintf(stderr, "Failed to load timestep %d\n", i);
			exit(1);
		}
		
		if (!datRaw_writeTimestep(&info2, (char*)ibuffer + offset,
		                          DR_FORMAT_RAW, 1, i)) {
			fprintf(stderr, "Writing step %d failed\n", i);
			exit(1);
		}
	}
	
	datRaw_close(&info2);	
	datRaw_freeInfo(&info);
	datRaw_freeInfo(&info2);

	free(ibuffer);
	return 0;
}
