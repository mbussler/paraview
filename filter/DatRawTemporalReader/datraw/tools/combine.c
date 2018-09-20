#define _GNU_SOURCE
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

void usage(const char *pname)
{
	fprintf(stderr, "Usage: %s  infile_1.dat infile_2.dat ... outfileName\n", 
			pname);
	exit(1);
}

int main(int argc, char *argv[])
{
    DatRawFileInfo *info, info2;
    int i, result, t, bytePerElement;
	char **fd, **input, *output, *fo;
	unsigned long numElements, j;

    if (argc < 4) {
		usage(argv[0]);
    }

	if (strchr(argv[argc - 1], '.')) {
		fprintf(stderr, "Output file name <%s> already has an extension!\n"
						"Continue? (y/n)", argv[2]);
		if (toupper(getchar()) != 'Y') {
			exit(0);
		}
	}

	if (!(info = (DatRawFileInfo*)malloc(sizeof(DatRawFileInfo)*(argc - 2)))) {
		fprintf(stderr, "Failed to allocate input info\n");
		exit(1);
	}
	if (!(input = (char**)malloc(sizeof(char*)*(argc - 2)))) {
		fprintf(stderr, "Failed to allocate input data array\n");
		exit(1);
	}
	if (!(fd = (char**)malloc(sizeof(char*)*(argc - 2)))) {
		fprintf(stderr, "Failed to allocate input pointer array\n");
		exit(1);
	}
	
	for (i = 0; i < argc - 2; i++) {
		result = datRaw_readHeader(argv[i + 1], &info[i], NULL);
	
		if (!result) {
			fprintf(stderr, "loading file %d failed\n", i);
			exit(1);
		}

		datRaw_printInfo(&info[i]);
		

		if (info[i].numComponents != 1 || info[i].gridType != DR_GRID_CARTESIAN) {
			fprintf(stderr, "'%s': Not scalar data on cartesian grid\n", argv[i + 1]);
			exit(1);
		}

		if (!(input[i] = malloc(datRaw_getBufferSize(&info[i], DR_FORMAT_RAW)))) {
			fprintf(stderr, "Failed to allocate memory for input data\n");
			exit(1);
		}

		if (!datRaw_checkInfo(&info[i], NULL, NULL,
		                      info[0].dimensions,
							  info[0].timeSteps,
							  info[0].gridType,
							  info[0].numComponents,
							  info[0].dataFormat,
							  info[0].sliceDist,
							  info[0].resolution)) {
			fprintf(stderr, "Input data does not match!\n");
			exit(1);
		}
		
	}
	
	datRaw_copyInfo(&info2, &info[0]);

	replaceFileName(&info2, argv[argc - 1]);
	
	info2.numComponents = argc - 2;

	datRaw_printInfo(&info2);

	if (!datRaw_writeHeaderFile(&info2, NULL)) {
		fprintf(stderr, "Writing of header file failed\n");	
		exit(1);
	}

	if (!(output = malloc(datRaw_getBufferSize(&info2, DR_FORMAT_RAW)))) {
		fprintf(stderr, "Failed to allocate memory for output data\n");
		exit(1);
	}

	bytePerElement = datRaw_getRecordSize(&info[0], DR_FORMAT_RAW);
	numElements = datRaw_getElementCount(&info[0]);

	for (t = 0; t < info[0].timeSteps; t++) {

		fo = output;
		
		for (i = 0; i < info2.numComponents; i++) {
			if (!datRaw_getNext(&info[i], (void*)&input[i], DR_FORMAT_RAW)) {
				fprintf(stderr, "Error reading timestep %d\n", t);
				exit(1);
			}
			fd[i] = input[i];
		}
		
		for (j = 0; j < numElements; j++) {
			for (i = 0; i < info2.numComponents; i++) {
				memcpy(fo, fd[i], bytePerElement);
				fo += bytePerElement;
				fd[i] += bytePerElement;
			}	
		}

		if (!datRaw_writeTimestep(&info2, output, DR_FORMAT_RAW, 1, t)) {
			fprintf(stderr, "Writing of timestep %d failed\n", t);
			exit(1);
		}
	}

	for (i = 0; i < info2.numComponents; i++) {
		datRaw_close(&info[i]);
		datRaw_freeInfo(&info[i]);
		free(input[i]);
	}	
	datRaw_freeInfo(&info2);

	free(output);

	return 0;
}

