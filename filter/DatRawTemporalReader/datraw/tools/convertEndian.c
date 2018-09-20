#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <datRaw.h>
#include <malloc.h>
#include <sys/time.h>
#include <stdlib.h>
#include <unistd.h>

#include "helpers.h"

void usage(const char *pname)
{
	fprintf(stderr, "Usage: %s  <infile.dat> <outfile>\n", 
			pname);
	exit(1);
}

int main(int argc, char *argv[])
{
    DatRawFileInfo info, info2;
    int t;
	char *input = NULL;

    if (argc < 3) {
		usage(argv[0]);
    }


	if (!datRaw_readHeader(argv[1], &info, NULL)) {
		fprintf(stderr, "Loading header file ""%s"" failed\n", argv[1]);
		exit(1);
	}

	datRaw_printInfo(&info);

	if (!(input = malloc(datRaw_getBufferSize(&info, DR_FORMAT_RAW)))) {
		fprintf(stderr, "Failed to allocate memory for input data\n");
		exit(1);
	}

	datRaw_copyInfo(&info2, &info);
	free(info2.descFileName);
	free(info2.dataFileName);

	if (!(info2.descFileName = malloc(strlen(argv[2]) + 5))) {
		fprintf(stderr, "Failed to allocate info2.descFileName\n");
		exit(1);
	}
	sprintf(info2.descFileName, "%s.dat", argv[2]);

	if (info.multiDataFiles) { 
		char *fileEnum = getMultifileEnumeration(info.dataFileName); 
		if (!(info2.dataFileName = malloc(strlen(argv[2]) + strlen(fileEnum) + 5))) {
			fprintf(stderr, "Failed to allocate info2.dataFileName\n");
			exit(1);
		}
		sprintf(info2.dataFileName, "%s%s.raw", argv[2], fileEnum);
		free(fileEnum);
	} else {
		if (!(info2.dataFileName = malloc(strlen(argv[2]) + 5))) {
			fprintf(stderr, "Failed to allocate info2.dataFileName\n");
			exit(1);
		}
		sprintf(info2.dataFileName, "%s.raw", argv[2]);
	}
	
	if (info.byteOrder == DR_LITTLE_ENDIAN) {
		info2.byteOrder = DR_BIG_ENDIAN;
	} else {
		info2.byteOrder = DR_LITTLE_ENDIAN;
	}

	if (!datRaw_writeHeaderFile(&info2, NULL)) {
		fprintf(stderr, "Writing of header file failed\n");	
		exit(1);
	}

	for (t = 0; t < info.timeSteps; t++) {

		if (datRaw_getNext(&info, (void*)&input, DR_FORMAT_RAW) <= 0) {
			fprintf(stderr, "Error reading timestep %d\n", t);
			exit(1);
		}

		if (!datRaw_writeTimestep(&info2, input, info.dataFormat, 1, t)) {
			fprintf(stderr, "Writing of timestep %d failed\n", t);
			exit(1);
		}
	}

	datRaw_close(&info);
	datRaw_freeInfo(&info);
	datRaw_freeInfo(&info2);

	free(input);

	return 0;
}
