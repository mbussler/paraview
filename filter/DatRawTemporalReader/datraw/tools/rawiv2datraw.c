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
#include <math.h>

void usage(const char *pname)
{
	fprintf(stderr, "Usage: %s  infile.rawiv outfileName\n", 
			pname);
	exit(1);
}

static void swapByteOrder32(void *data, unsigned long size) 
{
	unsigned long i;
	DR_INT v, sv;
	DR_INT *idata = (DR_INT*)data;
	for (i = 0; i < size; i++) {
		v = idata[i];
		sv = (v & 0x000000FF);
		sv = ((v & 0x0000FF00) >> 0x08) | (sv << 0x08);
		sv = ((v & 0x00FF0000) >> 0x10) | (sv << 0x08);
		sv = ((v & 0xFF000000) >> 0x18) | (sv << 0x08);
		idata[i] = sv;
	}
}

int main(int argc, char *argv[])
{
	float sliceDists[3];
    DatRawFileInfo info;
	float *data;
	int sizes[3];
	FILE *inFile;
	
    if (argc < 3) {
		usage(argv[0]);
    }

	if (strchr(argv[2], '.')) {
		fprintf(stderr, "Output file name <%s> already has an extension!\n"
						"Continue? (y/n)", argv[2]);
		if (toupper(getchar()) != 'Y') {
			exit(0);
		}
	}

	if (!(inFile = fopen(argv[1], "r"))) {
		fprintf(stderr, "Error opening input file\n");
		exit(1);
	}

	/* read header */
	if (fseek(inFile, 6*sizeof(DR_FLOAT) + 2*sizeof(DR_INT), SEEK_SET)) {
		perror("Invalid rawiv file:");
		exit(1);
	}
	if (fread(sizes, sizeof(DR_INT), 3, inFile) != 3) {
		fprintf(stderr, "Invalid rawiv file: header to short\n");
		exit(1);
	}
	if (fseek(inFile, 3*sizeof(DR_FLOAT), SEEK_CUR)) {
		perror("Invalid rawiv file:");
		exit(1);
	}
	if (fread(sliceDists, sizeof(DR_FLOAT), 3, inFile) != 3) {
		fprintf(stderr, "Invalid rawiv file: header to short\n");
		exit(1);
	}

	swapByteOrder32(sliceDists, 3);
	swapByteOrder32(sizes, 3);

	fprintf(stderr, "Input data size: %dx%dx%d\n", sizes[0], sizes[1], sizes[2]);

	if (!(data = malloc(sizes[0]*sizes[1]*sizes[2]*sizeof(DR_FLOAT)))) {
		fprintf(stderr, "Allocating data failed\n");
		exit(1);
	}

	if (fread(data, sizes[0]*sizes[1]*sizes[2]*sizeof(DR_FLOAT), 1, inFile) != 1) {
		fprintf(stderr, "Error reading data\n");
		exit(1);
	}
	swapByteOrder32(data, sizes[0]*sizes[1]*sizes[2]);
	
	fclose(inFile);

	if (!datRaw_createInfo(&info, "","", 3, 1,
							DR_GRID_CARTESIAN,
							1,
							DR_FORMAT_FLOAT,
							sliceDists, sizes)) {
		fprintf(stderr, "Failed to create file info\n");
		exit(1);
	}

	if (asprintf(&info.descFileName, "%s.dat", argv[2]) == -1) {
			fprintf(stderr, "out of memory\n");
			exit(1);
	}
	if (asprintf(&info.dataFileName, "%s.raw", argv[2]) == -1) {
			fprintf(stderr, "out of memory\n");
			exit(1);
	}

	if (!datRaw_write(&info, NULL, data, DR_FORMAT_FLOAT, 1)) {
		fprintf(stderr, "Writing data failed\n");
		exit(1);
	}

	datRaw_freeInfo(&info);

	free(data);

	return 0;
}
