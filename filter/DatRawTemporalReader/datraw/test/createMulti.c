#include <stdio.h>
#include <string.h>
#include <datRaw.h>
#include <malloc.h>
#include <stdlib.h>
#include <ctype.h>

#include "../datRaw_half.h"

void usage(const char *pname)
{
	fprintf(stderr, "Usage: %s <numTimeSteps>\n",
            pname);
	exit(1);
}

void writeTest(int format, int n)
{
    DatRawFileInfo info;
    int i;
	void *buffer;

	int resolution[] = {50, 60, 70};
	float sliceDist[] = {0.2f, 0.3f, 0.4f};
	
	if (!datRaw_createInfo(&info,
						  "test.dat",
						  "test%05+2*5d.raw",
						  3,
						  n,
						  DR_GRID_CARTESIAN,
						  3,
						  format,
						  resolution,
                          sliceDist)) {
		fprintf(stderr, "Create file info failed\n");
		exit(1);
	}
	
	/*datRaw_printInfo(&info);*/
	
	if (!datRaw_writeHeaderFile(&info, NULL)) {
		fprintf(stderr, "Writing header file failed\n");	
		exit(1);
	}

	if (!(buffer = malloc(datRaw_getBufferSize(&info, DR_FORMAT_UCHAR)))) {
		fprintf(stderr, "Failed to allocate buffer for input data\n");
		exit(1);
	}

	for (i = 0; i < n; i++) {
		
		memset(buffer, i%256, datRaw_getBufferSize(&info, DR_FORMAT_UCHAR));	
		/* write data and compress every other timestep */
		if (!datRaw_writeTimestep(&info, buffer, DR_FORMAT_UCHAR, i%2, i)) {
			fprintf(stderr, "Writing timestep %d failed\n", i);
			exit(1);
		}
	}

	datRaw_freeInfo(&info);
	free(buffer);
}

#define CHECK_VALUE_I(FORMAT, VALUE) \
case DR_FORMAT_##FORMAT: \
	fprintf(stderr, "%i: %i\n", i, VALUE); \
	if ((DR_##FORMAT)i != VALUE) { \
		fprintf(stderr, "Error!!\n"); \
		exit(1); \
	} \
	break;

#define CHECK_VALUE_LI(FORMAT, VALUE) \
case DR_FORMAT_##FORMAT: \
	fprintf(stderr, "%i: %li\n", i, VALUE); \
	if ((DR_##FORMAT)i != VALUE) { \
		fprintf(stderr, "Error!!\n"); \
		exit(1); \
	} \
	break;

#define CHECK_VALUE_F(FORMAT, VALUE) \
case DR_FORMAT_##FORMAT: \
	fprintf(stderr, "%i: %f\n", i, VALUE); \
	if ((DR_##FORMAT)i != VALUE) { \
		fprintf(stderr, "Error!!\n"); \
		exit(1); \
	} \
	break;
#define CHECK_VALUE_H(FORMAT, VALUE) \
case DR_FORMAT_HALF: \
	fprintf(stderr, "%i: %f\n", i, halfToFloat(VALUE)); \
	if ((DR_FLOAT)i != halfToFloat(VALUE)) { \
		fprintf(stderr, "Error!!\n"); \
		exit(1); \
	} \
	break;

void readTest(int format)
{
    DatRawFileInfo info;
    int i;
	void *buffer;
	if (!datRaw_readHeader("test.dat", &info, NULL)) {
		fprintf(stderr, "Loading header file failed\n");
		exit(1);
	}
	
	/*datRaw_printInfo(&info);*/

	if (!datRaw_checkInfo(&info,
						  NULL,
						  NULL,
						  3,
						  0,
						  DR_GRID_CARTESIAN,
						  3,
						  0,
						  NULL,
						  NULL)) {
		fprintf(stderr, "Check file info failed, hm?\n");
		exit(1);
	}

	if (!(buffer = malloc(datRaw_getBufferSize(&info, format)))) {
		fprintf(stderr, "Failed to allocate buffer for input data\n");
		exit(1);
	}

	for (i = 0; i < info.timeSteps; i++) {
		if (datRaw_getNext(&info, (void*)&buffer, format) <= 0) {
			fprintf(stderr, "Loading timestep %d failed\n", i);
			exit(1);
		}
		switch (format) {
			CHECK_VALUE_I(CHAR, ((DR_CHAR*)buffer)[0])
			CHECK_VALUE_I(UCHAR, ((DR_UCHAR*)buffer)[0])
			CHECK_VALUE_I(SHORT, ((DR_SHORT*)buffer)[0])
			CHECK_VALUE_I(USHORT, ((DR_USHORT*)buffer)[0])
			CHECK_VALUE_I(INT, ((DR_INT*)buffer)[0])
			CHECK_VALUE_I(UINT, ((DR_UINT*)buffer)[0])
			CHECK_VALUE_LI(LONG, ((DR_LONG*)buffer)[0])
			CHECK_VALUE_LI(ULONG, ((DR_ULONG*)buffer)[0])
			CHECK_VALUE_H(HALF, ((DR_HALF*)buffer)[0])
			CHECK_VALUE_F(FLOAT, ((DR_FLOAT*)buffer)[0])
			CHECK_VALUE_F(DOUBLE, ((DR_DOUBLE*)buffer)[0])
		}	
	}

	datRaw_close(&info);
	datRaw_freeInfo(&info);

	free(buffer);
}

int main(int argc, char *argv[])
{
	int n, i, j;
	
    if (argc != 2) {
		usage(argv[0]);
    }

	n = atoi(argv[1]);

	for (i = 1; i < DR_FORMAT_RAW; i++) {
		fprintf(stderr, "Write as %s\n",datRaw_getDataFormatName(i));
		writeTest(i, n);
		for (j = 1; j < DR_FORMAT_RAW; j++) {
			fprintf(stderr, "Read as %s\n",datRaw_getDataFormatName(j));
			readTest(j);
		}
	}
	
	return 0;
}

