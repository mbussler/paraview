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

enum {CONV_CLAMP, CONV_AUTO_CLAMP, CONV_RESCALE, CONV_NORMALIZE, CONV_NORMALIZE_POSITIVE}; 

void usage(const char *pname)
{
	fprintf(stderr, "Usage: %s  <infile.dat> <outfile> <outformat> <conversion>\n"
	                "       where <outformat> is one of:\n"
	                "          UCHAR\n"
	                "          CHAR\n"
	                "          USHORT\n"
	                "          SHORT\n"
	                "          UINT\n"
	                "          INT\n"
	                "          ULONG\n"
	                "          LONG\n"
	                "          FLOAT\n"
	                "          DOUBLE\n",
			pname);
	fprintf(stderr, "       and where <conversion> is one of:\n"
	                "          CLAMP:<min>,<max>    clamps the data elements at <min> and <max>\n"
	                "          AUTO_CLAMP           clamps the data elements at MIN(<outformat>) and MAX(<outformat>)\n"
	                "          RESCALE:<min>,<max>  rescales the data elements to [<min>,<max>]\n"
	                "          NORMALIZE            rescales the data elements [MIN(<outformat>),MAX(<outformat>)]\n"
	                "          NORMALIZE_POSITIVE   rescales the data elements [0, MAX(<outformat>)]\n");
	exit(1);
}

#define GET_DATA_MIN_MAX_BODY(TYPE) \
	case DR_FORMAT_##TYPE: \
		{ \
			unsigned long i; \
			DR_##TYPE *rdata = (DR_##TYPE *)srcData; \
			for (i = 0; i < size; i++) { \
				double f = (double)*rdata++; \
				if (f > *max) { \
					*max  = f; \
				} \
				if (f < *min) { \
					*min = f; \
				} \
			} \
		} \
		break;

void getDataMinMax(void *srcData, int srcFormat, unsigned long size,
                   double *min, double *max)
{
	*min = DBL_MAX;
	*max = -DBL_MAX;
	switch (srcFormat) {
		GET_DATA_MIN_MAX_BODY(CHAR)
		GET_DATA_MIN_MAX_BODY(UCHAR)
		GET_DATA_MIN_MAX_BODY(SHORT)
		GET_DATA_MIN_MAX_BODY(USHORT)
		GET_DATA_MIN_MAX_BODY(INT)
		GET_DATA_MIN_MAX_BODY(UINT)
		GET_DATA_MIN_MAX_BODY(LONG)
		GET_DATA_MIN_MAX_BODY(ULONG)
		GET_DATA_MIN_MAX_BODY(FLOAT)
		GET_DATA_MIN_MAX_BODY(DOUBLE)
		default:
			fprintf(stderr, "Error unsupported data format '%s'",
			        datRaw_getDataFormatName(srcFormat));
			break;
	}
}

#define RESCALE_BLOCK(SRC_FORMAT, DST_FORMAT) \
	{ \
		unsigned long i; \
		DR_##DST_FORMAT CMIN = (DR_##DST_FORMAT)convMin; \
		DR_##DST_FORMAT CMAX = (DR_##DST_FORMAT)convMax; \
		DR_##SRC_FORMAT *idata = (DR_##SRC_FORMAT *)srcData; \
		DR_##DST_FORMAT *odata = (DR_##DST_FORMAT *)dstData; \
		printf("Normalizing to %g - %g ...", (double)CMIN, (double)CMAX); \
		for (i = 0; i < size; i++) { \
			*odata = CMIN + (DR_##DST_FORMAT)(CMAX*(*idata - min)/(max - min)); \
			odata++; idata++; \
		} \
		printf("done\n"); \
	} \
	break;

#define RESCALE_DST_BLOCK(DST_FORMAT, SRC_FORMAT) \
case (DR_FORMAT_##DST_FORMAT): \
	RESCALE_BLOCK(SRC_FORMAT, DST_FORMAT) \
	break;

#define RESCALE_SRC_BLOCK(SRC_FORMAT) \
case (DR_FORMAT_##SRC_FORMAT): \
	switch(dstFormat) { \
		RESCALE_DST_BLOCK(CHAR, SRC_FORMAT) \
		RESCALE_DST_BLOCK(UCHAR, SRC_FORMAT) \
		RESCALE_DST_BLOCK(SHORT, SRC_FORMAT) \
		RESCALE_DST_BLOCK(USHORT, SRC_FORMAT) \
		RESCALE_DST_BLOCK(INT, SRC_FORMAT) \
		RESCALE_DST_BLOCK(UINT, SRC_FORMAT) \
		RESCALE_DST_BLOCK(LONG, SRC_FORMAT) \
		RESCALE_DST_BLOCK(ULONG, SRC_FORMAT) \
		RESCALE_DST_BLOCK(FLOAT, SRC_FORMAT) \
		RESCALE_DST_BLOCK(DOUBLE, SRC_FORMAT) \
		/* ADD NEW BASIC TYPES BEFORE THIS LINE !!!! */ \
		default: \
			fprintf(stderr, "Error unsupported data format '%s'\n", \
			        datRaw_getDataFormatName(dstFormat)); \
			break; \
	} \
	break;

void rescaleData(void *srcData, int srcFormat, unsigned long size, double min, double max,
                 void *dstData, int dstFormat, double convMin, double convMax)
{
	switch (srcFormat) {
		RESCALE_SRC_BLOCK(CHAR)
		RESCALE_SRC_BLOCK(UCHAR)
		RESCALE_SRC_BLOCK(SHORT)
		RESCALE_SRC_BLOCK(USHORT)
		RESCALE_SRC_BLOCK(INT)
		RESCALE_SRC_BLOCK(UINT)
		RESCALE_SRC_BLOCK(LONG)
		RESCALE_SRC_BLOCK(ULONG)
		RESCALE_SRC_BLOCK(FLOAT)
		RESCALE_SRC_BLOCK(DOUBLE)
		default:
			fprintf(stderr,
			        "Error unsupported data format '%s'",
			        datRaw_getDataFormatName(srcFormat));
			break;
	}
}

#define CLAMP_BLOCK(SRC_FORMAT, DST_FORMAT) \
	{ \
		unsigned long i; \
		DR_##DST_FORMAT CMIN = (DR_##DST_FORMAT)convMin; \
		DR_##DST_FORMAT CMAX = (DR_##DST_FORMAT)convMax; \
		DR_##SRC_FORMAT *idata = (DR_##SRC_FORMAT *)srcData; \
		DR_##DST_FORMAT *odata = (DR_##DST_FORMAT *)dstData; \
		printf("Clamping to %g - %g ...", (double)CMIN, (double)CMAX); \
		for (i = 0; i < size; i++) { \
			if (*idata < CMIN) { \
				*odata = CMIN; \
			} else if (*idata > CMAX) { \
				*odata = CMAX; \
			} else { \
				*odata = (DR_##DST_FORMAT)(*idata); \
			} \
			odata++; idata++; \
		} \
		printf("done\n"); \
	} \
	break;

#define CLAMP_DST_BLOCK(DST_FORMAT, SRC_FORMAT) \
case (DR_FORMAT_##DST_FORMAT): \
	CLAMP_BLOCK(SRC_FORMAT, DST_FORMAT) \
	break;

#define CLAMP_SRC_BLOCK(SRC_FORMAT) \
case (DR_FORMAT_##SRC_FORMAT): \
	switch(dstFormat) { \
		CLAMP_DST_BLOCK(CHAR, SRC_FORMAT) \
		CLAMP_DST_BLOCK(UCHAR, SRC_FORMAT) \
		CLAMP_DST_BLOCK(SHORT, SRC_FORMAT) \
		CLAMP_DST_BLOCK(USHORT, SRC_FORMAT) \
		CLAMP_DST_BLOCK(INT, SRC_FORMAT) \
		CLAMP_DST_BLOCK(UINT, SRC_FORMAT) \
		CLAMP_DST_BLOCK(LONG, SRC_FORMAT) \
		CLAMP_DST_BLOCK(ULONG, SRC_FORMAT) \
		CLAMP_DST_BLOCK(FLOAT, SRC_FORMAT) \
		CLAMP_DST_BLOCK(DOUBLE, SRC_FORMAT) \
		/* ADD NEW BASIC TYPES BEFORE THIS LINE !!!! */ \
		default: \
			fprintf(stderr, "Error unsupported data format '%s'\n", \
			        datRaw_getDataFormatName(dstFormat)); \
			break; \
	} \
	break;

void clampData(void *srcData, int srcFormat, unsigned long size,
               void *dstData, int dstFormat, double convMin, double convMax)
{
	switch (srcFormat) {
		CLAMP_SRC_BLOCK(CHAR)
		CLAMP_SRC_BLOCK(UCHAR)
		CLAMP_SRC_BLOCK(SHORT)
		CLAMP_SRC_BLOCK(USHORT)
		CLAMP_SRC_BLOCK(INT)
		CLAMP_SRC_BLOCK(UINT)
		CLAMP_SRC_BLOCK(LONG)
		CLAMP_SRC_BLOCK(ULONG)
		CLAMP_SRC_BLOCK(FLOAT)
		CLAMP_SRC_BLOCK(DOUBLE)
		default:
			fprintf(stderr,
			        "Error unsupported data format '%s'",
			        datRaw_getDataFormatName(srcFormat));
			break;
	}
}

int main(int argc, char *argv[])
{
    DatRawFileInfo info, info2;
	void *input;
	void *output;
	double min, max;
	int compress = 0;
	int outFormat = 0;
	int conversion = 0;
	double convMin = 0.0, convMax = 0.0;
	int t;

    if (argc < 5) {
		usage(argv[0]);
    }

	if (strchr(argv[2], '.')) {
		fprintf(stderr, "Output file name <%s> already has an extension!\n"
						"Continue? (y/n)", argv[2]);
		if (toupper(getchar()) != 'Y') {
			exit(0);
		}
	}

	if (!strcmp(argv[3], "UCHAR")) {
		outFormat = DR_FORMAT_UCHAR;
	} else if (!strcmp(argv[3], "CHAR")) {
		outFormat = DR_FORMAT_CHAR;
	} else if (!strcmp(argv[3], "USHORT")) {
		outFormat = DR_FORMAT_USHORT;
	} else if (!strcmp(argv[3], "SHORT")) {
		outFormat = DR_FORMAT_SHORT;
	} else if (!strcmp(argv[3], "UINT")) {
		outFormat = DR_FORMAT_UINT;
	} else if (!strcmp(argv[3], "INT")) {
		outFormat = DR_FORMAT_INT;
	} else if (!strcmp(argv[3], "ULONG")) {
		outFormat = DR_FORMAT_ULONG;
	} else if (!strcmp(argv[3], "LONG")) {
		outFormat = DR_FORMAT_LONG;
	} else if (!strcmp(argv[3], "FLOAT")) {
		outFormat = DR_FORMAT_FLOAT;
	} else if (!strcmp(argv[3], "DOUBLE")) {
		outFormat = DR_FORMAT_DOUBLE;
	} else {
		fprintf(stderr, "Unknown output format: ""%s""\n", argv[3]);
		usage(argv[0]);
	}

	if (!strncmp(argv[4], "CLAMP", 5)) {
		char *s = argv[4] + 5;
		char *s2 = NULL;
		if (*s != ':') {
			goto conv_error;
		}
		convMin = strtod(++s, &s2);
		if (*s2 != ',') {
			goto conv_error;
		}
		convMax = strtod(++s2, &s);
		if (*s != '\0') {
			goto conv_error;
		}
		conversion = CONV_CLAMP;
	} else if (!strcmp(argv[4], "AUTO_CLAMP")) {
		conversion = CONV_AUTO_CLAMP;
		switch (outFormat) {
			case DR_FORMAT_CHAR:
				convMin = CHAR_MIN;
				convMax = CHAR_MAX;
				break;
			case DR_FORMAT_UCHAR:
				convMin = 0;
				convMax = UCHAR_MAX;
				break;
			case DR_FORMAT_SHORT:
				convMin = SHRT_MIN;
				convMax = SHRT_MAX;
				break;
			case DR_FORMAT_USHORT:
				convMin = 0;
				convMax = USHRT_MAX;
				break;
			case DR_FORMAT_INT:
				convMin = INT_MIN;
				convMax = INT_MAX;
				break;
			case DR_FORMAT_UINT:
				convMin = 0;
				convMax = UINT_MAX;
				break;
			case DR_FORMAT_LONG:
				convMin = LONG_MIN;
				convMax = LONG_MAX;
				break;
			case DR_FORMAT_ULONG:
				convMin = 0;
				convMax = ULONG_MAX;
				break;
			case DR_FORMAT_FLOAT:
				convMin = FLT_MIN;
				convMax = FLT_MAX;
				break;
			case DR_FORMAT_DOUBLE:
				convMin = DBL_MIN;
				convMax = DBL_MAX;
				break;
			default:
			fprintf(stderr, "Invalid output format\n");
			usage(argv[0]);
		}
	} else if (!strncmp(argv[4], "RESCALE", 7)) {
		char *s = argv[4] + 7;
		char *s2 = NULL;
		if (*s != ':') {
			goto conv_error;
		}
		conversion = CONV_CLAMP;
		convMin = strtod(++s, &s2);
		if (*s2 != ',') {
			goto conv_error;
		}
		convMax = strtod(++s2, &s);
		if (*s != '\0') {
			goto conv_error;
		}
		conversion = CONV_RESCALE;
	} else if (!strcmp(argv[4], "NORMALIZE")) {
		conversion = CONV_NORMALIZE;
		switch (outFormat) {
			case DR_FORMAT_CHAR:
				convMin = CHAR_MIN;
				convMax = CHAR_MAX;
				break;
			case DR_FORMAT_UCHAR:
				convMin = 0;
				convMax = UCHAR_MAX;
				break;
			case DR_FORMAT_SHORT:
				convMin = SHRT_MIN;
				convMax = SHRT_MAX;
				break;
			case DR_FORMAT_USHORT:
				convMin = 0;
				convMax = USHRT_MAX;
				break;
			case DR_FORMAT_INT:
				convMin = INT_MIN;
				convMax = INT_MAX;
				break;
			case DR_FORMAT_UINT:
				convMin = 0;
				convMax = UINT_MAX;
				break;
			case DR_FORMAT_LONG:
				convMin = LONG_MIN;
				convMax = LONG_MAX;
				break;
			case DR_FORMAT_ULONG:
				convMin = 0;
				convMax = ULONG_MAX;
				break;
			case DR_FORMAT_FLOAT:
				convMin = FLT_MIN;
				convMax = FLT_MAX;
				break;
			case DR_FORMAT_DOUBLE:
				convMin = DBL_MIN;
				convMax = DBL_MAX;
				break;
			default:
			fprintf(stderr, "Invalid output format\n");
			usage(argv[0]);
		}
	} else if (!strcmp(argv[4], "NORMALIZE_POSITIVE")) {
		conversion = CONV_NORMALIZE_POSITIVE;
		switch (outFormat) {
			case DR_FORMAT_CHAR:
				convMin = 0;
				convMax = CHAR_MAX;
				break;
			case DR_FORMAT_UCHAR:
				convMin = 0;
				convMax = UCHAR_MAX;
				break;
			case DR_FORMAT_SHORT:
				convMin = 0;
				convMax = SHRT_MAX;
				break;
			case DR_FORMAT_USHORT:
				convMin = 0;
				convMax = USHRT_MAX;
				break;
			case DR_FORMAT_INT:
				convMin = 0;
				convMax = INT_MAX;
				break;
			case DR_FORMAT_UINT:
				convMin = 0;
				convMax = UINT_MAX;
				break;
			case DR_FORMAT_LONG:
				convMin = 0;
				convMax = LONG_MAX;
				break;
			case DR_FORMAT_ULONG:
				convMin = 0;
				convMax = ULONG_MAX;
				break;
			case DR_FORMAT_FLOAT:
				convMin = 0.0f;
				convMax = FLT_MAX;
				break;
			case DR_FORMAT_DOUBLE:
				convMin = 0.0;
				convMax = DBL_MAX;
				break;
			default:
			fprintf(stderr, "Invalid output format\n");
			usage(argv[0]);
		}
	} else {
conv_error:
		fprintf(stderr, "Invalid conversion: ""%s""\n", argv[4]);
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

	if (!(output = malloc(datRaw_getBufferSize(&info, outFormat)))) {
		fprintf(stderr, "Failed to allocate memory for output data\n");
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
	
	info2.dataFormat = outFormat;

	if (!datRaw_writeHeaderFile(&info2, NULL)) {
		fprintf(stderr, "Writing of header file failed\n");	
		exit(1);
	}

	for (t = 0; t < info.timeSteps; t++) {

		if (datRaw_getNext(&info, (void*)&input, DR_FORMAT_RAW) <= 0) {
			fprintf(stderr, "Error reading timestep %d\n", t);
			exit(1);
		}

		switch (conversion) {
			case CONV_CLAMP:
			case CONV_AUTO_CLAMP:
				clampData(input, info.dataFormat,
				          datRaw_getElementCount(&info)*info.numComponents,
				          output, outFormat, convMin, convMax);
				break;
			case CONV_RESCALE:
			case CONV_NORMALIZE:
			case CONV_NORMALIZE_POSITIVE:
				getDataMinMax(input, info.dataFormat,
				              datRaw_getElementCount(&info)*info.numComponents,
				              &min, &max);
				rescaleData(input, info.dataFormat, datRaw_getElementCount(&info)*info.numComponents,
				            min, max, output, outFormat, convMin, convMax);
				break;
			default:
				fprintf(stderr, "WTF: invalid conversion\n");
				exit(1);
		}

		if (!datRaw_writeTimestep(&info2, output, outFormat, compress, t)) {
			fprintf(stderr, "Writing of timestep %d failed\n", t);
			exit(1);
		}
	}

	datRaw_close(&info);
	datRaw_freeInfo(&info);
	datRaw_freeInfo(&info2);

	free(input);
	free(output);

	return 0;
}

