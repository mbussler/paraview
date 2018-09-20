#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <datRaw.h>
#include <malloc.h>
#include <sys/time.h>
#include <stdlib.h>

void usage(const char *pname)
{
	fprintf(stderr, "Usage: %s <datfileIn >\n",
			pname);
	exit(1);
}

int main(int argc, char *argv[])
{
    DatRawFileInfo info, info2;
    void *buffer = NULL, *obuffer = NULL;
    int i, j, result;
	unsigned long nsize, osize, offset;
	char *ext, *filename, *di, *dp, *base;
	
    if (argc < 2) {
		usage(argv[0]);
    }

	if (!(result = datRaw_load(argv[1], &info, NULL, &buffer, DR_FORMAT_RAW))) {
		fprintf(stderr, "loading file %s failed\n", argv[2]);
		exit(1);
	}
	datRaw_printInfo(&info);
	datRaw_close(&info);

	if (info.gridType != DR_GRID_CARTESIAN) {
		fprintf(stderr, "No cartesian grid!\n"); 
		exit(1);
	}

	for (j = 0; j < info.resolution[info.dimensions - 1]; j++) {
		
		datRaw_copyInfo(&info2, &info);

		if ((ext = strrchr(info2.descFileName, '.'))) {
			*ext ='\0';
			ext++;
		} else {
			ext = "";
		}
		if ((base = strrchr(info2.descFileName, '/'))) {
			base++;
		} else {
			base = info2.descFileName;
		}
		if (asprintf(&filename, "%s_%04d.%s", base, j, ext) == -1) {
			fprintf(stderr, "out of memory\n");
			exit(1);
		}
		free(info2.descFileName);
		info2.descFileName = filename;
		if ((ext = strrchr(info2.dataFileName, '.'))) {
			*ext ='\0';
			ext++;
		} else {
			ext = "";
		}
		if ((base = strrchr(info2.dataFileName, '/'))) {
			base++;
		} else {
			base = info2.dataFileName;
		}
		if (asprintf(&filename, "%s_%04d.%s", base, j, ext) == -1) {
			fprintf(stderr, "out of memory\n");
			exit(1);
		}
		free(info2.dataFileName);
		info2.dataFileName = filename;

		info2.dimensions--;  

		nsize = datRaw_getBufferSize(&info2, DR_FORMAT_RAW);
		osize = datRaw_getBufferSize(&info, DR_FORMAT_RAW);
		
		if (!(obuffer = malloc(nsize * info2.timeSteps))) {
			fprintf(stderr, "Failed to allocate buffer for selected data\n");
			exit(1);
		}

		offset = 1;
		for (i = 0; i < info.dimensions - 1; i++) {
			offset *= info.resolution[i];
		}
		offset *= j * datRaw_getRecordSize(&info, DR_FORMAT_RAW);

		di = (char*)buffer;
		dp = (char*)obuffer;
		for (i = 0; i < info2.timeSteps; i++) {
			memcpy(dp, di + offset, nsize);
			dp += nsize;
			dp += osize;
		}

		if (!datRaw_write(&info2, NULL, obuffer, DR_FORMAT_RAW, 0)) {
			fprintf(stderr, "Writing data failed\n");
			exit(1);
		}
		fprintf(stdout, "New dat-file: %s\nNew raw-file: %s\n", 
				info2.descFileName, info2.dataFileName);
		datRaw_close(&info2);	
	}
	datRaw_freeInfo(&info);
	datRaw_freeInfo(&info2);

	free(buffer);
	free(obuffer);
	return 0;
}
