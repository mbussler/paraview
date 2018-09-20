#include <stdio.h>
#include <string.h>
#include <datRaw.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    DatRawFileInfo info;
    float *fbuffer;
    int i, j, result;
	unsigned long center;

    if (argc < 2) {
        fprintf(stderr, "Usage: %s <file.dat>\n", argv[0]);
        exit(1);;
    }

	if (!datRaw_readHeader(argv[1], &info, NULL)) {
		fprintf(stderr, "Reading header file failed!\n");
		exit(1);
	}

	datRaw_printInfo(&info);
	
#if 0	
	size = datRaw_getBufferSize(&info, DR_FORMAT_FLOAT);
	if (!(fbuffer = (float*)malloc(size))) {
		fprintf(stderr, "Failed to allocate data buffer\n");
		exit(1);
	}
#else
	fbuffer = NULL;
#endif

	center = info.numComponents;	
	for (i = 0; i < info.dimensions; i++) {
		center *= info.resolution[i]/2;
	}

	for (i = 0; i < info.timeSteps; i++) {

		result = datRaw_getNext(&info, (void*)&fbuffer, DR_FORMAT_FLOAT);
	
		if (result > 0) {
			fprintf(stderr, "time step %d loaded as float\n value at center = ", i);
			for (j = 0; j < info.numComponents; j++) {
				fprintf(stderr, "%f ", fbuffer[center + j]);
			}	
			fprintf(stderr, "\n");
		} else {
			fprintf(stderr, " loading time step %d as float failed\n", i);
			exit(1);
		}
	}
	
	for (i = info.timeSteps - 2; i >= 0; i--) {

		result = datRaw_getPrevious(&info, (void*)&fbuffer, DR_FORMAT_FLOAT);
	
		if (result) {
			fprintf(stderr, "time step %d loaded as float\n", i);
		} else {
			fprintf(stderr, " loading time step %d as float failed\n", i);
			exit(1);
		}
	}
	
	datRaw_close(&info);
	datRaw_freeInfo(&info);

	free(fbuffer);

	return 0;
}

