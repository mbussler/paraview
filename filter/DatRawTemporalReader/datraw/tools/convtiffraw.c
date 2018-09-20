#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <tiffio.h>
#include <stdlib.h>
#include <ctype.h>

int main(int argc, char *argv[])
{
	
	TIFF* tif;
	
	if (argc != 3) {
		fprintf(stderr, "Usage: %s <in.tif> <out.raw>\n"
		                "       NOTE: testet only for 16bit gray tiffs!\n", 
		        argv[0]);
		exit(1);
	}
	
	tif = TIFFOpen(argv[1], "r");

	
	if (tif) {
		FILE *of;
		uint32 imagelength;
		tdata_t buf;
		uint32 row;

		of = fopen(argv[2], "w");

		if (of) {
			uint32 scanlinesize = TIFFScanlineSize(tif);
			TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &imagelength);
			
			buf = _TIFFmalloc(scanlinesize);
			for (row = 0; row < imagelength; row++) {
				if (TIFFReadScanline(tif, buf, row, 0) < 0) {
					fprintf(stderr, "reading scanline %i failed :(\n", row);
					exit(1);	
				}
				if (fwrite(buf, scanlinesize, 1, of) != 1) {
					fprintf(stderr, "writing scanline %i failed :(\n", row);
					exit(1);	
				}
			}
			_TIFFfree(buf);
		} else {
			fprintf(stderr, "open output file failed :(\n");
			exit(1);	
		}
		fclose(of);
		TIFFClose(tif);
		return 0;
	} else {
		fprintf(stderr, "open input file failed :(\n");
		exit(1);	
	}
	return 0;
}

