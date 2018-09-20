// test-rectilinear.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"


int create(char *datFile, char *rawFile) {
    DatRawFileInfo fi;
    void *data;
    float *cursor;

    int resolution[] = { 1, 2, 3 };
    float sliceDist[] = { 0.1f, 0.1f, 0.2f, 0.1f, 0.2f, 0.3f };

    if (!datRaw_createInfo(&fi, datFile, rawFile, 3, 1, DR_GRID_RECTILINEAR, 1, DR_FORMAT_FLOAT, resolution, sliceDist)) {
        std::cerr << "Creating dat file failed." << std::endl;
        return 1;
    }

    if (!datRaw_writeHeaderFile(&fi, NULL)) {
        std::cerr << "Saving dat file failed." << std::endl;
        return 2;
    }

    try {
        size_t size = datRaw_getBufferSize(&fi, DR_FORMAT_UCHAR);
        data = new unsigned char[size];
        ::memset(data, 0, size);
    } catch (std::bad_alloc) {
        std::cerr << "Could not allocate buffer for raw data." << std::endl;
        return 3;
    }

    // TODO: insert meaningful data.

    if (!datRaw_writeTimestep(&fi, data, DR_FORMAT_UCHAR, 0, 0)) {
        std::cerr << "Saving raw file failed." << std::endl;
        return 4;
    }

    datRaw_freeInfo(&fi);
    if (!datRaw_load(datFile, &fi, NULL, &data, DR_FORMAT_UCHAR)) {
        std::cerr << "Reading data failed." << std::endl;
        return 5;
    }

    datRaw_freeInfo(&fi);
    delete[] data;
    return 0;
}

int _tmain(int argc, TCHAR **argv) {
    ::create("rect.dat", "rect.raw");
    return 0;
}

