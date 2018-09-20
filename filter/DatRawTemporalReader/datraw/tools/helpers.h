#ifndef HELPERS_H
#define HELPERS_H

#include <datRaw.h>

int isMultifileDescription(const char *filename);

char *getMultifileEnumeration(char *filename);

void replaceFileName(DatRawFileInfo *info, const char *newName);

#endif

