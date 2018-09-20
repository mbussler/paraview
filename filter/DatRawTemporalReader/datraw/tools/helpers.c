#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "helpers.h"

int isMultifileDescription(const char *filename)
{
	const char *s = filename;
	while (*s) {
		if (s[0] == '%' && s[1] != '%') {
			return 1;
		}
		s++;
	}
	return 0;
}

char *getMultifileEnumeration(char *filename) 
{
	char *p, *s, *q;

	q = filename;
	s = NULL;
	while (*q) {
		if (q[0] == '%' && q[1] != '%') {
			if (s) {
				fprintf(stderr, "Multi file description contains more than "
				                "one varying\n");
				return NULL;
			} else {
				s = q;	
			}
		}
		q++;
	}

	if (!(s = p = strdup(s))) {
		fprintf(stderr, "Failed to allocate MultifileEnumeration string\n");
		exit(1);
	}

	if (!s || *s == '\0') {
		fprintf(stderr, "Strange input dataFileName: %s\n",
				filename);
		exit(1);
	}
	s++;	
	/* skip flags */
	if (*s == '0' || *s == ' ' || *s == '-') {
		s++;
	}	

	/* skip width */
	while (*s && isdigit(*s)) s++;

	/* skip offset */
	if (*s && *s == '+') {
		s++;
		while (*s && isdigit(*s)) s++;
	}
	
	/* skip stride */
	if (*s && *s == '*') {
		s++;
		while (*s && isdigit(*s)) s++;
	}
	
	/* skip conversion modifiers */
	if (*s == 'd') {
		s++;
	} else {
		fprintf(stderr, "invalid file enumeration: %s\n", s);
		exit(1);
	}

	*s = '\0';

	fprintf(stderr, "FileEnum: %s\n", p);

	return p;
}		

void replaceFileName(DatRawFileInfo *info, const char *newName) 
{
	char *fileEnum = NULL;
	
	if (info->multiDataFiles) {
		fileEnum = getMultifileEnumeration(info->dataFileName);
	}
	
	free(info->descFileName);
	free(info->dataFileName);

	if (!(info->descFileName = malloc(strlen(newName) + 5))) {
		fprintf(stderr, "Failed to replace descFileName\n");
		exit(1);
	}
	sprintf(info->descFileName, "%s.dat", newName);

	if (info->multiDataFiles) { 
		if (!(info->dataFileName = malloc(strlen(newName) + strlen(fileEnum) + 5))) {
			fprintf(stderr, "Failed to replace dataFileName\n");
			exit(1);
		}
		sprintf(info->dataFileName, "%s%s.raw", newName, fileEnum);
		free(fileEnum);
	} else {
		if (!(info->dataFileName = malloc(strlen(newName) + 5))) {
			fprintf(stderr, "Failed to replace dataFileName\n");
			exit(1);
		}
		sprintf(info->dataFileName, "%s.raw", newName);
	}
}

