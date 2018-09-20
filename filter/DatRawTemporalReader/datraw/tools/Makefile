include ../common.mk

INCLUDES  = -I..
LOADLIBES = -ldatRaw -lm
LDFLAGS   += $(CFLAGS) -L../$(TARGET_DIR) $(LIBS)

all: select slice convertEndian slice-all combine selectComponents \
     rawiv2datraw convtiffraw convertFormat datrawInfo

select: select.o helpers.o

selectComponents: selectComponents.o helpers.o

slice: slice.o helpers.o

slice-all.o: slice-all.c

convertEndian: convertEndian.o helpers.o

convertFormat: convertFormat.o helpers.o

combine: combine.o helpers.o

rawiv2datraw: rawiv2datraw.o helpers.o

convtiffraw: convtiffraw.o
	$(CC) -o $@ $< $(LDFLAGS) -ltiff

datrawInfo: datrawInfo.o

%.o: %.c
	$(CC) -c $(CFLAGS) $(INCLUDES) $<

tags:
	ctags --recurse=yes *.c *.h $(SUBDIRS)

depend:
	makedepend -- $(CFLAGS) -- *.[ch]

cleandepend:
	makedepend
	rm -f Makefile.bak

clean: cleandepend
	rm -f slice select convertEndian slice-all combine rawiv2datraw \
	      selectComponents convertFormat convtiffraw datrawInfo *.o core
	
.PHONY: clean cleandepend lib depend

# DO NOT DELETE
