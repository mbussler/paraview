include ../common.mk

CFLAGS += -fPIC

OBJS = adler32.o compress.o crc32.o gzio.o uncompr.o deflate.o trees.o \
       zutil.o inflate.o infback.o inftrees.o inffast.o

libz.a: $(OBJS)
	$(AR) $(ARFLAGS) $@ $^ 
	-@ ($(RANLIB) $@ || true) >/dev/null 2>&1

clean: cleandepend
	rm -f *.o *~ libz.* 

cleandepend:
	makedepend
	rm -f Makefile.bak

tags:
	ctags *.[ch]

depend:
	makedepend -- $(CFLAGS) -- *.[ch]

.PHONY: depend cleandepend clean tags

# DO NOT DELETE
