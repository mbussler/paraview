include common.mk

CFLAGS += -std=c99 -fPIC

TARGET = libdatRaw.a

OBJECTS = datRaw.o datRaw_log.o datRaw_half.o

# this is ugly, but since there is no "merge" provided by ar this is the easiest way
ZLIBOBJS = adler32.o compress.o crc32.o gzio.o uncompr.o deflate.o trees.o \
           zutil.o inflate.o infback.o inftrees.o inffast.o

SUBDIRS = zlib test tools

all: zlib $(TARGET) test tools

zlib:
	$(MAKE) -C zlib

test:
	$(MAKE) -C test

tools:
	$(MAKE) -C tools

$(TARGET): $(OBJECTS) 
	@echo $(DGGFDFG)
	$(AR) $(ARFLAGS) $@ $^
	for obj in $(ZLIBOBJS) ; do $(AR) $(ARFLAGS) $@ zlib/$$obj; done
	-$(RANLIB) $@
	mkdir -p $(TARGET_DIR)
	cp $(TARGET) $(TARGET_DIR)/$(TARGET)
	
# create shared library
#$(CC) -shared -o $@ -Wl,-soname,$@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) -c $(CFLAGS) $<

tags:
	ctags --recurse=yes *.c *.h $(SUBDIRS)

depend:
	makedepend -- $(CFLAGS) -- *.[ch]*
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done

clean: cleandepend
	rm -f $(TARGET) *.o tags
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done

cleandepend:
	makedepend
	rm -f Makefile.bak
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done

realclean: clean
	rm -rf lib*

.PHONY: clean cleandepend realclean tags zlib test tools depend


# DO NOT DELETE
