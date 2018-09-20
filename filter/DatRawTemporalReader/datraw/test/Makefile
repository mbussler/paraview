include ../common.mk
CFLAGS += -I..
LDFLAGS  += -L../$(TARGET_DIR) 
LIBS = -ldatRaw -lm

all: testMulti createMulti createMultiWithOptional testHalf

testMulti: testMulti.o
	 $(CC) -o $@ $(LDFLAGS) $^ $(LIBS)

createMulti: createMulti.o 
	 $(CC) -o $@ $(LDFLAGS) $^ $(LIBS)

createMultiWithOptional: createMultiWithOptional.o 
	 $(CC) -o $@ $(LDFLAGS) $^ $(LIBS)

testHalf: testHalf.o 
	 $(CC) -o $@ $(LDFLAGS) $^ $(LIBS)

.c.o:
	$(CC) -c $(CFLAGS) $<

depend:
	makedepend -- $(CFLAGS) -- *.[ch]*

cleandepend:
	makedepend
	rm -f Makefile.bak
	
clean: cleandepend
	rm -f createMulti createMultiWithOptional testMulti testHalf *.o core test*.dat test*.raw

.PHONY: clean cleandepend depend


# DO NOT DELETE
