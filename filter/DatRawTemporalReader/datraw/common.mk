
OPTIMIZATIONS = ON
SHELL = /bin/bash
#OPTIMIZATIONS = OFF

# determine presence of Intel C/C++ compiler
ifeq ($(shell which icc >& /dev/null && echo exists),exists)
	CC = icc
else
	CC = gcc
endif

AR = ar
ARFLAGS = rcs
RANLIB = ranlib

CFLAGS = -Wall

ifeq ($(CC),icc)
	# icc options
	LDFLAGS = -static-intel
	DBGCFLAGS = -ggdb -O0
	OPTCFLAGS = -O3 -no-prec-div -finline-functions -unroll -parallel \
	            -par-runtime-control
	#OPTCFLAGS += -par-report3
	OPTLDFLAGS = -parallel
else
	# gcc options
	CFLAGS += -pedantic
	DBGCFLAGS = -ggdb -O0
	OPTCFLAGS = -O3 -finline-functions -funroll-loops -fprefetch-loop-arrays -ffast-math
	OPTLDFLAGS =
endif

# architecture specific options
ARCH := $(shell uname -m | sed 's/i.86/i386/')
ifeq ($(ARCH),x86_64)
	TARGET_DIR = lib64
	ifeq ($(CC),icc)
		OPTCFLAGS  += -xP
	else
		OPTCFLAGS  += -march=k8
	endif
else 
	ifeq ($(ARCH),i386)
		TARGET_DIR = lib32
		ifeq ($(CC),icc)
			OPTCFLAGS  += -xW
		else
			OPTCFLAGS  += -march=pentium4
			#OPTCFLAGS  += -march=athlon
		endif
	else
		TARGET_DIR = lib_$(ARCH)
		OPTCFLAGS  += -O2
	endif
endif

ifeq ($(OPTIMIZATIONS),ON)
	CFLAGS += $(OPTCFLAGS)
	LDFLAGS += $(OPTLDFLAGS)
else
	CFLAGS += $(DBGCFLAGS)
endif

