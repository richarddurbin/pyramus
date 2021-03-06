# makefile for mover, developed on Richard's Mac.

#CFLAGS= -O3
CFLAGS= -g -Wall			# for debugging
#CFLAGS= -03 -DOMP -fopenmp		# for OMP parallelisation - doesn't compile on Mac

all: pyramus composition

clean:
	$(RM) *.o *~ pyramus
	$(RM) -r *.dSYM

### object files

UTILS_OBJS=hash.o dict.o array.o utils.o
UTILS_HEADERS=utils.h array.h dict.h hash.h
$(UTILS_OBJS): utils.h $(UTILS_HEADERS)

ONE_DIR = ../vgp-tools/Core
HTS_DIR = $(PWD)/../htslib
SEQIO_OPTS = -DONEIO -I$(ONE_DIR) -DBAMIO -I$(HTS_DIR)/htslib/
SEQIO_LIBS = -L$(ONE_DIR) -lONE -L$(HTS_DIR) -Wl,-rpath $(HTS_DIR) -lhts -lm -lbz2 -llzma -lcurl -lz 
# the "-Wl,-rpath $(HTS_DIR)" incantation is needed for local dynamic linking if htslib is not installed centrally

seqio.o: seqio.c seqio.h 
	$(CC) $(CFLAGS) $(SEQIO_OPTS) -c $^

### programs

pyramus: pyramus.c seqio.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(SEQIO_LIBS)

composition: composition.c seqio.o utils.o
	$(CC) $(CFLAGS) $^ -o $@ $(SEQIO_LIBS)


### end of file
