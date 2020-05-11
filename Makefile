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
BAM_DIR = ../htslib
SEQIO_OPTS = -DONEIO -I$(ONE_DIR) -DBAMIO -I$(BAM_DIR)/htslib/
SEQIO_LIBS = -L$(ONE_DIR) -lONE $(BAM_DIR)/libhts.a -lm -lbz2 -llzma -lcurl -lz 

seqio.o: seqio.c seqio.h 
	$(CC) $(SEQIO_OPTS) -c $^

### programs

pyramus: pyramus.c seqio.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(SEQIO_LIBS)

composition: composition.c seqio.o utils.o
	$(CC) $(CFLAGS) $^ -o $@ $(SEQIO_LIBS)


### end of file
