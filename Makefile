# makefile for mover, developed on Richard's Mac.

#CFLAGS= -O3
CFLAGS= -g -Wall			# for debugging
#CFLAGS= -03 -DOMP -fopenmp		# for OMP parallelisation - doesn't compile on Mac

all: pyramus composition

clean:
	$(RM) *.o *~ pyramus
	$(RM) -r *.dSYM

### object files

UTILS_OBJS=hash.o dict.o array.o # utils.o 	# must include utils.o if exclude VGP below
UTILS_HEADERS=utils.h array.h dict.h hash.h
$(UTILS_OBJS): utils.h $(UTILS_HEADERS)

VGP_DIR = ../vgp-tools
BAM_DIR = ../htslib
SEQIO_OPTS = -DVGPIO -I$(VGP_DIR)/include -DBAMIO -I$(BAM_DIR)/htslib/
SEQIO_LIBS = -L$(VGP_DIR)/lib -lVGP $(BAM_DIR)/libhts.a -lz -lm -lbz2 -llzma -lcurl -lz 

seqio.o: seqio.c seqio.h 
	$(CC) $(SEQIO_OPTS) -c $^

### programs

pyramus: pyramus.c seqio.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(SEQIO_LIBS)

composition: composition.c seqio.o
	$(CC) $(CFLAGS) $^ -o $@ $(SEQIO_LIBS)


### end of file
