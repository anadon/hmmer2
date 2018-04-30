CC        = gcc
CFLAGS    = -march=native -O1 -g -pipe --param=ssp-buffer-size=4 -D_FORTIFY_SOURCE=2
LIBS      =  -lm -pthread
MYLIBS    = -lhmmer$(PROGSUFFIX) -lsquid

SHELL     = /bin/bash

SHIVA = alignalign_test\
	evd_test\
	masks_test\
	parsingviterbi_test\
	tophits_test\
	trace_test\
	viterbi_exercise\
	weeviterbi_test

#######
## Targets defining how to make Shiva executables.
#######

.c.o:
	$(CC) $(CFLAGS) $(PTHREAD_CFLAGS) $(DEFS) -I../squid -I../src -c $<

all: 	$(SHIVA)

$(SHIVA): %: %.o
	$(CC) $(CFLAGS) $(PTHREAD_CFLAGS) $(DEFS) -o $@ -I../squid -I../src -L../squid -L../src $@.o $(MYLIBS) $(PTHREAD_LIBS) $(LIBS)

#######
## `make check` actually runs the tests.
#######

check:
	@echo
	@echo Running test suite exercises.
	@echo Warning: some tests may take several minutes to complete.
	@echo
	LD_LIBRARY_PATH=LD_LIBRARY_PATH:../src/ ./sqc 2 exercises.sqc . ../src


#######
## Miscellaneous
#######

clean:
	-rm -f *.o *~ Makefile.bak core $(SHIVA) TAGS gmon.out

TAGS:
	etags -t *.c *.h Makefile.in
