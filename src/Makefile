CC     = gcc
CFLAGS = -pipe -Wall -Wextra -fPIC -D_FORTIFY_SOURCE=2 -O1 -g
LIBS   = -lm -pthread

## archiving command, and ranlib command if you need one.
#	In general, you shouldn't need to change these, and they're
#	only used for building the testsuite anyway... e.g. we
#	make a "libhmmer.a" library for building the testsuite.
#
#AR		 = /usr/bin/ar rcv
#RANLIB = ranlib


MYLIBS	 = -lsquid
MYLIBDIR = -L../squid
MYINCDIR = -I../squid

HMMLIB = libhmmer$(PROGSUFFIX).so


OBJS =	alphabet.o\
	core_algorithms.o\
	debug.o\
	display.o\
	emit.o\
	emulation.o\
	fast_algorithms.o\
	histogram.o\
	hmmio.o\
	infocontent.o\
	lsj_eweight.o\
	mathsupport.o\
	masks.o\
	matrices.o\
	misc.o\
	modelmakers.o\
	plan7.o\
	plan9.o\
	postprob.o\
	prior.o\
	threads.o\
	tophits.o\
	trace.o

HDRS =	config.h\
	funcs.h\
	globals.h\
	lsjfuncs.h\
	postprob.h\
	structs.h

.c.o:
	${CC} ${CFLAGS} ${PTHREAD_CFLAGS} ${MYINCDIR} ${LIBS} -c $<

#################################################################
## Targets defining how to make HMMER executables.
##
all:	 $(PROGS)

${PROGS}: %: %.o ${OBJS}
	${CC} ${CFLAGS} ${MYLIBDIR} -o $@ $@.o ${OBJS} ${MYLIBS} ${PTHREAD_LIBS} ${LIBS}


#################################################################
## Targets used in making HMMER module for testsuite compilation.
##
module: $(HMMLIB)

$(HMMLIB): $(OBJS)
#	$(AR) libhmmer.a $(OBJS)
#	$(RANLIB) libhmmer.a
#	chmod 644 libhmmer.a
	gcc -shared $(OBJS) -o $@


clean:
	rm -f *.o *~ Makefile.bak core $(PROGS) TAGS gmon.out *.so
	rm -f hmmalign hmmbuild hmmcalibrate hmmconvert hmmemit hmmfetch hmmindex hmmpfam hmmsearch
