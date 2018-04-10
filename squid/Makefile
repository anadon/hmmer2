###############################################################
# Makefile for SQUID library
# SVN $Id: Makefile.in 1721 2006-11-29 21:51:21Z eddys $
#
# Note: The autoconf variables in this file must be coordinated
#       with HMMER. HMMER creates a Makefile from this 
#       Makefile.in using its own configure script, not SQUID's.
#################################################################
# HMMER - Biological sequence analysis with profile HMMs
# Copyright (C) 1992-2006 HHMI Janelia Farm
# All Rights Reserved
# 
#     This source code is distributed under the terms of the
#     GNU General Public License. See the files COPYING and LICENSE
#     for details.
#################################################################


### Installation points
###
prefix      = /usr
exec_prefix = ${prefix}
BINDIR      = ${exec_prefix}/bin
MANDIR      = ${datarootdir}/man
INCLUDEDIR  = ${prefix}/include
LIBDIR      = ${exec_prefix}/lib
SCRIPTDIR   = ${exec_prefix}/bin
DATAROOTDIR = ${prefix}/share

## your compiler and compiler flags
#
CC       = gcc
CFLAGS   = -march=native -Ofast -pipe --param=ssp-buffer-size=4
CPPFLAGS = -D_FORTIFY_SOURCE=2
LDFLAGS  = -Wl,-O2,--sort-common,--as-needed,-z,relro
DEFS     = -DHAVE_CONFIG_H
LIBS     =  -lm

## Archiver command
#
AR     = /usr/bin/ar rcv
RANLIB = ranlib

## instructions for installing man pages
#
INSTMAN   = cp
MANSUFFIX = 1

#######
## You should not need to modify below this line
#######
SHELL       = /bin/sh
BASENAME    = squid
PACKAGE     = SQUID
RELEASE     = 1.9g

PROGS =	afetch\
	alistat\
	compalign\
	compstruct\
	revcomp\
	seqsplit\
	seqstat\
	sfetch\
	shuffle\
	sindex\
	sreformat\
	translate\
	weight

MANS =  afetch\
	alistat\
	seqstat\
	sfetch\
	shuffle\
	sreformat\

READMES = README INSTALL Makefile.in 

SCRIPTS = 

PRECONFHDRS = \
	squid.h.in\
	squidconf.h.in

POSTCONFHDRS = \
	squid.h\
	squidconf.h

HDRS =	dirichlet.h\
	rk.h\
	sqfuncs.h\
	gki.h\
	msa.h\
	sre_random.h\
	sre_stack.h\
	ssi.h\
	stopwatch.h\
	vectorops.h

OBJS =	a2m.o\
	aligneval.o\
	alignio.o\
	clustal.o\
	cluster.o\
	dayhoff.o\
	dirichlet.o\
	eps.o\
	file.o\
	getopt.o\
	gki.o\
	hsregex.o\
	iupac.o\
	msa.o\
	msf.o\
	phylip.o\
	revcomp.o\
	rk.o\
	selex.o\
	seqencode.o\
	shuffle.o\
	sqerror.o\
	sqio.o\
	squidcore.o\
	sre_ctype.o\
	sre_math.o\
	sre_random.o\
	sre_stack.o\
	sre_string.o\
	ssi.o\
	stockholm.o\
	stopwatch.o\
	translate.o\
	types.o\
	vectorops.o\
	weight.o\
	wuss.o


################################################################
# Targets that actually build the squid executables
all: $(PROGS) libsquid.a

$(PROGS): %: %_main.o ${OBJS} 
	${CC} ${CFLAGS} ${LDFLAGS} ${PVMLIBDIR} ${DEFS} -o $@ $@_main.o ${OBJS} ${PVMLIBS} ${LIBS}

.c.o:
	${CC} ${CFLAGS} ${CPPFLAGS} ${PVMINCDIR} ${DEFS} -c $<		

################################################################


################################################################
# Targets expected by packages (e.g. HMMER) that
# include SQUID as a module. 
#
module: libsquid.a

libsquid.a: $(OBJS)
	$(AR) libsquid.a $(OBJS)
	$(RANLIB) libsquid.a
	chmod 644 libsquid.a
#################################################################


install:  $(PROGS) libsquid.a
	test -d $(LIBDIR)    || mkdir -p $(LIBDIR)
	test -d $(BINDIR)    || mkdir -p $(BINDIR)
	test -d $(SCRIPTDIR) || mkdir -p $(SCRIPTDIR)
	test -d $(INCLUDEDIR)|| mkdir -p $(INCLUDEDIR)
	test -d $(MANDIR)/man$(MANSUFFIX) || mkdir -p $(MANDIR)/man$(MANSUFFIX)
	cp libsquid.a $(LIBDIR)/
	cp $(HDRS)  $(INCLUDEDIR)/
	cp $(POSTCONFHDRS)  $(INCLUDEDIR)/
	cp $(PROGS) $(BINDIR)/
#	for scriptfile in $(SCRIPTS); do\
#	  cp Scripts/$$scriptfile $(SCRIPTDIR)/;\
#	done
	@for manpage in $(MANS); do\
	  $(INSTMAN) Man/$$manpage.man $(MANDIR)/man$(MANSUFFIX)/$$manpage.$(MANSUFFIX);\
	done


check:	libsquid.a 
	(cd testsuite; make CC="$(CC)" CFLAGS="$(CFLAGS)")
	(cd testsuite; make check)


clean:
	-rm -f *.o *~ core TAGS llib-lsquid.ln ccmalloc.log $(PROGS) 
	-rm -f libsquid.a
	(cd testsuite; make clean)


TAGS:
	etags -t *.h *.c Makefile.in


