# --- COMPILER ----------------------------------------
CC = mpicc -std=gnu99 -Wall -pedantic
CPP = cpp
MAKEDEP = $(CPP) -MM

# --- DO NOT CHANGE -----------------------------------
SRCDIR = src
BUILDDIR = build
GSRCDIR = $(BUILDDIR)/gsrc
SRC = $(patsubst $(SRCDIR)/%,%,$(filter-out %_generic.c,$(wildcard $(SRCDIR)/*.c)))
SRCGEN = $(patsubst $(SRCDIR)/%,%,$(wildcard $(SRCDIR)/*_generic.c))
GSRCFLT = $(patsubst %_generic.c,$(GSRCDIR)/%_float.c,$(SRCGEN))
GSRCDBL = $(patsubst %_generic.c,$(GSRCDIR)/%_double.c,$(SRCGEN))
GSRC = $(patsubst %,$(GSRCDIR)/%,$(SRC)) $(GSRCFLT) $(GSRCDBL)
HEA = $(patsubst $(SRCDIR)/%,%,$(filter-out %_generic.h,$(wildcard $(SRCDIR)/*.h)))
HEAGEN = $(patsubst $(SRCDIR)/%,%,$(wildcard $(SRCDIR)/*_generic.h))
GHEAFLT = $(patsubst %_generic.h,$(GSRCDIR)/%_float.h,$(HEAGEN))
GHEADBL = $(patsubst %_generic.h,$(GSRCDIR)/%_double.h,$(HEAGEN))
GHEA = $(patsubst %,$(GSRCDIR)/%,$(HEA)) $(GHEAFLT) $(GHEADBL)
OBJ = $(patsubst $(GSRCDIR)/%.c,$(BUILDDIR)/%.o,$(GSRC))
OBJDB = $(patsubst %.o,%_db.o,$(OBJ))
DEP = $(patsubst %.c,%.dep,$(GSRC))

# --- FLAGS -------------------------------------------
OPT_FLAGS = -fopenmp -DOPENMP -DSSE -msse4.2 
CFLAGS = -DPARAMOUTPUT -DTRACK_RES -DFGMRES_RESTEST -DPROFILING
# -DSINGLE_ALLREDUCE_ARNOLDI
# -DCOARSE_RES -DSCHWARZ_RES -DTESTVECTOR_ANALYSIS
OPT_VERSION_FLAGS =$(OPT_FLAGS) -O3 -ffast-math
DEBUG_VERSION_FLAGS = $(OPT_FLAGS)

# --- FLAGS FOR HDF5 ---------------------------------
# H5HEADERS=-DHAVE_HDF5 /usr/include
# H5LIB=-lhdf5 -lz

# --- FLAGS FOR LIME ---------------------------------
# LIMEH=-DHAVE_LIME -I$(LIMEDIR)/include
# LIMELIB= -L$(LIMEDIR)/lib -llime

all: wilson library documentation
wilson: dd_alpha_amg dd_alpha_amg_db
library: lib/libdd_alpha_amg.a include/dd_alpha_amg_parameters.h include/dd_alpha_amg.h
documentation: doc/user_doc.pdf

.PHONY: all wilson library
.SUFFIXES:
.SECONDARY:

dd_alpha_amg : $(OBJ)
	$(CC) $(OPT_VERSION_FLAGS) $(LIMEH) -o $@ $(OBJ) $(H5LIB) $(LIMELIB) -lm

dd_alpha_amg_db : $(OBJDB)
	$(CC) -g $(DEBUG_VERSION_FLAGS) $(LIMEH) -o $@ $(OBJDB) $(H5LIB) $(LIMELIB) -lm

lib/libdd_alpha_amg.a: $(OBJ)
	ar rc $@ $(OBJ)
	ar d $@ main.o
	ranlib $@

doc/user_doc.pdf: doc/user_doc.tex doc/user_doc.bib
	( cd doc; pdflatex user_doc; bibtex user_doc; pdflatex user_doc; pdflatex user_doc; )

include/dd_alpha_amg.h: src/dd_alpha_amg.h
	cp src/dd_alpha_amg.h $@

include/dd_alpha_amg_parameters.h: src/dd_alpha_amg_parameters.h
	cp src/dd_alpha_amg_parameters.h $@

$(BUILDDIR)/%.o: $(GSRCDIR)/%.c $(SRCDIR)/*.h
	$(CC) $(CFLAGS) $(OPT_VERSION_FLAGS) $(H5HEADERS) $(LIMEH) -c $< -o $@

$(BUILDDIR)/%_db.o: $(GSRCDIR)/%.c $(SRCDIR)/*.h
	$(CC) -g $(CFLAGS) $(DEBUG_VERSION_FLAGS) $(H5HEADERS) $(LIMEH) -DDEBUG -c $< -o $@

$(GSRCDIR)/%.h: $(SRCDIR)/%.h $(firstword $(MAKEFILE_LIST))
	cp $< $@

$(GSRCDIR)/%_float.h: $(SRCDIR)/%_generic.h $(firstword $(MAKEFILE_LIST))
	sed -f float.sed $< > $@

$(GSRCDIR)/%_double.h: $(SRCDIR)/%_generic.h $(firstword $(MAKEFILE_LIST))
	sed -f double.sed $< > $@

$(GSRCDIR)/%.c: $(SRCDIR)/%.c $(firstword $(MAKEFILE_LIST))
	cp $< $@

$(GSRCDIR)/%_float.c: $(SRCDIR)/%_generic.c $(firstword $(MAKEFILE_LIST))
	sed -f float.sed $< > $@

$(GSRCDIR)/%_double.c: $(SRCDIR)/%_generic.c $(firstword $(MAKEFILE_LIST))
	sed -f double.sed $< > $@

%.dep: %.c $(GHEA)
	$(MAKEDEP) $< | sed 's,\(.*\)\.o[ :]*,$(BUILDDIR)/\1.o $@ : ,g' > $@
	$(MAKEDEP) $< | sed 's,\(.*\)\.o[ :]*,$(BUILDDIR)/\1_db.o $@ : ,g' >> $@
clean:
	rm -f $(BUILDDIR)/*.o
	rm -f $(GSRCDIR)/*
	rm -f dd_alpha_amg
	rm -f dd_alpha_amg_db

-include $(DEP)
