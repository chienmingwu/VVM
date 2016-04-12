#-----------------------------------------------------------------------------
# Generic Makefile for models with Fortran source
#
# Requires:
#
#    GNU make		(works with version 3.79 or greater)
#    Makefile.host	(processor-dependent information)
#    Makefile.model	(model-dependent information)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Using this makefile:
#
# To make the model, type "make"
#
# To override a makefile variable, specify it on the command line, e.g.,
# make FFLAGS='-fpic -i8'
#
# To make a debug version of the model, type "make DEBUG=true"
# What a debug version is depends on the options set in Makefile.host; this 
# might include compiler options for profiling, array bounds checking, etc.
#
# To remove all generated files except the executable code, type "make clean"
#
# To remove all generated files, type "make clobber"
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Comments on the construction of this Makefile:
#
# It is named "GNUmakefile" so that it is read first by GNU make.  This will 
# be ignored by other versions of make, so you can write a different makefile 
# (or Makefile) if you're not using GNU make, # and still keep this one around.
# To get GNU make to use "makefile" instead, use "make -f makefile".
#
# GNU make is used since it is widely available (probably for any Unix system 
# and possibly for others).  It's the closest thing to a standard for "make" 
# as is available, and is standard for Linux.  Most variable names, etc. are 
# close to GNU standards.  Some features used here which might not work with 
# other versions of make include:
#
#	VPATH (allows easy access to source in multiple directories)
#	target-specific variable values (allows free/fixed source)
#
# Two other (partial) makefiles are read in:
#
#	MAKEFILE.HOST	for processor-dependent information
#	MAKEFILE.MODEL	for model-dependent information
#
# Their default names are Makefile.host and Makefile.model (in the current
# directory), but you can change these definitions below if you want to.
#
# It is assumed here that all source files are .F files (Fortran source
# with possibly some directives for the C preprocessor).  They will be
# passed through the C preprocessor to turn them into pure Fortran files;
# the suffix used (SUF) is defined in MAKEFILE_HOST since it depends on
# the Fortran compiler being used (some use .f, others use .f90, etc.).
#
# Future improvements (someday?):
#
#	support for C source
#	use automake/autoconfig to generate processor-dependent info?
#
# Author:  Scott R. Fulton
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Initializations
#-----------------------------------------------------------------------------
# specify the shell to avoid surprises in weird environments
SHELL	= /bin/csh

PETSC_DIR        = /opt/petsc360/lib/petsc
include ${PETSC_DIR}/conf/variables
#-----------------------------------------------------------------------------
# Read the universal makefile to define the processor-dependent information
#-----------------------------------------------------------------------------
MAKEFILE_HOST = ../Makefile.host
include $(MAKEFILE_HOST)

#-----------------------------------------------------------------------------
# Read the makefile to define the model-dependent information.
# Must define SRC_FREE and SRC_FIXED as lists of Fortran source (.F) files
# (free source and fixed source, respectively; either may be empty).
# Also must define the default target "all" depending on $(MAIN),
# and should list all nontrivial dependencies (e.g., on modules).
#-----------------------------------------------------------------------------
MAKEFILE_MODEL	= ../Makefile.model
include $(MAKEFILE_MODEL)

#-----------------------------------------------------------------------------
# Create the corresponding lists of preprocessed source files
#-----------------------------------------------------------------------------
SRC_FREE_F	= $(filter %.F,$(SRC_FREE))
SRC_FREE_F90  = $(filter %.F90,$(SRC_FREE))
SRC_FREE_f90  = $(filter %.f90,$(SRC_FREE))
#SRC_FREE_f = $(SRC_FREE_F) $(SRC_FREE_F90) $(SRC_FREE_f90)
#SRC_FIXED_F	= $(SRC_FIXED:.F=.$(SUF))

#-----------------------------------------------------------------------------
# Create the corresponding lists of objects
#-----------------------------------------------------------------------------
OBJS_FREE_F	= $(SRC_FREE_F:.F=.o)
OBJS_FREE_F90    = $(SRC_FREE_F90:.F90=.o)
OBJS_FREE_f90    = $(SRC_FREE_f90:.f90=.o)
OBJS_FIXED	= $(SRC_FIXED:.F=.o)
OBJS		= $(OBJ_FREE*) $(OBJ_FIXED)

#-----------------------------------------------------------------------------
# Target-specific variables:  flags for free source and fixed source
#-----------------------------------------------------------------------------

$(OBJS_FREE_F)	: FORM = $(FREE)
$(OBJS_FREE_F90)	: FORM = $(FREE)
$(OBJS_FREE_f90)	: FORM = $(FREE)
$(OBJS_FIXED)	: FORM = $(FIXED)

#-----------------------------------------------------------------------------
# Replace the built-in implicit rules with GNU-make "static pattern rules".
# This allows us to treat free source and fixed source separately,
# and to save the preprocessed files with suffix $(SUF) for debugging.
# It also overcomes pgf90's apparent refusal to invoke cpp.
#-----------------------------------------------------------------------------
.SUFFIXES:	# delete the default list of suffixes for safety
#.SUFFIXES:	.o .f .F 

#$(SRC_F) : %.$(SUF):  %.F 
#	$(CPP) $(CPPFLAGS) $(CPP_OPTS) $< > $@
#$(SRC_F90) : %.$(SUF):  %.F90 
#	$(CPP) $(CPPFLAGS) $(CPP_OPTS) $< > $@
#$(SRC_f90) : %.$(SUF):  %.f90
#	$(CPP) $(CPPFLAGS) $(CPP_OPTS) $< > $@

$(OBJS_FIXED) : %.o: %.F 
	$(CPP) $(CPPFLAGS) $(CPP_OPTS) $< > $(BINDIR)/$*.$(SUF)
	$(FC) -c $(FFLAGS) $(FORM) $(PRECISION) $(TARGET_ARCH) $(BINDIR)/$*.$(SUF)
$(OBJS_FREE_F) : %.o: %.F 
	$(CPP) $(CPPFLAGS) $(CPP_OPTS) $< > $(BINDIR)/$*.$(SUF)
	$(FC) -c $(FFLAGS) $(FORM) $(PRECISION) $(TARGET_ARCH) $(BINDIR)/$*.$(SUF) 
$(OBJS_FREE_F90) : %.o: %.F90 
	$(CPP) $(CPPFLAGS) $(CPP_OPTS) $< > $(BINDIR)/$*.$(SUF)
	$(FC) -c $(FFLAGS) $(FORM) $(PRECISION) $(TARGET_ARCH) $(BINDIR)/$*.$(SUF) 
$(OBJS_FREE_f90) : %.o: %.f90
	$(CPP) $(CPPFLAGS) $(CPP_OPTS) $< > $(BINDIR)/$*.$(SUF)
	$(FC) -c $(FFLAGS) $(FORM) $(PRECISION) $(TARGET_ARCH) $(BINDIR)/$*.$(SUF) 

#-----------------------------------------------------------------------------
# Define the rule(s) for building the executable code
#-----------------------------------------------------------------------------
$(MAIN): $(OBJS_FREE_F90) $(OBJS_FREE_f90) $(OBJS_FREE_F) $(OBJS_FIXED)
	$(LD) $(LDFLAGS) $^ -o $@ $(LDLIBS) 

#-----------------------------------------------------------------------------
# Target to check that the code still produces the same results.
# Results from standard test are stored in the file $(CHECK_OUTPUT)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# standard targets for miscellaneous chores
#-----------------------------------------------------------------------------

clean :
	rm -f *.o $(BINDIR)/*.$(SUF) *.h *.mod

clobber :
	rm -f *.o $(BINDIR)/*.$(SUF) *.h *.mod $(MAIN)
#	echo $(OBJS)
topography :
	ifort -c ../../CODE/RRTMG/kinds.F -free 
	${CPP} -p $(BINDIR)/PARMSLD.f90 > PARMSLD.f	
	ifort -c PARMSLD.f -free
	${CPP} -p nc2fortranBinary.F > nc2fortranBinary.f
	ifort -c nc2fortranBinary.f  -I/opt/intel-libs/include
	ifort -o $(BINDIR)/nc2fortranBinary nc2fortranBinary.o PARMSLD.o kinds.o -L/opt/intel-libs/lib -lnetcdf -lnetcdff 
#	$(CPP) -p nc2fortranBinary.F > nc2fortranBinary.f
#	ifort -c nc2fortranBinary.f -L/opt/netcdf4-intel/lib -lnetcdf -lnetcdff


