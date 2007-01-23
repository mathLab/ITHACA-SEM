ARCH	 = $(shell uname)
VPATH	 = ../../StdRegions 
LIB	 = libStdRegions.a

                       
# Must provide (relative) location or root director of libraries,i.e
# typiecally Nektar++/libs
NEK_LIB_ROOT = ../../

# Compilers and options

CPPFLAGS = 
include ../Flags/${ARCH}.inc

# These are the module names:

ALL   = StdExpansion StdExpansion1D StdSegExp

# ALL  =  SStdExpansion2D StdExpansion3D StdExpMap StdHexExp    \
# 	  StdNodalTriExp  StdNodalTetExp  StdPrismExp StdPyrExp \
#	  StdSegExp StdTetExp  StdTriExp StdQuadExp 

# ---------------------------------------------------------------- #
# The following expand to give the names of all sections (ALL),    #
# all library member dependencies (LMOD), all associated           #
# object files (OBJS), and all source code files (SRC).            #
# ---------------------------------------------------------------- #

LMOD  = $(foreach routine, $(ALL), $(LIB)($(routine).o))
OBJS  = $(foreach routine, $(ALL), $(routine).o)

# ----------------------------------------------------------------------------
# Make library, default action.
#
$(LIB): $(LMOD)

ifeq ($(ARCH),Darwin)
	ranlib $(LIB)
endif

headers:
	make -f ../LibUtilities/Makefile.defheaders

# ----------------------------------------------------------------------------
# Install library & associated header files.
#
install:
	cp -f libStdRegions.a            ../../; 

ifeq ($(ARCH),Darwin)
	ranlib ../../libStdRegions.a;
endif

help:
	@echo ''
	@echo ' Run make :'	
	@echo ''

clean:
	-rm -f *.o *~ *.a


