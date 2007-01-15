ARCH	 = $(shell uname)
VPATH	 = ../../MultiRegions ../../Polylib
LIB	 = libMultiRegions.a

# Must provide (relative) location or root director of libraries,i.e
# typiecally Nektar++/libs
NEK_LIB_ROOT = ../../

# Compilers and options
CPPFLAGS = 
include ../Flags/${ARCH}.inc

# These are the module names:

ALL   =  ContExpList1D ContExpList2D  DisContExpList1D ExpList ExpList1D ExpList2D ExpList3D Field Field1D Field2D Field3D LocalToGlobalMap  LocalToGlobalMap1D LocalToGlobalMap2D LocalToGlobalMap3D

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

# ----------------------------------------------------------------------------
# Install library & associated header files.
#
install:
	cp  -f libMultiRegions.a           ../../; 
ifeq ($(ARCH),Darwin)
	ranlib ../../libMultiRegions.a;
endif

help:
	@echo ''
	@echo ' Run make :'	
	@echo ' To compiled debuged version type "Make DEBUG=1"'
	@echo ''

clean:
	-rm -f *.o *~ *.a


