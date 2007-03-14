ARCH	 = $(shell uname)
VPATH	 = ../../MultiRegions ../../Polylib
LIB	 = libMultiRegions.a

# Must provide (relative) location or root director of libraries,i.e
# typiecally Nektar++/libs
NEK_LIB_ROOT = ../../

# Compilers and options
CPPFLAGS = -DNEKMAT
include ../Flags/${ARCH}.inc

# These are the module names:

ALL   =  ExpList ExpList1D 

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


