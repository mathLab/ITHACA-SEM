ARCH	 = $(shell uname)
VPATH	 = ../../LibUtilities/BasicUtils ../../LibUtilities/Foundations ../../LibUtilities/Polylib ../../LibUtilities/LinearAlgebra
LIB	 = libNekUtilities.a

# Must provide (relative) location or root director of libraries,i.e
# typiecally Nektar++/libs
NEK_LIB_ROOT = ../../

# Compilers and options

CPPFLAGS = 
include ../Flags/${ARCH}.inc

# These are the module names:

ALL = Basis FourierPoints GaussPoints Graph ManagerAccess NodalTetElec NodalTriElec NodalTriFekete  Points PolyEPoints Polylib lapack

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
	cp -f libNekUtilities.a            ../../; 
ifeq ($(ARCH),Darwin)
	ranlib ../../libNekUtilities.a;
endif

help:
	@echo ${ALL}
	@echo ${LMOD}

	@echo ''
	@echo ' Run make :'	
	@echo ''

clean:
	-rm -f *.o *~ *.a


