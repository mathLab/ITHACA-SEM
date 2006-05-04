ARCH	 = $(shell uname)
VPATH	 = ../../LibUtilities/Polylib/
LIB	 = libPolylib.a
INCDIR	 = ../../LibUtilities/Polylib

# Compilers and options

CPPFLAGS = -I$(INCDIR) 

ifeq ($(ARCH),Linux)			# Linux ==  Generic GCC platform
include ../Flags/Linux.inc
else                                    # Default look for file with uname
include ../Flags/${ARCH}.inc
endif

# These are the module names:

ALL   = Polylib

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
	ranlib $(LIB)

# ----------------------------------------------------------------------------
# Install library & associated header files.
#
install:
	mv -f libPolylib.a   ../../; 
help:
	@echo ''
	@echo ' Run make :'	
	@echo ''

clean:
	-rm -f *.o *~


