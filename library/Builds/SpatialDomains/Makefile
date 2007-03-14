ARCH	 = $(shell uname)
VPATH	 = ../../SpatialDomains ../../LibUtilities/Polylib
LIB	 = libSpatialDomains.a

# Must provide (relative) location or root director of libraries,i.e
# typiecally Nektar++/libs
NEK_LIB_ROOT = ../../

# Compilers and options

CPPFLAGS = 
include ../Flags/${ARCH}.inc

# These are the module names:

ALL = GeoFac Geometry Geometry1D EdgeComponent	\
MeshComponents MeshGraph MeshGraph1D 		\
InterfaceComponent SegGeom 

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
	cp -f libSpatialDomains.a           ../../; 
ifeq ($(ARCH),Darwin)
	ranlib ../../libSpatialDomains.a;
endif

help:
	@echo ''
	@echo ' Run make :'	
	@echo ' To compiled debuged version type "Make DEBUG=1"'
	@echo ''

clean:
	-rm -f *.o *~ *.a


