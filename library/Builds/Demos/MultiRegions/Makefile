ARCH	 = $(shell uname)
VPATH	 = ../../../Demos/MultiRegions

# Must provide (relative) location or root director of libraries,i.e
# typiecally Nektar++/libs

NEK_LIB_ROOT = ../../../

# Compilers and options
CPPFLAGS = 
include ../../Flags/${ARCH}.inc

LDLIBS	 = -L../../../  -lMultiRegions -lLocalRegions -lSpatialDomains -lStdRegions -lNekUtilities -lPolylib -lm  $(LIBS)


default:
# These are the module names:
	PROG = $<
	$(CXX) $(CPPFLAGS) -c $(VPATH)/$(PROG)
	$(CXX) $(CPPFLAGS) -o $(PROG)  $(PROG).o -L../../../ -LMultiRegions -lSpatialDomains  -lLocalRegions  -lStdRegions -lPolylib -lm  $(LIBS)




help:
	@echo ''
	@echo ' To run: make [DEBUG=1] "Program name"'	
	@echo ''

clean:
	-rm -f *.o *~


