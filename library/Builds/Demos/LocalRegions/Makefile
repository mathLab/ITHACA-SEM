ARCH	 = $(shell uname)
VPATH	 = ../../../Demos/LocalRegions

# Compilers and options

CPPFLAGS = -I../../../

include ../../Flags/${ARCH}.inc

LDLIBS	 = -L../../../  -lLocalRegions -lSpatialDomains -lStdRegions -lNekUtilities -lPolylib -lm  $(LIBS)


default:
# These are the module names:
	PROG = $<
	$(CXX) $(CPPFLAGS) -c $(VPATH)/$(PROG)
	$(CXX) $(CPPFLAGS) -o $(PROG)  $(PROG).o -L../../../ -lSpatialDomains  -lLocalRegions  -lStdRegions -lPolylib -lm  $(LIBS)




help:
	@echo ''
	@echo ' To run: make [DEBUG=1] "Program name"'	
	@echo ''

clean:
	-rm -f *.o *~


