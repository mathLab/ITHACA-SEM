ARCH	 = $(shell uname)
VPATH	 = ../../../Demos/LibUtilities

CPPFLAGS = -I../../../ \
           -I../../../../../ThirdParty/met \
           -I../../../../../ThirdParty/met/common \
           -I../../../../../ThirdParty/boost_1_33_0


include ../../Flags/${ARCH}.inc

LDLIBS	 =  -lm  -L../../../ -lNekUtilities $(LIBS) -L../../../../ThirdParty/boost_1_33_0/stage/lib -lboost_thread-gcc-1_33


default:
# These are the module names:
	PROG = $<
	$(CXX) $(CPPFLAGS) -c $(VPATH)/$(PROG)
	$(CXX) $(CPPFLAGS) -o $(PROG)  $(PROG).o  -lm  $(LIBS)


help:
	@echo ''
	@echo ' To run: make [DEBUG=1] "Program name"'	
	@echo ''

clean:
	-rm -f *.o *~


