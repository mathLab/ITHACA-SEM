# define default compilers 
COMPILER_PATH	= /usr/bin
CXX	= $(COMPILER_PATH)/CC
CC	= $(COMPILER_PATH)/cc


# define MPI compilers if flag MPI_PROG is activated.
ifdef MPI_PROG
	MPICOMPILER_PATH = 
	CXX	= $(MPICOMPILER_PATH)/mpiCC
	CC	= $(MPICOMPILER_PATH)/mpicc
	FC	= $(MPICOMPILER_PATH)/mpif77
endif

LD	= $(CXX)

#define Third party include directories
BOOST_VERSION = boost_1_33_0
LOKI_VERSION = loki-0.1.3

  CPPFLAGS += -I${NEK_LIB_ROOT}					\
  -I${NEK_LIB_ROOT}/../ThirdParty/${LOKI_VERSION}/include	\
  -I${NEK_LIB_ROOT}/../ThirdParty/${BOOST_VERSION}              \
  -I${NEK_LIB_ROOT}/../ThirdParty/

ifdef FULLDEBUG
	CPPFLAGS += -g -trapuv -fullwarn -LANG:std -DDEBUG -DFULLDEBUG	
else
ifdef DEBUG
	CPPFLAGS += -g -trapuv -fullwarn -LANG:std -DDEBUG
else
	CPPFLAGS += -O3 -LANG:std
endif
endif

LIBS	= -lcomplib.sgimath -lftn -lm
