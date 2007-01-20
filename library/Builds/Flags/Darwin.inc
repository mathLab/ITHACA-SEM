# define default compilers 
COMPILER_PATH	= /usr/bin
CXX	= $(COMPILER_PATH)/g++
CC	= $(COMPILER_PATH)/gcc 
FC	= $(COMPILER_PATH)/g77


# define MPI compilers if flag MPI_PROG is activated.
ifdef MPI_PROG
	MPICOMPILER_PATH = /usr/local/bin
	CXX	= $(MPICOMPILER_PATH)/mpiCC
	CC	= $(MPICOMPILER_PATH)/mpicc
	FC	= $(MPICOMPILER_PATH)/mpif77
endif

LD	= $(CXX)

#turn of long-double warning associated with singleton
  CPPFLAGS += -Wno-long-double
  CPPFLAGS += -Wno-deprecated
  CPPFLAGS += -DNEKTAR_USE_EXPRESSION_TEMPLATES  -DNEKTAR_USING_LAPACK

#define Third party include directories
BOOST_VERSION = boost_1_33_1
LOKI_VERSION = loki-0.1.3

  CPPFLAGS += -I${NEK_LIB_ROOT}					\
  -I${NEK_LIB_ROOT}/../ThirdParty/${LOKI_VERSION}/include	\
  -I${NEK_LIB_ROOT}/../ThirdParty/${BOOST_VERSION}              \
  -I${NEK_LIB_ROOT}/../ThirdParty/

ifdef FULLDEBUG
	CPPFLAGS += -g -gstabs+ -DDEBUG -DFULLDEBUG
else
ifdef DEBUG
	CPPFLAGS += -g -gstabs+ -DEBUG
else
	CPPFLAGS += -O3 $(MCPU) -ffast-math 
endif
endif

#define location and link for lapack and blas 
LIBS	= -framework Accelerate  
