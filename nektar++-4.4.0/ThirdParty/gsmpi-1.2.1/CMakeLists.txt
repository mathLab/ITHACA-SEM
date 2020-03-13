CMAKE_MINIMUM_REQUIRED(VERSION 2.8.7)

# compile in Release mode by default
SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")

PROJECT(GSMPI)

# First check to see if our compiler has MPI built in to avoid linking libraries etc.
INCLUDE (CheckIncludeFiles)
INCLUDE (CheckFunctionExists)
CHECK_INCLUDE_FILES  (mpi.h    HAVE_MPI_H)
CHECK_FUNCTION_EXISTS(MPI_Send HAVE_MPI_SEND)

IF (NOT "${HAVE_MPI_H}" OR NOT "${HAVE_MPI_SEND}")
  FIND_PACKAGE(MPI REQUIRED)
  IF (MPI_CXX_FOUND)
    INCLUDE_DIRECTORIES( ${MPI_CXX_INCLUDE_PATH} )
    MESSAGE(STATUS "Found MPI: ${MPI_CXX_LIBRARIES}")
  ELSE()
    MESSAGE(FATAL_ERROR "Unable to find MPI installation")
  ENDIF()
ENDIF()

ADD_DEFINITIONS(-DMPI -DPREFIX=nektar_ -DNO_NEK_EXITT -DGLOBAL_LONG)
ADD_DEFINITIONS(-DPRINT_MALLOCS=0 -DUSE_NAIVE_BLAS)

SET(GS_SOURCES
    src/gs.c
    src/sort.c
    src/sarray_transfer.c
    src/sarray_sort.c
    src/gs_local.c
    src/fail.c
    src/crystal.c
    src/comm.c
    src/tensor.c
)

SET(GS_HEADERS
    src/gs.h
)

SET(XXT_SOURCES
    src/xxt.c
    src/sparse_cholesky.c
)

SET(XXT_HEADERS
    src/sparse_cholesky.c
)

ADD_LIBRARY(gsmpi STATIC ${GS_SOURCES})
ADD_LIBRARY(xxt   STATIC ${XXT_SOURCES})

IF (NOT "${HAVE_MPI_H}" OR NOT "${HAVE_MPI_SEND}")
  TARGET_LINK_LIBRARIES(gsmpi ${MPI_CXX_LIBRARIES})
  TARGET_LINK_LIBRARIES(xxt   ${MPI_CXX_LIBRARIES})
  SET_PROPERTY(TARGET gsmpi APPEND PROPERTY LINK_FLAGS
      ${MPI_CXX_LINK_FLAGS})
  SET_PROPERTY(TARGET xxt APPEND PROPERTY LINK_FLAGS
      ${MPI_CXX_LINK_FLAGS})
ENDIF ()

SET_PROPERTY(TARGET gsmpi PROPERTY POSITION_INDEPENDENT_CODE TRUE)
SET_PROPERTY(TARGET xxt PROPERTY POSITION_INDEPENDENT_CODE TRUE)

# Fix for CMake < 2.8.9 which does not support PIC property
IF (CMAKE_VERSION VERSION_LESS 2.8.9)
    SET(picflag "-fPIC")
ENDIF()

SET_PROPERTY(TARGET gsmpi APPEND PROPERTY COMPILE_FLAGS "${picflag} ${MPI_CXX_COMPILE_FLAGS}")
SET_PROPERTY(TARGET xxt APPEND PROPERTY COMPILE_FLAGS "${picflag} ${MPI_CXX_COMPILE_FLAGS}")

INSTALL(TARGETS gsmpi DESTINATION lib)
INSTALL(TARGETS xxt DESTINATION lib)
