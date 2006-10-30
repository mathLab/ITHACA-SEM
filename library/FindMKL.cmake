
FIND_PATH(MKL_INCLUDE_DIR mkl_cblas.h /usr/include /usr/local/include
			/opt/intel/mk1/8.1.1/include )

#FIND_LIBRARY( BOOST_UNIT_TEST_LIB NAMES boost_unit_test_framework
#	      PATHS /usr/lib /usr/local/lib C:\\Boost\\lib )
#FIND_LIBRARY( BOOST_PROGRAM_OPTIONS_LIB NAMES boost_program_options
#	      PATHS /usr/lib /usr/local/lib C:\\Boost\\lib )
#FIND_LIBRARY( BOOST_FILESYSTEM_LIB NAMES boost_filesystem
#	      PATHS /usr/lib /usr/local/lib C:\\Boost\\lib )

SET(MKL_LIB_PATH /opt/intel/mkl/8.1.1/lib/32)

FIND_LIBRARY( MKL_LAPACK NAMES mkl_lapack PATHS ${MKL_LIB_PATH} )
FIND_LIBRARY( MKL NAMES mkl_ia32 PATHS ${MKL_LIB_PATH} )
FIND_LIBRARY( MKL_GUIDE NAMES GUIDE PATHS ${MKL_LIB_PATH} )

SET( MKL_BLAS_INCLUDE_FILE ${MKL_INCLUDE_DIR}/mkl_cblas.h )
SET( MKL_LAPACK_INCLUDE_FILE ${MKL_INCLUDE_DIR}/mkl_lapack.h )

IF (MKL_INCLUDE_DIR)
  SET(MKL_FOUND ON)
ENDIF (MKL_INCLUDE_DIR)

IF (MKL_FOUND)
  IF (NOT MKL_FIND_QUIETLY)
     MESSAGE(STATUS "Found MKL: ${MKL_INCLUDE_DIR}")
  ENDIF (NOT MKL_FIND_QUIETLY)
ELSE(MKL_FOUND)
  IF (MKL_FIND_REQUIRED)
     MESSAGE(FATAL_ERROR "Could not find MKL")
  ENDIF (MKL_FIND_REQUIRED)
ENDIF (MKL_FOUND)
