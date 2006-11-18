
FIND_PATH(ACML_INCLUDE_PATH acml.h /usr/include /usr/local/include
			/opt/acml3.5.0/gnu64/include
                        /opt/acml3.5.0/gfortran64_int64/include
                        /opt/acml3.5.0/gfortran64_mp_int64/include )

#FIND_LIBRARY( BOOST_UNIT_TEST_LIB NAMES boost_unit_test_framework
#	      PATHS /usr/lib /usr/local/lib C:\\Boost\\lib )
#FIND_LIBRARY( BOOST_PROGRAM_OPTIONS_LIB NAMES boost_program_options
#	      PATHS /usr/lib /usr/local/lib C:\\Boost\\lib )
#FIND_LIBRARY( BOOST_FILESYSTEM_LIB NAMES boost_filesystem
#	      PATHS /usr/lib /usr/local/lib C:\\Boost\\lib )

SET(ACML_LIB_PATH ${ACML_INCLUDE_PATH}/../lib)

FIND_LIBRARY( ACML NAMES acml PATHS ${ACML_LIB_PATH} )

SET( ACML_BLAS_INCLUDE_FILE ${ACML_INCLUDE_PATH}/acml.h )
SET( ACML_LAPACK_INCLUDE_FILE ${ACML_INCLUDE_PATH}/acml.h )

IF (ACML_INCLUDE_PATH)
  SET(ACML_FOUND ON)
ENDIF (ACML_INCLUDE_PATH)

IF (ACML_FOUND)
  IF (NOT ACML_FIND_QUIETLY)
     MESSAGE(STATUS "Found ACML: ${ACML_INCLUDE_PATH}")
  ENDIF (NOT ACML_FIND_QUIETLY)
ELSE(ACML_FOUND)
  IF (ACML_FIND_REQUIRED)
     MESSAGE(FATAL_ERROR "Could not find ACML")
  ENDIF (ACML_FIND_REQUIRED)
ENDIF (ACML_FOUND)
