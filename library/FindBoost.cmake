
FIND_PATH(BOOST_INCLUDE_DIR boost/weak_ptr.hpp /usr/include /usr/local/include 
	"C:\\Program Files\\Microsoft Visual Studio .NET 2003\\Vc7\\include"
	c:\\Boost\\include )

#FIND_LIBRARY( BOOST_UNIT_TEST_LIB NAMES boost_unit_test_framework
#	      PATHS /usr/lib /usr/local/lib C:\\Boost\\lib )
#FIND_LIBRARY( BOOST_PROGRAM_OPTIONS_LIB NAMES boost_program_options
#	      PATHS /usr/lib /usr/local/lib C:\\Boost\\lib )
#FIND_LIBRARY( BOOST_FILESYSTEM_LIB NAMES boost_filesystem
#	      PATHS /usr/lib /usr/local/lib C:\\Boost\\lib )

IF (BOOST_INCLUDE_DIR)
  SET(BOOST_FOUND TRUE)
ENDIF (BOOST_INCLUDE_DIR)

IF (BOOST_FOUND)
  IF (NOT Boost_FIND_QUIETLY)
     MESSAGE(STATUS "Found Boost: ${BOOST_INCLUDE_DIR}")
  ENDIF (NOT Boost_FIND_QUIETLY)
ELSE(BOOST_FOUND)
  IF (Boost_FIND_REQUIRED)
     MESSAGE(FATAL_ERROR "Could not find Boost")
  ENDIF (Boost_FIND_REQUIRED)
ENDIF (BOOST_FOUND)
