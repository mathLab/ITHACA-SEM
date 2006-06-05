
FIND_PATH(LOKI_INCLUDE_DIR loki/Typelist.h /usr/include /usr/local/include 
	"C:\\Program Files\\Microsoft Visual Studio .NET 2003\\Vc7\\include" )

#FIND_LIBRARY( BOOST_UNIT_TEST_LIB NAMES boost_unit_test_framework
#	      PATHS /usr/lib /usr/local/lib C:\\Boost\\lib )
#FIND_LIBRARY( BOOST_PROGRAM_OPTIONS_LIB NAMES boost_program_options
#	      PATHS /usr/lib /usr/local/lib C:\\Boost\\lib )
#FIND_LIBRARY( BOOST_FILESYSTEM_LIB NAMES boost_filesystem
#	      PATHS /usr/lib /usr/local/lib C:\\Boost\\lib )

IF (LOKI_INCLUDE_DIR)
  SET(LOKI_FOUND TRUE)
ENDIF (LOKI_INCLUDE_DIR)

IF (LOKI_FOUND)
  IF (NOT Loki_FIND_QUIETLY)
     MESSAGE(STATUS "Found Loki: ${LOKI_INCLUDE_DIR}")
  ENDIF (NOT Loki_FIND_QUIETLY)
ELSE(LOKI_FOUND)
  IF (Loki_FIND_REQUIRED)
     MESSAGE(FATAL_ERROR "Could not find Loki")
  ENDIF (Loki_FIND_REQUIRED)
ENDIF (LOKI_FOUND)
