SET(BOOST_INCLUDE_SEARCH_PATH /usr/include /usr/local/include)

IF( ${CMAKE_GENERATOR} STREQUAL "Visual Studio 7 .NET 2003" )
	SET(BOOST_INCLUDE_SEARCH_PATH ${BOOST_INCLUDE_SEARCH_PATH} 
	"C:\\Program Files\\Microsoft Visual Studio .NET 2003\\Vc7\\include"
	c:\\Boost\\include )
ENDIF( ${CMAKE_GENERATOR} STREQUAL "Visual Studio 7 .NET 2003" )

IF( ${CMAKE_GENERATOR} STREQUAL "Visual Studio 8 2005" )
	SET(BOOST_INCLUDE_SEARCH_PATH ${BOOST_INCLUDE_SEARCH_PATH} 
		"C:\\Program Files (x86)\\Microsoft Visual Studio 8\\VC\\include"
		"C:\\Program Files\\Microsoft Visual Studio 8\\VC\\include"
		c:\\Boost\\include)
ENDIF( ${CMAKE_GENERATOR} STREQUAL "Visual Studio 8 2005" )

FIND_PATH(BOOST_INCLUDE_DIR boost/weak_ptr.hpp ${BOOST_INCLUDE_SEARCH_PATH} )

SET(BoostFileSystemName "")
SET(BoostFileSystemDebugName "")

IF( ${CMAKE_GENERATOR} STREQUAL "Visual Studio 7 .NET 2003" )
	SET(BoostFileSystemName "libboost_filesystem-vc71-mt")
	SET(BoostFileSystemDebugName "libboost_filesystem-vc71-mt-gd")	
	SET(BoostThreadName "boost_thread-vc71-mt")
	SET(BoostThreadDebugName "boost_thread-vc71-mt-gd")
ENDIF( ${CMAKE_GENERATOR} STREQUAL "Visual Studio 7 .NET 2003" )

IF( ${CMAKE_GENERATOR} STREQUAL "Visual Studio 8 2005" )
	SET(BoostFileSystemName "libboost_filesystem-vc80-mt")
	SET(BoostFileSystemDebugName "libboost_filesystem-vc80-mt-gd")	
ENDIF( ${CMAKE_GENERATOR} STREQUAL "Visual Studio 8 2005" )

IF( ${CMAKE_GENERATOR} STREQUAL "Unix Makefiles" )
	IF( ${CMAKE_COMPILER_IS_GNUCXX} )
		SET(BoostFileSystemName "boost_filesystem-gcc-mt")
		SET(BoostFileSystemDebugName "boost_filesystem-gcc-mt-d")
	ENDIF( ${CMAKE_COMPILER_IS_GNUCXX} )
ENDIF( ${CMAKE_GENERATOR} STREQUAL "Unix Makefiles" )


FIND_LIBRARY( BOOST_FILESYSTEM_LIB NAMES ${BoostFileSystemName}
	      PATHS 
	      ${BOOST_INCLUDE_DIR}/../lib
	      /usr/local/lib 
	      /usr/lib 
	      C:\\Boost\\lib )

FIND_LIBRARY( BOOST_FILESYSTEM_DEBUG_LIB NAMES ${BoostFileSystemDebugName}
	      PATHS 
	      ${BOOST_INCLUDE_DIR}/../lib
	      /usr/local/lib 
	      /usr/lib 
	      C:\\Boost\\lib )

FIND_LIBRARY( BOOST_THREAD_LIB NAMES ${BoostThreadName}
	      PATHS 
	      ${BOOST_INCLUDE_DIR}/../lib
	      /usr/local/lib 
	      /usr/lib 
	      C:\\Boost\\lib )

FIND_LIBRARY( BOOST_THREAD_DEBUG_LIB NAMES ${BoostThreadDebugName}
	      PATHS 
	      ${BOOST_INCLUDE_DIR}/../lib
	      /usr/local/lib 
	      /usr/lib 
	      C:\\Boost\\lib )
	      
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
