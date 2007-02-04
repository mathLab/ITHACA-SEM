SET(BOOST_INCLUDE_SEARCH_PATH /usr/include /usr/local/include ${CMAKE_SOURCE_DIR}/../ThirdParty/include )

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

FIND_PATH(BOOST_BASE_DIR boost/weak_ptr.hpp ${BOOST_INCLUDE_SEARCH_PATH} )
SET(BOOST_INCLUDE_DIR ${BOOST_BASE_DIR})
    
SET(BoostFileSystemName "")
SET(BoostFileSystemDebugName "")

IF( ${CMAKE_GENERATOR} STREQUAL "Visual Studio 7 .NET 2003" )
    SET(BoostFileSystemName "libboost_filesystem-vc71-mt")
    SET(BoostFileSystemDebugName "libboost_filesystem-vc71-mt-gd")
    SET(BoostThreadName "boost_thread-vc71-mt")
    SET(BoostThreadDebugName "boost_thread-vc71-mt-gd")
ELSEIF( ${CMAKE_GENERATOR} STREQUAL "Visual Studio 8 2005" )
    SET(BoostFileSystemName "libboost_filesystem-vc80-mt")
    SET(BoostFileSystemDebugName "libboost_filesystem-vc80-mt-gd")
    SET(BoostThreadName "boost_thread-vc80-mt")
    SET(BoostThreadDebugName "boost_thread-vc80-mt-gd")
ELSE( ${CMAKE_GENERATOR} STREQUAL "Visual Studio 7 .NET 2003" )
    IF( ${CMAKE_COMPILER_IS_GNUCXX} )
        SET(BoostFileSystemName boost_filesystem-gcc boost_filesystem )
        SET(BoostFileSystemDebugName boost_filesystem-gcc-d boost_filesystem-d )
        SET(BoostThreadName boost_thread-gcc-mt boost_thread-gcc boost_thread-mt boost_thread  )
        SET(BoostThreadDebugName boost_thread-gcc-mt-d boost_thread-gcc-d  boost_thread-mt-d boost_thread-d )
    ENDIF( ${CMAKE_COMPILER_IS_GNUCXX} )
ENDIF( ${CMAKE_GENERATOR} STREQUAL "Visual Studio 7 .NET 2003" )

FIND_LIBRARY( BOOST_FILESYSTEM_LIB NAMES ${BoostFileSystemName}
          PATHS
          ${BOOST_BASE_DIR}/lib
	    ${BOOST_BASE_DIR}/../lib
          /usr/local/lib
          /usr/lib
          C:\\Boost\\lib )

FIND_LIBRARY( BOOST_FILESYSTEM_DEBUG_LIB NAMES ${BoostFileSystemDebugName}
          PATHS
          ${BOOST_BASE_DIR}/lib
          ${BOOST_BASE_DIR}/../lib
          /usr/local/lib
          /usr/lib
          C:\\Boost\\lib )

FIND_LIBRARY( BOOST_THREAD_LIB NAMES ${BoostThreadName}
          PATHS
          ${BOOST_BASE_DIR}/lib
	    ${BOOST_BASE_DIR}/../lib
          /usr/local/lib
          /usr/lib
          C:\\Boost\\lib )

FIND_LIBRARY( BOOST_THREAD_DEBUG_LIB NAMES ${BoostThreadDebugName}
          PATHS
          ${BOOST_BASE_DIR}/lib
          ${BOOST_BASE_DIR}/../lib
          /usr/local/lib
          /usr/lib
          C:\\Boost\\lib )

IF (BOOST_INCLUDE_DIR)
  SET(BOOST_FOUND TRUE)
ENDIF (BOOST_INCLUDE_DIR)

#SET (BOOST_LIB_DIR ${BOOST_BASE_DIR}/lib )
GET_FILENAME_COMPONENT(BOOST_LIB_DIR ${BOOST_THREAD_LIB} PATH CACHE)
LINK_DIRECTORIES(${BOOST_LIB_DIR})
    
IF (BOOST_FOUND)
  IF (NOT Boost_FIND_QUIETLY)
     MESSAGE(STATUS "Found Boost: ${BOOST_INCLUDE_DIR}")
  ENDIF (NOT Boost_FIND_QUIETLY)
ELSE(BOOST_FOUND)
  IF (Boost_FIND_REQUIRED)
     MESSAGE(FATAL_ERROR "Could not find Boost")
  ENDIF (Boost_FIND_REQUIRED)
ENDIF (BOOST_FOUND)
