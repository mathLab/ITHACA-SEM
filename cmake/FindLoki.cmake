
# Use LOKI_ADDITIONAL_INCLUDE_DIRS to add additional search directories.

FIND_PATH(LOKI_INCLUDE_DIR loki/Typelist.h 
	${LOKI_ADDITIONAL_INCLUDE_DIRS}
	/usr/include 
	/usr/local/include 
	${CMAKE_SOURCE_DIR}/ThirdParty/loki/include/ 
	${CMAKE_SOURCE_DIR}/ThirdParty/loki-0.1.3/include 
        ${CMAKE_SOURCE_DIR}/ThirdParty/loki-0.1.3/loki-0.1.3/include 
	${CMAKE_SOURCE_DIR}/../ThirdParty/loki/include/
	${CMAKE_SOURCE_DIR}/../ThirdParty/loki-0.1.3/include/
        ${CMAKE_SOURCE_DIR}/../ThirdParty/loki-0.1.3/loki-0.1.3/include 
	"C:\\Program Files\\Microsoft Visual Studio .NET 2003\\Vc7\\include" 
	"C:\\Program Files\\Microsoft Visual Studio 8\\VC\\include" 
)


IF (LOKI_INCLUDE_DIR)
  SET(LOKI_FOUND TRUE)
ENDIF (LOKI_INCLUDE_DIR)

SET (LOKI_LIB_DIR ${LOKI_INCLUDE_DIR}../lib )

IF (LOKI_FOUND)
  IF (NOT Loki_FIND_QUIETLY)
     MESSAGE(STATUS "Found Loki: ${LOKI_INCLUDE_DIR}")
  ENDIF (NOT Loki_FIND_QUIETLY)
ELSE(LOKI_FOUND)
  IF (Loki_FIND_REQUIRED)
     MESSAGE(FATAL_ERROR "Could not find Loki")
  ENDIF (Loki_FIND_REQUIRED)
ENDIF (LOKI_FOUND)

MARK_AS_ADVANCED(LOKI_INCLUDE_DIR)
MARK_AS_ADVANCED(LOKI_LIB_DIR)
