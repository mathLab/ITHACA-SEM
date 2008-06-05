# Find an installed Nektar++.  Will not find a Nektar++ development environment.
SET(NEKTAR_SEARCH_PATH 
	${CMAKE_SOURCE_DIR}/../../Nektar++/dist
	${CMAKE_BINARY_DIR}/../Nektar++/dist/
	/usr/include 
	/usr/local/include
)  

FIND_PATH(NEKTAR++_ROOT_DIR Nektar++Config.cmake  ${NEKTAR_SEARCH_PATH} )  	       

IF (NEKTAR++_ROOT_DIR)   
	SET(NEKTAR++_FOUND TRUE) 
ELSE (NEKTAR++_ROOT_DIR)  
	SET(NEKTAR++_FOUND FALSE)
ENDIF (NEKTAR++_ROOT_DIR)  

IF( NEKTAR++_FOUND )
	INCLUDE(${NEKTAR++_ROOT_DIR}/Nektar++Config.cmake)

	SET(NEKTAR++_INCLUDE_DIRS ${NEKTAR++_ROOT_DIR}/include)
	SET(NEKTAR++_LIBRARY_DIRS ${NEKTAR++_ROOT_DIR}/lib)

	SET(NEKTAR++_LIBRARIES LibUtilities LocalRegions MultiRegions SpatialDomains StdRegions)
	# TODO - This is a little awkward because Nektar++ already knows what 
	# preprocessor definitions are defined, but I can't find a way to 
	# access them.  This is a maintenance issue.
	IF( NEKTAR_USE_MEMORY_POOLS )
		ADD_DEFINITIONS(-DNEKTAR_MEMORY_POOL_ENABLED)
	ELSE( NEKTAR_USE_MEMORY_POOLS )
		REMOVE_DEFINITIONS(-DNEKTAR_MEMORY_POOL_ENABLED)
	ENDIF( NEKTAR_USE_MEMORY_POOLS )

ELSE( NEKTAR++_FOUND )
	SET(NEKTAR++_DEFINITIONS "")
	SET(NEKTAR++_INCLUDE_DIRS "")
ENDIF( NEKTAR++_FOUND )

