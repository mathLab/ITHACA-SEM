CMAKE_MINIMUM_REQUIRED(VERSION 2.4.2)
PROJECT(Nektar++)

# Find the modules included with Nektar
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR} ${CMAKE_MODULE_PATH})
SET(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/dist CACHE FILEPATH "" FORCE)

IF( ${CMAKE_COMPILER_IS_GNUCXX} )
	SET(USE_GCOV OFF CACHE BOOL "Use gcov coverage testing." )
ENDIF( ${CMAKE_COMPILER_IS_GNUCXX} )

SET(USE_MKL OFF CACHE BOOL "Use the Intel Math Kernel Library (MKL) for BLAS and Lapack support." )
SET(USE_ACML OFF CACHE BOOL "Use the AMD Core Math Library (ACML) for BLAS and Lapack support." )
SET(USE_ATLAS OFF CACHE BOOL "Use ATLAS for BLAS and Lapack support." )

SET(USE_EXPRESSION_TEMPLATES ON CACHE BOOL "")

INCLUDE (FindBoost)
INCLUDE (FindLoki)

IF( USE_MKL )
    INCLUDE (FindMKL)
ENDIF( USE_MKL )

IF( USE_ATLAS )
    INCLUDE( FindAtlas )
ENDIF( USE_ATLAS )

IF( USE_ACML )
    INCLUDE (FindACML)
ENDIF( USE_ACML )

INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR} ${BOOST_INCLUDE_DIR} ${LOKI_INCLUDE_DIR} ${PROJECT_BINARY_DIR})

SET(NEKTAR_VERSION 0.1.1)
SET(NEKTAR_BIN_DIR /bin)
SET(NEKTAR_LIB_DIR /lib)
SET(NEKTAR_INCLUDE_DIR /include)

MACRO(ADD_NEKTAR_LIBRARY name type)
    ADD_LIBRARY(${name} ${type} ${ARGN})

    SET_TARGET_PROPERTIES(${name} PROPERTIES VERSION ${NEKTAR_VERSION})

    SET_TARGET_PROPERTIES(${name} PROPERTIES DEBUG_POSTFIX -g)
    SET_TARGET_PROPERTIES(${name} PROPERTIES MINSIZEREL_POSTFIX -ms)
    SET_TARGET_PROPERTIES(${name} PROPERTIES RELWITHDEBINFO_POSTFIX -rg)

    IF( ${CMAKE_GENERATOR} STREQUAL "Visual Studio 8 2005" )	
    	ADD_DEFINITIONS(-D_USE_MATH_DEFINES)
    ENDIF( ${CMAKE_GENERATOR} STREQUAL "Visual Studio 8 2005" )	

    INSTALL_TARGETS(${NEKTAR_LIB_DIR} ${name})
    
ENDMACRO(ADD_NEKTAR_LIBRARY name)

IF( ${CMAKE_COMPILER_IS_GNUCXX} )	
	IF( USE_GCOV )
		ADD_DEFINITIONS(-fprofile-arcs -ftest-coverage)
	ENDIF( USE_GCOV )
ENDIF( ${CMAKE_COMPILER_IS_GNUCXX} )

IF( USE_EXPRESSION_TEMPLATES )
    ADD_DEFINITIONS(-DNEKTAR_USE_EXPRESSION_TEMPLATES -DNEKTAR_USING_CMAKE)
ENDIF( USE_EXPRESSION_TEMPLATES )

IF( USE_MKL )
    SET(NEKTAR_USING_BLAS TRUE)
    SET(NEKTAR_USING_LAPACK TRUE)
    SET(NEKTAR_USING_MKL TRUE)
ENDIF( USE_MKL )

IF( USE_ATLAS )
    SET(NEKTAR_USING_BLAS TRUE)
    SET(NEKTAR_USING_ATLAS TRUE)
ENDIF( USE_ATLAS )

IF( USE_ACML )
    SET(NEKTAR_USING_BLAS TRUE)
    SET(NEKTAR_USING_LAPACK TRUE)
    SET(NEKTAR_USING_ACML TRUE)
ENDIF( USE_ACML )

MACRO(ADD_NEKTAR_EXECUTABLE name sources)
    ADD_EXECUTABLE(${name} ${${sources}} )

    IF( ${CMAKE_GENERATOR} STREQUAL "Visual Studio 8 2005" )	
    	ADD_DEFINITIONS(-D_USE_MATH_DEFINES)
    ENDIF( ${CMAKE_GENERATOR} STREQUAL "Visual Studio 8 2005" )	

    IF( USE_MKL AND MKL_FOUND )
    #        MESSAGE("Using MKL")
            TARGET_LINK_LIBRARIES(${name}
                ${MKL_LAPACK}
                optimized ${MKL} debug ${MKL}
                ${MKL_GUIDE}
    
    #optimized mkl_ia32 debug mkl_ia32
    #            optimized mkl debug mkl
    #            optimized vml debug vml
    #            optimized guide debug guide
            )
    ENDIF( USE_MKL AND MKL_FOUND )
    
    IF( USE_ACML AND ACML_FOUND )
        TARGET_LINK_LIBRARIES(${name} ${ACML} gfortran )
    ENDIF( USE_ACML AND ACML_FOUND )
    
    IF( USE_ATLAS AND ATLAS_FOUND )
            TARGET_LINK_LIBRARIES(${name}
                ${ATLAS_LAPACK}
                ${ATLAS_CBLAS}
                ${ATLAS} )
    ENDIF( USE_ATLAS AND ATLAS_FOUND )
    
    IF( ${CMAKE_COMPILER_IS_GNUCXX} )
            IF( USE_GCOV )
            TARGET_LINK_LIBRARIES(${name} gcov)
        ENDIF( USE_GCOV )
    ENDIF( ${CMAKE_COMPILER_IS_GNUCXX} )
    
    TARGET_LINK_LIBRARIES(${name}
        optimized LibUtilities debug LibUtilities-g
        optimized StdRegions debug StdRegions-g
        optimized ${BoostThreadName} debug ${BoostThreadDebugName})
        
    INSTALL_TARGETS(${NEKTAR_BIN_DIR} ${name})
ENDMACRO(ADD_NEKTAR_EXECUTABLE name)

SUBDIRS(LibUtilities LocalRegions MultiRegions SpatialDomains StdRegions Stochastic UnitTests Demos)

