OPTION(NEKTAR_USE_SCOTCH
    "Use Scotch library for performing mesh partitioning." OFF)

CMAKE_DEPENDENT_OPTION(THIRDPARTY_BUILD_SCOTCH
    "Build Scotch library from ThirdParty" OFF
    "NEKTAR_USE_SCOTCH" OFF)

IF( NEKTAR_USE_SCOTCH )
    IF (THIRDPARTY_BUILD_SCOTCH)
        MESSAGE("Building Scotch is not yet supported.")
    
        INCLUDE(ExternalProject)
        EXTERNALPROJECT_ADD(
            scotch-6.0.0
            PREFIX ${TPSRC}
            URL ${TPURL}/scotch_6.0.0.tar.gz
            URL_MD5 "6c6816aea0f53db6c71b1d98ed4ad42b"
            DOWNLOAD_DIR ${TPSRC}
            CONFIGURE_COMMAND ${CMAKE_COMMAND}
                -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                -DCMAKE_INSTALL_PREFIX:PATH=${TPSRC}/dist
                -DCMAKE_C_FLAGS:STRING=-fPIC
                ${TPSRC}/src/scotch_6.0.0
        )
        SET(SCOTCH_LIB scotch CACHE FILEPATH
            "Scotch library" FORCE)
        SET(SCOTCHMETIS_LIB scotchmetis CACHE FILEPATH
            "Scotch Metis interface library" FORCE)
        MARK_AS_ADVANCED(SCOTCH_LIB)
        MARK_AS_ADVANCED(SCOTCHMETIS_LIB)
        LINK_DIRECTORIES(${TPSRC}/dist/lib)
        INCLUDE_DIRECTORIES(${TPSRC}/dist/include)
        MESSAGE(STATUS "Build Scotch: ${TPSRC}/dist/lib/lib${SCOTCH_LIB}.a")
    ELSE (THIRDPARTY_BUILD_SCOTCH)
        INCLUDE (FindScotch)
    ENDIF (THIRDPARTY_BUILD_SCOTCH)
ENDIF( NEKTAR_USE_SCOTCH )
