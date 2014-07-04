OPTION(NEKTAR_USE_SCOTCH
    "Use Scotch library for performing mesh partitioning." OFF)

CMAKE_DEPENDENT_OPTION(THIRDPARTY_BUILD_SCOTCH
    "Build Scotch library from ThirdParty" OFF
    "NEKTAR_USE_SCOTCH" OFF)

IF( NEKTAR_USE_SCOTCH )
    IF (THIRDPARTY_BUILD_SCOTCH)
        UNSET(FLEX CACHE)
        FIND_PROGRAM(FLEX flex)
        IF(NOT FLEX)
            MESSAGE(FATAL_ERROR
                "'flex' lexical parser not found. Cannot build scotch.")
        ENDIF(NOT FLEX)

        SET(SCOTCH_SRC ${TPSRC}/src/scotch-6.0.0/src)

        IF (APPLE)
            SET(SCOTCH_MAKE Makefile.inc.i686_mac_darwin8)
            SET(SCOTCH_LDFLAGS "")
        ELSE ()
            IF (CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
                SET(SCOTCH_MAKE Makefile.inc.x86-64_pc_linux2)
            ELSE ()
                SET(SCOTCH_MAKE Makefile.inc.i686_pc_linux2)
            ENDIF ()
            SET(SCOTCH_LDFLAGS "-lz -lm -lrt -lpthread")
        ENDIF ()

        INCLUDE(ExternalProject)
        EXTERNALPROJECT_ADD(
            scotch-6.0.0
            PREFIX ${TPSRC}
            URL ${TPURL}/scotch_6.0.0.tar.gz
            URL_MD5 "ba117428c0a6cd97d0c93e8b872bb3fe"
            DOWNLOAD_DIR ${TPSRC}
            CONFIGURE_COMMAND rm -f ${SCOTCH_SRC}/Makefile.inc
                COMMAND ln -s
                ${SCOTCH_SRC}/Make.inc/${SCOTCH_MAKE}
                ${SCOTCH_SRC}/Makefile.inc
            BUILD_COMMAND $(MAKE) -C ${TPSRC}/src/scotch-6.0.0/src
                "LDFLAGS=${SCOTCH_LDFLAGS}"
                "CLIBFLAGS=-fPIC" scotch
            INSTALL_COMMAND $(MAKE) -C ${TPSRC}/src/scotch-6.0.0/src
                prefix=${TPSRC}/dist install
        )
        SET(SCOTCH_LIB scotch CACHE FILEPATH
            "Scotch library" FORCE)
        SET(SCOTCHERR_LIB scotcherr CACHE FILEPATH
            "Scotch error library" FORCE)
        SET(SCOTCHMETIS_LIB scotchmetis CACHE FILEPATH
            "Scotch Metis interface library" FORCE)
        MARK_AS_ADVANCED(SCOTCH_LIB)
        MARK_AS_ADVANCED(SCOTCHERR_LIB)
        MARK_AS_ADVANCED(SCOTCHMETIS_LIB)
        LINK_DIRECTORIES(${TPSRC}/dist/lib)
        INCLUDE_DIRECTORIES(${TPSRC}/dist/include)
        MESSAGE(STATUS "Build Scotch: ${TPSRC}/dist/lib/lib${SCOTCH_LIB}.a")
    ELSE (THIRDPARTY_BUILD_SCOTCH)
        INCLUDE (FindScotch)
    ENDIF (THIRDPARTY_BUILD_SCOTCH)
ENDIF( NEKTAR_USE_SCOTCH )
