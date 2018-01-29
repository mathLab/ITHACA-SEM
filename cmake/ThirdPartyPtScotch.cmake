########################################################################
#
# ThirdParty configuration for Nektar++
#
# Scotch partitioner
#
########################################################################

IF (NOT WIN32)
    OPTION(NEKTAR_USE_PTSCOTCH
        "Use PtScotch library for parallel mesh partitioning." OFF)
ENDIF(NOT WIN32)

IF (NEKTAR_USE_PTSCOTCH)
    # First search for system TinyXML installs. Hint /opt/local for MacPorts.
    FIND_LIBRARY(PTSCOTCH_LIBRARY    NAMES ptscotch PATHS ${MACPORTS_PREFIX}/lib)
    FIND_LIBRARY(PTSCOTCHERR_LIBRARY NAMES ptscotcherr PATHS ${MACPORTS_PREFIX}/lib)
    FIND_PATH   (PTSCOTCH_INCLUDE_DIR ptscotch.h PATHS ${MACPORTS_PREFIX}/include /usr/include/scotch)

    MESSAGE(STATUS ${PTSCOTCH_LIBRARY} ${PTSCOTCHERR_LIBRARY} ${PTSCOTCH_INCLUDE_DIR})
    
    IF (PTSCOTCH_LIBRARY AND PTSCOTCHERR_LIBRARY AND PTSCOTCH_INCLUDE_DIR)
        SET(BUILD_PTSCOTCH OFF)
    ELSE()
        SET(BUILD_PTSCOTCH ON)
    ENDIF ()

    CMAKE_DEPENDENT_OPTION(THIRDPARTY_BUILD_PTSCOTCH
        "Build Scotch library from ThirdParty" ${BUILD_PTSCOTCH}
        "NEKTAR_USE_PTSCOTCH" OFF)

    IF (THIRDPARTY_BUILD_PTSCOTCH)
        UNSET(FLEX CACHE)
        FIND_PROGRAM(FLEX flex)
        IF(NOT FLEX)
            MESSAGE(FATAL_ERROR
                "'flex' lexical parser not found. Cannot build scotch.")
        ENDIF(NOT FLEX)
        MARK_AS_ADVANCED(FLEX)

        # Note that scotch is compiled in the source-tree, so we unpack the
        # source code in the ThirdParty builds directory.
        SET(PTSCOTCH_SRC ${TPBUILD}/ptscotch-6.0.4/src)

        IF (APPLE)
            SET(PTSCOTCH_MAKE Makefile.inc.i686_mac_darwin8)
            SET(PTSCOTCH_LDFLAGS "")
            SET(PTSCOTCH_CFLAGS "-O3 -Drestrict=__restrict -DCOMMON_PTHREAD -DCOMMON_RANDOM_FIXED_SEED -DCOMMON_TIMING_OLD -DSCOTCH_RENAME -DCOMMON_PTHREAD_BARRIER")
        ELSE ()
            IF (CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
                SET(PTSCOTCH_MAKE Makefile.inc.x86-64_pc_linux2)
                SET(PTSCOTCH_CFLAGS "-O3 -DCOMMON_FILE_COMPRESS_GZ -DCOMMON_PTHREAD -DCOMMON_RANDOM_FIXED_SEED -DSCOTCH_RENAME -Drestrict=__restrict -DIDXSIZE64")
            ELSE ()
                SET(PTSCOTCH_MAKE Makefile.inc.i686_pc_linux2)
                SET(PTSCOTCH_CFLAGS "-O3 -DCOMMON_FILE_COMPRESS_GZ -DCOMMON_PTHREAD -DCOMMON_RANDOM_FIXED_SEED -DSCOTCH_RENAME -Drestrict=__restrict")
            ENDIF ()
            SET(PTSCOTCH_LDFLAGS "-lz -lm -lrt -lpthread")
        ENDIF ()

        INCLUDE(ExternalProject)
        EXTERNALPROJECT_ADD(
            ptscotch-6.0.4
            PREFIX ${TPSRC}
            URL ${TPURL}/scotch_6.0.4.tar.gz
            URL_MD5 "d58b825eb95e1db77efe8c6ff42d329f"
            STAMP_DIR ${TPBUILD}/stamp
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPBUILD}/ptscotch-6.0.4
            BINARY_DIR ${TPBUILD}/ptscotch-6.0.4
            TMP_DIR ${TPBUILD}/ptscotch-6.0.4-tmp
            INSTALL_DIR ${TPDIST}
            CONFIGURE_COMMAND rm -f ${PTSCOTCH_SRC}/Makefile.inc
                COMMAND ln -s
                ${PTSCOTCH_SRC}/Make.inc/${PTSCOTCH_MAKE}
                ${PTSCOTCH_SRC}/Makefile.inc
            BUILD_COMMAND $(MAKE) -C ${PTSCOTCH_SRC}
                "CFLAGS=${PTSCOTCH_CFLAGS}"
                "LDFLAGS=${PTSCOTCH_LDFLAGS}"
                "CLIBFLAGS=-fPIC"
                "CCP=${MPI_C_COMPILER}"
                "CCD=${MPI_C_COMPILER}"
                ptscotch
            INSTALL_COMMAND $(MAKE) -C ${PTSCOTCH_SRC}
                prefix=${TPDIST} install
        )

        THIRDPARTY_LIBRARY(PTSCOTCH_LIBRARY STATIC ptscotch
            DESCRIPTION "PT-Scotch library")
        THIRDPARTY_LIBRARY(PTSCOTCHERR_LIBRARY STATIC ptscotcherr
            DESCRIPTION "PT-Scotch error library")
        SET(PTSCOTCH_INCLUDE_DIR ${TPDIST}/include CACHE FILEPATH
            "PT-Scotch include directory" FORCE)
        MESSAGE(STATUS "Build PT-Scotch: ${PTSCOTCH_LIBRARY}")
        SET(PTSCOTCH_CONFIG_INCLUDE_DIR ${TPINC})
    ELSE (THIRDPARTY_BUILD_PTSCOTCH)
        ADD_CUSTOM_TARGET(ptscotch-6.0.4 ALL)
        MESSAGE(STATUS "Found PT-Scotch: ${PTSCOTCH_LIBRARY}")
        SET(PTSCOTCH_CONFIG_INCLUDE_DIR ${PTSCOTCH_INCLUDE_DIR})
    ENDIF (THIRDPARTY_BUILD_PTSCOTCH)

    INCLUDE_DIRECTORIES(${PTSCOTCH_INCLUDE_DIR})

    MARK_AS_ADVANCED(PTSCOTCH_LIBRARY)
    MARK_AS_ADVANCED(PTSCOTCHERR_LIBRARY)
    MARK_AS_ADVANCED(PTSCOTCH_INCLUDE_DIR)
ENDIF()
