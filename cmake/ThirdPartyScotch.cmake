########################################################################
#
# ThirdParty configuration for Nektar++
#
# Scotch partitioner
#
########################################################################

IF (NOT WIN32)
    OPTION(NEKTAR_USE_SCOTCH
        "Use Scotch library for performing mesh partitioning." OFF)
ENDIF(NOT WIN32)

IF (NEKTAR_USE_SCOTCH)
    # First search for system TinyXML installs. Hint /opt/local for MacPorts.
    FIND_LIBRARY(SCOTCH_LIBRARY    NAMES scotch PATHS /opt/local/lib)
    FIND_LIBRARY(SCOTCHERR_LIBRARY NAMES scotcherr PATHS /opt/local/lib)
    FIND_PATH   (SCOTCH_INCLUDE_DIR scotch.h PATHS /opt/local/include)

    IF (SCOTCH_LIBRARY AND SCOTCHERR_LIBRARY AND SCOTCH_INCLUDE_DIR)
        SET(BUILD_SCOTCH OFF)
    ELSE()
        SET(BUILD_SCOTCH ON)
    ENDIF ()

    CMAKE_DEPENDENT_OPTION(THIRDPARTY_BUILD_SCOTCH
        "Build Scotch library from ThirdParty" ${BUILD_SCOTCH}
        "NEKTAR_USE_SCOTCH" OFF)

    IF (THIRDPARTY_BUILD_SCOTCH)
        UNSET(FLEX CACHE)
        FIND_PROGRAM(FLEX flex)
        IF(NOT FLEX)
            MESSAGE(FATAL_ERROR
                "'flex' lexical parser not found. Cannot build scotch.")
        ENDIF(NOT FLEX)
        MARK_AS_ADVANCED(FLEX)

        # Note that scotch is compiled in the source-tree, so we unpack the
        # source code in the ThirdParty builds directory.
        SET(SCOTCH_SRC ${TPBUILD}/scotch-6.0.0/src)

        IF (APPLE)
            SET(SCOTCH_MAKE Makefile.inc.i686_mac_darwin8)
            SET(SCOTCH_LDFLAGS "")
            SET(SCOTCH_CFLAGS "-O3 -Drestrict=__restrict -DCOMMON_PTHREAD -DCOMMON_RANDOM_FIXED_SEED -DCOMMON_TIMING_OLD -DSCOTCH_PTHREAD -DSCOTCH_RENAME -DCOMMON_PTHREAD_BARRIER")
        ELSE ()
            IF (CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
                SET(SCOTCH_MAKE Makefile.inc.x86-64_pc_linux2)
                SET(SCOTCH_CFLAGS "-O3 -DCOMMON_FILE_COMPRESS_GZ -DCOMMON_PTHREAD -DCOMMON_RANDOM_FIXED_SEED -DSCOTCH_RENAME -DSCOTCH_PTHREAD -Drestrict=__restrict -DIDXSIZE64")
            ELSE ()
                SET(SCOTCH_MAKE Makefile.inc.i686_pc_linux2)
                SET(SCOTCH_CFLAGS "-O3 -DCOMMON_FILE_COMPRESS_GZ -DCOMMON_PTHREAD -DCOMMON_RANDOM_FIXED_SEED -DSCOTCH_RENAME -DSCOTCH_PTHREAD -Drestrict=__restrict")
            ENDIF ()
            SET(SCOTCH_LDFLAGS "-lz -lm -lrt -lpthread")
        ENDIF ()

        INCLUDE(ExternalProject)
        EXTERNALPROJECT_ADD(
            scotch-6.0.0
            PREFIX ${TPSRC}
            URL ${TPURL}/scotch_6.0.0.tar.gz
            URL_MD5 "ba117428c0a6cd97d0c93e8b872bb3fe"
            STAMP_DIR ${TPBUILD}/stamp
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPBUILD}/scotch-6.0.0
            BINARY_DIR ${TPBUILD}/scotch-6.0.0
            TMP_DIR ${TPBUILD}/scotch-6.0.0-tmp
            INSTALL_DIR ${TPDIST}
            CONFIGURE_COMMAND rm -f ${SCOTCH_SRC}/Makefile.inc
                COMMAND ln -s
                ${SCOTCH_SRC}/Make.inc/${SCOTCH_MAKE}
                ${SCOTCH_SRC}/Makefile.inc
            BUILD_COMMAND $(MAKE) -C ${SCOTCH_SRC}
                "CFLAGS=${SCOTCH_CFLAGS}"
                "LDFLAGS=${SCOTCH_LDFLAGS}"
                "CLIBFLAGS=-fPIC" scotch
            INSTALL_COMMAND $(MAKE) -C ${SCOTCH_SRC}
                prefix=${TPDIST} install
        )

        SET(SCOTCH_LIBRARY scotch CACHE FILEPATH
            "Scotch library" FORCE)
        SET(SCOTCHERR_LIBRARY scotcherr CACHE FILEPATH
            "Scotch error library" FORCE)
        SET(SCOTCH_INCLUDE_DIR ${TPDIST}/include CACHE FILEPATH
            "Scotch include directory" FORCE)

        LINK_DIRECTORIES(${TPDIST}/lib)

        MESSAGE(STATUS "Build Scotch: ${TPDIST}/lib/lib${SCOTCH_LIBRARY}.a")
        SET(SCOTCH_CONFIG_INCLUDE_DIR ${TPINC})
    ELSE (THIRDPARTY_BUILD_SCOTCH)
        ADD_CUSTOM_TARGET(scotch-6.0.0 ALL)
        MESSAGE(STATUS "Found Scotch: ${SCOTCH_LIBRARY}")
        SET(SCOTCH_CONFIG_INCLUDE_DIR ${SCOTCH_INCLUDE_DIR})
    ENDIF (THIRDPARTY_BUILD_SCOTCH)

    INCLUDE_DIRECTORIES(${SCOTCH_INCLUDE_DIR})

    MARK_AS_ADVANCED(SCOTCH_LIBRARY)
    MARK_AS_ADVANCED(SCOTCHERR_LIBRARY)
    MARK_AS_ADVANCED(SCOTCH_INCLUDE_DIR)
ENDIF()
