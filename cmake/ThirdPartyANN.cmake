########################################################################
#
# ThirdParty configuration for Nektar++
#
# libann partitioner
#
########################################################################

IF (NOT WIN32)
   OPTION(NEKTAR_USE_ANN
       "Use ANN routines for performing Approximate Nearest Neighbour searches." OFF)
ENDIF(NOT WIN32)

IF( NEKTAR_USE_MESHGEN )
    SET(NEKTAR_USE_ANN ON CACHE BOOL "" FORCE)
ENDIF()

IF (NEKTAR_USE_ANN)
    # First search for system ANN installs. Hint /opt/local for MacPorts and
    # /usr/local/opt/ann for Homebrew.
    FIND_LIBRARY(ANN_LIBRARY NAMES ANN
        PATHS /opt/local/lib /usr/local/opt/ann/lib $ENV{ANN_ROOT}/lib)
    FIND_PATH   (ANN_INCLUDE_DIR ANN.h
        PATHS /opt/local/include /usr/local/opt/ann/include $ENV{ANN_ROOT}/include
        PATH_SUFFIXES ANN)
    GET_FILENAME_COMPONENT(ANN_LIBRARY_PATH ${ANN_LIBRARY} PATH)

    IF (ANN_LIBRARY AND ANN_INCLUDE_DIR)
        SET(BUILD_ANN OFF)
    ELSE()
        SET(BUILD_ANN ON)
    ENDIF ()

    OPTION(THIRDPARTY_BUILD_ANN "Build ANN library from ThirdParty" ${BUILD_ANN})

    IF (THIRDPARTY_BUILD_ANN)
        # Note that ANN is compiled in the source-tree, so we unpack the
        # source code in the ThirdParty builds directory.
        SET(ANN_DIR ${TPBUILD}/ann-1.1.2)
        SET(ANN_SRC ${ANN_DIR}/src)

        IF (APPLE)
            SET(ANN_CFLAGS "-O3 -fPIC")
            SET(ANN_MAKELIB "libtool -static -o")
        ELSE ()
            SET(ANN_CFLAGS "-O3 -fPIC")
            SET(ANN_MAKELIB "ar ruv")
        ENDIF ()

        INCLUDE(ExternalProject)
        EXTERNALPROJECT_ADD(
            ann-1.1.2
            PREFIX ${TPSRC}
            URL ${TPURL}/ann_1.1.2.tar.gz
            URL_MD5 "9f99653b76798ecb1cfadc88950c4707"
            STAMP_DIR ${TPBUILD}/stamp
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPBUILD}/ann-1.1.2
            BINARY_DIR ${TPBUILD}/ann-1.1.2
            TMP_DIR ${TPBUILD}/ann-1.1.2-tmp
            INSTALL_DIR ${TPDIST}
            CONFIGURE_COMMAND ${CMAKE_COMMAND} -E remove -f ${ANN_DIR}/Makefile
            BUILD_COMMAND cd src
                 COMMAND $(MAKE) -C ${ANN_SRC} targets
                "ANNLIB  = libANN.a"
                "C++     = ${CMAKE_CXX_COMPILER}"
                "CFLAGS  = ${ANN_CFLAGS}"
                "MAKELIB = ${ANN_MAKELIB}"
                "RANLIB  = true"
            INSTALL_COMMAND ${CMAKE_COMMAND} -E make_directory ${TPDIST}/lib
                COMMAND ${CMAKE_COMMAND} -E copy ${ANN_DIR}/lib/libANN.a
                                                 ${TPDIST}/lib
                COMMAND ${CMAKE_COMMAND} -E copy_directory ${ANN_DIR}/include
                                                 ${TPDIST}/include
        )

        SET(ANN_LIBRARY ANN CACHE FILEPATH
            "ANN library" FORCE)
        SET(ANN_INCLUDE_DIR ${TPDIST}/include CACHE FILEPATH
            "ANN include directory" FORCE)

        LINK_DIRECTORIES(${TPDIST}/lib)

        MESSAGE(STATUS "Build ANN: ${TPDIST}/lib/lib${ANN_LIBRARY}.a")
        SET(ANN_CONFIG_INCLUDE_DIR ${TPINC})
    ELSE (THIRDPARTY_BUILD_ANN)
        ADD_CUSTOM_TARGET(ann-1.1.2 ALL)
        MESSAGE(STATUS "Found ANN: ${ANN_LIBRARY}")
        SET(ANN_CONFIG_INCLUDE_DIR ${ANN_INCLUDE_DIR})
    ENDIF (THIRDPARTY_BUILD_ANN)

    INCLUDE_DIRECTORIES(${ANN_INCLUDE_DIR})

    MARK_AS_ADVANCED(ANN_LIBRARY)
    MARK_AS_ADVANCED(ANN_INCLUDE_DIR)
ENDIF()
