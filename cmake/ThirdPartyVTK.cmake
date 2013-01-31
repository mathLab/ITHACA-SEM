OPTION(NEKTAR_USE_VTK "Use VTK library for utilities." OFF)

CMAKE_DEPENDENT_OPTION(THIRDPARTY_BUILD_VTK
    "Build VTK library from ThirdParty" OFF
    "NEKTAR_USE_VTK; FALSE" OFF)

IF( NEKTAR_USE_VTK )
    IF( THIRDPARTY_BUILD_VTK )
        INCLUDE( ExternalProject )

        # The cmake package has been modified due to a bug in the CMake files
        # which causes it to produce an error when the path includes a '+'.
        # Obviously this is inconvenient for us.
        EXTERNALPROJECT_ADD(
            vtk-5.10.1
            PREFIX ${TPSRC}
            URL ${TPURL}/vtk-5.10.1-nek.tar.bz2
            URL_MD5 "f4e2c6b848d3873d44479baa9e7e4d35"
            DOWNLOAD_DIR ${TPSRC}
            CONFIGURE_COMMAND ${CMAKE_COMMAND} -DCMAKE_INSTALL_PREFIX:PATH=${TPSRC}/dist -DBUILD_SHARED_LIBS:BOOL=ON -DCMAKE_BUILD_TYPE:STRING=Release ${TPSRC}/src/vtk-5.10.1
        )
        SET(VTK_DIR ${TPSRC}/dist/lib/vtk-5.10)
        SET(VTK_FOUND 1)
        SET(VTK_USE_FILE ${VTK_DIR}/UseVTK.cmake)
        INCLUDE (${VTK_DIR}/VTKConfig.cmake)
    ELSE()
        INCLUDE (FindVTK)
        IF (NOT VTK_FOUND)
            MESSAGE(FATAL_ERROR "VTK not found")
        ENDIF (NOT VTK_FOUND)
    ENDIF()
    
    INCLUDE (${VTK_USE_FILE})
    MARK_AS_ADVANCED(VTK_DIR)
    ADD_DEFINITIONS(-DNEKTAR_USING_VTK)
    INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS})
ENDIF( NEKTAR_USE_VTK )

