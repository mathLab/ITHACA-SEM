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
            URL ${TPURL}/vtk-5.10.1-nek.tar.bz2
            URL_MD5 "f4e2c6b848d3873d44479baa9e7e4d35"
            STAMP_DIR ${TPBUILD}/stamp
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPSRC}/vtk-5.10.1
            BINARY_DIR ${TPBUILD}/vtk-5.10.1
            TMP_DIR ${TPBUILD}/vtk-5.10.1-tmp
            INSTALL_DIR ${TPDIST}
            CONFIGURE_COMMAND ${CMAKE_COMMAND} 
                -G ${CMAKE_GENERATOR}
                -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST} 
                -DBUILD_SHARED_LIBS:BOOL=ON 
                -DCMAKE_BUILD_TYPE:STRING=Release 
                ${TPSRC}/vtk-5.10.1
        )
        SET(VTK_DIR ${TPDIST}/lib/vtk-5.10)
        SET(VTK_FOUND 1)
        SET(VTK_USE_FILE ${VTK_DIR}/UseVTK.cmake)
        INCLUDE (${VTK_DIR}/VTKConfig.cmake)
    ELSE()
        FIND_PACKAGE(VTK COMPONENTS
            vtkFiltersGeometry vtkIOLegacy vtkIOXML vtkIOImage vtkRenderingCore)
        IF (VTK_FOUND)
            MESSAGE(STATUS "Found VTK: ${VTK_USE_FILE}")
            IF (VTK_MAJOR_VERSION EQUAL 6 AND VTK_MINOR_VERSION EQUAL 0 AND VTK_BUILD_VERSION EQUAL 0)
                ADD_DEFINITIONS(-DNEKTAR_HAS_VTK_6_0_0)
            ENDIF()
        ELSE (VTK_FOUND)
            MESSAGE(FATAL_ERROR "VTK not found")
        ENDIF (VTK_FOUND)
    ENDIF()
    
    INCLUDE (${VTK_USE_FILE})

    # Force VTK headers to be treated as system headers.
    INCLUDE_DIRECTORIES(SYSTEM ${VTK_INCLUDE_DIRS})

    MARK_AS_ADVANCED(VTK_DIR)
    ADD_DEFINITIONS(-DNEKTAR_USING_VTK)
ENDIF( NEKTAR_USE_VTK )

