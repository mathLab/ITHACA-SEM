# Doxygen support
# adapted from http://www.bluequartz.net/projects/EIM_Segmentation/             SoftwareDocumentation/html/usewithcmakeproject.html
OPTION(NEKTAR_BUILD_DOC "Build source code documentation using doxygen" OFF)

CMAKE_DEPENDENT_OPTION(NEKTAR_BUILD_DOC_QHP 
    "Use Doxygen to create documentation for Qt Creator" OFF
    "NEKTAR_BUILD_DOC" OFF)
CMAKE_DEPENDENT_OPTION(NEKTAR_BUILD_DOC_XCODE 
    "Use Doxygen to create documentation for Apple Xcode" OFF
    "NEKTAR_BUILD_DOC" OFF)
CMAKE_DEPENDENT_OPTION(NEKTAR_BUILD_DOC_ECLIPSE
    "Use Doxygen to create documentation for Eclipse" OFF
    "NEKTAR_BUILD_DOC" OFF)

IF (NEKTAR_BUILD_DOC)
    FIND_PACKAGE(Doxygen)
    IF(NOT DOXYGEN_FOUND)
        MESSAGE(WARNING 
                "Doxygen not found. Building the documentation will fail.")
    ENDIF()
    IF(NEKTAR_BUILD_DOC_QHP)
        SET(DOXYGEN_GENERATE_QHP "YES")
    ENDIF()
    IF(NEKTAR_BUILD_DOC_XCODE)
        SET(DOXYGEN_GENERATE_DOCSET "YES")
    ENDIF()
    IF(NEKTAR_BUILD_DOC_ECLIPSE)
        SET(DOXYGEN_GENERATE_ECLIPSEHELP "YES")
    ENDIF()
    INSTALL(DIRECTORY ${PROJECT_BINARY_DIR}/doxygen/ 
        DESTINATION ${NEKTAR_DOC_DIR})

    CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/docs/doxygen/Doxyfile.in     
                   ${PROJECT_BINARY_DIR}/Doxyfile @ONLY IMMEDIATE)
    ADD_CUSTOM_TARGET(doc ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
        SOURCES ${PROJECT_BINARY_DIR}/Doxyfile)

ENDIF (NEKTAR_BUILD_DOC)

