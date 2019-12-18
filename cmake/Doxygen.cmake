# Doxygen support
# adapted from http://www.bluequartz.net/projects/EIM_Segmentation/
#    SoftwareDocumentation/html/usewithcmakeproject.html
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
CMAKE_DEPENDENT_OPTION(NEKTAR_BUILD_DOC_FIXEDWIDTH
    "Use a fixed-width style sheet for doxygen output" OFF
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

    SET(DOXYGEN_EXTRA_CSS "")
    IF(NEKTAR_BUILD_DOC_FIXEDWIDTH)
        SET(DOXYGEN_EXTRA_CSS "docs/doxygen/doxygen-fixed-width.css")
    ENDIF()

    INSTALL(DIRECTORY ${PROJECT_BINARY_DIR}/docs/doxygen/html/
        DESTINATION ${NEKTAR_DOC_DIR}/doxygen)

    CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/docs/doxygen/Doxyfile.in
                   ${PROJECT_BINARY_DIR}/docs/doxygen/Doxyfile @ONLY IMMEDIATE)

    ADD_CUSTOM_TARGET(doc
        COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_SOURCE_DIR}/docs/doxygen/doxygen-fixed-width.css
            ${PROJECT_BINARY_DIR}/docs/doxygen/doxygen-fixed-width.css
        COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/docs/doxygen/Doxyfile
        SOURCES ${PROJECT_BINARY_DIR}/docs/doxygen/Doxyfile)

ENDIF (NEKTAR_BUILD_DOC)

