##
## NektarCommon.cmake
##
## Frequently used Nektar++ CMake configuration macros and functions
##

#
# THIRDPARTY_LIBRARY(varname DESCRIPTION <description> [STATIC|SHARED] lib1 [lib2]...)
#
# Updates a variable containing the name of a third-party shared or static
# library to point to an absolute path defining its location instead of adding
# `-llibname` to the linker flags. This avoids the issue of e.g. linking against
# an outdated system zlib installation.
#
# Arguments:
#   - `varname`: variable name containing the third-party library name. On
#     output will be altered to update the correct path.
#   - `DESCRIPTION`: a brief description of the variable (used in the SET
#     command).
#   - `SHARED`: if the library will be built as a shared library
#   - `STATIC`: if the library will be built as a static library
#
MACRO(THIRDPARTY_LIBRARY varname)
    CMAKE_PARSE_ARGUMENTS(TPLIB "" "DESCRIPTION" "STATIC;SHARED" ${ARGN})

    IF(TPLIB_SHARED)
        SET(LIBTYPE "SHARED")
        SET(TPLIBS ${TPLIB_SHARED})
    ELSEIF(TPLIB_STATIC)
        SET(LIBTYPE "STATIC")
        SET(TPLIBS ${TPLIB_STATIC})
    ENDIF()

    FOREACH (lib ${TPLIBS})
        LIST(APPEND tmplist "${TPDIST}/lib/${CMAKE_${LIBTYPE}_LIBRARY_PREFIX}${lib}${CMAKE_${LIBTYPE}_LIBRARY_SUFFIX}")
    ENDFOREACH()

    SET(${varname} ${tmplist} CACHE FILEPATH ${TPLIB_DESCRIPTION} FORCE)
    UNSET(tmplist)
    UNSET(LIBTYPE)
    UNSET(TPLIBS)
    UNSET(TPLIB_SHARED)
    UNSET(TPLIB_STATIC)
    UNSET(lib)
ENDMACRO()

#
# SET_COMMON_PROPERTIES(target)
#
# Sets properties that are common to either library or executable targets. This
# includes:
#
# - Name suffixes: -g for debug, -ms for minsize, -rg for release w/debug.
# - Disable some MSVC compiler warnings
# - Add -pg flag if NEKTAR_ENABLE_PROFILE is switched on and we're using gcc
# - Add compiler definitions and appropriate warning levels to gcc-like
#   compilers (e.g. clang)
# - Define versions for the target
# - Make sure that -fPIC is enabled for library code if building shared
#   libraries.
#
# Arguments:
#   - `target`: target name
#
MACRO(SET_COMMON_PROPERTIES name)
    SET_TARGET_PROPERTIES(${name} PROPERTIES DEBUG_POSTFIX -g)
    SET_TARGET_PROPERTIES(${name} PROPERTIES MINSIZEREL_POSTFIX -ms)
    SET_TARGET_PROPERTIES(${name} PROPERTIES RELWITHDEBINFO_POSTFIX -rg)
    
    IF (MSVC)
        # Enable production-level warnings
        TARGET_COMPILE_OPTIONS(${name} PRIVATE /W4)
        # Temporarily disable signed/unsigned comparison warning
        TARGET_COMPILE_OPTIONS(${name} PRIVATE /wd4018)
        # Temporarily disable narrowing warnings
        TARGET_COMPILE_OPTIONS(${name} PRIVATE /wd4244 /wd4267)
        # Enable source-level parallel builds
        TARGET_COMPILE_OPTIONS(${name} PRIVATE /MP)
        # Specify minimum Windows version (501=WinXP, 601=Windows 7)
        TARGET_COMPILE_DEFINITIONS(${name} PRIVATE _WIN32_WINNT=0x0601)
    ELSE ()
        # Enable all warnings
        TARGET_COMPILE_OPTIONS(${name} PRIVATE -Wall -Wextra)
        IF (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
            # For GNU compilers add pedantic warnings
            TARGET_COMPILE_OPTIONS(${name} PRIVATE -Wpedantic)
        ENDIF()
        # Temporarily disable warnings about comparing signed and unsigned
        TARGET_COMPILE_OPTIONS(${name} PRIVATE -Wno-sign-compare)
        # Temporarily disable warnings about narrowing of data types
        TARGET_COMPILE_OPTIONS(${name} PRIVATE -Wno-narrowing -Wno-conversion)
        IF (NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
            TARGET_COMPILE_OPTIONS(${name} PRIVATE -fpermissive)
        ENDIF()
        
        # Disable dignostic about partially overloaded virtual functions
        IF (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
            TARGET_COMPILE_OPTIONS(${name} PRIVATE -diag-disable 654)
        ENDIF()

        IF ( NEKTAR_ERROR_ON_WARNINGS )
            TARGET_COMPILE_OPTIONS(${name} PRIVATE -Werror)
        ENDIF()
    ENDIF()

    IF (${CMAKE_COMPILER_IS_GNUCXX})
        IF(NEKTAR_ENABLE_PROFILE)
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
            SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pg")
            SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -pg")
            SET(LINK_FLAGS "${LINK_FLAGS} -pg")
        ENDIF(NEKTAR_ENABLE_PROFILE)
    ENDIF()

    # Prevent including these common flags multiple times.
    IF (NOT ${CMAKE_CXX_FLAGS_DEBUG} MATCHES ".*DNEKTAR_DEBUG.*")
        SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DNEKTAR_DEBUG")

        IF ( NEKTAR_FULL_DEBUG )
            SET(CMAKE_CXX_FLAGS_DEBUG 
                    "${CMAKE_CXX_FLAGS_DEBUG} -DNEKTAR_FULLDEBUG")
        ENDIF( NEKTAR_FULL_DEBUG)
   
        # Define version
        SET_PROPERTY(TARGET ${name}
            APPEND PROPERTY COMPILE_DEFINITIONS
            NEKTAR_VERSION=\"${NEKTAR_VERSION}\")

        SET(CMAKE_CXX_FLAGS_RELEASE 
                "${CMAKE_CXX_FLAGS_RELEASE} -DNEKTAR_RELEASE")
    ENDIF(NOT ${CMAKE_CXX_FLAGS_DEBUG} MATCHES ".*DNEKTAR_DEBUG.*")
        
    IF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
        # The static libraries must be compiled with position independent
        # code on 64 bit Linux.
        SET_PROPERTY(TARGET ${name} APPEND PROPERTY COMPILE_FLAGS "-fPIC")
    ENDIF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
ENDMACRO(SET_COMMON_PROPERTIES name)

#
# ADD_NEKTAR_EXECUTABLE(name COMPONENT <component> [DEPENDS dep1 ...] [SOURCES src1 ...])
#
# Adds a new executable to a component with the supplied component dependencies
# and sources files.
#
# Arguments:
#   - `name`: target name to construct
#   - `COMPONENT`: component name in which this target will live (e.g. demos)
#   - `DEPENDS`: a list of components on which this target depends on
#   - `SOURCES`: a list of source files for this target
#
MACRO(ADD_NEKTAR_EXECUTABLE name)
    CMAKE_PARSE_ARGUMENTS(NEKEXE "" "COMPONENT" "DEPENDS;SOURCES" ${ARGN})
    ADD_EXECUTABLE(${name} ${NEKEXE_SOURCES})
    SET_COMMON_PROPERTIES(${name})

    IF (${CMAKE_SYSTEM} MATCHES "Linux.*")
        SET_PROPERTY(TARGET ${name} APPEND_STRING PROPERTY COMPILE_FLAGS " -pthread")
        SET_PROPERTY(TARGET ${name} APPEND_STRING PROPERTY LINK_FLAGS " -pthread")
    ENDIF()

    STRING(TOLOWER ${NEKEXE_COMPONENT} NEKEXE_COMPONENT)
    STRING(TOUPPER ${NEKEXE_COMPONENT} NEKEXE_COMPVAR)

    SET_PROPERTY(TARGET ${name} PROPERTY FOLDER ${NEKEXE_COMPONENT})
    INSTALL(TARGETS ${name}
        RUNTIME DESTINATION ${NEKTAR_BIN_DIR} COMPONENT ${NEKEXE_COMPONENT} OPTIONAL
        ARCHIVE DESTINATION ${NEKTAR_LIB_DIR} COMPONENT ${NEKEXE_COMPONENT} OPTIONAL
        LIBRARY DESTINATION ${NEKTAR_LIB_DIR} COMPONENT ${NEKEXE_COMPONENT} OPTIONAL)

    # Add dependencies for executable.
    TARGET_LINK_LIBRARIES(${name} LINK_PUBLIC ${NEKEXE_DEPENDS})
ENDMACRO()

#
# ADD_NEKTAR_LIBRARY(name
#                    DESCRIPTION <description>
#                    SUMMARY <summary>
#                    DEPENDS dep1 dep2 ...
#                    SOURCES src1 src2 ...
#                    HEADERS head1 head2 ...)
#
# Adds a new library to a component with the supplied component dependencies and
# sources files. A new component will be set up automatically with a lower-case
# name: e.g. if the supplied library name is `LibUtilities` the corresponding
# component is `libutilities`.
#
# Arguments:
#   - `name`: target name to construct
#   - `SUMMARY`: a brief summary of the library
#   - `DESCRIPTION`: a more detailed description of the library
#   - `DEPENDS`: a list of components on which this target depends on
#   - `SOURCES`: a list of source files for this target
#   - `HEADERS`: a list of header files for this target. These will be
#     automatically put into a `dev` package.
#
MACRO(ADD_NEKTAR_LIBRARY name)
    CMAKE_PARSE_ARGUMENTS(NEKLIB "" "DESCRIPTION;SUMMARY" "DEPENDS;SOURCES;HEADERS" ${ARGN})

    ADD_LIBRARY(${name} ${NEKTAR_LIBRARY_TYPE} ${NEKLIB_SOURCES} ${NEKLIB_HEADERS})

    # Infer component name from lower-case library name, variables should use
    # upper-case.
    STRING(TOLOWER ${name} NEKLIB_COMPONENT)
    STRING(TOUPPER ${name} NEKLIB_COMPVAR)

    # Add name to a list so that we know for constructing dependencies.
    SET(NEKTAR++_LIBRARIES ${NEKTAR++_LIBRARIES} ${name} CACHE INTERNAL "")

    SET_PROPERTY(TARGET ${name} PROPERTY FOLDER ${NEKLIB_COMPONENT})
    SET_PROPERTY(TARGET ${name} PROPERTY VERSION ${NEKTAR_VERSION})

    SET_COMMON_PROPERTIES(${name})

    INSTALL(TARGETS ${name}
        EXPORT Nektar++Libraries
        RUNTIME DESTINATION ${NEKTAR_BIN_DIR} COMPONENT ${NEKLIB_COMPONENT} OPTIONAL
        ARCHIVE DESTINATION ${NEKTAR_LIB_DIR} COMPONENT ${NEKLIB_COMPONENT} OPTIONAL
        LIBRARY DESTINATION ${NEKTAR_LIB_DIR} COMPONENT ${NEKLIB_COMPONENT} OPTIONAL)

    FOREACH(HEADER ${NEKLIB_HEADERS})
        STRING(REGEX MATCH "(.*)[/\\]" DIR ${HEADER})
        INSTALL(FILES ${HEADER}
            DESTINATION ${NEKTAR_INCLUDE_DIR}/${name}/${DIR}
            COMPONENT dev)
    ENDFOREACH()

    # If we have dependencies then link against them.
    IF(NEKLIB_DEPENDS)
        TARGET_LINK_LIBRARIES(${name} LINK_PUBLIC ${NEKLIB_DEPENDS})
    ENDIF()
ENDMACRO()

#
# ADD_NEKTAR_TEST(name [LENGTHY])
#
# Adds a test with a given name.  The Test Definition File should be in a
# subdirectory called Tests relative to the CMakeLists.txt file calling this
# macros. The test file should be called NAME.tst, where NAME is given as a
# parameter to this macro. If the LENGTHY flag is given, the test will only be
# run if `NEKTAR_TEST_ALL` is enabled.
#
# Arguments:
#   - `name`: name of the test file
#   - `LENGTHY`: denotes a test that requires extended runtime.
#
MACRO(ADD_NEKTAR_TEST name)
    CMAKE_PARSE_ARGUMENTS(NEKTEST "LENGTHY" "" "" ${ARGN})

    IF (NOT NEKTEST_LENGTHY OR NEKTAR_TEST_ALL)
        GET_FILENAME_COMPONENT(dir ${CMAKE_CURRENT_SOURCE_DIR} NAME)
        ADD_TEST(NAME ${dir}_${name}
            COMMAND Tester ${CMAKE_CURRENT_SOURCE_DIR}/Tests/${name}.tst)
    ENDIF()
ENDMACRO(ADD_NEKTAR_TEST)

#
# ADD_NEKPY_LIBRARY(name SOURCES src1 src2 ...)
#
# Adds a new NekPy library with the given sources.
#
MACRO(ADD_NEKPY_LIBRARY name)
    CMAKE_PARSE_ARGUMENTS(NEKPY "" "DEPENDS" "SOURCES" ${ARGN})

    # Create library.
    ADD_LIBRARY(_${name} SHARED ${NEKPY_SOURCES})

    # Python requires a .so extension, even on OS X.
    SET_TARGET_PROPERTIES(_${name} PROPERTIES PREFIX "")
    SET_TARGET_PROPERTIES(_${name} PROPERTIES SUFFIX ".so")

    ADD_DEPENDENCIES(_${name} boost-numpy)

    # Add target link libraries.
    TARGET_LINK_LIBRARIES(_${name}
        ${Boost_SYSTEM_LIBRARY}
        ${BOOST_PYTHON_LIB}
        ${BOOST_NUMPY_LIB}
        ${PYTHON_LIBRARIES}
        ${name})

    # Install __init__.py files.
    SET(TMPOUT "")
    IF (NEKPY_DEPENDS)
        SET(TMPOUT "from ..${NEKPY_DEPENDS} import _${NEKPY_DEPENDS}\n")
    ENDIF()
    SET(TMPOUT "${TMPOUT}from ._${name} import *")

    FILE(WRITE ${CMAKE_BINARY_DIR}/NekPy/${name}/__init__.py ${TMPOUT})
    INSTALL(TARGETS _${name} DESTINATION ${CMAKE_BINARY_DIR}/NekPy/${name})
ENDMACRO()

MACRO(ADD_NEKPY_EXECUTABLE name source)
    # Copy the files into binary directory.
    INSTALL(FILES ${source} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
    FILE(COPY ${source} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    ADD_CUSTOM_TARGET(${name} SOURCES ${source})
ENDMACRO()
