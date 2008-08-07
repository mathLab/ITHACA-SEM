######################################################
# This CMakeLists.txt is used to generate an utilities
# project.  Customization points for your
# projects are marked below.
######################################################

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(utilities)

# INCLUDE(CMakeExportBuildSettings)

# Use the CMake 2.4 policy for CMake 2.6
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 OLD)
endif(COMMAND cmake_policy)


SET(TIXML_USE_STL ON CACHE BOOL "Add support for STL std::string")
MARK_AS_ADVANCED(TIXML_USE_STL)


# Setup where CMake will look for modules.  You probably 
# won't need to modify this.
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR} ${CMAKE_MODULE_PATH})


IF( TIXML_USE_STL )
    ADD_DEFINITIONS(-DTIXML_USE_STL)
ENDIF( TIXML_USE_STL )


# Finds Nektar++.  This may introduce other dependencies, such 
# as Boost and Loki.
# FIND_PACKAGE(Nektar++ REQUIRED)

# INCLUDE (FindBoost)
# INCLUDE (FindLoki)
INCLUDE (FindTinyXml)

INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR} ${Boost_INCLUDE_DIRS} ${LOKI_INCLUDE_DIR} ${TINYXML_BASE})

MACRO(SET_COMMON_PROPERTIES name)

    SET_TARGET_PROPERTIES(${name} PROPERTIES DEBUG_POSTFIX -g)
    SET_TARGET_PROPERTIES(${name} PROPERTIES MINSIZEREL_POSTFIX -ms)
    SET_TARGET_PROPERTIES(${name} PROPERTIES RELWITHDEBINFO_POSTFIX -rg)


ENDMACRO(SET_COMMON_PROPERTIES name)

MACRO(ADD_UTILITIES_EXECUTABLE name source)
    ADD_EXECUTABLE(${name} ${source})

    TARGET_LINK_LIBRARIES(${name}
        optimized tinyxml debug tinyxml)

    SET_COMMON_PROPERTIES(${name})

ENDMACRO(ADD_UTILITIES_EXECUTABLE name source)


MACRO(ADD_UTILITIES_LIBRARY name type)
    ADD_LIBRARY(${name} ${type} ${ARGN})

#     MESSAGE(${name})
#     MESSAGE(${type})
#     MESSAGE(${ARGN})
#     MESSAGE(${ARGC})
#     MESSAGE(${ARGV1})
#     MESSAGE(${ARGV2})
#     MESSAGE(${ARGV3})

    SET_COMMON_PROPERTIES(${name})

ENDMACRO(ADD_UTILITIES_LIBRARY name type)


SUBDIRS( PreProcessing )

