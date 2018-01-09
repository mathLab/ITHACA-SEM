MACRO(SOLVER_COMPONENT_NAME name var)
    # Lower-case solver name
    STRING(TOLOWER ${name} ${var})
    # Then replace ending "solver" with "-solver"
    STRING(REPLACE "solver" "-solver" ${var} ${${var}})
ENDMACRO()

MACRO(ADD_SOLVER_EXECUTABLE name)
    CMAKE_PARSE_ARGUMENTS(SOLVER "" "COMPONENT" "SOURCES;LIBRARY_SOURCES" ${ARGN})

    IF (SOLVER_LIBRARY_SOURCES)
        ADD_LIBRARY(${name}-obj OBJECT ${SOLVER_LIBRARY_SOURCES})
        ADD_DEPENDENCIES(${name}-obj SolverUtils)
        LIST(APPEND SOLVER_SOURCES $<TARGET_OBJECTS:${name}-obj>)
    ENDIF()

    IF (NOT SOLVER_COMPONENT)
        SOLVER_COMPONENT_NAME(${name} SOLVER_COMPONENT)
    ENDIF()

    ADD_NEKTAR_EXECUTABLE(${name}
        COMPONENT ${SOLVER_COMPONENT} DEPENDS SolverUtils SOURCES ${SOLVER_SOURCES})
ENDMACRO()

MACRO(ADD_SOLVER_UTILITY name solver)
    CMAKE_PARSE_ARGUMENTS(UTIL "USELIBRARY" "COMPONENT" "SOURCES" ${ARGN})

    IF (NOT UTIL_COMPONENT)
        SOLVER_COMPONENT_NAME(${solver} UTIL_COMPONENT)
    ENDIF()

    IF (UTIL_USELIBRARY)
        ADD_NEKTAR_EXECUTABLE(${name}
            COMPONENT ${UTIL_COMPONENT} DEPENDS SolverUtils
            SOURCES ${UTIL_SOURCES} $<TARGET_OBJECTS:${solver}-obj>)
    ELSE()
        ADD_NEKTAR_EXECUTABLE(${name}
            COMPONENT ${UTIL_COMPONENT} DEPENDS SolverUtils
            SOURCES ${UTIL_SOURCES})
    ENDIF()
ENDMACRO()

# Generate list of available subdirectories
FILE(GLOB dir_list "*")
FOREACH(dir ${dir_list})
    IF(IS_DIRECTORY ${dir} AND EXISTS ${dir}/CMakeLists.txt)    
        ADD_SUBDIRECTORY(${dir})
    ENDIF(IS_DIRECTORY ${dir} AND EXISTS ${dir}/CMakeLists.txt)
ENDFOREACH(dir ${dir_list})
