MACRO(CONSTRUCT_DEBIAN_DEPS depends outvar)
    SET(${outvar} "")
    FOREACH (pkg ${depends})
        STRING(TOLOWER ${pkg} pkg_lower)
        SET(${outvar} "${DEB_DEPS}, nektar++-${pkg_lower} (>= ${NEKTAR_VERSION})")
    ENDFOREACH()

    # Remove starting ", "
    STRING(SUBSTRING ${${outvar}} 2 -1 ${outvar})
ENDMACRO()

MACRO(FINALISE_CPACK_COMPONENT name)
    CMAKE_PARSE_ARGUMENTS(COMP "" "DESCRIPTION" "" ${ARGN})
    SET(CPACK_COMPONENT_${name}_DISPLAY_NAME nektar++-${name})
    SET(CPACK_COMPONENT_${name}_DESCRIPTION ${COMP_DESCRIPTION})

    LIST(REMOVE_DUPLICATES "CPACK_COMPONENT_${name}_DEPENDS")
    CONSTRUCT_DEBIAN_DEPS(${CPACK_COMPONENT_${name}_DEPENDS} "CPACK_DEBIAN_${name}_DEPENDS")
ENDMACRO()

MACRO(SET_COMMON_PROPERTIES name)
    SET_TARGET_PROPERTIES(${name} PROPERTIES DEBUG_POSTFIX -g)
    SET_TARGET_PROPERTIES(${name} PROPERTIES MINSIZEREL_POSTFIX -ms)
    SET_TARGET_PROPERTIES(${name} PROPERTIES RELWITHDEBINFO_POSTFIX -rg)
    
    IF( MSVC )
        # Disable the warnings about duplicate copy/assignment methods 
        #   (4521, 4522)
        # Disable the warning that arrays are default intialized (4351)	
        # Disable "forcing value to bool 'true' or 'false' (performance
        #   warning)" warning (4800)
        # 4250 - Inheritance via dominance.  Nektar appears to be handling the 
        # diamond correctly.
        # 4373 - Overriding a virtual method with parameters that differ by const
        #        or volatile conforms to the standard.
        # /Za is necessary to prevent temporaries being bound to reference
        #   parameters.
        SET_TARGET_PROPERTIES(${name} PROPERTIES COMPILE_FLAGS 
                                "/wd4521 /wd4522 /wd4351 /wd4018 /wd4800 /wd4250 /wd4373")

        # Enable source level parallel builds.
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
    ENDIF( MSVC )	
    
    IF( ${CMAKE_COMPILER_IS_GNUCXX} )
        IF(NEKTAR_ENABLE_PROFILE)
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
            SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pg")
            SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -pg")
            SET(LINK_FLAGS "${LINK_FLAGS} -pg")
        ENDIF(NEKTAR_ENABLE_PROFILE)
    ENDIF( ${CMAKE_COMPILER_IS_GNUCXX} )

    # Prevent including these common flags multiple times.
    IF (NOT ${CMAKE_CXX_FLAGS_DEBUG} MATCHES ".*DNEKTAR_DEBUG.*")
        SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DNEKTAR_DEBUG")

        IF ( NEKTAR_FULL_DEBUG )
            SET(CMAKE_CXX_FLAGS_DEBUG 
                    "${CMAKE_CXX_FLAGS_DEBUG} -DNEKTAR_FULLDEBUG")
        ENDIF( NEKTAR_FULL_DEBUG)
   
        IF( NOT MSVC )
            SET(CMAKE_CXX_FLAGS_DEBUG 
                "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wno-deprecated -Wno-sign-compare")
            SET(CMAKE_CXX_FLAGS_RELEASE 
                    "${CMAKE_CXX_FLAGS_RELEASE} -Wall -Wno-deprecated -Wno-sign-compare")
            SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO
                    "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -Wall -Wno-deprecated -Wno-sign-compare")
            IF (NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
                SET(CMAKE_CXX_FLAGS_DEBUG 
                    "${CMAKE_CXX_FLAGS_DEBUG} -fpermissive")
            ENDIF()
        ENDIF( NOT MSVC)

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

MACRO(ADD_NEKTAR_EXECUTABLE name)
    CMAKE_PARSE_ARGUMENTS(NEKEXE "" "COMPONENT" "DEPENDS;SOURCES" ${ARGN})
    ADD_EXECUTABLE(${name} ${NEKEXE_SOURCES})
    SET_COMMON_PROPERTIES(${name})

    IF (${CMAKE_SYSTEM} MATCHES "Linux.*")
        # The boost thread library needs pthread on linux.
        GET_TARGET_PROPERTY(THE_COMPILE_FLAGS ${name} COMPILE_FLAGS)
        GET_TARGET_PROPERTY(THE_LINK_FLAGS ${name} LINK_FLAGS)

        # It is possible these flags have not been set yet.
        IF(NOT THE_COMPILE_FLAGS)
            SET(THE_COMPILE_FLAGS "")
        ENDIF(NOT THE_COMPILE_FLAGS)

        IF(NOT THE_LINK_FLAGS )
            SET(THE_LINK_FLAGS "")
        ENDIF(NOT THE_LINK_FLAGS)

        SET_TARGET_PROPERTIES(${name}
                PROPERTIES COMPILE_FLAGS "${THE_COMPILE_FLAGS} -pthread")
        SET_TARGET_PROPERTIES(${name}
                PROPERTIES LINK_FLAGS "${THE_LINK_FLAGS} -pthread")
    ENDIF()

    SET_PROPERTY(TARGET ${name} PROPERTY FOLDER ${NEKEXE_COMPONENT})
    INSTALL(TARGETS ${name}
        RUNTIME DESTINATION ${NEKTAR_BIN_DIR} COMPONENT ${NEKEXE_COMPONENT} OPTIONAL
        ARCHIVE DESTINATION ${NEKTAR_LIB_DIR} COMPONENT ${NEKEXE_COMPONENT} OPTIONAL
        LIBRARY DESTINATION ${NEKTAR_LIB_DIR} COMPONENT ${NEKEXE_COMPONENT} OPTIONAL)

    # Add dependencies for executable. We append the dependencies to the CPack
    # component list so that we can resolve these later.
    TARGET_LINK_LIBRARIES(${name} LINK_PUBLIC ${NEKEXE_DEPENDS})
    FOREACH(dep ${NEKEXE_DEPENDS})
        LIST(APPEND CPACK_COMPONENT_${NEKEXE_COMPONENT}_DEPENDS ${dep})
    ENDFOREACH()
ENDMACRO()

MACRO(ADD_NEKTAR_LIBRARY name)
    CMAKE_PARSE_ARGUMENTS(NEKLIB "" "DESCRIPTION;GROUP" "DEPENDS;SOURCES;HEADERS" ${ARGN})

    ADD_LIBRARY(${name} ${NEKTAR_LIBRARY_TYPE} ${NEKLIB_SOURCES} ${NEKLIB_HEADERS})

    # Infer component name from lower-case library name
    STRING(TOLOWER ${name} NEKLIB_COMPONENT)

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

    # Add CPack information
    SET(CPACK_COMPONENT_${NEKLIB_COMPONENT}_DISPLAY_NAME nektar++-${NEKLIB_COMPONENT})
    SET(CPACK_COMPONENT_${NEKLIB_COMPONENT}_DISPLAY_GROUP lib)
    SET(CPACK_COMPONENT_${NEKLIB_COMPONENT}_DESCRIPTION ${NEKLIB_DESCRIPTION})

    # If we have dependencies then link against them, and also configure CPack
    # Debian dependencies, which are a special case for some reason. Then set up
    # standard CPack components.
    IF(NEKLIB_DEPENDS)
        TARGET_LINK_LIBRARIES(${name} LINK_PUBLIC ${NEKLIB_DEPENDS})
        CONSTRUCT_DEBIAN_DEPS(${NEKLIB_DEPENDS} "CPACK_DEBIAN_${NEKLIB_NAME}_DEPENDS")
        SET(CPACK_COMPONENT_${NEKLIB_COMPONENT}_DEPENDS ${NEKLIB_DEPENDS})
    ENDIF()
ENDMACRO()

# Adds a test with a given name.
# The Test Definition File should be in a subdirectory called Tests relative
# to the CMakeLists.txt file calling this macros. The test file should be called
# NAME.tst, where NAME is given as a parameter to this macro.
MACRO(ADD_NEKTAR_TEST name)
    GET_FILENAME_COMPONENT(dir ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    ADD_TEST(NAME ${dir}_${name}
         COMMAND Tester ${CMAKE_CURRENT_SOURCE_DIR}/Tests/${name}.tst)
ENDMACRO(ADD_NEKTAR_TEST)

MACRO(ADD_NEKTAR_TEST_LENGTHY name)
    IF (NEKTAR_TEST_ALL)
        GET_FILENAME_COMPONENT(dir ${CMAKE_CURRENT_SOURCE_DIR} NAME)
        ADD_TEST(NAME ${dir}_${name}
             COMMAND Tester ${CMAKE_CURRENT_SOURCE_DIR}/Tests/${name}.tst)
    ENDIF(NEKTAR_TEST_ALL)
ENDMACRO(ADD_NEKTAR_TEST_LENGTHY)
