MACRO(CHANGE_EXTENSION output var new_ext)
    GET_FILENAME_COMPONENT(FileName ${var} NAME_WE)
    GET_FILENAME_COMPONENT(Path ${var} PATH)
    SET(${output} ${Path}/${FileName}.${new_ext})
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

MACRO(SETUP_PRECOMPILED_HEADERS sourceFiles precompiledHeader)
    IF( NEKTAR_USE_PRECOMPILED_HEADERS )
        IF( MSVC )	
            # /Yu"stdafx.h" 
            #MESSAGE(${${precompiledHeader}})
    	    #MESSAGE(${${sourceFiles}})
            SET_SOURCE_FILES_PROPERTIES(${${sourceFiles}} 
                PROPERTIES COMPILE_FLAGS "/Yu\"${${precompiledHeader}}\"")
            LIST(GET ${sourceFiles} 0 OUTVAR)
            #MESSAGE(${OUTVAR})
            SET_SOURCE_FILES_PROPERTIES(${OUTVAR} 
                PROPERTIES COMPILE_FLAGS "/Yc\"${${precompiledHeader}}\"")
            
        ENDIF()	
    ENDIF()
ENDMACRO()

MACRO(ADD_NEKTAR_EXECUTABLE name component sources)
    IF( ${ARGC} LESS 4 )
        ADD_EXECUTABLE(${name} ${${sources}})
    ELSE( ${ARGC} LESS 4 )
        ADD_EXECUTABLE(${name} ${${sources}} ${${ARGV3}})
    ENDIF( ${ARGC} LESS 4)
	
    SET_COMMON_PROPERTIES(${name})
    
    IF( ${CMAKE_SYSTEM} MATCHES "Linux.*" )
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
	
    ENDIF( ${CMAKE_SYSTEM} MATCHES "Linux.*" )

    SET_PROPERTY(TARGET ${name} PROPERTY FOLDER ${component})
	INSTALL(TARGETS ${name} 
		RUNTIME DESTINATION ${NEKTAR_BIN_DIR} COMPONENT ${component} OPTIONAL
		ARCHIVE DESTINATION ${NEKTAR_LIB_DIR} COMPONENT ${component} OPTIONAL
		LIBRARY DESTINATION ${NEKTAR_LIB_DIR} COMPONENT ${component} OPTIONAL)

ENDMACRO(ADD_NEKTAR_EXECUTABLE name component sources)

MACRO(ADD_NEKTAR_LIBRARY name component type)
    ADD_LIBRARY(${name} ${type} ${ARGN})

    SET_PROPERTY(TARGET ${name} PROPERTY FOLDER ${component})
    SET_PROPERTY(TARGET ${name} PROPERTY VERSION ${NEKTAR_VERSION})

    SET_COMMON_PROPERTIES(${name})

    # Set properties for building shared libraries
    IF( ${type} STREQUAL "SHARED" )
        # Properties specific to Mac OSX
        IF( ${CMAKE_SYSTEM} MATCHES "Darwin-*")
            # We allow undefined symbols to be looked up dynamically at runtime
            # from the boost libraries linked by the executables.
            SET_TARGET_PROPERTIES(${name} 
                PROPERTIES LINK_FLAGS "-Wl,-undefined,dynamic_lookup")
        ENDIF( ${CMAKE_SYSTEM} MATCHES "Darwin-*")
    ENDIF( ${type} STREQUAL "SHARED" )

    INSTALL(TARGETS ${name} 
        EXPORT Nektar++Libraries 
		RUNTIME DESTINATION ${NEKTAR_BIN_DIR} COMPONENT ${component} OPTIONAL
		ARCHIVE DESTINATION ${NEKTAR_LIB_DIR} COMPONENT ${component} OPTIONAL
		LIBRARY DESTINATION ${NEKTAR_LIB_DIR} COMPONENT ${component} OPTIONAL) 

ENDMACRO(ADD_NEKTAR_LIBRARY name component type)

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
