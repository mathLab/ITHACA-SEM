######################################################
# This CMakeLists.txt is used to generate an executable
# that uses Nektar++.  Customization points for your 
# projects are marked below.
######################################################

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# TODO - Put the name of your project here.
SET(ProjectName Sample)
PROJECT(${ProjectName})


# Setup where CMake will look for modules.  You probably 
# won't need to modify this.
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR} ${CMAKE_MODULE_PATH})

# Finds Nektar++.  This may introduce other dependencies, such 
# as Boost and Loki.
FIND_PACKAGE(Nektar++ REQUIRED)

# TODO - Add an entry here for every source file in your project. Must not be empty!
SET(SourceFileNames
	main.cpp
)

# TODO - Add an entry here for every header file in your project. Must not be empty!
SET(HeaderFileNames
        main.cpp
)

IF( NEKTAR++_FOUND )
    INCLUDE_DIRECTORIES(${NEKTAR++_INCLUDE_DIRS})
    LINK_DIRECTORIES(${NEKTAR++_LIBRARY_DIRS})
    ADD_DEFINITIONS(${NEKTAR++_DEFINITIONS})

    # temporary workarounds:
    ADD_DEFINITIONS("-DNEKTAR_USING_LAPACK -DNEKTAR_USING_BLAS")
    # workaround for libLibUtilities: http://stackoverflow.com/questions/17150075/undefined-reference-to-clock-gettime-although-lrt-is-given
    # to fix, add rt to TARGET_LINK_LIBRARIES(LibUtilities ... in library/LibUtilities/CMakeLists.txt:341, but i have no idea how non-Linux
    # platforms handle this
    SET(CMAKE_EXE_LINKER_FLAGS "-Wl,--no-as-needed")

    ADD_NEKTAR_EXECUTABLE(${ProjectName} SourceFileNames HeaderFileNames)

    # TODO - If your executable needs to link to any other libraries,  
    # add them here.
    TARGET_LINK_LIBRARIES(${ProjectName}
	${NEKTAR++_LIBRARIES}
	${Boost_THREAD_LIBRARY} 
	)

    SET_LAPACK_LINK_LIBRARIES(${ProjectName})
ENDIF( NEKTAR++_FOUND )

