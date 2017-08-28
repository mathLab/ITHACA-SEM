ADD_NEKTAR_EXECUTABLE(TimingCGHelmSolve2D 
                      COMPONENT timing 
                      DEPENDS MultiRegions 
                      SOURCES TimingCGHelmSolve2D.cpp)

SET_TARGET_PROPERTIES(TimingCGHelmSolve2D PROPERTIES COMPILE_FLAGS -DBASE_PATH="\\"${CMAKE_CURRENT_SOURCE_DIR}\\"")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Scripts/do_TimingCGHelmSolve2D       
    ${CMAKE_BINARY_DIR}/dist/bin/do_TimingCGHelmSolve2D COPYONLY)

ADD_NEKTAR_EXECUTABLE(TimingCGHelmSolve3D 
                      COMPONENT timing 
                      DEPENDS MultiRegions
                      SOURCES TimingCGHelmSolve3D.cpp)

SET_TARGET_PROPERTIES(TimingCGHelmSolve3D PROPERTIES COMPILE_FLAGS -DBASE_PATH="\\"${CMAKE_CURRENT_SOURCE_DIR}\\"")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Scripts/do_TimingCGHelmSolve3D       
    ${CMAKE_BINARY_DIR}/dist/bin/do_TimingCGHelmSolve3D COPYONLY)

ADD_NEKTAR_EXECUTABLE(TimingHDGHelmSolve2D 
                      COMPONENT timing 
                      DEPENDS MultiRegions
                      SOURCES TimingHDGHelmSolve2D.cpp)

SET_TARGET_PROPERTIES(TimingHDGHelmSolve2D PROPERTIES COMPILE_FLAGS -DBASE_PATH="\\"${CMAKE_CURRENT_SOURCE_DIR}\\"")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Scripts/do_TimingHDGHelmSolve2D       
    ${CMAKE_BINARY_DIR}/dist/bin/do_TimingHDGHelmSolve2D COPYONLY)

ADD_NEKTAR_EXECUTABLE(TimingHDGHelmSolve3D 
                      COMPONENT timing 
                      DEPENDS MultiRegions
                      SOURCES TimingHDGHelmSolve3D.cpp)

SET_TARGET_PROPERTIES(TimingHDGHelmSolve3D PROPERTIES COMPILE_FLAGS -DBASE_PATH="\\"${CMAKE_CURRENT_SOURCE_DIR}\\"")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Scripts/do_TimingHDGHelmSolve3D       
    ${CMAKE_BINARY_DIR}/dist/bin/do_TimingHDGHelmSolve3D COPYONLY)

ADD_NEKTAR_EXECUTABLE(TimingGeneralMatrixOp2D
                      COMPONENT timing 
                      DEPENDS MultiRegions
                      SOURCES TimingGeneralMatrixOp2D.cpp)

SET_TARGET_PROPERTIES(TimingGeneralMatrixOp2D PROPERTIES COMPILE_FLAGS -DBASE_PATH="\\"${CMAKE_CURRENT_SOURCE_DIR}\\"")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Scripts/do_TimingGeneralMatrixOp2D       
    ${CMAKE_BINARY_DIR}/dist/bin/do_TimingGeneralMatrixOp2D COPYONLY)

ADD_NEKTAR_EXECUTABLE(TimingGeneralMatrixOp3D
                      COMPONENT timing 
                      DEPENDS MultiRegions
                      SOURCES TimingGeneralMatrixOp3D.cpp)

SET_TARGET_PROPERTIES(TimingGeneralMatrixOp3D PROPERTIES COMPILE_FLAGS -DBASE_PATH="\\"${CMAKE_CURRENT_SOURCE_DIR}\\"")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Scripts/do_TimingGeneralMatrixOp3D   
    ${CMAKE_BINARY_DIR}/dist/bin/do_TimingGeneralMatrixOp3D COPYONLY)

ADD_NEKTAR_EXECUTABLE(TimingCGGeneralMatrixOp3D
                       COMPONENT timing 
                       DEPENDS MultiRegions
                       SOURCES TimingCGGeneralMatrixOp3D.cpp)

SET_TARGET_PROPERTIES(TimingCGGeneralMatrixOp3D PROPERTIES COMPILE_FLAGS -DBASE_PATH="\\"${CMAKE_CURRENT_SOURCE_DIR}\\"")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Scripts/do_TimingCGGeneralMatrixOp3D
    ${CMAKE_BINARY_DIR}/dist/bin/do_TimingCGGeneralMatrixOp3D COPYONLY)

    
