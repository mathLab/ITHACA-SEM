SET(LinkLibraries MultiRegions)

SET(TimingCGHelmSolve2DSource TimingCGHelmSolve2D.cpp)
ADD_NEKTAR_EXECUTABLE(TimingCGHelmSolve2D timing TimingCGHelmSolve2DSource)
TARGET_LINK_LIBRARIES(TimingCGHelmSolve2D ${LinkLibraries})
SET_TARGET_PROPERTIES(TimingCGHelmSolve2D PROPERTIES COMPILE_FLAGS -DBASE_PATH="\\"${CMAKE_CURRENT_SOURCE_DIR}\\"")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Scripts/do_TimingCGHelmSolve2D       
    ${CMAKE_BINARY_DIR}/dist/bin/do_TimingCGHelmSolve2D COPYONLY)

SET(TimingCGHelmSolve3DSource TimingCGHelmSolve3D.cpp)
ADD_NEKTAR_EXECUTABLE(TimingCGHelmSolve3D timing TimingCGHelmSolve3DSource)
TARGET_LINK_LIBRARIES(TimingCGHelmSolve3D ${LinkLibraries})
SET_TARGET_PROPERTIES(TimingCGHelmSolve3D PROPERTIES COMPILE_FLAGS -DBASE_PATH="\\"${CMAKE_CURRENT_SOURCE_DIR}\\"")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Scripts/do_TimingCGHelmSolve3D       
    ${CMAKE_BINARY_DIR}/dist/bin/do_TimingCGHelmSolve3D COPYONLY)

SET(TimingHDGHelmSolve2DSource TimingHDGHelmSolve2D.cpp)
ADD_NEKTAR_EXECUTABLE(TimingHDGHelmSolve2D timing TimingHDGHelmSolve2DSource)
TARGET_LINK_LIBRARIES(TimingHDGHelmSolve2D ${LinkLibraries})
SET_TARGET_PROPERTIES(TimingHDGHelmSolve2D PROPERTIES COMPILE_FLAGS -DBASE_PATH="\\"${CMAKE_CURRENT_SOURCE_DIR}\\"")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Scripts/do_TimingHDGHelmSolve2D       
    ${CMAKE_BINARY_DIR}/dist/bin/do_TimingHDGHelmSolve2D COPYONLY)

SET(TimingHDGHelmSolve3DSource TimingHDGHelmSolve3D.cpp)
ADD_NEKTAR_EXECUTABLE(TimingHDGHelmSolve3D timing TimingHDGHelmSolve3DSource)
TARGET_LINK_LIBRARIES(TimingHDGHelmSolve3D ${LinkLibraries})
SET_TARGET_PROPERTIES(TimingHDGHelmSolve3D PROPERTIES COMPILE_FLAGS -DBASE_PATH="\\"${CMAKE_CURRENT_SOURCE_DIR}\\"")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Scripts/do_TimingHDGHelmSolve3D       
    ${CMAKE_BINARY_DIR}/dist/bin/do_TimingHDGHelmSolve3D COPYONLY)
#----------------------------------------------------------------------
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Scripts/do_TimingHelmSolve3D       
    ${CMAKE_BINARY_DIR}/dist/bin/do_TimingHelmSolve3D COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Scripts/do_TimingHelmSolve2D       
    ${CMAKE_BINARY_DIR}/dist/bin/do_TimingHelmSolve2D COPYONLY)
#----------------------------------------------------------------------

SET(TimingGeneralMatrixOp2DSource TimingGeneralMatrixOp2D.cpp)
ADD_NEKTAR_EXECUTABLE(TimingGeneralMatrixOp2D
    timing TimingGeneralMatrixOp2DSource)
TARGET_LINK_LIBRARIES(TimingGeneralMatrixOp2D ${LinkLibraries})
SET_TARGET_PROPERTIES(TimingGeneralMatrixOp2D PROPERTIES COMPILE_FLAGS -DBASE_PATH="\\"${CMAKE_CURRENT_SOURCE_DIR}\\"")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Scripts/do_TimingGeneralMatrixOp2D       
    ${CMAKE_BINARY_DIR}/dist/bin/do_TimingGeneralMatrixOp2D COPYONLY)

SET(TimingGeneralMatrixOp3DSource TimingGeneralMatrixOp3D.cpp)
ADD_NEKTAR_EXECUTABLE(TimingGeneralMatrixOp3D
    timing TimingGeneralMatrixOp3DSource)
TARGET_LINK_LIBRARIES(TimingGeneralMatrixOp3D ${LinkLibraries})
SET_TARGET_PROPERTIES(TimingGeneralMatrixOp3D PROPERTIES COMPILE_FLAGS -DBASE_PATH="\\"${CMAKE_CURRENT_SOURCE_DIR}\\"")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Scripts/do_TimingGeneralMatrixOp3D   
    ${CMAKE_BINARY_DIR}/dist/bin/do_TimingGeneralMatrixOp3D COPYONLY)

SET(TimingCGGeneralMatrixOp3DSource TimingCGGeneralMatrixOp3D.cpp)
ADD_NEKTAR_EXECUTABLE(TimingCGGeneralMatrixOp3D
    timing TimingCGGeneralMatrixOp3DSource)
TARGET_LINK_LIBRARIES(TimingCGGeneralMatrixOp3D ${LinkLibraries})
SET_TARGET_PROPERTIES(TimingCGGeneralMatrixOp3D PROPERTIES COMPILE_FLAGS -DBASE_PATH="\\"${CMAKE_CURRENT_SOURCE_DIR}\\"")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Scripts/do_TimingCGGeneralMatrixOp3D
    ${CMAKE_BINARY_DIR}/dist/bin/do_TimingCGGeneralMatrixOp3D COPYONLY)

    
