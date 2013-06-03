SET(LinkLibraries 
    optimized MultiRegions debug MultiRegions-g    
    optimized LocalRegions debug LocalRegions-g    
    optimized SpatialDomains debug SpatialDomains-g    
    optimized StdRegions debug StdRegions-g
    optimized LibUtilities debug LibUtilities-g
    optimized ${Boost_THREAD_LIBRARY_RELEASE} debug ${Boost_THREAD_LIBRARY_DEBUG}
    optimized ${Boost_IOSTREAMS_LIBRARY_RELEASE} debug ${Boost_IOSTREAMS_LIBRARY_DEBUG}
    optimized ${ZLIB_LIBRARY_RELEASE} debug ${ZLIB_LIBRARY_DEBUG}	
    optimized ${Boost_DATE_TIME_LIBRARY_RELEASE} debug ${Boost_DATE_TIME_LIBRARY_DEBUG}
    optimized ${Boost_FILESYSTEM_LIBRARY_RELEASE} debug ${Boost_FILESYSTEM_LIBRARY_DEBUG}
    optimized ${Boost_SYSTEM_LIBRARY_RELEASE} debug ${Boost_SYSTEM_LIBRARY_DEBUG}
    optimized tinyxml debug tinyxml
    optimized sbtk    debug sbtk)

SET(TimingCGHelmSolve2DSource TimingCGHelmSolve2D.cpp)
ADD_NEKTAR_EXECUTABLE(TimingCGHelmSolve2D timing TimingCGHelmSolve2DSource)
TARGET_LINK_LIBRARIES(TimingCGHelmSolve2D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(TimingCGHelmSolve2D)
SET_TARGET_PROPERTIES(TimingCGHelmSolve2D PROPERTIES COMPILE_FLAGS -DBASE_PATH="\\"${CMAKE_CURRENT_SOURCE_DIR}\\"")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Scripts/do_TimingCGHelmSolve2D       
    ${CMAKE_BINARY_DIR}/dist/bin/do_TimingCGHelmSolve2D COPYONLY)

SET(TimingHDGHelmSolve2DSource TimingHDGHelmSolve2D.cpp)
ADD_NEKTAR_EXECUTABLE(TimingHDGHelmSolve2D timing TimingHDGHelmSolve2DSource)
TARGET_LINK_LIBRARIES(TimingHDGHelmSolve2D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(TimingHDGHelmSolve2D)
SET_TARGET_PROPERTIES(TimingHDGHelmSolve2D PROPERTIES COMPILE_FLAGS -DBASE_PATH="\\"${CMAKE_CURRENT_SOURCE_DIR}\\"")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Scripts/do_TimingHDGHelmSolve2D       
    ${CMAKE_BINARY_DIR}/dist/bin/do_TimingHDGHelmSolve2D COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Scripts/do_TimingHelmSolve2D       
    ${CMAKE_BINARY_DIR}/dist/bin/do_TimingHelmSolve2D COPYONLY)

SET(TimingGeneralMatrixOp2DSource TimingGeneralMatrixOp2D.cpp)
ADD_NEKTAR_EXECUTABLE(TimingGeneralMatrixOp2D
    timing TimingGeneralMatrixOp2DSource)
TARGET_LINK_LIBRARIES(TimingGeneralMatrixOp2D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(TimingGeneralMatrixOp2D)
SET_TARGET_PROPERTIES(TimingGeneralMatrixOp2D PROPERTIES COMPILE_FLAGS -DBASE_PATH="\\"${CMAKE_CURRENT_SOURCE_DIR}\\"")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Scripts/do_TimingGeneralMatrixOp2D       
    ${CMAKE_BINARY_DIR}/dist/bin/do_TimingGeneralMatrixOp2D COPYONLY)

SET(TimingGeneralMatrixOp3DSource TimingGeneralMatrixOp3D.cpp)
ADD_NEKTAR_EXECUTABLE(TimingGeneralMatrixOp3D
    timing TimingGeneralMatrixOp3DSource)
TARGET_LINK_LIBRARIES(TimingGeneralMatrixOp3D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(TimingGeneralMatrixOp3D)
SET_TARGET_PROPERTIES(TimingGeneralMatrixOp3D PROPERTIES COMPILE_FLAGS -DBASE_PATH="\\"${CMAKE_CURRENT_SOURCE_DIR}\\"")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Scripts/do_TimingGeneralMatrixOp3D   
    ${CMAKE_BINARY_DIR}/dist/bin/do_TimingGeneralMatrixOp3D COPYONLY)

SET(TimingCGGeneralMatrixOp3DSource TimingCGGeneralMatrixOp3D.cpp)
ADD_NEKTAR_EXECUTABLE(TimingCGGeneralMatrixOp3D
    timing TimingCGGeneralMatrixOp3DSource)
TARGET_LINK_LIBRARIES(TimingCGGeneralMatrixOp3D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(TimingCGGeneralMatrixOp3D)
SET_TARGET_PROPERTIES(TimingCGGeneralMatrixOp3D PROPERTIES COMPILE_FLAGS -DBASE_PATH="\\"${CMAKE_CURRENT_SOURCE_DIR}\\"")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Scripts/do_TimingCGGeneralMatrixOp3D
    ${CMAKE_BINARY_DIR}/dist/bin/do_TimingCGGeneralMatrixOp3D COPYONLY)

    
