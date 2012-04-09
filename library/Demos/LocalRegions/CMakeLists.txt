SET(LinkLibraries 
    optimized LocalRegions debug LocalRegions-g  
    optimized SpatialDomains debug SpatialDomains-g
    optimized StdRegions debug StdRegions-g
    optimized LibUtilities debug LibUtilities-g
    optimized ${Boost_THREAD_LIBRARY_RELEASE} debug ${Boost_THREAD_LIBRARY_DEBUG}
    optimized ${Boost_IOSTREAMS_LIBRARY_RELEASE} debug ${Boost_IOSTREAMS_LIBRARY_DEBUG}
    optimized ${Boost_ZLIB_LIBRARY_RELEASE} debug ${Boost_ZLIB_LIBRARY_DEBUG}
    optimized ${TINYXML_LIB} debug ${TINYXML_LIB}
)

SET(LocProject1DSource  LocProject1D.cpp)
ADD_NEKTAR_EXECUTABLE(LocProject1D demos LocProject1DSource)
TARGET_LINK_LIBRARIES(LocProject1D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(LocProject1D)


SET(LocProject2DSource  LocProject2D.cpp)
ADD_NEKTAR_EXECUTABLE(LocProject2D demos LocProject2DSource)
TARGET_LINK_LIBRARIES(LocProject2D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(LocProject2D)

SET(LocProject3DSource  LocProject3D.cpp)
ADD_NEKTAR_EXECUTABLE(LocProject3D demos LocProject3DSource)
TARGET_LINK_LIBRARIES(LocProject3D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(LocProject3D)


SET(LocProject_Diff1DSource  LocProject_Diff1D.cpp)
ADD_NEKTAR_EXECUTABLE(LocProject_Diff1D demos LocProject_Diff1DSource)
TARGET_LINK_LIBRARIES(LocProject_Diff1D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(LocProject_Diff1D)

SET(LocProject_Diff2DSource  LocProject_Diff2D.cpp)
ADD_NEKTAR_EXECUTABLE(LocProject_Diff2D demos LocProject_Diff2DSource)
TARGET_LINK_LIBRARIES(LocProject_Diff2D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(LocProject_Diff2D)

SET(LocProject_Diff3DSource  LocProject_Diff3D.cpp)
ADD_NEKTAR_EXECUTABLE(LocProject_Diff3D demos LocProject_Diff3DSource)
TARGET_LINK_LIBRARIES(LocProject_Diff3D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(LocProject_Diff3D)

SET(LocTetExpDemoSource  LocTetExpDemo.cpp)
ADD_NEKTAR_EXECUTABLE(LocTetExpDemo demos LocTetExpDemoSource)
TARGET_LINK_LIBRARIES(LocTetExpDemo ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(LocTetExpDemo)

SET(LocPyrExpDemoSource  LocPyrExpDemo.cpp)
ADD_NEKTAR_EXECUTABLE(LocPyrExpDemo demos LocPyrExpDemoSource)
TARGET_LINK_LIBRARIES(LocPyrExpDemo ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(LocPyrExpDemo)

SET(LocPrismExpDemoSource  LocPrismExpDemo.cpp)
ADD_NEKTAR_EXECUTABLE(LocPrismExpDemo demos LocPrismExpDemoSource)
TARGET_LINK_LIBRARIES(LocPrismExpDemo ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(LocPrismExpDemo)

SET(LocHexExpDemoSource  LocHexExpDemo.cpp)
ADD_NEKTAR_EXECUTABLE(LocHexExpDemo demos LocHexExpDemoSource)
TARGET_LINK_LIBRARIES(LocHexExpDemo ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(LocHexExpDemo)

SET(LocHexExpDemoSource  LocHexExpDemo_2.cpp)
ADD_NEKTAR_EXECUTABLE(LocHexExpDemo_2 demos LocHexExpDemoSource)
TARGET_LINK_LIBRARIES(LocHexExpDemo_2 ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(LocHexExpDemo_2)

SET(LocPrismExpDemoSource  LocPrismExpDemo_2.cpp)
ADD_NEKTAR_EXECUTABLE(LocPrismExpDemo_2 demos LocPrismExpDemoSource)
TARGET_LINK_LIBRARIES(LocPrismExpDemo_2 ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(LocPrismExpDemo_2)

SET(LocPyrExpDemoSource  LocPyrExpDemo_2.cpp)
ADD_NEKTAR_EXECUTABLE(LocPyrExpDemo_2 demos LocPyrExpDemoSource)
TARGET_LINK_LIBRARIES(LocPyrExpDemo_2 ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(LocPyrExpDemo_2)

SET(LocTetExpDemoSource  LocTetExpDemo_2.cpp)
ADD_NEKTAR_EXECUTABLE(LocTetExpDemo_2 demos LocTetExpDemoSource)
TARGET_LINK_LIBRARIES(LocTetExpDemo_2 ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(LocTetExpDemo_2)

# Generate list of available subdirectories
FILE(GLOB dir_list "*")
FOREACH(dir ${dir_list})
    IF(IS_DIRECTORY ${dir} AND EXISTS ${dir}/CMakeLists.txt)
        ADD_SUBDIRECTORY(${dir})
    ENDIF(IS_DIRECTORY ${dir} AND EXISTS ${dir}/CMakeLists.txt)
ENDFOREACH(dir ${dir_list})
