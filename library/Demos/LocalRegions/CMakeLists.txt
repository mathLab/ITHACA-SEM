SET(LinkLibraries 
    optimized LocalRegions debug LocalRegions-g  
    optimized SpatialDomains debug SpatialDomains-g
    optimized StdRegions debug StdRegions-g
    optimized LibUtilities debug LibUtilities-g
    optimized ${Boost_THREAD_LIBRARY_RELEASE} debug ${Boost_THREAD_LIBRARY_DEBUG}
    optimized ${Boost_IOSTREAMS_LIBRARY_RELEASE} debug ${Boost_IOSTREAMS_LIBRARY_DEBUG}
    optimized ${Boost_ZLIB_LIBRARY_RELEASE} debug ${Boost_ZLIB_LIBRARY_DEBUG}
    optimized tinyxml debug tinyxml
)

SET(LocProject1DSource  Project1D.cpp)
ADD_NEKTAR_EXECUTABLE(LocProject1D LocProject1DSource)
TARGET_LINK_LIBRARIES(LocProject1D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(LocProject1D)


SET(LocProject2DSource  Project2D.cpp)
ADD_NEKTAR_EXECUTABLE(LocProject2D LocProject2DSource)
TARGET_LINK_LIBRARIES(LocProject2D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(LocProject2D)


SET(LocProject_Diff1DSource  Project_Diff1D.cpp)
ADD_NEKTAR_EXECUTABLE(LocProject_Diff1D LocProject_Diff1DSource)
TARGET_LINK_LIBRARIES(LocProject_Diff1D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(LocProject_Diff1D)

SET(LocTetExpDemoSource  LocTetExpDemo.cpp)
ADD_NEKTAR_EXECUTABLE(LocTetExpDemo LocTetExpDemoSource)
TARGET_LINK_LIBRARIES(LocTetExpDemo ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(LocTetExpDemo)

SET(LocPyrExpDemoSource  LocPyrExpDemo.cpp)
ADD_NEKTAR_EXECUTABLE(LocPyrExpDemo LocPyrExpDemoSource)
TARGET_LINK_LIBRARIES(LocPyrExpDemo ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(LocPyrExpDemo)

SET(LocPrismExpDemoSource  LocPrismExpDemo.cpp)
ADD_NEKTAR_EXECUTABLE(LocPrismExpDemo LocPrismExpDemoSource)
TARGET_LINK_LIBRARIES(LocPrismExpDemo ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(LocPrismExpDemo)

SET(LocHexExpDemoSource  LocHexExpDemo.cpp)
ADD_NEKTAR_EXECUTABLE(LocHexExpDemo LocHexExpDemoSource)
TARGET_LINK_LIBRARIES(LocHexExpDemo ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(LocHexExpDemo)
