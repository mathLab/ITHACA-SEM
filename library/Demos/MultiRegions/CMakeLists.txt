SET(LinkLibraries 
    optimized MultiRegions debug MultiRegions-g    
    optimized LocalRegions debug LocalRegions-g    
    optimized SpatialDomains debug SpatialDomains-g    
    optimized StdRegions debug StdRegions-g 
    optimized LibUtilities debug LibUtilities-g
    optimized ${Boost_THREAD_LIBRARY_RELEASE} debug ${Boost_THREAD_LIBRARY_DEBUG}
    optimized ${Boost_IOSTREAMS_LIBRARY_RELEASE} debug ${Boost_IOSTREAMS_LIBRARY_DEBUG}
    optimized ${Boost_ZLIB_LIBRARY_RELEASE} debug ${Boost_ZLIB_LIBRARY_DEBUG}
    optimized tinyxml debug tinyxml
)

SET(HelmholtzCont2DSource  Helmholtz2D.cpp)
ADD_NEKTAR_EXECUTABLE(Helmholtz2D HelmholtzCont2DSource)
TARGET_LINK_LIBRARIES(Helmholtz2D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(Helmholtz2D)

SET(HDGHelmholtz2DSource  HDGHelmholtz2D.cpp)
ADD_NEKTAR_EXECUTABLE(HDGHelmholtz2D HDGHelmholtz2DSource)
TARGET_LINK_LIBRARIES(HDGHelmholtz2D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(HDGHelmholtz2D)

SET(LaplaceCont2DSource  Laplace2D.cpp)
ADD_NEKTAR_EXECUTABLE(Laplace2D LaplaceCont2DSource)
TARGET_LINK_LIBRARIES(Laplace2D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(Laplace2D)
