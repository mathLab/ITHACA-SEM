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
    optimized sbtk    debug sbtk
)

SET(HelmholtzCont1DSource  Helmholtz1D.cpp)
ADD_NEKTAR_EXECUTABLE(Helmholtz1D HelmholtzCont1DSource)
TARGET_LINK_LIBRARIES(Helmholtz1D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(Helmholtz1D)

SET(HelmholtzCont2DSource  Helmholtz2D.cpp)
ADD_NEKTAR_EXECUTABLE(Helmholtz2D HelmholtzCont2DSource)
TARGET_LINK_LIBRARIES(Helmholtz2D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(Helmholtz2D)

SET(HelmholtzCont3DSource  Helmholtz3D.cpp)
ADD_NEKTAR_EXECUTABLE(Helmholtz3D HelmholtzCont3DSource)
TARGET_LINK_LIBRARIES(Helmholtz3D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(Helmholtz3D)

SET(Helmholtz3DHomo1DSource  Helmholtz3DHomo1D.cpp)
ADD_NEKTAR_EXECUTABLE(Helmholtz3DHomo1D Helmholtz3DHomo1DSource)
TARGET_LINK_LIBRARIES(Helmholtz3DHomo1D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(Helmholtz3DHomo1D)

SET(HDGHelmholtz1DSource  HDGHelmholtz1D.cpp)
ADD_NEKTAR_EXECUTABLE(HDGHelmholtz1D HDGHelmholtz1DSource)
TARGET_LINK_LIBRARIES(HDGHelmholtz1D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(HDGHelmholtz1D)

SET(HDGHelmholtz2DSource  HDGHelmholtz2D.cpp)
ADD_NEKTAR_EXECUTABLE(HDGHelmholtz2D HDGHelmholtz2DSource)
TARGET_LINK_LIBRARIES(HDGHelmholtz2D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(HDGHelmholtz2D)

SET(SteadyAdvectionDiffusionReactionCont2DSource  SteadyAdvectionDiffusionReaction2D.cpp)
ADD_NEKTAR_EXECUTABLE(SteadyAdvectionDiffusionReaction2D SteadyAdvectionDiffusionReactionCont2DSource)
TARGET_LINK_LIBRARIES(SteadyAdvectionDiffusionReaction2D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(SteadyAdvectionDiffusionReaction2D)

# Uncomment next line to compile demos in subdirectory
#SUBDIRS(ExtraDemos)
