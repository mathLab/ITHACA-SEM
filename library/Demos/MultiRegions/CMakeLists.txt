SET(LinkLibraries 
    optimized MultiRegions debug MultiRegions-g    
    optimized LocalRegions debug LocalRegions-g    
    optimized SpatialDomains debug SpatialDomains-g    
    optimized StdRegions debug StdRegions-g 
    optimized LibUtilities debug LibUtilities-g
    optimized ${Boost_THREAD_LIBRARY_RELEASE} debug ${Boost_THREAD_LIBRARY_DEBUG}
    optimized ${Boost_IOSTREAMS_LIBRARY_RELEASE} debug ${Boost_IOSTREAMS_LIBRARY_DEBUG}
    optimized ${ZLIB_LIBRARY_RELEASE} debug ${ZLIB_LIBRARY_DEBUG}
    optimized ${TINYXML_LIB} debug ${TINYXML_LIB}
    optimized sbtk    debug sbtk
)

SET(HelmholtzCont1DSource  Helmholtz1D.cpp)
ADD_NEKTAR_EXECUTABLE(Helmholtz1D demos HelmholtzCont1DSource)
TARGET_LINK_LIBRARIES(Helmholtz1D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(Helmholtz1D)

SET(HelmholtzCont2DSource  Helmholtz2D.cpp)
ADD_NEKTAR_EXECUTABLE(Helmholtz2D demos HelmholtzCont2DSource)
TARGET_LINK_LIBRARIES(Helmholtz2D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(Helmholtz2D)

#SET(Helmholtz2DHomo1DSource  Helmholtz2DHomo1D.cpp)
#ADD_NEKTAR_EXECUTABLE(Helmholtz2DHomo1D Helmholtz2DHomo1DSource)
#TARGET_LINK_LIBRARIES(Helmholtz2DHomo1D ${LinkLibraries})
#SET_LAPACK_LINK_LIBRARIES(Helmholtz2DHomo1D)

SET(HelmholtzCont3DSource  Helmholtz3D.cpp)
ADD_NEKTAR_EXECUTABLE(Helmholtz3D demos HelmholtzCont3DSource)
TARGET_LINK_LIBRARIES(Helmholtz3D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(Helmholtz3D)

SET(Helmholtz3DHomo1DSource  Helmholtz3DHomo1D.cpp)
ADD_NEKTAR_EXECUTABLE(Helmholtz3DHomo1D demos Helmholtz3DHomo1DSource)
TARGET_LINK_LIBRARIES(Helmholtz3DHomo1D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(Helmholtz3DHomo1D)

SET(Helmholtz3DHomo2DSource  Helmholtz3DHomo2D.cpp)
ADD_NEKTAR_EXECUTABLE(Helmholtz3DHomo2D demos Helmholtz3DHomo2DSource)
TARGET_LINK_LIBRARIES(Helmholtz3DHomo2D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(Helmholtz3DHomo2D)

SET(HDGHelmholtz1DSource  HDGHelmholtz1D.cpp)
ADD_NEKTAR_EXECUTABLE(HDGHelmholtz1D demos HDGHelmholtz1DSource)
TARGET_LINK_LIBRARIES(HDGHelmholtz1D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(HDGHelmholtz1D)

SET(HDGHelmholtz2DSource  HDGHelmholtz2D.cpp)
ADD_NEKTAR_EXECUTABLE(HDGHelmholtz2D demos HDGHelmholtz2DSource)
TARGET_LINK_LIBRARIES(HDGHelmholtz2D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(HDGHelmholtz2D)

SET(HDGHelmholtz3DSource  HDGHelmholtz3D.cpp)
ADD_NEKTAR_EXECUTABLE(HDGHelmholtz3D demos HDGHelmholtz3DSource)
TARGET_LINK_LIBRARIES(HDGHelmholtz3D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(HDGHelmholtz3D)

SET(HDGHelmholtz3DHomo1DSource  HDGHelmholtz3DHomo1D.cpp)
ADD_NEKTAR_EXECUTABLE(HDGHelmholtz3DHomo1D demos HDGHelmholtz3DHomo1DSource)
TARGET_LINK_LIBRARIES(HDGHelmholtz3DHomo1D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(HDGHelmholtz3DHomo1D)

SET(Deriv3DHomo1D_SingleModeSource  Deriv3DHomo1D_SingleMode.cpp)
ADD_NEKTAR_EXECUTABLE(Deriv3DHomo1D_SingleMode demos Deriv3DHomo1D_SingleModeSource)
TARGET_LINK_LIBRARIES(Deriv3DHomo1D_SingleMode ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(Deriv3DHomo1D_SingleMode)

SET(Deriv3DHomo1DSource  Deriv3DHomo1D.cpp)
ADD_NEKTAR_EXECUTABLE(Deriv3DHomo1D demos Deriv3DHomo1DSource)
TARGET_LINK_LIBRARIES(Deriv3DHomo1D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(Deriv3DHomo1D)

SET(Deriv3DHomo2DSource  Deriv3DHomo2D.cpp)
ADD_NEKTAR_EXECUTABLE(Deriv3DHomo2D demos Deriv3DHomo2DSource)
TARGET_LINK_LIBRARIES(Deriv3DHomo2D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(Deriv3DHomo2D)

SET(SteadyAdvectionDiffusionReactionCont2DSource  SteadyAdvectionDiffusionReaction2D.cpp)
ADD_NEKTAR_EXECUTABLE(SteadyAdvectionDiffusionReaction2D demos SteadyAdvectionDiffusionReactionCont2DSource)
TARGET_LINK_LIBRARIES(SteadyAdvectionDiffusionReaction2D ${LinkLibraries})
SET_LAPACK_LINK_LIBRARIES(SteadyAdvectionDiffusionReaction2D)

# Generate list of available subdirectories
FILE(GLOB dir_list "*")
FOREACH(dir ${dir_list})
    IF(IS_DIRECTORY ${dir} AND EXISTS ${dir}/CMakeLists.txt)
        ADD_SUBDIRECTORY(${dir})
    ENDIF(IS_DIRECTORY ${dir} AND EXISTS ${dir}/CMakeLists.txt)
ENDFOREACH(dir ${dir_list})

ADD_NEKTAR_TEST(Helmholtz1D_CG_P8)
ADD_NEKTAR_TEST(Helmholtz1D_CG_P8_RBC)
ADD_NEKTAR_TEST(Helmholtz1D_HDG_P8)
ADD_NEKTAR_TEST(Helmholtz1D_HDG_P8_RBC)

ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes)
ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Nodes)
ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_AllBCs_mlsc)
ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_AllBCs_sc)
#ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_AllBCs_full)
ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_AllBCs_iter_sc)
ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_AllBCs_iter_ml)
ADD_NEKTAR_TEST(Helmholtz2D_CG_P9_Modes_varcoeff)
ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_curved_quad)
ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_curved_tri)
ADD_NEKTAR_TEST(Helmholtz2D_HDG_P7_Modes)
ADD_NEKTAR_TEST(Helmholtz2D_HDG_P7_Modes_AllBCs)

ADD_NEKTAR_TEST(Helmholtz3D_CG_Hex)
ADD_NEKTAR_TEST(Helmholtz3D_CG_Hex_AllBCs)
ADD_NEKTAR_TEST(Helmholtz3D_CG_Hex_AllBCs_iter_ml)
ADD_NEKTAR_TEST(Helmholtz3D_CG_Hex_AllBCs_iter_sc_cont)
ADD_NEKTAR_TEST(Helmholtz3D_CG_Tet)
ADD_NEKTAR_TEST(Helmholtz3D_CG_Prism)
ADD_NEKTAR_TEST(Helmholtz3D_CG_Prism_iter_ml)
ADD_NEKTAR_TEST(Helmholtz3D_CG_Prism_Deformed)
ADD_NEKTAR_TEST(Helmholtz3D_CG_Homo1D)
ADD_NEKTAR_TEST(Helmholtz3D_HDG_Homo1D)
ADD_NEKTAR_TEST(Helmholtz3D_HDG_Prism)
ADD_NEKTAR_TEST(Helmholtz3D_HDG_Hex_AllBCs)
ADD_NEKTAR_TEST(Helmholtz3D_HDG_Tet)

ADD_NEKTAR_TEST(LinearAdvDiffReact2D_P7_Modes)
ADD_NEKTAR_TEST(Deriv3D_Homo1D)
ADD_NEKTAR_TEST(Deriv3D_Homo2D)

IF (NEKTAR_USE_MPI)
    ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_AllBCs_xxt_full)
    ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_AllBCs_xxt_sc)
    ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_AllBCs_iter_sc_par3)
    ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_AllBCs_iter_ml_par3)
    ADD_NEKTAR_TEST(Helmholtz3D_CG_Hex_AllBCs_iter_ml_par3)
    ADD_NEKTAR_TEST(Helmholtz3D_CG_Prism_iter_ml_par3)
    ADD_NEKTAR_TEST(Helmholtz3D_CG_Hex_AllBCs_xxt_sc_par3)
    ADD_NEKTAR_TEST(Helmholtz3D_CG_Tet_iter_global_par3)
    ADD_NEKTAR_TEST(Helmholtz2D_HDG_P7_Modes_AllBCs_par2)
    ADD_NEKTAR_TEST(Helmholtz3D_HDG_Prism_par2)
    ADD_NEKTAR_TEST(Helmholtz3D_HDG_Hex_AllBCs_par2)
ENDIF (NEKTAR_USE_MPI)
