SET(NekMeshSources
    InputModules/InputGmsh.cpp
    InputModules/InputNek.cpp
    InputModules/InputNek5000.cpp
    InputModules/InputNekpp.cpp
    InputModules/InputPly.cpp
    InputModules/InputSem.cpp
    InputModules/InputSwan.cpp
    InputModules/InputStarTec.cpp
    NekMesh.cpp
    OutputModules/OutputGmsh.cpp
    OutputModules/OutputNekpp.cpp
    OutputModules/OutputSTL.cpp
    ProcessModules/ProcessBL.cpp
    ProcessModules/ProcessCurve.cpp
    ProcessModules/ProcessCurvedEdges.cpp
    ProcessModules/ProcessCyl.cpp
    ProcessModules/ProcessDetectSurf.cpp
    ProcessModules/ProcessExtractSurf.cpp
    ProcessModules/ProcessExtractTetPrismInterface.cpp
    ProcessModules/ProcessJac.cpp
    ProcessModules/ProcessLinkCheck.cpp
    ProcessModules/ProcessLinear.cpp
    ProcessModules/ProcessPerAlign.cpp
    ProcessModules/ProcessScalar.cpp
    ProcessModules/ProcessSpherigon.cpp
    ProcessModules/ProcessTetSplit.cpp
    ProcessModules/ProcessOptiExtract.cpp
    ProcessModules/ProcessInsertSurface.cpp
    ProcessModules/ProcessExtrude.cpp
    ProcessModules/ProcessVarOpti/ProcessVarOpti.cpp
    ProcessModules/ProcessVarOpti/PreProcessing.cpp
    ProcessModules/ProcessVarOpti/NodeOpti.cpp
    ProcessModules/ProcessVarOpti/ElUtil.cpp
)

IF (NEKTAR_USE_CFI)
    ADD_DEFINITIONS(-DNEKTAR_USE_CFI)
    SET(NekMeshSources ${NekMeshSources} InputModules/InputCADfix.cpp)
    SET(NekMeshSources ${NekMeshSources} OutputModules/OutputCADfix.cpp)
ENDIF (NEKTAR_USE_CFI)

IF (NEKTAR_USE_CCM)
    SET(NekMeshSources ${NekMeshSources} InputModules/InputStar.cpp)
ENDIF (NEKTAR_USE_CCM)

IF (NEKTAR_USE_VTK)
    SET(NekMeshSources ${NekMeshSources} InputModules/InputVtk.cpp OutputModules/OutputVtk.cpp)
ENDIF (NEKTAR_USE_VTK)

IF (NEKTAR_USE_MESHGEN)
    SET(NekMeshSources ${NekMeshSources}
            ProcessModules/ProcessVarOpti/NodeOptiCAD.cpp
            ProcessModules/ProcessProjectCAD.cpp
            InputModules/InputMCF.cpp)
ENDIF (NEKTAR_USE_MESHGEN)

# Don't use macro -- this way we can just link against NekMeshUtils and
# LocalRegions.
ADD_UTILITIES_EXECUTABLE(NekMesh COMPONENT nekmesh
    SOURCES ${NekMeshSources}
    DEPENDS NekMeshUtils LocalRegions)

IF (NEKTAR_USE_CCM)
   TARGET_LINK_LIBRARIES(NekMesh LINK_PRIVATE ${CCMIO_LIBRARIES})
   ADD_DEPENDENCIES(NekMesh libccmio-2.6.1)
ENDIF (NEKTAR_USE_CCM)

IF (NEKTAR_USE_VTK)
    IF (VTK_MAJOR_VERSION LESS 6)
        TARGET_LINK_LIBRARIES(NekMesh LINK_PUBLIC vtkCommon vtkIO vtkGraphics)
    ELSE ()
        TARGET_LINK_LIBRARIES(NekMesh LINK_PUBLIC vtkCommonCore vtkIOLegacy)
    ENDIF ()
ENDIF (NEKTAR_USE_VTK)

# Nektar++
ADD_NEKTAR_TEST(Nektar++/bl_quad)
ADD_NEKTAR_TEST(Nektar++/bl_tube)
ADD_NEKTAR_TEST(Nektar++/extract_curved_edge)
ADD_NEKTAR_TEST(Nektar++/extract_curved_face)
ADD_NEKTAR_TEST(Nektar++/extract_detectbnd_tube)
ADD_NEKTAR_TEST(Nektar++/extrude)
ADD_NEKTAR_TEST(Nektar++/jac_list_tet_face)
ADD_NEKTAR_TEST(Nektar++/linearise_invalid_quad)

# MeshGen test
IF(NEKTAR_USE_MESHGEN)
    ADD_NEKTAR_TEST(MeshGen/GEO/2d_adapt)
    ADD_NEKTAR_TEST(MeshGen/GEO/2d_bl)
    ADD_NEKTAR_TEST(MeshGen/GEO/2d_t106c)
    ADD_NEKTAR_TEST(MeshGen/GEO/annulus)

    ADD_NEKTAR_TEST(MeshGen/NACA/2d_bl)

    ADD_NEKTAR_TEST(MeshGen/STEP/2d_bl_aerofoil)
    ADD_NEKTAR_TEST(MeshGen/STEP/2d_bl_t106)
    ADD_NEKTAR_TEST(MeshGen/STEP/3d_bl_cyl)
    ADD_NEKTAR_TEST(MeshGen/STEP/3d_bl_wing LENGTHY)
    ADD_NEKTAR_TEST(MeshGen/STEP/3d_sphere LENGTHY)

    IF(NEKTAR_USE_CCM)
        ADD_NEKTAR_TEST(MeshGen/StarCCM/projectcad_ahmed)
    ENDIF()
ENDIF()

# Nektar tests
ADD_NEKTAR_TEST(Nektar/bl_hex_lin)
ADD_NEKTAR_TEST(Nektar/gz_uk)
ADD_NEKTAR_TEST(Nektar/tetsplit_tube)

# Gmsh tests
ADD_NEKTAR_TEST(Gmsh/v2.2/cube_all)
ADD_NEKTAR_TEST(Gmsh/v2.2/cube_hex_lin)
ADD_NEKTAR_TEST(Gmsh/v2.2/cube_hex)
ADD_NEKTAR_TEST(Gmsh/v2.2/cube_prism_lin)
ADD_NEKTAR_TEST(Gmsh/v2.2/cube_prism)
ADD_NEKTAR_TEST(Gmsh/v2.2/cube_tet_lin LENGTHY)
ADD_NEKTAR_TEST(Gmsh/v2.2/cube_tet)
ADD_NEKTAR_TEST(Gmsh/v2.2/square_quad_lin)
ADD_NEKTAR_TEST(Gmsh/v2.2/square_quad)
ADD_NEKTAR_TEST(Gmsh/v2.2/square_tri_lin)
ADD_NEKTAR_TEST(Gmsh/v2.2/square_tri)

ADD_NEKTAR_TEST(Gmsh/v4.0/cube_all)
ADD_NEKTAR_TEST(Gmsh/v4.0/cube_hex_lin)
ADD_NEKTAR_TEST(Gmsh/v4.0/cube_hex)
ADD_NEKTAR_TEST(Gmsh/v4.0/cube_prism_lin)
ADD_NEKTAR_TEST(Gmsh/v4.0/cube_prism)
ADD_NEKTAR_TEST(Gmsh/v4.0/cube_tet_lin LENGTHY)
ADD_NEKTAR_TEST(Gmsh/v4.0/cube_tet)
ADD_NEKTAR_TEST(Gmsh/v4.0/extract_square_circ)
ADD_NEKTAR_TEST(Gmsh/v4.0/square_quad_lin)
ADD_NEKTAR_TEST(Gmsh/v4.0/square_quad)
ADD_NEKTAR_TEST(Gmsh/v4.0/square_tri_lin)
ADD_NEKTAR_TEST(Gmsh/v4.0/square_tri)

ADD_NEKTAR_TEST(Gmsh/v4.1/cube_all)
ADD_NEKTAR_TEST(Gmsh/v4.1/cube_hex_lin)
ADD_NEKTAR_TEST(Gmsh/v4.1/cube_hex)
ADD_NEKTAR_TEST(Gmsh/v4.1/cube_prism_lin)
ADD_NEKTAR_TEST(Gmsh/v4.1/cube_prism)
ADD_NEKTAR_TEST(Gmsh/v4.1/cube_tet_lin LENGTHY)
ADD_NEKTAR_TEST(Gmsh/v4.1/cube_tet)
ADD_NEKTAR_TEST(Gmsh/v4.1/extract_square_circ)
ADD_NEKTAR_TEST(Gmsh/v4.1/square_quad_lin)
ADD_NEKTAR_TEST(Gmsh/v4.1/square_quad)
ADD_NEKTAR_TEST(Gmsh/v4.1/square_tri_lin)
ADD_NEKTAR_TEST(Gmsh/v4.1/square_tri)

# Rot peralign tests
ADD_NEKTAR_TEST(Gmsh/v2.2/peralign_rot_cyl)

# VarOpti tests
ADD_NEKTAR_TEST(Gmsh/v2.2/varopti_cube_sphere)
ADD_NEKTAR_TEST(Gmsh/v2.2/varopti_naca)

# HDF5 tests
IF (NEKTAR_USE_HDF5)
    ADD_NEKTAR_TEST(Gmsh/v4.0/hdf5_cube_tet)
ENDIF()

# StarCCM tests
ADD_NEKTAR_TEST(StarTec/peralign_bl_cube)
ADD_NEKTAR_TEST(StarTec/spherigon_bl_straight_rw LENGTHY)

# Nek5000 tests
ADD_NEKTAR_TEST(Nek5000/annulus)
ADD_NEKTAR_TEST(Nek5000/box)
