SET(FieldConvertSources FieldConvert.cpp)

ADD_UTILITIES_EXECUTABLE(FieldConvert util ${FieldConvertSources})
TARGET_LINK_LIBRARIES(FieldConvert FieldUtils)

#ADD_NEKTAR_TEST(chan3D_tec)  # This test is covered in bfs_tec
#ADD_NEKTAR_TEST(chan3D_vtu)  # This does not currently work -- hash metric
#ADD_NEKTAR_TEST(chan3D_vort) # This is covered by bfs_vort
#ADD_NEKTAR_TEST(chan3D_pts)  # This does not currently work -- hash metric

ADD_NEKTAR_TEST(chan3D_tec_n10)
ADD_NEKTAR_TEST(chan3D_interppointsplane)
ADD_NEKTAR_TEST(chan3D_interppointsbox)
ADD_NEKTAR_TEST(compositeid)
ADD_NEKTAR_TEST(bfs_probe)
ADD_NEKTAR_TEST(bfs_tec)
ADD_NEKTAR_TEST(bfs_tec_rng)
ADD_NEKTAR_TEST(bfs_vort)
ADD_NEKTAR_TEST(bfs_vort_rng)
ADD_NEKTAR_TEST(pointdatatofld)
# ADD_NEKTAR_TEST(chan3D_pts)
ADD_NEKTAR_TEST(chan3DH1D_meanmode)
ADD_NEKTAR_TEST(chan3DH1D_plane)
ADD_NEKTAR_TEST(chan3DH1D_stretch)
ADD_NEKTAR_TEST(chan3D_probe)
ADD_NEKTAR_TEST(cube_prismhex)
ADD_NEKTAR_TEST(outflow_pointdatatofld)
ADD_NEKTAR_TEST(chan3D_equispacedoutput)
ADD_NEKTAR_TEST(chan3D_isocontour)
ADD_NEKTAR_TEST(interpfield)
ADD_NEKTAR_TEST(naca0012_bnd_equispacedoutput)
ADD_NEKTAR_TEST(smallmesh_isocontour)

IF (NEKTAR_USE_MPI)
    ADD_NEKTAR_TEST(chan3D_tec_par)
    #ADD_NEKTAR_TEST(chan3D_vtu_par)
    ADD_NEKTAR_TEST(chan3D_interppointsplane_par)
    ADD_NEKTAR_TEST(chan3D_interppointsbox_par)
    ADD_NEKTAR_TEST(chan3D_vort_par)
ENDIF (NEKTAR_USE_MPI)
