CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_IMAGE_WARPING
    "Build the Image Warping solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF( NEKTAR_SOLVER_IMAGE_WARPING )
    SET(ImageWarpingSolverSource 
       ./ImageWarpingSolver.cpp 
       ./EquationSystems/ImageWarpingSystem.cpp)
    ADD_SOLVER_EXECUTABLE(ImageWarpingSolver SOURCES ${ImageWarpingSolverSource})
    FINALISE_CPACK_COMPONENT(imagewarping-solver
        SUMMARY "Nektar++ Image Warping Solver"
        DESCRIPTION "Solves model image-warping equations.")
ENDIF( NEKTAR_SOLVER_IMAGE_WARPING )
