CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_IMAGE_WARPING
    "Build the Image Warping solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF( NEKTAR_SOLVER_IMAGE_WARPING )
    SET(ImageWarpingSolverSource 
       ./ImageWarpingSolver.cpp 
       ./EquationSystems/ImageWarpingSystem.cpp)

    ADD_SOLVER_EXECUTABLE(ImageWarpingSolver solvers 
			${ImageWarpingSolverSource})
ENDIF( NEKTAR_SOLVER_IMAGE_WARPING )
