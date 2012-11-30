SET(NEKTAR_SOLVER_IMAGE_WARPING ON CACHE BOOL
    "Build the Image Warping solver.")

IF( NEKTAR_SOLVER_IMAGE_WARPING )
    SET(ImageWarpingSolverSource 
       ./ImageWarpingSolver.cpp 
       ./EquationSystems/ImageWarpingSystem.cpp)

    ADD_SOLVER_EXECUTABLE(ImageWarpingSolver solvers 
			${ImageWarpingSolverSource})
ENDIF( NEKTAR_SOLVER_IMAGE_WARPING )
