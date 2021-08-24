CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_IMAGE_WARPING
    "Build the Image Warping solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF (NOT NEKTAR_BUILD_SOLVERS)
    SET(NEKTAR_SOLVER_IMAGE_WARPING OFF CACHE INTERNAL "")
ENDIF()

IF( NEKTAR_SOLVER_IMAGE_WARPING )
    SET(ImageWarpingSolverSource 
       ./ImageWarpingSolver.cpp 
       ./EquationSystems/ImageWarpingSystem.cpp)
    ADD_SOLVER_EXECUTABLE(ImageWarpingSolver SOURCES ${ImageWarpingSolverSource})
ENDIF( NEKTAR_SOLVER_IMAGE_WARPING )
