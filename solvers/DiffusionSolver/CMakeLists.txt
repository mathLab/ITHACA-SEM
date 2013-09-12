CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_DIFFUSION
    "Build the Diffusion Solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF( NEKTAR_SOLVER_DIFFUSION )
    SET(DiffusionTestSource 
       ./DiffusionTest.cpp)
    SET(DiffusionTestTISource
        ./DiffusionTestTI.cpp)
        
    ADD_SOLVER_EXECUTABLE(DiffusionTest solvers 
			${DiffusionTestSource})
    ADD_SOLVER_EXECUTABLE(DiffusionTestTI solvers 
            ${DiffusionTestTISource})

ENDIF( NEKTAR_SOLVER_DIFFUSION )
