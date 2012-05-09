SET(ImageWarpingSolverSource 
       ./ImageWarpingSolver.cpp 
       ./EquationSystems/ImageWarpingSystem.cpp)

ADD_SOLVER_EXECUTABLE(ImageWarpingSolver solvers 
			${ImageWarpingSolverSource})
