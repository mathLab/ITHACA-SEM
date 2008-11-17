SET(BoussinesqSolverSource 
       ./BoussinesqSolver.cpp 
       ./BoussinesqEquations.cpp
       ./StreamFunctionWaves.cpp
       ../ShallowWaterSolver/ShallowWaterEquations.cpp
       ../Auxiliary/ADRBase.cpp)

ADD_SOLVER_EXECUTABLE(BoussinesqSolver 
			${BoussinesqSolverSource})
