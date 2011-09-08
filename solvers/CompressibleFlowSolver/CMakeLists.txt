SET(CompressibleFlowSolverSource 
       ./CompressibleFlowSolver.cpp 
       ./EquationSystems/CompressibleFlowSystem.cpp
       ./EquationSystems/EulerCFE.cpp
       ./EquationSystems/EulerArtificialDiffusionCFE.cpp
       ./EquationSystems/NavierStokesCFE.cpp
       ../Auxiliary/EquationSystem.cpp
       ./EquationSystems/UnsteadySystem.cpp)

ADD_SOLVER_EXECUTABLE(CompressibleFlowSolver solvers 
			${CompressibleFlowSolverSource})
