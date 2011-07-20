SET(CompressibleFlowSolverSource 
       ./CompressibleFlowSolver.cpp 
       ./EquationSystems/CompressibleFlowSystem.cpp
       ./EquationSystems/EulerCFE.cpp
       ./EquationSystems/EulerArtificialDiffusionCFE.cpp
       ./EquationSystems/NavierStokesCFE.cpp
       ../Auxiliary/EquationSystem.cpp
       ./EquationSystems/UnsteadySystem.cpp
       ../Auxiliary/ADRBase.cpp)

ADD_SOLVER_EXECUTABLE(CompressibleFlowSolver 
			${CompressibleFlowSolverSource})
