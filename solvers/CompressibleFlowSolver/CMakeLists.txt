SET(CompressibleFlowSolverSource 
       ./CompressibleFlowSolver.cpp 
       ../ADRSolver/EquationSystem.cpp
       ./EquationSystems/UnsteadySystem.cpp
       ./EquationSystems/CompressibleFlowSystem.cpp
       ./EquationSystems/EulerCFE.cpp
       ./EquationSystems/EulerArtificialDiffusionCFE.cpp
       ./EquationSystems/NavierStokesCFE.cpp
       ../Auxiliary/ADRBase.cpp)

ADD_SOLVER_EXECUTABLE(CompressibleFlowSolver 
			${CompressibleFlowSolverSource})
