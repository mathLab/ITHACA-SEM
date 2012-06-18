SET(CompressibleFlowSolverSource 
       ./CompressibleFlowSolver.cpp 
       ./EquationSystems/CompressibleFlowSystem.cpp
       ./EquationSystems/EulerCFE.cpp
       ./EquationSystems/EulerArtificialDiffusionCFE.cpp
       ./EquationSystems/NavierStokesCFE.cpp
       ./EquationSystems/UnsteadySystem.cpp
       ./RiemannSolvers/CompressibleSolver.cpp
       ./RiemannSolvers/AverageSolver.cpp
       ./RiemannSolvers/ExactSolver.cpp)

ADD_SOLVER_EXECUTABLE(CompressibleFlowSolver solvers 
			${CompressibleFlowSolverSource})
