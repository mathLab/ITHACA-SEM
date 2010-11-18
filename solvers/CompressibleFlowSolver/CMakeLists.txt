SET(CompressibleFlowSolverSource 
       ./CompressibleFlowSolver.cpp 
       ../ADRSolver/EquationSystem.cpp
       ../ADRSolver/SessionReader.cpp
       ./EquationSystems/CompressibleFlowSystem.cpp
       	./EquationSystems/EulerCFE.cpp	
       ../Auxiliary/ADRBase.cpp)

ADD_SOLVER_EXECUTABLE(CompressibleFlowSolver 
			${CompressibleFlowSolverSource})
