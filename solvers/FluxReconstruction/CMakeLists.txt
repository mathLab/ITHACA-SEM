SET(UnsteadyAdvection1DFRSources
        ./UnsteadyAdvection1DFRSolver.cpp
	./Advection1DFR.cpp)

ADD_SOLVER_EXECUTABLE(UnsteadyAdvection1DFRSolver solvers 
			${UnsteadyAdvection1DFRSources})


