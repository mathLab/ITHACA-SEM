SET(CardiacEPSolverSource 
       ./CardiacEPSolver.cpp
       ./EquationSystems/Monodomain.cpp
       ./CellModels/CellModel.cpp
       ./CellModels/FitzhughNagumo.cpp
       ./CellModels/AlievPanfilov.cpp
       ../Auxiliary/UnsteadySystem.cpp
       ../Auxiliary/EquationSystem.cpp
       ../Auxiliary/Driver.cpp
       ../Auxiliary/DriverStandard.cpp)

ADD_SOLVER_EXECUTABLE(CardiacEPSolver solvers-extra 
			${CardiacEPSolverSource})
