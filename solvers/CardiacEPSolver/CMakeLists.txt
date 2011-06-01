SET(CardiacEPSolverSource 
       ./CardiacEPSolver.cpp
       ../ADRSolver/EquationSystems/UnsteadySystem.cpp
       ../ADRSolver/EquationSystem.cpp
       ./EquationSystems/Monodomain.cpp
       ./CellModels/CellModel.cpp
       ./CellModels/FitzhughNagumo.cpp
       ./CellModels/AlievPanfilov.cpp
       ../Auxiliary/ADRBase.cpp)

ADD_SOLVER_EXECUTABLE(CardiacEPSolver 
			${CardiacEPSolverSource})
