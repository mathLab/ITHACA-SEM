SET(CardiacEPSolverSource 
       ./CardiacEPSolver.cpp
       ./EquationSystems/Monodomain.cpp
       ./CellModels/CellModel.cpp
       ./CellModels/FitzhughNagumo.cpp
       ./CellModels/AlievPanfilov.cpp
       ../Auxiliary/UnsteadySystem.cpp
       ../Auxiliary/EquationSystem.cpp
       ../Auxiliary/ADRBase.cpp)

ADD_SOLVER_EXECUTABLE(CardiacEPSolver 
			${CardiacEPSolverSource})
