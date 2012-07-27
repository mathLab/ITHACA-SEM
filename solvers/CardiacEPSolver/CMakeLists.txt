SET(CardiacEPSolverSource 
       ./CardiacEPSolver.cpp
       ./EquationSystems/Monodomain.cpp
       ./EquationSystems/Bidomain.cpp
       ./CellModels/CellModel.cpp
       ./CellModels/FitzhughNagumo.cpp
       ./CellModels/AlievPanfilov.cpp
       ./CellModels/CourtemancheRamirezNattel98.cpp
       ./CellModels/LuoRudy91.cpp
       ./CellModels/Fox02.cpp
       ./CellModels/Winslow99.cpp
       ./CellModels/TenTusscher06Epi.cpp
       ./CellModels/TenTusscher06M.cpp
       ./CellModels/TenTusscher06Endo.cpp
       ./Filters/FilterCheckpointCellModel.cpp
)

ADD_SOLVER_EXECUTABLE(CardiacEPSolver solvers-extra 
			${CardiacEPSolverSource})
