SET(ADRSolverSource 
       ./ADRSolver.cpp 
       ./EquationSystem.cpp
       ./EquationSystems/Laplace.cpp
       ./EquationSystems/Poisson.cpp
       ./EquationSystems/Helmholtz.cpp
       ./EquationSystems/UnsteadyDiffusion.cpp
       ./EquationSystems/UnsteadyAdvection.cpp
       ./EquationSystems/UnsteadyAdvectionDiffusion.cpp
       ./EquationSystems/UnsteadyInviscidBurger.cpp
       ./EquationSystems/AlievPanfilov.cpp
       ./SessionReader.cpp
       ../Auxiliary/ADRBase.cpp)

ADD_SOLVER_EXECUTABLE(ADRSolver 
			${ADRSolverSource})
