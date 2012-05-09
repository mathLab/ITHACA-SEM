SET(ADRSolverSource 
       ./ADRSolver.cpp 
       ./EquationSystems/Laplace.cpp
       ./EquationSystems/Poisson.cpp       
       ./EquationSystems/Helmholtz.cpp
       ./EquationSystems/SteadyAdvectionDiffusion.cpp
       ./EquationSystems/SteadyAdvectionDiffusionReaction.cpp
       ./EquationSystems/UnsteadyAdvection.cpp
       ./EquationSystems/UnsteadyDiffusion.cpp
       ./EquationSystems/UnsteadyAdvectionDiffusion.cpp
       ./EquationSystems/UnsteadyInviscidBurger.cpp
       ./EquationSystems/CFLtester.cpp
       ./EquationSystems/EigenValuesAdvection.cpp)

ADD_SOLVER_EXECUTABLE(ADRSolver solvers 
			${ADRSolverSource})
