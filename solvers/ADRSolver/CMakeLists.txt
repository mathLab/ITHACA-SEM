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
       ./EquationSystems/EigenValuesAdvection.cpp
       ../Auxiliary/Driver.cpp
       ../Auxiliary/DriverStandard.cpp
       ../Auxiliary/DriverArnoldi.cpp
       ../Auxiliary/DriverModifiedArnoldi.cpp
       ../Auxiliary/Filters/Filter.cpp
       ../Auxiliary/Filters/FilterThresholdMax.cpp
       ../Auxiliary/Filters/FilterCheckpoint.cpp
       ../Auxiliary/UnsteadySystem.cpp
       ../Auxiliary/EquationSystem.cpp)

IF (NEKTAR_USE_ARPACK)
    SET(ADRSolverSource ${ADRSolverSource}
               ../Auxiliary/DriverArpack.cpp)
ENDIF (NEKTAR_USE_ARPACK)

ADD_SOLVER_EXECUTABLE(ADRSolver solvers 
			${ADRSolverSource})
