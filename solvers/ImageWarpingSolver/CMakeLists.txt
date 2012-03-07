SET(ImageWarpingSolverSource 
       ./ImageWarpingSolver.cpp 
       ./EquationSystems/ImageWarpingSystem.cpp
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
    SET(ImageWarpingSolverSource ${ImageWarpingSolverSource}
               ../Auxiliary/DriverArpack.cpp)
ENDIF (NEKTAR_USE_ARPACK)

ADD_SOLVER_EXECUTABLE(ImageWarpingSolver solvers 
			${ImageWarpingSolverSource})
