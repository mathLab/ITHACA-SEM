SET(SOLVER_UTILS_SOURCES
  Core/Misc.cpp
  AdvectionSystem.cpp
  Advection/Advection.cpp
  Advection/Advection3DHomogeneous1D.cpp
  Advection/AdvectionFR.cpp
  Advection/AdvectionNonConservative.cpp
  Advection/AdvectionWeakDG.cpp
  Diffusion/Diffusion.cpp
  Diffusion/Diffusion3DHomogeneous1D.cpp
  Diffusion/DiffusionLDG.cpp
  Diffusion/DiffusionLDGNS.cpp
  Diffusion/DiffusionLFR.cpp
  Diffusion/DiffusionLFRNS.cpp
  Driver.cpp
  DriverArnoldi.cpp
  DriverModifiedArnoldi.cpp
  DriverStandard.cpp
  DriverSteadyState.cpp	
  EquationSystem.cpp
  Filters/Filter.cpp
  Filters/FilterAeroForces.cpp
  Filters/FilterAverageFields.cpp
  Filters/FilterCheckpoint.cpp
  Filters/FilterEnergyBase.cpp
  Filters/FilterHistoryPoints.cpp
  Filters/FilterModalEnergy.cpp
  Filters/FilterThresholdMax.cpp
  Filters/FilterThresholdMin.cpp
  RiemannSolvers/RiemannSolver.cpp
  RiemannSolvers/UpwindSolver.cpp
  RiemannSolvers/UpwindLDGSolver.cpp
  UnsteadySystem.cpp
  Forcing/Forcing.cpp
  Forcing/ForcingAbsorption.cpp
  Forcing/ForcingBody.cpp
  Forcing/ForcingNoise.cpp
  Forcing/ForcingProgrammatic.cpp
)

SET(SOLVER_UTILS_HEADERS
  Core/Misc.h
  AdvectionSystem.h
  Advection/Advection.h
  Advection/AdvectionFR.h
  Advection/Advection3DHomogeneous1D.h
  Advection/AdvectionNonConservative.h
  Advection/AdvectionWeakDG.h
  Diffusion/Diffusion.h
  Diffusion/Diffusion3DHomogeneous1D.h
  Diffusion/DiffusionLDG.h
  Diffusion/DiffusionLDGNS.h
  Diffusion/DiffusionLFR.h
  Diffusion/DiffusionLFRNS.h
  Driver.h
  DriverArnoldi.h
  DriverModifiedArnoldi.h
  DriverStandard.h
  DriverSteadyState.h
  EquationSystem.h
  Filters/Filter.h
  Filters/FilterAeroForces.h
  Filters/FilterAverageFields.h
  Filters/FilterCheckpoint.h
  Filters/FilterEnergyBase.h
  Filters/FilterHistoryPoints.h
  Filters/FilterModalEnergy.h
  Filters/FilterThresholdMax.h
  Filters/FilterThresholdMin.h
  RiemannSolvers/RiemannSolver.h
  RiemannSolvers/UpwindSolver.h
  RiemannSolvers/UpwindLDGSolver.h
  SolverUtils.hpp
  SolverUtilsDeclspec.h
  UnsteadySystem.h
  Forcing/Forcing.h
  Forcing/ForcingAbsorption.h
  Forcing/ForcingBody.h
  Forcing/ForcingNoise.h
  Forcing/ForcingProgrammatic.h
)

IF (NEKTAR_USE_ARPACK)
  SET(SOLVER_UTILS_SOURCES ${SOLVER_UTILS_SOURCES}
    DriverArpack.cpp)
  SET(SOLVER_UTILS_HEADERS ${SOLVER_UTILS_HEADERS}
    DriverArpack.h)
ENDIF (NEKTAR_USE_ARPACK)

ADD_DEFINITIONS(-DSOLVER_UTILS_EXPORTS)

ADD_NEKTAR_LIBRARY(SolverUtils lib ${NEKTAR_LIBRARY_TYPE} ${SOLVER_UTILS_SOURCES} ${SOLVER_UTILS_HEADERS})
TARGET_LINK_LIBRARIES(SolverUtils MultiRegions)

INSTALL(DIRECTORY ./ DESTINATION ${NEKTAR_INCLUDE_DIR}/SolverUtils COMPONENT dev FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
