SET(SOLVER_UTILS_SOURCES
  Core/Misc.cpp
  Advection/Advection.cpp
  Advection/AdvectionFR.cpp
  Advection/AdvectionNonConservative.cpp
  Advection/AdvectionWeakDG.cpp
  Diffusion/Diffusion.cpp
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
  Filters/FilterCheckpoint.cpp
  Filters/FilterHistoryPoints.cpp
  Filters/FilterModalEnergy.cpp
  Filters/FilterThresholdMax.cpp
  RiemannSolvers/RiemannSolver.cpp
  RiemannSolvers/UpwindSolver.cpp
  RiemannSolvers/UpwindLDGSolver.cpp
  UnsteadySystem.cpp
  Forcing/Forcing.cpp
  Forcing/ForcingSponge.cpp
  Forcing/ForcingBody.cpp
  Forcing/ForcingProgrammatic.cpp
)

SET(SOLVER_UTILS_HEADERS
  Core/Misc.h
  Advection/Advection.h
  Advection/AdvectionFR.h
  Advection/AdvectionNonConservative.h
  Advection/AdvectionWeakDG.h
  Diffusion/Diffusion.h
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
  Filters/FilterCheckpoint.h
  Filters/FilterHistoryPoints.h
  Filters/FilterModalEnergy.h
  Filters/FilterThresholdMax.h
  RiemannSolvers/RiemannSolver.h
  RiemannSolvers/UpwindSolver.h
  RiemannSolvers/UpwindLDGSolver.h
  SolverUtils.hpp
  SolverUtilsDeclspec.h
  UnsteadySystem.h
  Forcing/Forcing.h
  Forcing/ForcingSponge.h
  Forcing/ForcingBody.h
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
TARGET_LINK_LIBRARIES(SolverUtils
  optimized LibUtilities   debug LibUtilities-g
  optimized SpatialDomains debug SpatialDomains-g
  optimized StdRegions     debug StdRegions-g
  optimized MultiRegions   debug MultiRegions-g
  optimized LocalRegions   debug LocalRegions-g
  optimized ${TINYXML_LIB} debug ${TINYXML_LIB}
  ${Boost_THREAD_LIBRARY}
  ${Boost_IOSTREAMS_LIBRARY}
  ${Boost_DATE_TIME_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${ZLIB_LIBRARY}
  )
SET_LAPACK_LINK_LIBRARIES(SolverUtils)

INSTALL(FILES ${SOLVER_UTILS_HEADERS} DESTINATION ${NEKTAR_INCLUDE_DIR}/SolverUtils COMPONENT dev)
