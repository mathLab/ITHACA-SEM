SUBDIRS(ADRSolver IncNavierStokesSolver ShallowWaterSolver CompressibleFlowSolver)

IF(NEKTAR_USE_ARPACK)
SUBDIRS(IncNavierStokesStability)
ENDIF(NEKTAR_USE_ARPACK)

