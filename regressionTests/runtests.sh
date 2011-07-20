#!/bin/bash
Demos/StdRegions/StdRegionsDemoTest && \
Demos/LocalRegions/LocalRegionsDemoTest && \
Demos/MultiRegions/MultiRegionsDemoTest && \
Solvers/ADRSolver/ADRSolverTest && \
Solvers/IncNavierStokesSolver/IncNavierStokesSolverTest && \
Solvers/IncNavierStokesStability/IncNavierStokesStabilityTest && \
Solvers/CompressibleFlowSolver/CompressibleFlowSolverTest
rm -f *.dat *.err *.pos
