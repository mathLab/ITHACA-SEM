#!/bin/bash
Demos/StdRegions/StdRegionsDemoTest && \
Demos/LocalRegions/LocalRegionsDemoTest && \
Demos/MultiRegions/MultiRegionsDemoTest && \
Solvers/ADRSolver/ADRSolverTest && \
Solvers/IncNavierStokesSolver/IncNavierStokesSolverTest && \
Solvers/CompressibleFlowSolver/CompressibleFlowSolverTest && \
Solvers/ShallowWaterSolver/ShallowWaterSolverTest
rm -f *.dat *.err *.pos
