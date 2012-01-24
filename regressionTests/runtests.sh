#!/bin/bash
Demos/StdRegions/StdRegionsDemoTest && \
Demos/LocalRegions/LocalRegionsDemoTest && \
Demos/MultiRegions/MultiRegionsDemoTest && \
Solvers/ADRSolver/ADRSolverTest && \
Solvers/IncNavierStokesSolver/IncNavierStokesSolverTest && \
Solvers/CompressibleFlowSolver/CompressibleFlowSolverTest && \
Solvers/ShallowWaterSolver/ShallowWaterSolverTest && \
Solvers/APESolver/APESolverTest && \
Solvers/PulseWaveSolver/PulseWaveSolverTest
rm -f *.dat *.err *.pos *.0 *.1 *.bse *_eig_? *.evl
