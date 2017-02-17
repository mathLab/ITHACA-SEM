ADD_SOLVER_UTILITY(ExtractSurface2DCFS CompressibleFlowSolver SOURCES ExtractSurface2DCFS.cpp)
ADD_SOLVER_UTILITY(ExtractSurface3DCFS CompressibleFlowSolver SOURCES ExtractSurface3DCFS.cpp)
ADD_SOLVER_UTILITY(CompressibleBL CompressibleFlowSolver SOURCES CompressibleBL.cpp)

FINALISE_CPACK_COMPONENT(compressibleflow-solver
    SUMMARY "Nektar++ compressible Navier-Stokes flow solver"
    DESCRIPTION "Solves the compressible Euler and compressible Navier-Stokes equations.")
