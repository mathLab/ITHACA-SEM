ADD_LIBRARY(CompRiemann STATIC
    AverageSolver.cpp
    AUSM0Solver.cpp
    AUSM1Solver.cpp
    AUSM2Solver.cpp
    AUSM3Solver.cpp
    CompressibleSolver.cpp
    ExactSolverToro.cpp
    HLLSolver.cpp
    HLLCSolver.cpp
    LaxFriedrichsSolver.cpp
    RoeSolver.cpp)

ADD_DEPENDENCIES(CompRiemann SolverUtils)
# TARGET_LINK_LIBRARIES(CompRiemann SolverUtils)

IF(NEKTAR_BUILD_UNIT_TESTS)
    ADD_SUBDIRECTORY(UnitTests)
ENDIF()