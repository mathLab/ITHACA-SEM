SET(SOLVERS_AUXILIARY_HEADERS
AdvectionDiffusionReaction.h
ShallowWaterEquations.h
EulerEquations.h
)

SET(SOLVERS_AUXILIARY_SOURCES
AdvectionDiffusionReaction.cpp
ShallowWaterEquations.cpp
EulerEquations.cpp
)

ADD_SOLVER_LIBRARY(SolverAuxiliary STATIC ${SOLVERS_AUXILIARY_HEADERS} ${SOLVERS_AUXILIARY_SOURCES})
