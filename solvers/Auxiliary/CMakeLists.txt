SET(SOLVERS_AUXILIARY_HEADERS
AdvectionDiffusionReaction.h
ShallowWaterEquations.h
BoussinesqEquations.h
EulerEquations.h
)

SET(SOLVERS_AUXILIARY_SOURCES
AdvectionDiffusionReaction.cpp
ShallowWaterEquations.cpp
BoussinesqEquations.cpp
EulerEquations.cpp
)

ADD_SOLVER_LIBRARY(SolverAuxiliary STATIC ${SOLVERS_AUXILIARY_HEADERS} ${SOLVERS_AUXILIARY_SOURCES})
