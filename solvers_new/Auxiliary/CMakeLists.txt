SET(SOLVERS_AUXILIARY_HEADERS
ADRBase.h
AdvectionDiffusionReaction.h
#ShallowWaterEquations.h
#EulerEquations.h
)

SET(SOLVERS_AUXILIARY_SOURCES
ADRBase.cpp
AdvectionDiffusionReaction.cpp
#ShallowWaterEquations.cpp
#EulerEquations.cpp
)

ADD_SOLVER_LIBRARY(SolverAuxiliary STATIC ${SOLVERS_AUXILIARY_HEADERS} ${SOLVERS_AUXILIARY_SOURCES})
