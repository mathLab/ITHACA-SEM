SET(AdvectionDiffusionReactionSolverSource 
       ./AdvectionDiffusionReactionSolver.cpp 
       ./AdvectionDiffusionReaction.cpp
       ../Auxiliary/ADRBase.cpp)

ADD_SOLVER_EXECUTABLE(AdvectionDiffusionReactionSolver 
			${AdvectionDiffusionReactionSolverSource})
