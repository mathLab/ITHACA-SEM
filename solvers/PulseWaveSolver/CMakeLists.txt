SET(PulseWaveSolverSource 
       ./PulseWaveSolver.cpp 
       ./EquationSystems/PulseWaveSystem.h
       ./EquationSystems/PulseWaveSystem.cpp
       ./EquationSystems/PulseWavePropagation.h
       ./EquationSystems/PulseWavePropagation.cpp)

ADD_SOLVER_EXECUTABLE(PulseWaveSolver solvers 
			${PulseWaveSolverSource})

#ADD_NEKTAR_TEST(Bifurcation)
#ADD_NEKTAR_TEST(Junction)
#ADD_NEKTAR_TEST(Paper1_DG)
#ADD_NEKTAR_TEST(Paper2_DG)
#ADD_NEKTAR_TEST(Paper3_DG)
#ADD_NEKTAR_TEST(Paper4_55Network)
#ADD_NEKTAR_TEST(Q_inflow)
#ADD_NEKTAR_TEST(TwoBifurcations)
