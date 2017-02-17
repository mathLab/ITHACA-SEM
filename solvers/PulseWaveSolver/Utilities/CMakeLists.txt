ADD_SOLVER_UTILITY(Fld2Tecplot PulseWaveSolver SOURCES Fld2Tecplot.cpp USELIBRARY)

FINALISE_CPACK_COMPONENT(pulsewave-solver
    SUMMARY "Nektar++ pulse wave solver"
    DESCRIPTION "Models pulse wave propagation through a 1-dimensional network of arteries.")

