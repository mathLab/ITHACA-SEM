SET(LinkLibraries MultiRegions)

SET(PP_SOURCES ./Prepacing.cpp)
ADD_SOLVER_UTILITY(PrePacing CardiacEPSolver SOURCES Prepacing.cpp USELIBRARY)

ADD_NEKTAR_TEST(Courtemanche)
ADD_NEKTAR_TEST(CourtemancheAF)
ADD_NEKTAR_TEST(FentonKarma)

FINALISE_CPACK_COMPONENT(cardiacep-solver
    SUMMARY "Nektar++ cardiac electrophysiology solver"
    DESCRIPTION "Solves the cardiac electrophysiology monodomain equations.")
