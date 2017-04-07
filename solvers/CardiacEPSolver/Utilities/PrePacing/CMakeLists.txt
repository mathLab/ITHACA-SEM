ADD_SOLVER_UTILITY(PrePacing CardiacEPSolver SOURCES Prepacing.cpp USELIBRARY)
TARGET_LINK_LIBRARIES(PrePacing MultiRegions)

ADD_NEKTAR_TEST(Courtemanche)
ADD_NEKTAR_TEST(CourtemancheAF)
ADD_NEKTAR_TEST(FentonKarma)

FINALISE_CPACK_COMPONENT(cardiacep-solver
    SUMMARY "Nektar++ cardiac electrophysiology solver"
    DESCRIPTION "Solves the cardiac electrophysiology monodomain equations.")
