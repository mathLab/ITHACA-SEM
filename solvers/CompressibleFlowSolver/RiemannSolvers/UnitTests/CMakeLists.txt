SET(RiemannCompUnitTestsSources
    main.cpp
    TestRiemann.cpp
)

ADD_SOLVER_UTILITY(RiemannCompUnitTests CompressibleFlowSolver
    SOURCES ${RiemannCompUnitTestsSources} USELIBRARY)
ADD_TEST(NAME RiemannCompUnitTests COMMAND RiemannCompUnitTests --detect_memory_leaks=0)