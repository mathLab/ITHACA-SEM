SET(ExpressionTemplatesUnitTestsSources
    ExpressionTemplateObjects.cpp
    main.cpp 
    TestAliasing.cpp
    TestAssociativeTransform.cpp
    TestBackwardsAssociativeTransform.cpp
    TestCommutativeTransform.cpp
    TestComplicatedMatrixExpressions.cpp
    TestCoupledLinearNSError.cpp
    TestFindACNodesToSwap.cpp
    TestForwardTransform.cpp
    TestInverseAssociative.cpp
    TestInvertNode.cpp
    TestTemporaryCount.cpp
    TestConstantConstantTree.cpp
    TestMatrixBlasOptimizations.cpp
    TestMatrixBufferSize.cpp 
    TestMatrixVectorOperations.cpp
    TestPushDownUnaryNodes.cpp
    TestRemoveUnecessaryTemporaries.cpp
    TestSortACClusters.cpp
    TestUnaryOperations.cpp
    TestUnrolledVectors.cpp
    ../../util.cpp
)

ADD_DEFINITIONS(-DENABLE_NEKTAR_EXCEPTIONS -DNEKTAR_UNIT_TESTS)
ADD_NEKTAR_EXECUTABLE(ExpressionTemplatesUnitTests
    COMPONENT unit-test DEPENDS LibUtilities SOURCES ${ExpressionTemplatesUnitTestsSources})
ADD_TEST(NAME ExpressionTemplatesUnitTests COMMAND ExpressionTemplatesUnitTests --detect_memory_leaks=0)
