
SET(TARGET_NAME ExpressionTemplateUnitTests)

SET(Sources 
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

SET(Headers
    CountedObjectExpression.h
    ExpressionTemplateObjects.h
	../../util.h
)

ADD_DEFINITIONS(-DENABLE_NEKTAR_EXCEPTIONS -DNEKTAR_UNIT_TESTS)
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

ADD_NEKTAR_EXECUTABLE(${TARGET_NAME} unit-test Sources Headers)


TARGET_LINK_LIBRARIES(${TARGET_NAME}
    optimized LibUtilities debug LibUtilities-g
    optimized ${Boost_THREAD_LIBRARY_RELEASE} debug ${Boost_THREAD_LIBRARY_DEBUG}
)

SET_LAPACK_LINK_LIBRARIES(${TARGET_NAME})
ADD_TEST(NAME ExpressionTemplates COMMAND ${TARGET_NAME} --detect_memory_leaks=0)
