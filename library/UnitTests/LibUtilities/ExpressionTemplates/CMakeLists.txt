
SET(TARGET_NAME ExpressionTemplateUnitTests)

SET(Sources 
    main.cpp 
    TestAliasing.cpp
    TestAssociativeTransform.cpp
    TestAssociativeTreeTransform.cpp
    TestComplicatedMatrixExpressions.cpp
    TestCoupledLinearNSError.cpp
    TestDotProduct.cpp
    TestTemporaryCount.cpp
    TestTwoParameters.cpp 
    TestThreeParameters.cpp
    TestFourParameters.cpp 
    TestConstantConstantTree.cpp
    TestMatrixBlasOptimizations.cpp
    TestMatrixBufferSize.cpp 
    TestNode.cpp
    TestRemoveUnecessaryTemporaries.cpp
    TestTraits.cpp
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
ADD_TEST(NAME ExpressionTemplates WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}/${NEKTAR_BIN_DIR} COMMAND ${TARGET_NAME} --detect_memory_leaks=0)
