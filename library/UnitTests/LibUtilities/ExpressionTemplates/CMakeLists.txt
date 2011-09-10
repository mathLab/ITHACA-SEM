
SET(TARGET_NAME ExpressionTemplateUnitTests)

SET(Sources 
    main.cpp 
    TestComplicatedMatrixExpressions.cpp
    TestTwoParameters.cpp 
    TestThreeParameters.cpp
    TestFourParameters.cpp 
    TestConstantConstantTree.cpp
    TestAliasing.cpp
    TestMatrixBlasOptimizations.cpp
    TestMatrixBufferSize.cpp 
    TestTraits.cpp
    TestUnrolledVectors.cpp
    ../../util.cpp
)

SET(Headers
    CountedObjectExpression.h
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

