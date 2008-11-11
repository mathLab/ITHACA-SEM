
SET(TARGET_NAME ExpressionTemplateUnitTests)

SET(Sources 
    main.cpp 
    TestTwoParameters.cpp 
    TestThreeParameters.cpp
    TestFourParameters.cpp 
    TestConstantBinarySpecialization0.cpp
    TestConstantBinarySpecialization1.cpp
    TestConstantBinarySpecialization2.cpp
    TestConstantBinarySpecialization3.cpp
    TestConstantBinarySpecialization4.cpp
    TestConstantBinarySpecialization5.cpp
    TestConstantBinarySpecialization0_ParentOpType.cpp
    TestConstantBinarySpecialization1_ParentOpType.cpp
    TestBinaryExpressionEvaluator_ConstantConstantSpecialization0.cpp
    TestBinaryExpressionEvaluator_ConstantConstantSpecialization1.cpp
    TestBinaryExpressionEvaluator_GeneralSpecialization0.cpp
    TestBinaryExpressionEvaluator_GeneralSpecialization1.cpp
    TestBinaryExpressionEvaluator_GeneralSpecialization2.cpp	
    TestBinaryExpressionEvaluator_GeneralSpecialization3.cpp
    TestBinaryExpressionEvaluator_GeneralSpecialization_ParentsOpType1.cpp
    TestBinaryConstantSpecialization0.cpp
    TestBinaryConstantSpecialization1.cpp
    TestBinaryConstantSpecialization2.cpp
    TestBinaryExpressionEvaluatorWithConstantConstantParameters.cpp
    TestBinaryExpressionEvaluatorWithBinaryConstantParameters.cpp
    TestAliasing.cpp
    TestMatrixBufferSize.cpp 
    ../../util.cpp
)

SET(Headers
	../../util.h
)

ADD_DEFINITIONS(-DENABLE_NEKTAR_EXCEPTIONS -DNEKTAR_UNIT_TESTS)
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

ADD_NEKTAR_EXECUTABLE(${TARGET_NAME} Sources Headers)


TARGET_LINK_LIBRARIES(${TARGET_NAME}
    optimized LibUtilities debug LibUtilities-g
    optimized ${Boost_THREAD_LIBRARY_RELEASE} debug ${Boost_THREAD_LIBRARY_DEBUG}
)

SET_LAPACK_LINK_LIBRARIES(${TARGET_NAME})

