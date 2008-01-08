
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
    TestBinaryConstantSpecialization0.cpp
    TestBinaryConstantSpecialization1.cpp
    TestBinaryConstantSpecialization2.cpp
    TestBinaryExpressionEvaluatorWithConstantConstantParameters.cpp
    TestBinaryExpressionEvaluatorWithBinaryConstantParameters.cpp
    TestAliasing.cpp
)

SET(Headers
)

ADD_DEFINITIONS(-DENABLE_NEKTAR_EXCEPTIONS -DNEKTAR_UNIT_TESTS)

ADD_NEKTAR_EXECUTABLE(${TARGET_NAME} Sources Headers)


TARGET_LINK_LIBRARIES(${TARGET_NAME}
    optimized LibUtilities debug LibUtilities-g
    optimized ${BOOST_THREAD_LIB} debug ${BOOST_THREAD_DEBUG_LIB}
)

SET_LAPACK_LINK_LIBRARIES(${TARGET_NAME})

