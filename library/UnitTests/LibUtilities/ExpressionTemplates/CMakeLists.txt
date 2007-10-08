
SET(TARGET_NAME ExpressionTemplateUnitTests)

SET(Sources
    main.cpp   
    TestTwoParameters.cpp
    TestThreeParameters.cpp
    TestFourParameters.cpp
)

SET(Headers
)

ADD_DEFINITIONS(-DENABLE_NEKTAR_EXCEPTIONS)

ADD_NEKTAR_EXECUTABLE(${TARGET_NAME} Sources Headers)


TARGET_LINK_LIBRARIES(${TARGET_NAME}
    optimized LibUtilities debug LibUtilities-g
    optimized ${BOOST_THREAD_LIB} debug ${BOOST_THREAD_DEBUG_LIB}
)

SET_LAPACK_LINK_LIBRARIES(${TARGET_NAME})

INSTALL(TARGETS ${TARGET_NAME} RUNTIME DESTINATION ${NEKTAR_BIN_DIR})
