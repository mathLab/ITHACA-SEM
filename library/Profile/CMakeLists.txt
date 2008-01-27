
SET(Sources
    main.cpp
    StringConcat.cpp
    StringConcatExprTemp.cpp
    IntWrapper.cpp
    IntWrapperExprTemp.cpp
    VectorOpsExprTemp.cpp
    MatrixOps.cpp
)

SET(Headers
    StringConcat.h
    StringConcatExprTemp.h
    IntWrapper.h
    IntWrapperExprTemp.h
    VectorOps.h
    VectorOpsExprTemp.h
    MatrixOps.h
    MatrixOpsExprTemp.h
)

ADD_NEKTAR_EXECUTABLE(Profile Sources Headers)

TARGET_LINK_LIBRARIES(Profile
    optimized LibUtilities debug LibUtilities-g
    optimized ${BOOST_THREAD_LIB} debug ${BOOST_THREAD_DEBUG_LIB}
)

SET_LAPACK_LINK_LIBRARIES(Profile)

SET(GeneratorSources
	Generator.cpp)

ADD_NEKTAR_EXECUTABLE(ProfileGenerator GeneratorSources)

