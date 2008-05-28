
SET(Sources
    main.cpp
    StringConcat.cpp
    StringConcatExprTemp.cpp
    IntWrapper.cpp
    IntWrapperExprTemp.cpp
    VectorOpsExprTemp.cpp
    MatrixOps.cpp
    VariableCostObject.cpp
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
	VariableCostObject.h
)

ADD_NEKTAR_EXECUTABLE(Profile Sources Headers)

TARGET_LINK_LIBRARIES(Profile
    optimized LibUtilities debug LibUtilities-g
    optimized ${Boost_THREAD_LIBRARY_RELEASE} debug ${Boost_THREAD_LIBRARY_DEBUG}
)

SET_LAPACK_LINK_LIBRARIES(Profile)

SET(GeneratorSources
	Generator.cpp)

ADD_NEKTAR_EXECUTABLE(ProfileGenerator GeneratorSources)

