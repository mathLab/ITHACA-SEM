
SET(Sources
    main.cpp
    StringConcat.cpp
    StringConcatExprTemp.cpp
)

SET(Headers
    StringConcat.h
    StringConcatExprTemp.h
)

ADD_NEKTAR_EXECUTABLE(Profile Sources Headers)

TARGET_LINK_LIBRARIES(Profile
    optimized LibUtilities debug LibUtilities-g
    optimized ${BOOST_THREAD_LIB} debug ${BOOST_THREAD_DEBUG_LIB}
)

SET_LAPACK_LINK_LIBRARIES(Profile)

