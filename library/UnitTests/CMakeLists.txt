
SET(UnitTestSources
CountedObject.cpp
main.cpp
testNekManager.cpp
testNekMatrix.cpp
testNekVector.cpp
testNekPoint.cpp
testExpressionTemplates.cpp
testLinearSystem.cpp
)

SET(UnitTestHeaders
CountedObject.h
testNekManager.h
testNekMatrix.h
testNekPoint.h
testNekVector.h
testExpressionTemplates.h
testLinearSystem.h
)

ADD_DEFINITIONS(-DENABLE_NEKTAR_EXCEPTIONS)

ADD_EXECUTABLE(UnitTests ${UnitTestSources} ${UnitTestHeaders})

IF( ${CMAKE_COMPILER_IS_GNUCXX} )
        IF( USE_GCOV )
		TARGET_LINK_LIBRARIES(UnitTests gcov)
	ENDIF( USE_GCOV )
ENDIF( ${CMAKE_COMPILER_IS_GNUCXX} )

TARGET_LINK_LIBRARIES(UnitTests
    optimized LibUtilities debug LibUtilities-g
    optimized ${BoostThreadName} debug ${BoostThreadDebugName}
)

#MESSAGE(${USE_MKL})
#MESSAGE(${MKL_FOUND})
#IF( USE_MKL )
#    IF( MKL_FOUND )
#        MESSAGE("Using MKL>>>>>")
#    ENDIF( MKL_FOUND )
#ENDIF( USE_MKL )

IF( USE_MKL AND MKL_FOUND )
#        MESSAGE("Using MKL")
        TARGET_LINK_LIBRARIES(UnitTests
            ${MKL_LAPACK}
            optimized ${MKL} debug ${MKL}
            ${MKL_GUIDE}

#optimized mkl_ia32 debug mkl_ia32
#            optimized mkl debug mkl
#            optimized vml debug vml
#            optimized guide debug guide
        )
ENDIF( USE_MKL AND MKL_FOUND )

IF( USE_ACML AND ACML_FOUND )
    TARGET_LINK_LIBRARIES(UnitTests ${ACML} gfortran )
ENDIF( USE_ACML AND ACML_FOUND )

IF( USE_ATLAS AND ATLAS_FOUND )
        TARGET_LINK_LIBRARIES(UnitTests
            ${ATLAS_LAPACK}
            ${ATLAS_CBLAS}
            ${ATLAS} )
ENDIF( USE_ATLAS AND ATLAS_FOUND )

IF( ${CMAKE_COMPILER_IS_GNUCXX} )
        IF( USE_GCOV )
        TARGET_LINK_LIBRARIES(UnitTests gcov)
    ENDIF( USE_GCOV )
ENDIF( ${CMAKE_COMPILER_IS_GNUCXX} )

INSTALL_TARGETS(${NEKTAR_BIN_DIR} UnitTests)
