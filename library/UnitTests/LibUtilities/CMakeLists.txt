
SET(UnitTestSources
        TestNekMatrixOperations.cpp
	TestMatrixStoragePolicies.cpp
	TestUpperTriangularMatrix.cpp
    main.cpp   
)

SET(UnitTestHeaders
    TestNekMatrixOperations.h
    TestMatrixStoragePolicies.h
    TestUpperTriangularMatrix.h
)

ADD_DEFINITIONS(-DENABLE_NEKTAR_EXCEPTIONS)

ADD_NEKTAR_EXECUTABLE(LibUtilitiesUnitTests UnitTestSources UnitTestHeaders)


IF( ${CMAKE_COMPILER_IS_GNUCXX} )
        IF( USE_GCOV )
		TARGET_LINK_LIBRARIES(LibUtilitiesUnitTests gcov)
	ENDIF( USE_GCOV )
ENDIF( ${CMAKE_COMPILER_IS_GNUCXX} )

TARGET_LINK_LIBRARIES(LibUtilitiesUnitTests
    optimized LibUtilities debug LibUtilities-g
    optimized StdRegions debug StdRegions-g
    optimized ${BOOST_THREAD_LIB} debug ${BOOST_THREAD_DEBUG_LIB}
)


#IF( USE_MKL AND MKL_FOUND )
#        TARGET_LINK_LIBRARIES(LibUtilitiesUnitTests
#            ${MKL_LAPACK}
#            optimized ${MKL} debug ${MKL}
#            ${MKL_GUIDE}
#        )
#ENDIF( USE_MKL AND MKL_FOUND )

#IF( USE_ACML AND ACML_FOUND )
#    TARGET_LINK_LIBRARIES(LibUtilitiesUnitTests ${ACML_TARGET_LINK_LIBRARIES}  )
#ENDIF( USE_ACML AND ACML_FOUND )

#IF( USE_ATLAS AND ATLAS_FOUND )
#        TARGET_LINK_LIBRARIES(LibUtilitiesUnitTests
#            ${ATLAS_LAPACK}
#            ${ATLAS_CBLAS}
#            ${ATLAS} )
#ENDIF( USE_ATLAS AND ATLAS_FOUND )

#IF( ${CMAKE_COMPILER_IS_GNUCXX} )
#        IF( USE_GCOV )
#        TARGET_LINK_LIBRARIES(LibUtilitiesUnitTests gcov)
#    ENDIF( USE_GCOV )
#ENDIF( ${CMAKE_COMPILER_IS_GNUCXX} )

#IF( USE_FRAMEWORK )
#    TARGET_LINK_LIBRARIES(LibUtilitiesUnitTests ${FRAMEWORK_LINK_FLAGS})
#ENDIF ( USE_FRAMEWORK )

#IF( USE_WIN32_LAPACK )
#	TARGET_LINK_LIBRARIES(LibUtilitiesUnitTests ${WIN32_LAPACK} ${WIN32_BLAS})	
#ENDIF( USE_WIN32_LAPACK )

INSTALL(TARGETS LibUtilitiesUnitTests RUNTIME DESTINATION ${NEKTAR_BIN_DIR})
