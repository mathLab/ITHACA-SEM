SET(CompressibleFlowSolverTestSource
	   ./CompressibleFlowSolverTests.cpp
         ../../Auxiliary/RegressBase.cpp)

SET(CompressibleFlowSolverTestHeaders
            ../../Auxiliary/RegressBase.h)
       
ADD_REGRESSION_EXECUTABLE(CompressibleFlowSolverTest
			${CompressibleFlowSolverTestHeaders} ${CompressibleFlowSolverTestSource})

IF(REGRESSION_MAKE_OK_FILES_EXEC)
  SET(CompressibleFlowSolverMakeOkFileSource
	   ./CompressibleFlowSolverMakeOkFile.cpp
            ../../Auxiliary/RegressBase.cpp)

  SET(CompressibleFlowSolverMakeOkFileHeaders
            ../../Auxiliary/RegressBase.h)
       
   ADD_REGRESSION_EXECUTABLE(CompressibleFlowSolverMakeOkFile
			${CompressibleFlowSolverMakeOkFileHeaders} 
                        ${CompressibleFlowSolverMakeOkFileSource})
ENDIF(REGRESSION_MAKE_OK_FILES_EXEC)


