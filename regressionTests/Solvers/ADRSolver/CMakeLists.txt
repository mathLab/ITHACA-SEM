SET(ADRSolverTestSource
	   ./ADRSolverTests.cpp
         ../../Auxiliary/RegressBase.cpp)

SET(ADRSolverTestHeaders
            ../../Auxiliary/RegressBase.h)
       
ADD_REGRESSION_EXECUTABLE(ADRSolverTest
			${ADRSolverTestHeaders} ${ADRSolverTestSource})

IF(REGRESSION_MAKE_OK_FILES_EXEC)
  SET(ADRSolverMakeOkFileSource
	   ./ADRSolverMakeOkFile.cpp
            ../../Auxiliary/RegressBase.cpp)

  SET(ADRSolverMakeOkFileHeaders
            ../../Auxiliary/RegressBase.h)
       
   ADD_REGRESSION_EXECUTABLE(ADRSolverMakeOkFile
			${ADRSolverMakeOkFileHeaders} 
                        ${ADRSolverMakeOkFileSource})
ENDIF(REGRESSION_MAKE_OK_FILES_EXEC)


