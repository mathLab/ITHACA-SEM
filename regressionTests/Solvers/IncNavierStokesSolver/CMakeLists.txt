SET(IncNavierStokesSolverTestSource
	   ./IncNavierStokesSolverTests.cpp
         ../../Auxiliary/RegressBase.cpp)

SET(IncNavierStokesSolverTestHeaders
            ../../Auxiliary/RegressBase.h)
       
ADD_REGRESSION_EXECUTABLE(IncNavierStokesSolverTest
			${IncNavierStokesSolverTestHeaders} ${IncNavierStokesSolverTestSource})

IF(REGRESSION_MAKE_OK_FILES_EXEC)
  SET(IncNavierStokesSolverMakeOkFileSource
	   ./IncNavierStokesSolverMakeOkFile.cpp
            ../../Auxiliary/RegressBase.cpp)

  SET(IncNavierStokesSolverMakeOkFileHeaders
            ../../Auxiliary/RegressBase.h)
       
   ADD_REGRESSION_EXECUTABLE(IncNavierStokesSolverMakeOkFile
			${IncNavierStokesSolverMakeOkFileHeaders} 
                        ${IncNavierStokesSolverMakeOkFileSource})
ENDIF(REGRESSION_MAKE_OK_FILES_EXEC)


