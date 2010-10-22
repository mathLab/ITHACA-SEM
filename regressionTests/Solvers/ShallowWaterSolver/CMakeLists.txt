SET(ShallowWaterSolverTestSource
	   ./ShallowWaterSolverTests.cpp
         ../../Auxiliary/RegressBase.cpp)

SET(ShallowWaterSolverTestHeaders
            ../../Auxiliary/RegressBase.h)
       
ADD_REGRESSION_EXECUTABLE(ShallowWaterSolverTest
			${ShallowWaterSolverTestHeaders} ${ShallowWaterSolverTestSource})

IF(REGRESSION_MAKE_OK_FILES_EXEC)
  SET(ShallowWaterSolverMakeOkFileSource
	   ./ShallowWaterSolverMakeOkFile.cpp
            ../../Auxiliary/RegressBase.cpp)

  SET(ShallowWaterSolverMakeOkFileHeaders
            ../../Auxiliary/RegressBase.h)
       
   ADD_REGRESSION_EXECUTABLE(ShallowWaterSolverMakeOkFile
			${ShallowWaterSolverMakeOkFileHeaders} 
                        ${ShallowWaterSolverMakeOkFileSource})
ENDIF(REGRESSION_MAKE_OK_FILES_EXEC)


