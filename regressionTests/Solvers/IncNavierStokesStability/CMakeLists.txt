
SET(IncNavierStokesStabilityTestSource
	   ./IncNavierStokesStabilityTests.cpp
         ../../Auxiliary/RegressBase.cpp)

SET(IncNavierStokesStabilityTestHeaders
            ../../Auxiliary/RegressBase.h)
       
ADD_REGRESSION_EXECUTABLE(IncNavierStokesStabilityTest
			${IncNavierStokesStabilityTestHeaders} ${IncNavierStokesStabilityTestSource})


IF(REGRESSION_MAKE_OK_FILES_EXEC)

  SET(IncNavierStokesStabilityMakeOkFileHeaders
            ../../Auxiliary/RegressBase.h)
       
   ADD_REGRESSION_EXECUTABLE(IncNavierStokesStabilityMakeOkFile
			${IncNavierStokesStabilityMakeOkFileHeaders} 
                        ${IncNavierStokesStabilityMakeOkFileSource})
ENDIF(REGRESSION_MAKE_OK_FILES_EXEC)

