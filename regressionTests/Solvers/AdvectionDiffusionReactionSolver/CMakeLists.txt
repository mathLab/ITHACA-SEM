SET(AdvectionDiffusionReactionSolverTestSource
	   ./AdvectionDiffusionReactionSolverTests.cpp
         ../../Auxiliary/RegressBase.cpp)

SET(AdvectionDiffusionReactionSolverTestHeaders
            ../../Auxiliary/RegressBase.h)
       
ADD_REGRESSION_EXECUTABLE(AdvectionDiffusionReactionSolverTest
			${AdvectionDiffusionReactionSolverTestHeaders} ${AdvectionDiffusionReactionSolverTestSource})

IF(REGRESSION_MAKE_OK_FILES_EXEC)
  SET(AdvectionDiffusionReactionSolverMakeOkFileSource
	   ./AdvectionDiffusionReactionSolverMakeOkFile.cpp
            ../../Auxiliary/RegressBase.cpp)

  SET(AdvectionDiffusionReactionSolverMakeOkFileHeaders
            ../../Auxiliary/RegressBase.h)
       
   ADD_REGRESSION_EXECUTABLE(AdvectionDiffusionReactionSolverMakeOkFile
			${AdvectionDiffusionReactionSolverMakeOkFileHeaders} 
                        ${AdvectionDiffusionReactionSolverMakeOkFileSource})
ENDIF(REGRESSION_MAKE_OK_FILES_EXEC)


