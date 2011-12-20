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

FILE(GLOB input_files "InputFiles/*.xml" "InputFiles/*.rst" "InputFiles/*.bse"
						"InputFiles/*.bc")
INSTALL(FILES ${input_files} 
		DESTINATION ${NEKTAR_SHARE_DIR}/test/Solvers/CompressibleFlowSolver/InputFiles COMPONENT test)
FILE(GLOB ok_files "OkFiles/*.ok")
INSTALL(FILES ${ok_files}
		DESTINATION ${NEKTAR_SHARE_DIR}/test/Solvers/CompressibleFlowSolver/OkFiles COMPONENT test)
ADD_TEST(NAME CompressibleFlowSolverTest COMMAND CompressibleFlowSolverTest)
