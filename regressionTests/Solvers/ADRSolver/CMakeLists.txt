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

FILE(GLOB input_files "InputFiles/*.xml")
INSTALL(FILES ${input_files} 
		DESTINATION ${NEKTAR_SHARE_DIR}/test/Solvers/ADRSolver/InputFiles COMPONENT test)
FILE(GLOB ok_files "OkFiles/*.ok")
INSTALL(FILES ${ok_files}
		DESTINATION ${NEKTAR_SHARE_DIR}/test/Solvers/ADRSolver/OkFiles COMPONENT test)


ADD_TEST(NAME ADRSolverTest COMMAND ADRSolverTest)
