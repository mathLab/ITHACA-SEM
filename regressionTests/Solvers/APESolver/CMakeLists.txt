SET(APESolverTestSource
	   ./APESolverTests.cpp
         ../../Auxiliary/RegressBase.cpp)

SET(APESolverTestHeaders
            ../../Auxiliary/RegressBase.h)
       
ADD_REGRESSION_EXECUTABLE(APESolverTest
			${APESolverTestHeaders} ${APESolverTestSource})

IF(REGRESSION_MAKE_OK_FILES_EXEC)
  SET(APESolverMakeOkFileSource
	   ./APESolverMakeOkFile.cpp
            ../../Auxiliary/RegressBase.cpp)

  SET(APESolverMakeOkFileHeaders
            ../../Auxiliary/RegressBase.h)
       
   ADD_REGRESSION_EXECUTABLE(APESolverMakeOkFile
			${APESolverMakeOkFileHeaders} 
                        ${APESolverMakeOkFileSource})
ENDIF(REGRESSION_MAKE_OK_FILES_EXEC)

FILE(GLOB input_files "InputFiles/*.xml" "InputFiles/*.rst" "InputFiles/*.bse"
						"InputFiles/*.bc")
INSTALL(FILES ${input_files} 
		DESTINATION ${NEKTAR_SHARE_DIR}/test/Solvers/APESolver/InputFiles COMPONENT test)
FILE(GLOB ok_files "OkFiles/*.ok")
INSTALL(FILES ${ok_files}
		DESTINATION ${NEKTAR_SHARE_DIR}/test/Solvers/APESolver/OkFiles COMPONENT test)

ADD_TEST(NAME APESolverTest COMMAND APESolverTest)
