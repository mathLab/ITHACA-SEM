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

FILE(GLOB input_files "InputFiles/*.xml" "InputFiles/*.rst" "InputFiles/*.bse"
						"InputFiles/*.bc")
INSTALL(FILES ${input_files} 
		DESTINATION ${NEKTAR_SHARE_DIR}/test/Solvers/IncNavierStokesSolver/InputFiles COMPONENT test)
FILE(GLOB ok_files "OkFiles/*.ok")
INSTALL(FILES ${ok_files}
		DESTINATION ${NEKTAR_SHARE_DIR}/test/Solvers/IncNavierStokesSolver/OkFiles COMPONENT test)

ADD_TEST(NAME IncNavierStokesSolverTest COMMAND IncNavierStokesSolverTest)
