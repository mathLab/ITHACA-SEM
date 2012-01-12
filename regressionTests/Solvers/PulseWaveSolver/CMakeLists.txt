SET(PulseWaveSolverTestSource
	   ./PulseWaveSolverTests.cpp
         ../../Auxiliary/RegressBase.cpp)

SET(PulseWaveSolverTestHeaders
            ../../Auxiliary/RegressBase.h)
       
ADD_REGRESSION_EXECUTABLE(PulseWaveSolverTest
			${PulseWaveSolverTestHeaders} ${PulseWaveSolverTestSource})

IF(REGRESSION_MAKE_OK_FILES_EXEC)
  SET(PulseWaveSolverMakeOkFileSource
	   ./PulseWaveSolverMakeOkFile.cpp
            ../../Auxiliary/RegressBase.cpp)

  SET(PulseWaveSolverMakeOkFileHeaders
            ../../Auxiliary/RegressBase.h)
       
   ADD_REGRESSION_EXECUTABLE(PulseWaveSolverMakeOkFile
			${PulseWaveSolverMakeOkFileHeaders} 
                        ${PulseWaveSolverMakeOkFileSource})
ENDIF(REGRESSION_MAKE_OK_FILES_EXEC)

FILE(GLOB input_files "InputFiles/*.xml" "InputFiles/*.rst" "InputFiles/*.bse"
						"InputFiles/*.bc")
INSTALL(FILES ${input_files} 
		DESTINATION ${NEKTAR_SHARE_DIR}/test/Solvers/PulseWaveSolver/InputFiles COMPONENT test)
FILE(GLOB ok_files "OkFiles/*.ok")
INSTALL(FILES ${ok_files}
		DESTINATION ${NEKTAR_SHARE_DIR}/test/Solvers/PulseWaveSolver/OkFiles COMPONENT test)

ADD_TEST(NAME PulseWaveSolverTest COMMAND PulseWaveSolverTest)
