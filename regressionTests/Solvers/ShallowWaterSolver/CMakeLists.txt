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

FILE(GLOB input_files "InputFiles/*.xml" "InputFiles/*.rst" "InputFiles/*.bse"
						"InputFiles/*.bc")
INSTALL(FILES ${input_files} 
		DESTINATION ${NEKTAR_SHARE_DIR}/test/Solvers/ShallowWaterSolver/InputFiles COMPONENT test)
FILE(GLOB ok_files "OkFiles/*.ok")
INSTALL(FILES ${ok_files}
		DESTINATION ${NEKTAR_SHARE_DIR}/test/Solvers/ShallowWaterSolver/OkFiles COMPONENT test)

ADD_TEST(NAME ShallowWaterSolverTest COMMAND ShallowWaterSolverTest)
