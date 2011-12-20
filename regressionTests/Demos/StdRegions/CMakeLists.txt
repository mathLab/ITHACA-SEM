SET(StdRegionsDemoTestSource
	   ./StdRegionsDemoTests.cpp
         ../../Auxiliary/RegressBase.cpp)

SET(StdRegionsDemoTestHeaders
            ../../Auxiliary/RegressBase.h)
       
ADD_REGRESSION_EXECUTABLE(StdRegionsDemoTest
			${StdRegionsDemoTestHeaders} ${StdRegionsDemoTestSource})

IF(REGRESSION_MAKE_OK_FILES_EXEC)
  SET(StdRegionsDemoMakeOkFileSource
	   ./StdRegionsDemoMakeOkFile.cpp
            ../../Auxiliary/RegressBase.cpp)

  SET(StdRegionsDemoMakeOkFileHeaders
            ../../Auxiliary/RegressBase.h)
       
   ADD_REGRESSION_EXECUTABLE(StdRegionsDemoMakeOkFile
			${StdRegionsDemoMakeOkFileHeaders} 
                        ${StdRegionsDemoMakeOkFileSource})
ENDIF(REGRESSION_MAKE_OK_FILES_EXEC)

FILE(GLOB input_files "InputFiles/*.xml")
INSTALL(FILES ${input_files} 
		DESTINATION ${NEKTAR_SHARE_DIR}/test/Demos/StdRegions/InputFiles COMPONENT test)
FILE(GLOB ok_files "OkFiles/*.ok")
INSTALL(FILES ${ok_files}
		DESTINATION ${NEKTAR_SHARE_DIR}/test/Demos/StdRegions/OkFiles COMPONENT test)

ADD_TEST(NAME StdRegionsDemoTest COMMAND StdRegionsDemoTest)
