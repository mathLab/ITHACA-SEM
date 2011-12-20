SET(LocalRegionsDemoTestSource
	   ./LocalRegionsDemoTests.cpp
         ../../Auxiliary/RegressBase.cpp)

SET(LocalRegionsDemoTestHeaders
            ../../Auxiliary/RegressBase.h)
       
ADD_REGRESSION_EXECUTABLE(LocalRegionsDemoTest
			${LocalRegionsDemoTestHeaders} ${LocalRegionsDemoTestSource})

IF(REGRESSION_MAKE_OK_FILES_EXEC)
  SET(LocalRegionsDemoMakeOkFileSource
	   ./LocalRegionsDemoMakeOkFile.cpp
            ../../Auxiliary/RegressBase.cpp)

  SET(LocalRegionsDemoMakeOkFileHeaders
            ../../Auxiliary/RegressBase.h)
       
   ADD_REGRESSION_EXECUTABLE(LocalRegionsDemoMakeOkFile
			${LocalRegionsDemoMakeOkFileHeaders} 
                        ${LocalRegionsDemoMakeOkFileSource})
ENDIF(REGRESSION_MAKE_OK_FILES_EXEC)

FILE(GLOB input_files "InputFiles/*.xml")
INSTALL(FILES ${input_files} 
		DESTINATION ${NEKTAR_SHARE_DIR}/test/Demos/LocalRegions/InputFiles COMPONENT test)
FILE(GLOB ok_files "OkFiles/*.ok")
INSTALL(FILES ${ok_files}
		DESTINATION ${NEKTAR_SHARE_DIR}/test/Demos/LocalRegions/OkFiles COMPONENT test)
ADD_TEST(NAME LocalRegionsDemoTest COMMAND LocalRegionsDemoTest)
