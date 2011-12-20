SET(MultiRegionsDemoTestSource
	   ./MultiRegionsDemoTests.cpp
         ../../Auxiliary/RegressBase.cpp)

SET(MultiRegionsDemoTestHeaders
            ../../Auxiliary/RegressBase.h)
       
ADD_REGRESSION_EXECUTABLE(MultiRegionsDemoTest
			${MultiRegionsDemoTestHeaders} ${MultiRegionsDemoTestSource})

IF(REGRESSION_MAKE_OK_FILES_EXEC)
  SET(MultiRegionsDemoMakeOkFileSource
	   ./MultiRegionsDemoMakeOkFile.cpp
            ../../Auxiliary/RegressBase.cpp)

  SET(MultiRegionsDemoMakeOkFileHeaders
            ../../Auxiliary/RegressBase.h)
       
   ADD_REGRESSION_EXECUTABLE(MultiRegionsDemoMakeOkFile
			${MultiRegionsDemoMakeOkFileHeaders} 
                        ${MultiRegionsDemoMakeOkFileSource})
ENDIF(REGRESSION_MAKE_OK_FILES_EXEC)

FILE(GLOB input_files "InputFiles/*.xml")
INSTALL(FILES ${input_files} 
		DESTINATION ${NEKTAR_SHARE_DIR}/test/Demos/MultiRegions/InputFiles COMPONENT test)
FILE(GLOB ok_files "OkFiles/*.ok")
INSTALL(FILES ${ok_files}
		DESTINATION ${NEKTAR_SHARE_DIR}/test/Demos/MultiRegions/OkFiles COMPONENT test)
ADD_TEST(NAME MultiRegionsDemoTest COMMAND MultiRegionsDemoTest)
