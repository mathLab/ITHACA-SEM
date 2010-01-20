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


