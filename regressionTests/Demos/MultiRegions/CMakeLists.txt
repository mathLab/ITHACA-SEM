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


