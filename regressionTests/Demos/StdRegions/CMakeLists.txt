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


