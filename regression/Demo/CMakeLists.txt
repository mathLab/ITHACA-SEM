SET(demoTestSource
	   ./demoTest.cpp
	   ./regDemo.cpp 
         ../Auxiliary/regBase.cpp)

SET(demoTestHeaders
	     ./regDemo.h 
            ../Auxiliary/regBase.h)
       
ADD_REGRESSION_EXECUTABLE(demoTest
			${demoTestHeaders} ${demoTestSource})

