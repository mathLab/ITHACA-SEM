SET(solverTestSource
	./solverTest.cpp
      ./regSolver.cpp 
	../Auxiliary/regBase.cpp 

)

SET(solverTestHeaders
	./regSolver.h
	../Auxiliary/regBase.h 	
)
       
ADD_REGRESSION_EXECUTABLE(solverTest
			${solverTestHeaders} ${solverTestSource})

