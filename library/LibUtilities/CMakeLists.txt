SET(LibUtilitySources
./LinearAlgebra/BlockMat.cpp
./Foundations/Graph.cpp
./Polylib/Polylib.cpp
)

SET(LibUtilityHeaders
./LinearAlgebra/NekLinSys.h
./LinearAlgebra/NekLinAlg.h
./LinearAlgebra/BlockMat.h
./Foundations/Graph.h
./Polylib/Polylib.h
./LinearAlgebra/NekMatrix.hpp
./LinearAlgebra/NekVector.hpp
./LinearAlgebra/NekPoint.hpp
./ExpressionTemplates/BinaryExpressionTraits.hpp
./ExpressionTemplates/BinaryExpression.hpp
./ExpressionTemplates/ConstantExpression.hpp
./ExpressionTemplates/ArithmeticConcepts.hpp
./ExpressionTemplates/ExpressionTemplateOperators.hpp
./ExpressionTemplates/ExpressionReturnType.hpp
./ExpressionTemplates/ExpressionTemplates.hpp
./ExpressionTemplates/ExpressionTemplateConcepts.hpp
./BasicUtils/BoostUtil.hpp
./BasicUtils/ErrorUtil.hpp
./BasicUtils/Lapack.hpp
./BasicUtils/Concepts.hpp
./BasicUtils/Vmath.hpp
./BasicUtils/Blas.hpp
./BasicUtils/TransF77.hpp
./Memory/ThreadSpecificPool.hpp
./Memory/NekMemoryManager.hpp
)


ADD_NEKTAR_LIBRARY(LibUtilities STATIC ${LibUtilitySources} ${LibUtilityHeaders})

IF( ${CMAKE_GENERATOR} STREQUAL "Visual Studio 7 .NET 2003" )
	# Needed for M_PI to be visible in visual studio.
	ADD_DEFINITIONS(-D_USE_MATH_DEFINES)
ENDIF( ${CMAKE_GENERATOR} STREQUAL "Visual Studio 7 .NET 2003" )



# This doesn't preserve the BasicUtils, Memory, etc. Find a way to fix this.
INSTALL_FILES(${NEKTAR_INCLUDE_DIR} FILES ${LibUtilityHeaders})




