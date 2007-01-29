SET(LibUtilitySources
./LinearAlgebra/BlockMat.cpp
./Foundations/Graph.cpp
./Polylib/Polylib.cpp
)

SET(LibUtilityHeaders
./LinearAlgebra/NekLinSys.h
./LinearAlgebra/NekLinSys.hpp
./LinearAlgebra/NekLinAlg.h
./LinearAlgebra/NekLinAlgAlgorithms.hpp
./LinearAlgebra/BlockMat.h
./Foundations/Graph.h
./Polylib/Polylib.h
./LinearAlgebra/NekMatrix.hpp
./LinearAlgebra/NekMatrixImpl.hpp
./LinearAlgebra/NekMatrixMetadata.hpp
./LinearAlgebra/NekVector.hpp
./LinearAlgebra/NekVectorMetadata.hpp
./LinearAlgebra/NekPoint.hpp
./LinearAlgebra/blas.cpp
./LinearAlgebra/lapack.cpp
./ExpressionTemplates/Accumulator.hpp
./ExpressionTemplates/ArithmeticConcepts.hpp
./ExpressionTemplates/ArithmeticTraits.hpp
./ExpressionTemplates/AssignableTraits.hpp
./ExpressionTemplates/AssociativeExpressionTraits.hpp
./ExpressionTemplates/AssociativeTraits.hpp
./ExpressionTemplates/BinaryExpression.hpp
./ExpressionTemplates/BinaryExpressionEvaluator.hpp
./ExpressionTemplates/BinaryExpressionOperators.hpp
./ExpressionTemplates/BinaryExpressionPolicy.hpp
./ExpressionTemplates/BinaryExpressionPolicyFwd.hpp
./ExpressionTemplates/BinaryExpressionTraits.hpp
./ExpressionTemplates/CommutativeExpressionTraits.hpp
./ExpressionTemplates/CommutativeTraits.hpp
./ExpressionTemplates/ConstantExpression.hpp
./ExpressionTemplates/ConstantExpressionTraits.hpp
./ExpressionTemplates/Expression.hpp
./ExpressionTemplates/ExpressionConcepts.hpp
./ExpressionTemplates/ExpressionMetadata.hpp
./ExpressionTemplates/ExpressionReturnType.hpp
./ExpressionTemplates/ExpressionTemplateConcepts.hpp
./ExpressionTemplates/ExpressionTemplateOperators.hpp
./ExpressionTemplates/ExpressionTemplates.hpp
./ExpressionTemplates/ExpressionTraits.hpp
./ExpressionTemplates/NegateOp.hpp
./ExpressionTemplates/NullOp.hpp
./ExpressionTemplates/UnaryExpression.hpp
./ExpressionTemplates/UnaryExpressionTraits.hpp
./BasicUtils/BoostUtil.hpp
./BasicUtils/ConsistentObjectAccess.hpp
./BasicUtils/ErrorUtil.hpp
./BasicUtils/Lapack.hpp
./BasicUtils/Concepts.hpp
./BasicUtils/Vmath.hpp
./BasicUtils/Blas.hpp
./BasicUtils/TransF77.hpp
./Memory/ThreadSpecificPool.hpp
./Memory/NekMemoryManager.hpp
)

CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/LibUtilities/LinearAlgebra/blas.h.in ${PROJECT_BINARY_DIR}/LibUtilities/LinearAlgebra/blas.h)
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/LibUtilities/LinearAlgebra/lapack.h.in ${PROJECT_BINARY_DIR}/LibUtilities/LinearAlgebra/lapack.h)

ADD_NEKTAR_LIBRARY(LibUtilities STATIC ${LibUtilitySources} ${LibUtilityHeaders})

IF( ${CMAKE_GENERATOR} STREQUAL "Visual Studio 7 .NET 2003" )
	# Needed for M_PI to be visible in visual studio.
	ADD_DEFINITIONS(-D_USE_MATH_DEFINES)
ENDIF( ${CMAKE_GENERATOR} STREQUAL "Visual Studio 7 .NET 2003" )

IF( ${CMAKE_GENERATOR} STREQUAL "Visual Studio 8 2005" )
	# Needed for M_PI to be visible in visual studio.
	ADD_DEFINITIONS(-D_USE_MATH_DEFINES)
ENDIF( ${CMAKE_GENERATOR} STREQUAL "Visual Studio 8 2005" )

#IF( ${USE_MKL} AND ${MKL_FOUND} )
#        TARGET_LINK_LIBRARIES(${name}
#            optimized ${MKL_LAPACK} debug ${MKL_LAPACK}
#        )
#ENDIF( ${USE_MKL} AND ${MKL_FOUND} )

IF( USE_ATLAS AND ATLAS_FOUND )
        TARGET_LINK_LIBRARIES(LibUtilities
            ${ATLAS_LAPACK}
            ${ATLAS_CBLAS}
            ${ATLAS} )
ENDIF( USE_ATLAS AND ATLAS_FOUND )

# This doesn't preserve the BasicUtils, Memory, etc. Find a way to fix this.
INSTALL_FILES(${NEKTAR_INCLUDE_DIR} FILES ${LibUtilityHeaders})




