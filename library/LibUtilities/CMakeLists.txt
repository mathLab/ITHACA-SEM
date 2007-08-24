
SET(BasicConstHeaders
	./BasicConst/NektarUnivConsts.hpp
	./BasicConst/NektarUnivTypeDefs.hpp
)

SET(LinearAlgebraHeaders
    ./LinearAlgebra/Blas.hpp
    ./LinearAlgebra/BlockMatrix.hpp
    ./LinearAlgebra/DiagonalMatrixStoragePolicy.hpp
    ./LinearAlgebra/FullMatrixStoragePolicy.hpp
    ./LinearAlgebra/Lapack.hpp
    ./LinearAlgebra/MatrixBase.hpp
    ./LinearAlgebra/MatrixOperations.cpp
    ./LinearAlgebra/MatrixOperations.hpp
    ./LinearAlgebra/MatrixStoragePolicy.hpp
    ./LinearAlgebra/MatrixStorageType.h
    ./LinearAlgebra/MatrixType.h
    ./LinearAlgebra/NekLinAlgAlgorithms.hpp
    ./LinearAlgebra/NekLinSys.hpp
    ./LinearAlgebra/NekMatrix.hpp
    ./LinearAlgebra/NekMatrixObject.hpp
    ./LinearAlgebra/NekMatrixFwd.hpp
    ./LinearAlgebra/NekMatrixMetadata.hpp
    ./LinearAlgebra/NekPoint.hpp
    ./LinearAlgebra/NekTypeDefs.hpp
    ./LinearAlgebra/NekVector.hpp
    ./LinearAlgebra/NekVectorObject.hpp
    ./LinearAlgebra/NekVectorCommon.hpp
    ./LinearAlgebra/NekVectorConstantSized.hpp
    ./LinearAlgebra/NekVectorFwd.hpp
    ./LinearAlgebra/NekVectorMetadata.hpp
    ./LinearAlgebra/NekVectorTypeTraits.hpp
    ./LinearAlgebra/NekVectorVariableSized.hpp
    ./LinearAlgebra/NormalMatrix.hpp
    ./LinearAlgebra/PointerWrapper.h
    ./LinearAlgebra/ScaledMatrix.hpp
    ./LinearAlgebra/TriangularMatrixStoragePolicy.hpp
    ./LinearAlgebra/TransF77.hpp
)

SET(FoundationHeaders
    ./Foundations/Graph.h
    ./Foundations/Basis.h
    ./Foundations/FourierPoints.h
    ./Foundations/GaussPoints.h
    ./Foundations/ManagerAccess.h
    ./Foundations/NodalTetElec.h
    ./Foundations/NodalTriElec.h
    ./Foundations/NodalTriFekete.h
    ./Foundations/Points.h
    ./Foundations/PolyEPoints.h
    ./Foundations/NodalUtil.h
)

SET(InterpreterHeaders
        ./Interpreter/Interpreter.hpp
        ./Interpreter/defaults.h       
)       

SET(PolyLibHeaders
	./Polylib/Polylib.h
)	

SET(ExpressionTemplates
    ./ExpressionTemplates/Accumulator.hpp
    ./ExpressionTemplates/ArithmeticConcepts.hpp
    ./ExpressionTemplates/ArithmeticTraits.hpp
    ./ExpressionTemplates/AssignableTraits.hpp
    ./ExpressionTemplates/AssociativeExpressionTraits.hpp
    ./ExpressionTemplates/AssociativeTraits.hpp
    ./ExpressionTemplates/BinaryExpression.hpp
    ./ExpressionTemplates/BinaryExpressionEvaluator.hpp
    ./ExpressionTemplates/BinaryExpressionEvaluatorFwd.hpp
    ./ExpressionTemplates/BinaryExpressionOperators.hpp
    ./ExpressionTemplates/BinaryExpressionPolicy.hpp
    ./ExpressionTemplates/BinaryExpressionPolicyFwd.hpp
    ./ExpressionTemplates/BinaryExpressionTraits.hpp
    ./ExpressionTemplates/BinaryOperators.hpp
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
)

SET(BasicUtils
    ./BasicUtils/ArrayPolicies.hpp
    ./BasicUtils/BoostUtil.hpp
    ./BasicUtils/Concepts.hpp
    ./BasicUtils/ConsistentObjectAccess.hpp
    ./BasicUtils/ErrorUtil.hpp
    ./BasicUtils/mojo.hpp
    ./BasicUtils/NekManager.hpp
    ./BasicUtils/OperatorGenerators.hpp
    ./BasicUtils/SharedArray.hpp
    ./BasicUtils/Vmath.hpp
)

SET(Memory
	./Memory/ThreadSpecificPool.hpp
	./Memory/NekMemoryManager.hpp	
	./Memory/DeleteNothing.hpp
)


SET(FoundationSources
	./Foundations/Basis.cpp
	./Foundations/Graph.cpp
	./Foundations/FourierPoints.cpp
	./Foundations/GaussPoints.cpp
	./Foundations/ManagerAccess.cpp
	./Foundations/NodalTetElec.cpp
	./Foundations/NodalTriElec.cpp
	./Foundations/NodalTriFekete.cpp
	./Foundations/Points.cpp
	./Foundations/PolyEPoints.cpp
	./Foundations/NodalUtil.cpp
)

SET(InterpreterSources
        ./Interpreter/initial.y
)       

SET(PolylibSources
	./Polylib/Polylib.cpp
)

SET(LinearAlgebraSources
	./LinearAlgebra/blas.cpp
#	./LinearAlgebra/lapack.cpp
)

SOURCE_GROUP("Basic Const" FILES ${BasicConstHeaders})
SOURCE_GROUP("Basic Utils" FILES ${BasicUtils})
SOURCE_GROUP("Expression Templates" FILES ${ExpressionTemplates})
SOURCE_GROUP("Foundation Headers" FILES ${FoundationHeaders})
SOURCE_GROUP("Foundation Sources" FILES ${FoundationSources})
SOURCE_GROUP("Linear Algebra Headers" FILES ${LinearAlgebraHeaders})
SOURCE_GROUP("Linear Algebra Sources" FILES ${LinearAlgebraSources})
SOURCE_GROUP("Memory" FILES ${Memory})

SOURCE_GROUP("Polylib Headers" FILES ${PolyLibHeaders})
SOURCE_GROUP("Polylib Sources" FILES ${PolylibSources})

SOURCE_GROUP("Interpreter Headers" FILES ${InterpreterHeaders})
SOURCE_GROUP("Interpreter Sources" FILES ${InterpreterSources})

SET(LibUtilityHeaders
	${Memory}
	${BasicUtils}
	${ExpressionTemplates}
	${PolyLibHeaders}
	${FoundationHeaders}
	${LinearAlgebraHeaders}
	${BasicConstHeaders}
	${InterpreterHeaders}
)

SET(LibUtilitySources
	${FoundationSources}
	${PolylibSources}
	${LinearAlgebraSources}
	${InterpreterSources}
)


ADD_NEKTAR_LIBRARY(LibUtilities STATIC ${LibUtilitySources} ${LibUtilityHeaders})

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
INSTALL(FILES ${BasicConstHeaders} DESTINATION ${NEKTAR_INCLUDE_DIR}/LibUtilities/BasicConst)
INSTALL(FILES ${LinearAlgebraHeaders} DESTINATION ${NEKTAR_INCLUDE_DIR}/LibUtilities/LinearAlgebra)
INSTALL(FILES ${FoundationHeaders} DESTINATION ${NEKTAR_INCLUDE_DIR}/LibUtilities/Foundations)
INSTALL(FILES ${PolyLibHeaders} DESTINATION ${NEKTAR_INCLUDE_DIR}/LibUtilities/Polylib)
INSTALL(FILES ${ExpressionTemplates} DESTINATION ${NEKTAR_INCLUDE_DIR}/LibUtilities/ExpressionTemplates)
INSTALL(FILES ${BasicUtils} DESTINATION ${NEKTAR_INCLUDE_DIR}/LibUtilities/BasicUtils)
INSTALL(FILES ${Memory} DESTINATION ${NEKTAR_INCLUDE_DIR}/LibUtilities/Memory)
	
#INSTALL(TARGETS LibUtilities ARCHIVE DESTINATION ${NEKTAR_LIB_DIR})



