
SET(BasicConstHeaders
	./BasicConst/NektarUnivConsts.hpp
	./BasicConst/NektarUnivTypeDefs.hpp
)

SET(LinearAlgebraHeaders
    ./LinearAlgebra/Blas.hpp
    ./LinearAlgebra/SparseBlas.hpp
    ./LinearAlgebra/BlasArray.hpp
    ./LinearAlgebra/BlockMatrix.hpp
    ./LinearAlgebra/CanGetRawPtr.hpp
    ./LinearAlgebra/Lapack.hpp
    ./LinearAlgebra/MatrixBase.hpp
    ./LinearAlgebra/MatrixFuncs.h
    ./LinearAlgebra/MatrixOperations.cpp
    ./LinearAlgebra/MatrixOperations.hpp
    ./LinearAlgebra/MatrixStorageType.h
    ./LinearAlgebra/MatrixType.h
    ./LinearAlgebra/MatrixVectorMultiplication.hpp
    ./LinearAlgebra/NekLinAlgAlgorithms.hpp
    ./LinearAlgebra/NekLinSys.hpp
    ./LinearAlgebra/NekMatrix.hpp
    ./LinearAlgebra/NekMatrixFwd.hpp
    ./LinearAlgebra/NekMatrixMetadata.hpp
    ./LinearAlgebra/NekPoint.hpp
    ./LinearAlgebra/NekTypeDefs.hpp
    ./LinearAlgebra/NekVector.hpp
    ./LinearAlgebra/NekVectorCommon.hpp
    ./LinearAlgebra/NekVectorConstantSized.hpp
    ./LinearAlgebra/NekVectorFwd.hpp
    ./LinearAlgebra/NekVectorMetadata.hpp
    ./LinearAlgebra/NekVectorTypeTraits.hpp
    ./LinearAlgebra/NekVectorVariableSized.hpp
    ./LinearAlgebra/PointerWrapper.h
    ./LinearAlgebra/ScaledMatrix.hpp
    ./LinearAlgebra/Space.h
    ./LinearAlgebra/StandardMatrix.hpp
    ./LinearAlgebra/TransF77.hpp
    ./LinearAlgebra/SparseStandardMatrix.hpp
)

SET(FoundationHeaders
    ./Foundations/Graph.h
    ./Foundations/Basis.h
    ./Foundations/Foundations.hpp
    ./Foundations/FourierPoints.h
    ./Foundations/GaussPoints.h
    ./Foundations/Interp.h
    ./Foundations/ManagerAccess.h
    ./Foundations/NodalTetElec.h
    ./Foundations/NodalTriElec.h
    ./Foundations/NodalTriFekete.h
    ./Foundations/Points.h
    ./Foundations/PolyEPoints.h
    ./Foundations/NodalUtil.h
    ./Foundations/NodalTriEvenlySpaced.h
    ./Foundations/NodalTetEvenlySpaced.h
)

SET(InterpreterHeaders
	  ./Interpreter/ExpressionEvaluator.h   
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
    ./ExpressionTemplates/BinaryExpressionEvaluator_BinaryConstant.hpp
    ./ExpressionTemplates/BinaryExpressionEvaluator_BinaryBinary.hpp
    ./ExpressionTemplates/BinaryExpressionEvaluatorFwd.hpp
    ./ExpressionTemplates/BinaryExpressionOperators.hpp
    ./ExpressionTemplates/BinaryExpressionPolicy.hpp
    ./ExpressionTemplates/BinaryExpressionPolicyFwd.hpp
    ./ExpressionTemplates/BinaryExpressionTraits.hpp
    ./ExpressionTemplates/BinaryOperators.hpp
    ./ExpressionTemplates/BinaryOperatorsFwd.hpp
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
    ./BasicUtils/ArrayEqualityComparison.cpp
    ./BasicUtils/BoostUtil.hpp
    ./BasicUtils/Concepts.hpp
    ./BasicUtils/ConsistentObjectAccess.hpp
    ./BasicUtils/ErrorUtil.hpp
    ./BasicUtils/ErrorUtil.cpp
    ./BasicUtils/mojo.hpp
    ./BasicUtils/NekManager.hpp
    ./BasicUtils/NekPtr.hpp
    ./BasicUtils/OperatorGenerators.hpp
    ./BasicUtils/RawType.hpp
    ./BasicUtils/SharedArray.hpp
    ./BasicUtils/Vmath.hpp
    ./BasicUtils/VmathArray.hpp
    ./BasicUtils/Metis.hpp
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
	./Foundations/Interp.cpp
	./Foundations/ManagerAccess.cpp
	./Foundations/NodalTetElec.cpp
	./Foundations/NodalTriElec.cpp
	./Foundations/NodalTriFekete.cpp
	./Foundations/Points.cpp
	./Foundations/PolyEPoints.cpp
	./Foundations/NodalUtil.cpp
        ./Foundations/NodalTriEvenlySpaced.cpp
        ./Foundations/NodalTetEvenlySpaced.cpp
)

SET(InterpreterSources
        ./Interpreter/ExpressionEvaluator.cpp
)       

SET(PolylibSources
	./Polylib/Polylib.cpp
)

SET(LinearAlgebraSources
	./LinearAlgebra/blas.cpp
	./LinearAlgebra/MatrixFuncs.cpp
#	./LinearAlgebra/lapack.cpp
)

SET(PrecompiledHeaders
	LibUtilities.h
)

SET(KernelHeaders
	./Kernel/kernel.h
)

SET(KernelSources
	./Kernel/kernel.cpp
)

SET(TimeIntegrationHeaders
        ./TimeIntegration/TimeIntegrationScheme.h
) 

SET(TimeIntegrationSources
        ./TimeIntegration/TimeIntegrationScheme.cpp
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

SOURCE_GROUP("Precompiled Headers" FILES ${PrecompiledHeaders})

SOURCE_GROUP("Kernel Headers" FILES ${KernelHeaders})
SOURCE_GROUP("Kernel Sources" FILES ${KernelSources})

SOURCE_GROUP("TimeIntegration Headers" FILES ${TimeIntegrationHeaders})
SOURCE_GROUP("TimeIntegration Sources" FILES ${TimeIntegrationSources})

SET(LibUtilityHeaders
	${Memory}
	${BasicUtils}
	${ExpressionTemplates}
	${PolyLibHeaders}
	${FoundationHeaders}
	${LinearAlgebraHeaders}
	${BasicConstHeaders}
	${InterpreterHeaders}
	${PrecompiledHeaders}
	${KernelHeaders}
        ${TimeIntegrationHeaders}
)

SET(LibUtilitySources
	${FoundationSources}
	${PolylibSources}
	${LinearAlgebraSources}
	${InterpreterSources}
	${KernelSources}
        ${TimeIntegrationSources}
)

ADD_NEKTAR_LIBRARY(LibUtilities STATIC ${LibUtilitySources} ${LibUtilityHeaders})

SET(PrecompiledHeaderName LibUtilities/LibUtilities.h)
SETUP_PRECOMPILED_HEADERS(LibUtilitySources PrecompiledHeaderName)

# This doesn't preserve the BasicUtils, Memory, etc. Find a way to fix this.
INSTALL(FILES ${BasicConstHeaders} DESTINATION ${NEKTAR_INCLUDE_DIR}/LibUtilities/BasicConst)
INSTALL(FILES ${LinearAlgebraHeaders} DESTINATION ${NEKTAR_INCLUDE_DIR}/LibUtilities/LinearAlgebra)
INSTALL(FILES ${FoundationHeaders} DESTINATION ${NEKTAR_INCLUDE_DIR}/LibUtilities/Foundations)
INSTALL(FILES ${InterpreterHeaders} DESTINATION ${NEKTAR_INCLUDE_DIR}/LibUtilities/Interpreter)
INSTALL(FILES ${PolyLibHeaders} DESTINATION ${NEKTAR_INCLUDE_DIR}/LibUtilities/Polylib)
INSTALL(FILES ${ExpressionTemplates} DESTINATION ${NEKTAR_INCLUDE_DIR}/LibUtilities/ExpressionTemplates)
INSTALL(FILES ${BasicUtils} DESTINATION ${NEKTAR_INCLUDE_DIR}/LibUtilities/BasicUtils)
INSTALL(FILES ${Memory} DESTINATION ${NEKTAR_INCLUDE_DIR}/LibUtilities/Memory)
INSTALL(FILES ${KernelHeaders} DESTINATION ${NEKTAR_INCLUDE_DIR}/LibUtilities/Kernel)
INSTALL(FILES ${TimeIntegrationHeaders} DESTINATION ${NEKTAR_INCLUDE_DIR}/LibUtilities/TimeIntegration)

