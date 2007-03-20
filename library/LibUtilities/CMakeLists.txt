SET(LibUtilitySources
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
./Polylib/Polylib.cpp
./LinearAlgebra/blas.cpp
./LinearAlgebra/lapack.cpp
)

SET(LibUtilityHeaders
./BasicConst/NektarUnivConsts.hpp
./BasicConst/NektarUnivTypeDefs.hpp
./LinearAlgebra/NekLinSys.h
./LinearAlgebra/NekLinSys.hpp
./LinearAlgebra/NekLinAlg.h
./LinearAlgebra/NekLinAlgAlgorithms.hpp
./LinearAlgebra/BlockMat.h
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
./Polylib/Polylib.h
./LinearAlgebra/NekMatrix.hpp
./LinearAlgebra/NekMatrixImpl.hpp
./LinearAlgebra/NekMatrixMetadata.hpp
./LinearAlgebra/NekVector.hpp
./LinearAlgebra/NekVectorFwd.hpp
./LinearAlgebra/NekVectorConstantSized.hpp
./LinearAlgebra/NekVectorVariableSized.hpp
./LinearAlgebra/NekVectorMetadata.hpp
./LinearAlgebra/NekVectorTypeTraits.hpp
./LinearAlgebra/NekVectorCommon.hpp
./LinearAlgebra/NekPoint.hpp
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
INSTALL(FILES ${LibUtilityHeaders} DESTINATION ${NEKTAR_INCLUDE_DIR})
#INSTALL(TARGETS LibUtilities ARCHIVE DESTINATION ${NEKTAR_LIB_DIR})



