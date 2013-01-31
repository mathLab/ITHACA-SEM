
SET(BasicConstHeaders
	./BasicConst/NektarUnivConsts.hpp
	./BasicConst/NektarUnivTypeDefs.hpp
)

SET(LinearAlgebraHeaders
    ./LinearAlgebra/Arpack.hpp
    ./LinearAlgebra/Blas.hpp
    ./LinearAlgebra/BlasOverrideUtil.hpp
    ./LinearAlgebra/SparseBlas.hpp
    ./LinearAlgebra/BlasArray.hpp
    ./LinearAlgebra/BlockMatrix.hpp
    ./LinearAlgebra/CanGetRawPtr.hpp
    ./LinearAlgebra/DgemmOverride.hpp
    ./LinearAlgebra/DgemvOverride.hpp
    ./LinearAlgebra/ExplicitInstantiation.h
    ./LinearAlgebra/Lapack.hpp
    ./LinearAlgebra/MatrixBase.hpp
    ./LinearAlgebra/MatrixFuncs.h
    ./LinearAlgebra/MatrixOperationsDeclarations.hpp
    ./LinearAlgebra/MatrixOperations.cpp
    ./LinearAlgebra/MatrixOperations.hpp
    ./LinearAlgebra/MatrixSize.hpp
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
    ./Foundations/BasisType.h
    ./Foundations/Foundations.hpp
    ./Foundations/FourierPoints.h
    ./Foundations/FourierSingleModePoints.h
    ./Foundations/BLPoints.h
    ./Foundations/GaussPoints.h
    ./Foundations/Interp.h
    ./Foundations/InterpCoeff.h
    ./Foundations/ManagerAccess.h
    ./Foundations/NodalTetElec.h
    ./Foundations/NodalTriElec.h
    ./Foundations/NodalTriFekete.h
    ./Foundations/Points.h
    ./Foundations/PointsType.h
    ./Foundations/PolyEPoints.h
    ./Foundations/NodalUtil.h
    ./Foundations/NodalTriEvenlySpaced.h
    ./Foundations/NodalTetEvenlySpaced.h
    ./Foundations/NodalPrismEvenlySpaced.h
    ./Foundations/PhysGalerkinProject.h
)

SET(InterpreterHeaders
	  ./Interpreter/AnalyticExpressionEvaluator.hpp   
#	  ./Interpreter/ExpressionEvaluator.h
)       

SET(PolyLibHeaders
	./Polylib/Polylib.h
)	

SET(CommunicationHeaders
    ./Communication/Comm.h
    ./Communication/CommSerial.h
    ./Communication/GsLib.hpp
    ./Communication/Transposition.h
    ./Communication/Xxt.hpp
)

SET(ExpressionTemplates
    ../ExpressionTemplates/AssociativeTraits.hpp
    ../ExpressionTemplates/AssociativeTransform.hpp
    ../ExpressionTemplates/BackwardInverseTransform.hpp
    ../ExpressionTemplates/CommutativeTraits.hpp
    ../ExpressionTemplates/CommutativeTransform.hpp
    ../ExpressionTemplates/CreateFromTree.hpp
    ../ExpressionTemplates/ExpressionEvaluator.hpp
    ../ExpressionTemplates/ExpressionTemplates.hpp
    ../ExpressionTemplates/ForwardInverseTransform.hpp
    ../ExpressionTemplates/InverseOperatorTypeTraits.hpp
    ../ExpressionTemplates/InvertNode.hpp
    ../ExpressionTemplates/Node.hpp
    ../ExpressionTemplates/Operators.hpp
    ../ExpressionTemplates/PerformCommutativeTransformIfNeeded.hpp
    ../ExpressionTemplates/PushDownUnaryNodes.hpp
    ../ExpressionTemplates/RemoveAllUnecessaryTemporaries.hpp
    ../ExpressionTemplates/SortAssociativeCommutativeClusters.hpp
)

SET(BasicUtilsHeaders
    ./BasicUtils/ArrayPolicies.hpp
    ./BasicUtils/ArrayEqualityComparison.cpp
    ./BasicUtils/BoostUtil.hpp
    ./BasicUtils/Concepts.hpp
    ./BasicUtils/ConsistentObjectAccess.hpp
    ./BasicUtils/Equation.h
    ./BasicUtils/ErrorUtil.hpp
    ./BasicUtils/MeshPartition.h
    ./BasicUtils/NekManager.hpp
    ./BasicUtils/NekFactory.hpp
    ./BasicUtils/NekPtr.hpp
    ./BasicUtils/OperatorGenerators.hpp
    ./BasicUtils/ParseUtils.hpp
    ./BasicUtils/Timer.h
    ./BasicUtils/Timer.cpp
    ./BasicUtils/RawType.hpp
    ./BasicUtils/SessionReader.h
    ./BasicUtils/SharedArray.hpp
    ./BasicUtils/Vmath.hpp
    ./BasicUtils/VmathArray.hpp
    ./BasicUtils/Metis.hpp
    ./BasicUtils/XmlUtil.h
    ./BasicUtils/XmlUtil.cpp
)

SET(BasicUtilsSources
    ./BasicUtils/Equation.cpp
    ./BasicUtils/MeshPartition.cpp
    ./BasicUtils/SessionReader.cpp
    ./BasicUtils/Vmath.cpp
)

SET(Memory
	./Memory/ThreadSpecificPool.cpp	
	./Memory/ThreadSpecificPool.hpp
	./Memory/NekMemoryManager.hpp	
)


SET(FoundationSources
	./Foundations/Basis.cpp
	./Foundations/Graph.cpp
	./Foundations/FourierPoints.cpp
	./Foundations/FourierSingleModePoints.cpp
	./Foundations/BLPoints.cpp
	./Foundations/GaussPoints.cpp
	./Foundations/Interp.cpp
	./Foundations/InterpCoeff.cpp
	./Foundations/ManagerAccess.cpp
	./Foundations/NodalTetElec.cpp
	./Foundations/NodalTriElec.cpp
	./Foundations/NodalTriFekete.cpp
	./Foundations/Points.cpp
	./Foundations/PolyEPoints.cpp
	./Foundations/NodalUtil.cpp
        ./Foundations/NodalTriEvenlySpaced.cpp
        ./Foundations/NodalTetEvenlySpaced.cpp
        ./Foundations/NodalPrismEvenlySpaced.cpp
        ./Foundations/PhysGalerkinProject.cpp       
)

SET(InterpreterSources
        ./Interpreter/AnalyticExpressionEvaluator.cpp
#        ./Interpreter/ExpressionEvaluator.cpp
)       

SET(PolylibSources
	./Polylib/Polylib.cpp
)

SET(CommunicationSources
    ./Communication/Comm.cpp
    ./Communication/CommSerial.cpp
    ./Communication/Transposition.cpp
)

IF(NEKTAR_USE_MPI)
    SET(CommunicationHeaders ${CommunicationHeaders} 
        ./Communication/CommMpi.h)
    SET(CommunicationSources ${CommunicationSources}
        ./Communication/CommMpi.cpp)
ENDIF(NEKTAR_USE_MPI)

SET(LinearAlgebraSources
	./LinearAlgebra/blas.cpp
        ./LinearAlgebra/BlockMatrix.cpp
        ./LinearAlgebra/MatrixBase.cpp
	./LinearAlgebra/MatrixFuncs.cpp
        ./LinearAlgebra/MatrixVectorMultiplication.cpp
        ./LinearAlgebra/NekVector.cpp
        ./LinearAlgebra/ScaledMatrix.cpp
        ./LinearAlgebra/StandardMatrix.cpp
#	./LinearAlgebra/lapack.cpp
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


SET(FFTHeaders
       ./FFT/NektarFFT.h
)	

SET(FFTSources
       ./FFT/NektarFFT.cpp
)

IF(NEKTAR_USE_FFTW)
    SET(FFTHeaders ${FFTHeaders} ./FFT/NekFFTW.h)
    SET(FFTSources ${FFTSources} ./FFT/NekFFTW.cpp)
ENDIF(NEKTAR_USE_FFTW)

SET(BaseHeaders
	LibUtilitiesDeclspec.h
)

SOURCE_GROUP("Basic Const" FILES ${BasicConstHeaders})
SOURCE_GROUP("Basic Utils Headers" FILES ${BasicUtilsHeaders})
SOURCE_GROUP("Basic Utils Sources" FILES ${BasicUtilsSources})
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

SOURCE_GROUP("Communication Headers" FILES ${CommunicationHeaders})
SOURCE_GROUP("Communication Sources" FILES ${CommunicationSources})

SOURCE_GROUP("Kernel Headers" FILES ${KernelHeaders})
SOURCE_GROUP("Kernel Sources" FILES ${KernelSources})

SOURCE_GROUP("TimeIntegration Headers" FILES ${TimeIntegrationHeaders})
SOURCE_GROUP("TimeIntegration Sources" FILES ${TimeIntegrationSources})

SOURCE_GROUP("FFT Headers" FILES ${FFTHeaders})
SOURCE_GROUP("FFT Sources" FILES ${FFTSources})  

SOURCE_GROUP("Base" FILES ${BaseHeaders})

SET(LibUtilityHeaders
	${Memory}
	${BasicUtilsHeaders}
	${ExpressionTemplates}
	${PolyLibHeaders}
    ${CommunicationHeaders}
	${FoundationHeaders}
	${LinearAlgebraHeaders}
	${BasicConstHeaders}
	${InterpreterHeaders}
	${KernelHeaders}
        ${TimeIntegrationHeaders}
        ${FFTHeaders}
        ${BaseHeaders}
)

SET(LibUtilitySources
    ${BasicUtilsSources}
    ${FoundationSources}
    ${PolylibSources}
    ${CommunicationSources}
	${LinearAlgebraSources}
	${InterpreterSources}
	${KernelSources}
    ${TimeIntegrationSources}
  	${FFTSources}
)


ADD_DEFINITIONS(-DLIB_UTILITIES_EXPORTS)
ADD_NEKTAR_LIBRARY(LibUtilities lib ${NEKTAR_LIBRARY_TYPE} ${LibUtilitySources} ${LibUtilityHeaders})

TARGET_LINK_LIBRARIES(LibUtilities
	${Boost_THREAD_LIBRARY}
	${Boost_DATE_TIME_LIBRARY}
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
	${TINYXML_LIB}
)
SET_LAPACK_LINK_LIBRARIES(LibUtilities)

# This doesn't preserve the BasicUtils, Memory, etc. Find a way to fix this.
INSTALL(FILES ${BasicConstHeaders} DESTINATION ${NEKTAR_INCLUDE_DIR}/LibUtilities/BasicConst COMPONENT dev)
INSTALL(FILES ${LinearAlgebraHeaders} DESTINATION ${NEKTAR_INCLUDE_DIR}/LibUtilities/LinearAlgebra COMPONENT dev)
INSTALL(FILES ${FoundationHeaders} DESTINATION ${NEKTAR_INCLUDE_DIR}/LibUtilities/Foundations COMPONENT dev)
INSTALL(FILES ${InterpreterHeaders} DESTINATION ${NEKTAR_INCLUDE_DIR}/LibUtilities/Interpreter COMPONENT dev)
INSTALL(FILES ${PolyLibHeaders} DESTINATION ${NEKTAR_INCLUDE_DIR}/LibUtilities/Polylib COMPONENT dev)
INSTALL(FILES ${CommunicationHeaders} DESTINATION ${NEKTAR_INCLUDE_DIR}/LibUtilities/Communication COMPONENT dev)
INSTALL(FILES ${ExpressionTemplates} DESTINATION ${NEKTAR_INCLUDE_DIR}/ExpressionTemplates COMPONENT dev)
INSTALL(FILES ${BasicUtilsHeaders} DESTINATION ${NEKTAR_INCLUDE_DIR}/LibUtilities/BasicUtils COMPONENT dev)
INSTALL(FILES ${Memory} DESTINATION ${NEKTAR_INCLUDE_DIR}/LibUtilities/Memory COMPONENT dev)
INSTALL(FILES ${KernelHeaders} DESTINATION ${NEKTAR_INCLUDE_DIR}/LibUtilities/Kernel COMPONENT dev)
INSTALL(FILES ${TimeIntegrationHeaders} DESTINATION ${NEKTAR_INCLUDE_DIR}/LibUtilities/TimeIntegration COMPONENT dev)
INSTALL(FILES ${BaseHeaders} DESTINATION ${NEKTAR_INCLUDE_DIR}/LibUtilities COMPONENT dev)
INSTALL(FILES ${FFTHeaders} DESTINATION ${NEKTAR_INCLUDE_DIR}/LibUtilities/FFT COMPONENT dev)
