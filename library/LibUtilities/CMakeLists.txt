
SET(BasicConstHeaders
	./BasicConst/NektarUnivConsts.hpp
	./BasicConst/NektarUnivTypeDefs.hpp
)

SET(BasicUtilsHeaders
    ./BasicUtils/ArrayPolicies.hpp
    ./BasicUtils/BoostUtil.hpp
    ./BasicUtils/Concepts.hpp
    ./BasicUtils/ConsistentObjectAccess.hpp
    ./BasicUtils/Equation.h
    ./BasicUtils/FieldIO.h
    ./BasicUtils/ErrorUtil.hpp
    ./BasicUtils/MeshPartition.h
    ./BasicUtils/NekManager.hpp
    ./BasicUtils/NekFactory.hpp
    ./BasicUtils/NekPtr.hpp
    ./BasicUtils/OperatorGenerators.hpp
    ./BasicUtils/ParseUtils.hpp
    ./BasicUtils/Timer.h
    ./BasicUtils/RawType.hpp
    ./BasicUtils/SessionReader.h
    ./BasicUtils/ShapeType.hpp
    ./BasicUtils/SharedArray.hpp
    ./BasicUtils/Vmath.hpp
    ./BasicUtils/VmathArray.hpp
    ./BasicUtils/Metis.hpp
    ./BasicUtils/XmlUtil.h
)

SET(BasicUtilsSources
    ./BasicUtils/ArrayEqualityComparison.cpp
    ./BasicUtils/Equation.cpp
    ./BasicUtils/FieldIO.cpp
    ./BasicUtils/MeshPartition.cpp
    ./BasicUtils/SessionReader.cpp
    ./BasicUtils/Timer.cpp
    ./BasicUtils/Vmath.cpp
    ./BasicUtils/XmlUtil.cpp
)


SET(CommunicationHeaders
    ./Communication/Comm.h
    ./Communication/CommSerial.h
    ./Communication/GsLib.hpp
    ./Communication/Transposition.h
    ./Communication/Xxt.hpp
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

SET(FoundationHeaders
    ./Foundations/Basis.h
    ./Foundations/BasisType.h
    ./Foundations/BLPoints.h
    ./Foundations/FoundationsFwd.hpp
    ./Foundations/Foundations.hpp
    ./Foundations/FourierPoints.h
    ./Foundations/FourierSingleModePoints.h
    ./Foundations/GaussPoints.h
    ./Foundations/Graph.h
    ./Foundations/InterpCoeff.h
    ./Foundations/Interp.h
    ./Foundations/ManagerAccess.h
    ./Foundations/NodalPrismEvenlySpaced.h
    ./Foundations/NodalTetElecData.h
    ./Foundations/NodalTetElec.h
    ./Foundations/NodalTetEvenlySpaced.h
    ./Foundations/NodalTriElecData.h
    ./Foundations/NodalTriElec.h
    ./Foundations/NodalTriEvenlySpaced.h
    ./Foundations/NodalTriFeketeData.h
    ./Foundations/NodalTriFekete.h
    ./Foundations/NodalUtil.h
    ./Foundations/PhysGalerkinProject.h
    ./Foundations/Points.h
    ./Foundations/PointsType.h
    ./Foundations/PolyEPoints.h
)

SET(FoundationSources
	./Foundations/Basis.cpp
	./Foundations/BLPoints.cpp
	./Foundations/FourierPoints.cpp
	./Foundations/FourierSingleModePoints.cpp
	./Foundations/GaussPoints.cpp
	./Foundations/Graph.cpp
	./Foundations/InterpCoeff.cpp
	./Foundations/Interp.cpp
	./Foundations/ManagerAccess.cpp
    ./Foundations/NodalPrismEvenlySpaced.cpp
	./Foundations/NodalTetElec.cpp
    ./Foundations/NodalTetEvenlySpaced.cpp
	./Foundations/NodalTriElec.cpp
    ./Foundations/NodalTriEvenlySpaced.cpp
	./Foundations/NodalTriFekete.cpp
	./Foundations/NodalUtil.cpp
    ./Foundations/PhysGalerkinProject.cpp       
	./Foundations/Points.cpp
	./Foundations/PolyEPoints.cpp
)

SET(InterpreterHeaders
	  ./Interpreter/AnalyticExpressionEvaluator.hpp   
)       

SET(InterpreterSources
        ./Interpreter/AnalyticExpressionEvaluator.cpp
)    

SET(KernelHeaders
	./Kernel/kernel.h
)

SET(KernelSources
	./Kernel/kernel.cpp
)   

SET(LinearAlgebraHeaders
    ./LinearAlgebra/Arpack.hpp
    ./LinearAlgebra/BlasArray.hpp
    ./LinearAlgebra/Blas.hpp
    ./LinearAlgebra/BlasOverrideUtil.hpp
    ./LinearAlgebra/BlockMatrix.hpp
    ./LinearAlgebra/CanGetRawPtr.hpp
    ./LinearAlgebra/DgemmOverride.hpp
    ./LinearAlgebra/DgemvOverride.hpp
    ./LinearAlgebra/ExplicitInstantiation.h
    ./LinearAlgebra/IsDgemmTraits.hpp
    ./LinearAlgebra/Lapack.hpp
    ./LinearAlgebra/MatrixBase.hpp
    ./LinearAlgebra/MatrixFuncs.h
    ./LinearAlgebra/MatrixOperationsDeclarations.hpp
    ./LinearAlgebra/MatrixOperations.hpp
    ./LinearAlgebra/MatrixSize.hpp
    ./LinearAlgebra/MatrixStorageType.h
    ./LinearAlgebra/MatrixType.h
    ./LinearAlgebra/MatrixVectorMultiplication.hpp
    ./LinearAlgebra/NekLinAlgAlgorithms.hpp
    ./LinearAlgebra/NekLinSys.hpp
    ./LinearAlgebra/NekMatrixFwd.hpp
    ./LinearAlgebra/NekMatrix.hpp
    ./LinearAlgebra/NekMatrixMetadata.hpp
    ./LinearAlgebra/NekPoint.hpp
    ./LinearAlgebra/NekTypeDefs.hpp
    ./LinearAlgebra/NekVectorCommon.hpp
    ./LinearAlgebra/NekVectorConstantSized.hpp
    ./LinearAlgebra/NekVectorFwd.hpp
    ./LinearAlgebra/NekVector.hpp
    ./LinearAlgebra/NekVectorMetadata.hpp
    ./LinearAlgebra/NekVectorTypeTraits.hpp
    ./LinearAlgebra/NekVectorVariableSized.hpp
    ./LinearAlgebra/PointerWrapper.h
    ./LinearAlgebra/ScaledMatrix.hpp
    ./LinearAlgebra/Space.h
    ./LinearAlgebra/SparseBlas.hpp
    ./LinearAlgebra/SparseStandardMatrix.hpp
    ./LinearAlgebra/StandardMatrix.hpp
    ./LinearAlgebra/TransF77.hpp
)

SET(LinearAlgebraSources
	./LinearAlgebra/blas.cpp
    ./LinearAlgebra/BlockMatrix.cpp
    ./LinearAlgebra/MatrixBase.cpp
    ./LinearAlgebra/MatrixFuncs.cpp
    ./LinearAlgebra/MatrixOperations.cpp
    ./LinearAlgebra/MatrixVectorMultiplication.cpp
    ./LinearAlgebra/NekVector.cpp
    ./LinearAlgebra/ScaledMatrix.cpp
    ./LinearAlgebra/StandardMatrix.cpp
)

SET(MemoryHeaders
	./Memory/ThreadSpecificPool.hpp
	./Memory/NekMemoryManager.hpp	
)

SET(MemorySources
	./Memory/ThreadSpecificPool.cpp
)    

SET(PolyLibHeaders
	./Polylib/Polylib.h
)	

SET(PolylibSources
	./Polylib/Polylib.cpp
)

SET(TimeIntegrationHeaders
        ./TimeIntegration/TimeIntegrationScheme.h
) 

SET(TimeIntegrationSources
        ./TimeIntegration/TimeIntegrationScheme.cpp
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
SOURCE_GROUP("Memory Headers" FILES ${MemoryHeaders})
SOURCE_GROUP("Memory Sources" FILES ${MemorySources})


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
    ${MemoryHeaders}
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
    ${MemorySources}
)


ADD_DEFINITIONS(-DLIB_UTILITIES_EXPORTS)
ADD_NEKTAR_LIBRARY(LibUtilities lib ${NEKTAR_LIBRARY_TYPE} ${LibUtilitySources} ${LibUtilityHeaders})

TARGET_LINK_LIBRARIES(LibUtilities
    ${Boost_THREAD_LIBRARY}
    ${Boost_IOSTREAMS_LIBRARY}
    ${Boost_DATE_TIME_LIBRARY}
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${Boost_ZLIB_LIBRARY}
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
INSTALL(FILES ${MemoryHeaders} DESTINATION ${NEKTAR_INCLUDE_DIR}/LibUtilities/Memory COMPONENT dev)
INSTALL(FILES ${KernelHeaders} DESTINATION ${NEKTAR_INCLUDE_DIR}/LibUtilities/Kernel COMPONENT dev)
INSTALL(FILES ${TimeIntegrationHeaders} DESTINATION ${NEKTAR_INCLUDE_DIR}/LibUtilities/TimeIntegration COMPONENT dev)
INSTALL(FILES ${BaseHeaders} DESTINATION ${NEKTAR_INCLUDE_DIR}/LibUtilities COMPONENT dev)
INSTALL(FILES ${FFTHeaders} DESTINATION ${NEKTAR_INCLUDE_DIR}/LibUtilities/FFT COMPONENT dev)
