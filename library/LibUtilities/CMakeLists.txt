
SET(BasicConstHeaders
    ./BasicConst/NektarUnivConsts.hpp
    ./BasicConst/NektarUnivTypeDefs.hpp
    ./BasicConst/GitRevision.h
)

SET(BasicUtilsHeaders
    ./BasicUtils/ArrayPolicies.hpp
    ./BasicUtils/BoostUtil.hpp
    ./BasicUtils/Concepts.hpp
    ./BasicUtils/ConsistentObjectAccess.hpp
    ./BasicUtils/Equation.h
    ./BasicUtils/FieldIO.h
    ./BasicUtils/FileSystem.h
    ./BasicUtils/ErrorUtil.hpp
    ./BasicUtils/MeshPartition.h
    ./BasicUtils/MeshPartitionMetis.h
    ./BasicUtils/NekManager.hpp
    ./BasicUtils/NekFactory.hpp
    ./BasicUtils/NekPtr.hpp
    ./BasicUtils/OperatorGenerators.hpp
    ./BasicUtils/ParseUtils.hpp
    ./BasicUtils/Progressbar.hpp
    ./BasicUtils/PtsField.h
    ./BasicUtils/PtsIO.h
    ./BasicUtils/Timer.h
    ./BasicUtils/RawType.hpp
    ./BasicUtils/SessionReader.h
    ./BasicUtils/ShapeType.hpp
    ./BasicUtils/SharedArray.hpp
    ./BasicUtils/Thread.h
    ./BasicUtils/ThreadBoost.h
    ./BasicUtils/Vmath.hpp
    ./BasicUtils/VmathArray.hpp
    ./BasicUtils/Metis.hpp
    ./BasicUtils/XmlUtil.h
)

SET(BasicUtilsSources
    ./BasicUtils/ArrayEqualityComparison.cpp
    ./BasicUtils/Equation.cpp
    ./BasicUtils/FieldIO.cpp
    ./BasicUtils/FileSystem.cpp
    ./BasicUtils/MeshPartition.cpp
    ./BasicUtils/MeshPartitionMetis.cpp
    ./BasicUtils/PtsField.cpp
    ./BasicUtils/PtsIO.cpp
    ./BasicUtils/SessionReader.cpp
    ./BasicUtils/Thread.cpp
    ./BasicUtils/ThreadBoost.cpp
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

IF(NEKTAR_USE_SCOTCH)
    SET(BasicUtilsHeaders ${BasicUtilsHeaders}
        ./BasicUtils/MeshPartitionScotch.h)
    SET(BasicUtilsSources ${BasicUtilsSources}
        ./BasicUtils/MeshPartitionScotch.cpp)
ENDIF(NEKTAR_USE_SCOTCH)

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
    ./LinearAlgebra/StandardMatrix.hpp
    ./LinearAlgebra/TransF77.hpp

    ./LinearAlgebra/StorageSmvBsr.hpp
    ./LinearAlgebra/NistSparseDescriptors.hpp
    ./LinearAlgebra/SparseDiagBlkMatrix.hpp
    ./LinearAlgebra/SparseMatrix.hpp
    ./LinearAlgebra/SparseMatrixFwd.hpp
    ./LinearAlgebra/SparseUtils.hpp
)

IF(NEKTAR_USING_SMV)
    SET(LinearAlgebraHeaders ${LinearAlgebraHeaders} ./LinearAlgebra/LibSMV.hpp)
ENDIF(NEKTAR_USING_SMV)

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
    ./LinearAlgebra/SparseUtils.cpp
    ./LinearAlgebra/StorageSmvBsr.cpp
    ./LinearAlgebra/SparseDiagBlkMatrix.cpp
    ./LinearAlgebra/SparseMatrix.cpp
)

IF(NEKTAR_USING_SMV)
    SET(LinearAlgebraSources ${LinearAlgebraSources} ./LinearAlgebra/LibSMV.cpp)
ENDIF(NEKTAR_USING_SMV)

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
        ./TimeIntegration/TimeIntegrationWrapper.h
)

SET(TimeIntegrationSources
        ./TimeIntegration/TimeIntegrationScheme.cpp
        ./TimeIntegration/TimeIntegrationWrapper.cpp
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

IF(NEKTAR_USE_OCC)
    SET(CADSystemSources
        ./CADSystem/CADSystem.cpp
        ./CADSystem/CADCurve.cpp
        ./CADSystem/CADSurf.cpp
        )

    SET(CADSystemHeaders
        ./CADSystem/CADSystem.h
        ./CADSystem/CADCurve.h
        ./CADSystem/CADSurf.h
        ./CADSystem/OpenCascade.h
        )

    SOURCE_GROUP("CAD System Headers" FILES ${CADSystemHeaders})
    SOURCE_GROUP("CAD System Sources" FILES ${CADSystemSources})

ENDIF(NEKTAR_USE_OCC)

SOURCE_GROUP("Base" FILES ${BaseHeaders})

SET(LibUtilityHeaders
    ${MemoryHeaders}
    ${BasicUtilsHeaders}
    ${ExpressionTemplates}
    ${PolyLibHeaders}
    ${CommunicationHeaders}
    ${FoundationHeaders}
    ${CADSystemHeaders}
    ${LinearAlgebraHeaders}
    ${BasicConstHeaders}
    ${InterpreterHeaders}
    ${KernelHeaders}
    ${TimeIntegrationHeaders}
    ${FFTHeaders}
    ${BaseHeaders}
)

# Define the git branch and SHA1 hash if we are in a git repository. This is
# stored in an automatically generated file so that the source tree is not built
# on commit/update.
SET(GitRevisionFile "${CMAKE_CURRENT_BINARY_DIR}/GitRevision.cpp")
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/BasicConst/GitRevision.cpp.in
               ${GitRevisionFile})

SET(LibUtilitySources
    ${BasicUtilsSources}
    ${FoundationSources}
    ${CADSystemSources}
    ${PolylibSources}
    ${CommunicationSources}
    ${LinearAlgebraSources}
    ${InterpreterSources}
    ${KernelSources}
    ${TimeIntegrationSources}
    ${FFTSources}
    ${MemorySources}
    ${GitRevisionFile}
)

ADD_DEFINITIONS(-DLIB_UTILITIES_EXPORTS)

ADD_NEKTAR_LIBRARY(LibUtilities lib ${NEKTAR_LIBRARY_TYPE} ${LibUtilitySources} ${LibUtilityHeaders})

ADD_DEPENDENCIES(LibUtilities modmetis-5.1.0 boost tinyxml-2.6.2 zlib-1.2.7)


### Libraries needed specifically by LibUtilities ###
### Static libraries do not need to be linked by
# Scotch
IF( NEKTAR_USE_SCOTCH )
    TARGET_LINK_LIBRARIES(LibUtilities LINK_PRIVATE
        ${SCOTCH_LIBRARY} ${SCOTCHERR_LIBRARY})
    ADD_DEPENDENCIES(LibUtilities scotch-6.0.0)
ENDIF ()

# Metis
TARGET_LINK_LIBRARIES(LibUtilities LINK_PRIVATE ${METIS_LIB})

# FFTW
IF( NEKTAR_USE_FFTW )
    TARGET_LINK_LIBRARIES(LibUtilities LINK_PUBLIC ${FFTW_LIBRARY})
    ADD_DEPENDENCIES(LibUtilities fftw-3.2.2)
    SET_SOURCE_FILES_PROPERTIES(./FFT/NekFFTW.cpp
        PROPERTY COMPILE_FLAGS "-I${FFTW_INCLUDE_DIR}")
ENDIF( NEKTAR_USE_FFTW )

# ARPACK
IF( NEKTAR_USE_ARPACK )
    TARGET_LINK_LIBRARIES(LibUtilities LINK_PUBLIC ${ARPACK_LIBRARY})
ENDIF( NEKTAR_USE_ARPACK )

# Boost
TARGET_LINK_LIBRARIES(LibUtilities LINK_PUBLIC
    ${Boost_THREAD_LIBRARY}
    ${Boost_IOSTREAMS_LIBRARY}
    ${Boost_DATE_TIME_LIBRARY}
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${Boost_FILESYSTEM_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
    ${Boost_TIMER_LIBRARY}
    debug ${ZLIB_LIBRARY_DEBUG} optimized ${ZLIB_LIBRARY_RELEASE}
)

# TinyXML
TARGET_LINK_LIBRARIES(LibUtilities LINK_PUBLIC ${TINYXML_LIBRARY})
IF( ${CMAKE_SYSTEM} MATCHES "Linux.*" )
    TARGET_LINK_LIBRARIES(LibUtilities LINK_PUBLIC rt)
ENDIF()

# MPI
IF( NEKTAR_USE_MPI )
    TARGET_LINK_LIBRARIES(LibUtilities LINK_PRIVATE ${GSMPI_LIBRARY} ${XXT_LIBRARY})
    IF (NOT MPI_BUILTIN)
        TARGET_LINK_LIBRARIES(LibUtilities LINK_PUBLIC ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
        SET_TARGET_PROPERTIES(LibUtilities
            PROPERTIES LINK_FLAGS "${THE_LINK_FLAGS} ${MPI_LINK_FLAGS}")
        SET_TARGET_PROPERTIES(LibUtilities
            PROPERTIES COMPILE_FLAGS "${THE_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
    ENDIF()
    ADD_DEPENDENCIES(LibUtilities gsmpi-1.2)
ENDIF( NEKTAR_USE_MPI )

# Lapack and Blas
IF( NEKTAR_USE_BLAS_LAPACK )
    IF( NEKTAR_USE_MKL AND MKL_FOUND )
        TARGET_LINK_LIBRARIES(LibUtilities LINK_PUBLIC ${MKL} )
        SET_TARGET_PROPERTIES(LibUtilities
            PROPERTIES COMPILE_FLAGS "${THE_COMPILE_FLAGS} -DMKL_ILP64")
    ENDIF( NEKTAR_USE_MKL AND MKL_FOUND )

    IF( NEKTAR_USE_ACML AND ACML_FOUND )
        TARGET_LINK_LIBRARIES(LibUtilities LINK_PUBLIC ${ACML_TARGET_LINK_LIBRARIES}  )
    ENDIF( NEKTAR_USE_ACML AND ACML_FOUND )

    IF( NEKTAR_USE_ACCELERATE_FRAMEWORK )
        TARGET_LINK_LIBRARIES(LibUtilities LINK_PUBLIC ${ACCELERATE_FRAMEWORK_LINK_FLAGS})
    ENDIF ( NEKTAR_USE_ACCELERATE_FRAMEWORK )

    IF( NEKTAR_USE_WIN32_LAPACK )
        TARGET_LINK_LIBRARIES(LibUtilities LINK_PUBLIC ${WIN32_LAPACK} ${WIN32_BLAS})
        INSTALL(FILES ${WIN32_LAPACK_DLL} ${WIN32_BLAS_DLL}
                DESTINATION ${NEKTAR_BIN_DIR})
    ENDIF( NEKTAR_USE_WIN32_LAPACK )

    IF( NEKTAR_USE_OPENBLAS AND OPENBLAS_FOUND )
        TARGET_LINK_LIBRARIES(LibUtilities LINK_PUBLIC ${OPENBLAS})
    ENDIF( NEKTAR_USE_OPENBLAS AND OPENBLAS_FOUND )

    IF( NEKTAR_USE_SMV AND SMV_FOUND )
        TARGET_LINK_LIBRARIES(LibUtilities LINK_PUBLIC ${SMV_LIBRARY})
    ENDIF( NEKTAR_USE_SMV AND SMV_FOUND )

    IF( NEKTAR_USE_SYSTEM_BLAS_LAPACK )
        TARGET_LINK_LIBRARIES(LibUtilities LINK_PUBLIC ${NATIVE_LAPACK} ${NATIVE_BLAS})
    ENDIF( NEKTAR_USE_SYSTEM_BLAS_LAPACK )

ENDIF( NEKTAR_USE_BLAS_LAPACK )

IF( NEKTAR_USE_PETSC )
    TARGET_LINK_LIBRARIES(LibUtilities LINK_PRIVATE ${PETSC_LIBRARIES})
    ADD_DEPENDENCIES(LibUtilities petsc-3.5.2)
ENDIF( NEKTAR_USE_PETSC )

IF(NEKTAR_USE_OCC)
    SET(OCC_DEF LIN LININTEL HAVE_WOK_CONFIG_H HAVE_CONFIG_H CSFDB)
    SET_TARGET_PROPERTIES(LibUtilities PROPERTIES COMPILE_DEFINITIONS "${OCC_DEF}")
    TARGET_LINK_LIBRARIES(LibUtilities LINK_PUBLIC ${OCC_LIBS})
    ADD_DEPENDENCIES(LibUtilities opencascade-6.8)
ENDIF(NEKTAR_USE_OCC)

INSTALL(FILES ${ExpressionTemplates} DESTINATION ${NEKTAR_INCLUDE_DIR}/ExpressionTemplates COMPONENT dev)
INSTALL(DIRECTORY ./ DESTINATION ${NEKTAR_INCLUDE_DIR}/LibUtilities COMPONENT dev FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
