
SET(BasicConstHeaders
    ./BasicConst/NektarUnivConsts.hpp
    ./BasicConst/NektarUnivTypeDefs.hpp
    ./BasicConst/GitRevision.h
)

SET(BasicUtilsHeaders
    ./BasicUtils/ArrayPolicies.hpp
    ./BasicUtils/CheckedCast.hpp
    ./BasicUtils/ConsistentObjectAccess.hpp
    ./BasicUtils/CompressData.h
    ./BasicUtils/CsvIO.h
    ./BasicUtils/Equation.h
    ./BasicUtils/FieldIO.h
    ./BasicUtils/FieldIOXml.h
    ./BasicUtils/FileSystem.h
    ./BasicUtils/ErrorUtil.hpp
    ./BasicUtils/Interpolator.h
    ./BasicUtils/HashUtils.hpp
    ./BasicUtils/NekManager.hpp
    ./BasicUtils/NekFactory.hpp
    ./BasicUtils/ParseUtils.h
    ./BasicUtils/Progressbar.hpp
    ./BasicUtils/PtsField.h
    ./BasicUtils/PtsIO.h
    ./BasicUtils/Timer.h
    ./BasicUtils/RawType.hpp
    ./BasicUtils/RealComparison.hpp
    ./BasicUtils/SessionReader.h
    ./BasicUtils/ShapeType.hpp
    ./BasicUtils/SharedArray.hpp
    ./BasicUtils/Thread.h
    ./BasicUtils/ThreadBoost.h
    ./BasicUtils/Vmath.hpp
    ./BasicUtils/VmathArray.hpp
    ./BasicUtils/Metis.hpp
    ./BasicUtils/VtkUtil.hpp
)

SET(BasicUtilsSources
    ./BasicUtils/CompressData.cpp
    ./BasicUtils/CsvIO.cpp
    ./BasicUtils/Equation.cpp
    ./BasicUtils/ErrorUtil.cpp
    ./BasicUtils/FieldIO.cpp
    ./BasicUtils/FieldIOXml.cpp
    ./BasicUtils/FileSystem.cpp
    ./BasicUtils/Interpolator.cpp
    ./BasicUtils/ParseUtils.cpp
    ./BasicUtils/PtsField.cpp
    ./BasicUtils/PtsIO.cpp
    ./BasicUtils/SessionReader.cpp
    ./BasicUtils/Thread.cpp
    ./BasicUtils/ThreadBoost.cpp
    ./BasicUtils/Timer.cpp
    ./BasicUtils/Vmath.cpp
)

# HDF5
IF( NEKTAR_USE_HDF5 )
    SET(BasicUtilsHeaders ${BasicUtilsHeaders}
	./BasicUtils/H5.h
        ./BasicUtils/FieldIOHdf5.h)
    SET(BasicUtilsSources ${BasicUtilsSources}
	./BasicUtils/H5.cpp
        ./BasicUtils/FieldIOHdf5.cpp)
ENDIF( NEKTAR_USE_HDF5 )


SET(CommunicationHeaders
    ./Communication/Comm.h
    ./Communication/CommSerial.h
    ./Communication/CommDataType.h
    ./Communication/GsLib.hpp
    ./Communication/Transposition.h
    ./Communication/Xxt.hpp
)

SET(CommunicationSources
    ./Communication/Comm.cpp
    ./Communication/CommDataType.cpp
    ./Communication/CommSerial.cpp
    ./Communication/Transposition.cpp
)

IF(NEKTAR_USE_MPI)
    SET(CommunicationHeaders ${CommunicationHeaders}
        ./Communication/CommMpi.h)
    SET(CommunicationSources ${CommunicationSources}
        ./Communication/CommMpi.cpp)
ENDIF(NEKTAR_USE_MPI)

IF ( NEKTAR_USE_CWIPI )
    SET(CommunicationHeaders ${CommunicationHeaders}
        ./Communication/CommCwipi.h)
    SET(CommunicationSources ${CommunicationSources}
        ./Communication/CommCwipi.cpp)
ENDIF ()

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
    ./Foundations/NodalPrismElec.h
    ./Foundations/NodalTetElecData.h
    ./Foundations/NodalTetElec.h
    ./Foundations/NodalTetEvenlySpaced.h
    ./Foundations/NodalTriElecData.h
    ./Foundations/NodalTriElec.h
    ./Foundations/NodalTriEvenlySpaced.h
    ./Foundations/NodalTriFeketeData.h
    ./Foundations/NodalTriFekete.h
    ./Foundations/NodalTriSPIData.h
    ./Foundations/NodalTriSPI.h
    ./Foundations/NodalTetSPI.h
    ./Foundations/NodalPrismSPI.h
    ./Foundations/NodalTetSPIData.h
    ./Foundations/NodalQuadElec.h
    ./Foundations/NodalHexElec.h
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
    ./Foundations/NodalPrismElec.cpp
    ./Foundations/NodalTetElec.cpp
    ./Foundations/NodalTetEvenlySpaced.cpp
    ./Foundations/NodalTriElec.cpp
    ./Foundations/NodalTriEvenlySpaced.cpp
    ./Foundations/NodalTriFekete.cpp
    ./Foundations/NodalTriSPI.cpp
    ./Foundations/NodalTetSPI.cpp
    ./Foundations/NodalQuadElec.cpp
    ./Foundations/NodalHexElec.cpp
    ./Foundations/NodalPrismSPI.cpp
    ./Foundations/NodalUtil.cpp
    ./Foundations/PhysGalerkinProject.cpp
    ./Foundations/Points.cpp
    ./Foundations/PolyEPoints.cpp
)

SET(InterpreterHeaders
    ./Interpreter/Interpreter.h
)

SET(InterpreterSources
    ./Interpreter/Interpreter.cpp
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
    ./LinearAlgebra/BlockMatrix.hpp
    ./LinearAlgebra/CanGetRawPtr.hpp
    ./LinearAlgebra/ExplicitInstantiation.h
    ./LinearAlgebra/Lapack.hpp
    ./LinearAlgebra/MatrixBase.hpp
    ./LinearAlgebra/MatrixFuncs.h
    ./LinearAlgebra/MatrixOperations.hpp
    ./LinearAlgebra/MatrixStorageType.h
    ./LinearAlgebra/MatrixType.h
    ./LinearAlgebra/NekLinAlgAlgorithms.hpp
    ./LinearAlgebra/NekLinSys.hpp
    ./LinearAlgebra/NekMatrixFwd.hpp
    ./LinearAlgebra/NekMatrix.hpp
    ./LinearAlgebra/NekPoint.hpp
    ./LinearAlgebra/NekTypeDefs.hpp
    ./LinearAlgebra/NekVectorFwd.hpp
    ./LinearAlgebra/NekVector.hpp
    ./LinearAlgebra/PointerWrapper.h
    ./LinearAlgebra/ScaledMatrix.hpp
    ./LinearAlgebra/StandardMatrix.hpp
    ./LinearAlgebra/TransF77.hpp
    ./LinearAlgebra/StorageSmvBsr.hpp
    ./LinearAlgebra/NistSparseDescriptors.hpp
    ./LinearAlgebra/SparseDiagBlkMatrix.hpp
    ./LinearAlgebra/SparseMatrix.hpp
    ./LinearAlgebra/SparseMatrixFwd.hpp
    ./LinearAlgebra/SparseUtils.hpp
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
    ./LinearAlgebra/SparseUtils.cpp
    ./LinearAlgebra/StorageSmvBsr.cpp
    ./LinearAlgebra/SparseDiagBlkMatrix.cpp
    ./LinearAlgebra/SparseMatrix.cpp
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
        ./TimeIntegration/AdamsBashforthTimeIntegrationSchemes.h
        ./TimeIntegration/AdamsMoultonTimeIntegrationSchemes.h
        ./TimeIntegration/BDFImplicitTimeIntegrationSchemes.h
        ./TimeIntegration/EulerTimeIntegrationSchemes.h
        ./TimeIntegration/CNABTimeIntegrationScheme.h
        ./TimeIntegration/DIRKTimeIntegrationSchemes.h
        ./TimeIntegration/IMEXGearTimeIntegrationScheme.h
        ./TimeIntegration/IMEXTimeIntegrationSchemes.h
        ./TimeIntegration/IMEXdirkTimeIntegrationSchemes.h
        ./TimeIntegration/MCNABTimeIntegrationScheme.h
        ./TimeIntegration/RungeKuttaTimeIntegrationSchemes.h
        ./TimeIntegration/TimeIntegrationScheme.h
        ./TimeIntegration/TimeIntegrationSchemeData.h
        ./TimeIntegration/TimeIntegrationSolution.h
)

SET(TimeIntegrationSources
        ./TimeIntegration/SchemeInitializer.cpp
        ./TimeIntegration/TimeIntegrationScheme.cpp
        ./TimeIntegration/TimeIntegrationSchemeData.cpp
        ./TimeIntegration/TimeIntegrationSolution.cpp
)


SET(BaseHeaders
    LibUtilitiesDeclspec.h
)

SOURCE_GROUP("Basic Const" FILES ${BasicConstHeaders})
SOURCE_GROUP("Basic Utils Headers" FILES ${BasicUtilsHeaders})
SOURCE_GROUP("Basic Utils Sources" FILES ${BasicUtilsSources})
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

ADD_NEKTAR_LIBRARY(LibUtilities
    SOURCES ${LibUtilitySources}
    HEADERS ${LibUtilityHeaders}
    SUMMARY "Nektar++ LibUtilities library"
    DESCRIPTION "This library provides core routines including linear algebra and integration with ThirdParty libraries.")

ADD_DEPENDENCIES(LibUtilities boost tinyxml-2.6.2 zlib-1.2.7)

### Libraries needed specifically by LibUtilities ###

# CWIPI
IF ( NEKTAR_USE_CWIPI )
    TARGET_LINK_LIBRARIES(LibUtilities LINK_PRIVATE ${CWIPI_LIBRARY} ${CWIPI_LIBRARY_FVMC} ${CWIPI_LIBRARY_BFTC})
    ADD_DEPENDENCIES(LibUtilities cwipi-0.8.2)
ENDIF ()

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
    ADD_DEPENDENCIES(LibUtilities arpack-ng-1.0)
ENDIF( NEKTAR_USE_ARPACK )

# Boost
TARGET_LINK_LIBRARIES(LibUtilities LINK_PUBLIC
    ${Boost_THREAD_LIBRARY}
    ${Boost_IOSTREAMS_LIBRARY}
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${Boost_FILESYSTEM_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
    optimized ${ZLIB_LIBRARIES} debug ${ZLIB_LIBRARIES_DEBUG}
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
        TARGET_LINK_LIBRARIES(LibUtilities LINK_PUBLIC ${MPI_CXX_LIBRARIES})
        SET_TARGET_PROPERTIES(LibUtilities
            PROPERTIES LINK_FLAGS "${THE_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS}")
        SET_TARGET_PROPERTIES(LibUtilities
            PROPERTIES COMPILE_FLAGS "${THE_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
    ENDIF()

    ADD_DEPENDENCIES(LibUtilities gsmpi-1.2.1)

    IF (THIRDPARTY_BUILD_MPI)
        ADD_DEPENDENCIES(LibUtilities openmpi-1.10.3)
    ENDIF()
ENDIF( NEKTAR_USE_MPI )

# Lapack and Blas
IF(NOT BLAS_LAPACK_BUILTIN)
    TARGET_LINK_LIBRARIES(LibUtilities LINK_PUBLIC ${BLAS_LAPACK})
ENDIF()

IF(THIRDPARTY_BUILD_BLAS_LAPACK)
    ADD_DEPENDENCIES(LibUtilities lapack-3.7.0)
ENDIF()

IF( NEKTAR_USE_PETSC )
    TARGET_LINK_LIBRARIES(LibUtilities LINK_PRIVATE ${PETSC_LIBRARIES})
    TARGET_LINK_LIBRARIES(LibUtilities LINK_PUBLIC ${CMAKE_DL_LIBS})
    ADD_DEPENDENCIES(LibUtilities petsc-3.11.4)
ENDIF( NEKTAR_USE_PETSC )

IF (NEKTAR_BUILD_PYTHON)
    SUBDIRS(Python)
ENDIF()

# HDF5
IF( NEKTAR_USE_HDF5 )
    TARGET_LINK_LIBRARIES(LibUtilities LINK_PUBLIC ${HDF5_LIBRARIES})
    ADD_DEPENDENCIES(LibUtilities hdf5-1.8.16)
ENDIF( NEKTAR_USE_HDF5 )

INSTALL(DIRECTORY ./
    DESTINATION ${NEKTAR_INCLUDE_DIR}/LibUtilities
    COMPONENT dev
    FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
