ADD_NEKTAR_EXECUTABLE(PartitionAnalyse
    COMPONENT demos DEPENDS LibUtilities SOURCES PartitionAnalyse.cpp)
ADD_NEKTAR_EXECUTABLE(FoundationDemo
    COMPONENT demos DEPENDS LibUtilities SOURCES FoundationDemo.cpp)
ADD_NEKTAR_EXECUTABLE(NodalTriFeketeDemo
    COMPONENT demos DEPENDS LibUtilities SOURCES NodalTriFeketeDemo.cpp)
ADD_NEKTAR_EXECUTABLE(NodalTriElecDemo
    COMPONENT demos DEPENDS LibUtilities SOURCES NodalTriElecDemo.cpp)
ADD_NEKTAR_EXECUTABLE(NodalTriEvenlySpacedDemo
    COMPONENT demos DEPENDS LibUtilities SOURCES NodalTriEvenlySpacedDemo.cpp)
ADD_NEKTAR_EXECUTABLE(NodalTetEvenlySpacedDemo
    COMPONENT demos DEPENDS LibUtilities SOURCES NodalTetEvenlySpacedDemo.cpp)
ADD_NEKTAR_EXECUTABLE(NodalTetElecDemo
    COMPONENT demos DEPENDS LibUtilities SOURCES NodalTetElecDemo.cpp)
ADD_NEKTAR_EXECUTABLE(TimeIntegrationDemo
    COMPONENT demos DEPENDS LibUtilities SOURCES TimeIntegrationDemo.cpp)
ADD_NEKTAR_EXECUTABLE(FieldIOBenchmarker
    COMPONENT demos DEPENDS LibUtilities SOURCES FieldIOBenchmarker.cpp)
