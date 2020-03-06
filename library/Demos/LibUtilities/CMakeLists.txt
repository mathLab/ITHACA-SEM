ADD_NEKTAR_EXECUTABLE(FoundationDemo
    COMPONENT demos DEPENDS LibUtilities SOURCES FoundationDemo.cpp)
ADD_NEKTAR_EXECUTABLE(NodalDemo
    COMPONENT demos DEPENDS LibUtilities SOURCES NodalDemo.cpp)
ADD_NEKTAR_EXECUTABLE(TimeIntegrationDemo
    COMPONENT demos DEPENDS LibUtilities SOURCES TimeIntegrationDemo.cpp)
ADD_NEKTAR_EXECUTABLE(ErrorStream
    COMPONENT demos DEPENDS LibUtilities SOURCES ErrorStream.cpp)

IF(NEKTAR_USE_MPI)
    ADD_NEKTAR_EXECUTABLE(FieldIOBenchmarker
        COMPONENT demos DEPENDS LibUtilities SOURCES FieldIOBenchmarker.cpp)
ENDIF()

ADD_NEKTAR_TEST(NodalDemo_Tri_Deriv_P8)
ADD_NEKTAR_TEST(NodalDemo_Tri_Integral_P6)
ADD_NEKTAR_TEST(NodalDemo_Tri_Interp_P7)
ADD_NEKTAR_TEST(NodalDemo_Prism_Deriv_P8)
ADD_NEKTAR_TEST(NodalDemo_Prism_Integral_P6)
ADD_NEKTAR_TEST(NodalDemo_Prism_Interp_P7)
ADD_NEKTAR_TEST(NodalDemo_Tet_Deriv_P8)
ADD_NEKTAR_TEST(NodalDemo_Tet_Integral_P6)
ADD_NEKTAR_TEST(NodalDemo_Tet_Interp_P7)
ADD_NEKTAR_TEST(TimeIntegrationDemoIMEXGear)
ADD_NEKTAR_TEST(TimeIntegrationDemoIMEXOrder1)
ADD_NEKTAR_TEST(TimeIntegrationDemoIMEXOrder2)
ADD_NEKTAR_TEST(TimeIntegrationDemoIMEXOrder3)
ADD_NEKTAR_TEST(TimeIntegrationDemoDIRKIMEXOrder2)
ADD_NEKTAR_TEST(TimeIntegrationDemoDIRKIMEXOrder3)
ADD_NEKTAR_TEST(TimeIntegrationDemoCNAB)
ADD_NEKTAR_TEST(TimeIntegrationDemoMCNAB)
ADD_NEKTAR_TEST(ErrorStream)
