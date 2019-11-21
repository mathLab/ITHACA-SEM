#include "AVXHelmholtzGlobal.h"

#if defined(__AVX2__)
std::string __register_HelmholtzGlobal_Quad_AVX = GetOperatorFactory().RegisterCreatorFunction(
    std::string("HelmholtzGlobal_Quad_Regular_AVX"), &AVXHelmholtzGlobalQuad<4>::Create);

std::string __register_HelmholtzGlobal_Quad_Deformed_AVX = GetOperatorFactory().RegisterCreatorFunction(
    std::string("HelmholtzGlobal_Quad_Deformed_AVX"), &AVXHelmholtzGlobalQuad<4, true>::Create);

std::string __register_HelmholtzGlobal_Tri_AVX = GetOperatorFactory().RegisterCreatorFunction(
    std::string("HelmholtzGlobal_Tri_Regular_AVX"), &AVXHelmholtzGlobalTri<4>::Create);

std::string __register_HelmholtzGlobal_Tri_Deformed_AVX = GetOperatorFactory().RegisterCreatorFunction(
    std::string("HelmholtzGlobal_Tri_Deformed_AVX"), &AVXHelmholtzGlobalTri<4, true>::Create);

std::string __register_HelmholtzGlobal_Prism_AVX = GetOperatorFactory().RegisterCreatorFunction(
    std::string("HelmholtzGlobal_Prism_Regular_AVX"), &AVXHelmholtzGlobalPrism<4>::Create);

std::string __register_HelmholtzGlobal_Prism_Deformed_AVX = GetOperatorFactory().RegisterCreatorFunction(
    std::string("HelmholtzGlobal_Prism_Deformed_AVX"), &AVXHelmholtzGlobalPrism<4, true>::Create);

std::string __register_HelmholtzGlobal_Tet_AVX = GetOperatorFactory().RegisterCreatorFunction(
    std::string("HelmholtzGlobal_Tet_Regular_AVX"), &AVXHelmholtzGlobalTet<4>::Create);

std::string __register_HelmholtzGlobal_Tet_Deformed_AVX = GetOperatorFactory().RegisterCreatorFunction(
    std::string("HelmholtzGlobal_Tet_Deformed_AVX"), &AVXHelmholtzGlobalTet<4, true>::Create);

std::string __register_HelmholtzGlobal_Hex_AVX = GetOperatorFactory().RegisterCreatorFunction(
    std::string("HelmholtzGlobal_Hex_Regular_AVX"), &AVXHelmholtzGlobalHex<4>::Create);

std::string __register_HelmholtzGlobal_Hex_Deformed_AVX = GetOperatorFactory().RegisterCreatorFunction(
    std::string("HelmholtzGlobal_Hex_Deformed_AVX"), &AVXHelmholtzGlobalHex<4, true>::Create);
#endif

#if defined(__AVX512F__)
std::string __register_HelmholtzGlobal_Quad_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
    std::string("HelmholtzGlobal_Quad_Regular_AVX512"), &AVXHelmholtzGlobalQuad<8>::Create);

std::string __register_HelmholtzGlobal_Quad_Deformed_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
    std::string("HelmholtzGlobal_Quad_Deformed_AVX512"), &AVXHelmholtzGlobalQuad<8, true>::Create);

std::string __register_HelmholtzGlobal_Tri_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
    std::string("HelmholtzGlobal_Tri_Regular_AVX512"), &AVXHelmholtzGlobalTri<8>::Create);

std::string __register_HelmholtzGlobal_Tri_Deformed_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
    std::string("HelmholtzGlobal_Tri_Deformed_AVX512"), &AVXHelmholtzGlobalTri<8, true>::Create);

std::string __register_HelmholtzGlobal_Prism_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
    std::string("HelmholtzGlobal_Prism_Regular_AVX512"), &AVXHelmholtzGlobalPrism<8>::Create);

std::string __register_HelmholtzGlobal_Prism_Deformed_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
    std::string("HelmholtzGlobal_Prism_Deformed_AVX512"), &AVXHelmholtzGlobalPrism<8, true>::Create);

std::string __register_HelmholtzGlobal_Tet_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
    std::string("HelmholtzGlobal_Tet_Regular_AVX512"), &AVXHelmholtzGlobalTet<8>::Create);

std::string __register_HelmholtzGlobal_Tet_Deformed_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
    std::string("HelmholtzGlobal_Tet_Deformed_AVX512"), &AVXHelmholtzGlobalTet<8, true>::Create);

std::string __register_HelmholtzGlobal_Hex_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
    std::string("HelmholtzGlobal_Hex_Regular_AVX512"), &AVXHelmholtzGlobalHex<8>::Create);

std::string __register_HelmholtzGlobal_Hex_Deformed_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
    std::string("HelmholtzGlobal_Hex_Deformed_AVX512"), &AVXHelmholtzGlobalHex<8, true>::Create);
#endif