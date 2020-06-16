#include "AVXHelmholtz.h"

#if defined(__AVX2__)
std::string __register_Helmholtz_Quad_AVX = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Quad_Regular_AVX"), &AVXHelmholtzQuad<4>::Create);

std::string __register_Helmholtz_Quad_Deformed_AVX = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Quad_Deformed_AVX"), &AVXHelmholtzQuad<4, true>::Create);

std::string __register_Helmholtz_Tri_AVX = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Tri_Regular_AVX"), &AVXHelmholtzTri<4>::Create);

std::string __register_Helmholtz_Tri_Deformed_AVX = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Tri_Deformed_AVX"), &AVXHelmholtzTri<4, true>::Create);

std::string __register_Helmholtz_Prism_AVX = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Prism_Regular_AVX"), &AVXHelmholtzPrism<4>::Create);

std::string __register_Helmholtz_Prism_Deformed_AVX = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Prism_Deformed_AVX"), &AVXHelmholtzPrism<4, true>::Create);

std::string __register_Helmholtz_Tet_AVX = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Tet_Regular_AVX"), &AVXHelmholtzTet<4>::Create);

std::string __register_Helmholtz_Tet_Deformed_AVX = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Tet_Deformed_AVX"), &AVXHelmholtzTet<4, true>::Create);

std::string __register_Helmholtz_Hex_AVX = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Hex_Regular_AVX"), &AVXHelmholtzHex<4>::Create);

std::string __register_Helmholtz_Hex_Deformed_AVX = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Hex_Deformed_AVX"), &AVXHelmholtzHex<4, true>::Create);
#endif

#if defined(__AVX512F__)
std::string __register_Helmholtz_Quad_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Quad_Regular_AVX512"), &AVXHelmholtzQuad<8>::Create);

std::string __register_Helmholtz_Quad_Deformed_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Quad_Deformed_AVX512"), &AVXHelmholtzQuad<8, true>::Create);

std::string __register_Helmholtz_Tri_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Tri_Regular_AVX512"), &AVXHelmholtzTri<8>::Create);

std::string __register_Helmholtz_Tri_Deformed_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Tri_Deformed_AVX512"), &AVXHelmholtzTri<8, true>::Create);

std::string __register_Helmholtz_Prism_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Prism_Regular_AVX512"), &AVXHelmholtzPrism<8>::Create);

std::string __register_Helmholtz_Prism_Deformed_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Prism_Deformed_AVX512"), &AVXHelmholtzPrism<8, true>::Create);

std::string __register_Helmholtz_Tet_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Tet_Regular_AVX512"), &AVXHelmholtzTet<8>::Create);

std::string __register_Helmholtz_Tet_Deformed_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Tet_Deformed_AVX512"), &AVXHelmholtzTet<8, true>::Create);

std::string __register_Helmholtz_Hex_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Hex_Regular_AVX512"), &AVXHelmholtzHex<8>::Create);

std::string __register_Helmholtz_Hex_Deformed_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Hex_Deformed_AVX512"), &AVXHelmholtzHex<8, true>::Create);
#endif

