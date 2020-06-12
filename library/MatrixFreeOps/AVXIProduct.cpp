#include "AVXIProduct.h"

namespace Nektar
{
namespace AVX
{

std::string __register_IProduct_Quad_AVX = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProduct_Quad_Regular_AVX"), &AVXIProductQuad<false>::Create);

std::string __register_IProduct_Quad_Deformed_AVX = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProduct_Quad_Deformed_AVX"), &AVXIProductQuad<true>::Create);

#if defined(__AVX2__)

// std::string __register_IProduct_Tri_AVX = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("IProduct_Tri_Regular_AVX"), &AVXIProductTri<4>::Create);

// std::string __register_IProduct_Tri_Deformed_AVX = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("IProduct_Tri_Deformed_AVX"), &AVXIProductTri<4,true>::Create);

// std::string __register_IProduct_Tet_AVX = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("IProduct_Tet_Regular_AVX"), &AVXIProductTet<4>::Create);

// std::string __register_IProduct_Tet_Deformed_AVX = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("IProduct_Tet_Deformed_AVX"), &AVXIProductTet<4, true>::Create);

// std::string __register_IProduct_Prism_AVX = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("IProduct_Prism_Regular_AVX"), &AVXIProductPrism<4>::Create);

// std::string __register_IProduct_Prism_Deformed_AVX = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("IProduct_Prism_Deformed_AVX"), &AVXIProductPrism<4, true>::Create);

// std::string __register_IProduct_Hex_AVX = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("IProduct_Hex_Regular_AVX"), &AVXIProductHex<4>::Create);

// std::string __register_IProduct_Hex_Deformed_AVX = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("IProduct_Hex_Deformed_AVX"), &AVXIProductHex<4,true>::Create);
#endif

#if defined(__AVX512F__)
// std::string __register_IProduct_Tri_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("IProduct_Tri_Regular_AVX512"), &AVXIProductTri<8>::Create);

// std::string __register_IProduct_Tri_Deformed_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("IProduct_Tri_Deformed_AVX512"), &AVXIProductTri<8,true>::Create);

// std::string __register_IProduct_Tet_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("IProduct_Tet_Regular_AVX512"), &AVXIProductTet<8>::Create);

// std::string __register_IProduct_Tet_Deformed_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("IProduct_Tet_Deformed_AVX512"), &AVXIProductTet<8, true>::Create);

// std::string __register_IProduct_Prism_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("IProduct_Prism_Regular_AVX512"), &AVXIProductPrism<8>::Create);

// std::string __register_IProduct_Prism_Deformed_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("IProduct_Prism_Deformed_AVX512"), &AVXIProductPrism<8, true>::Create);

// std::string __register_IProduct_Hex_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("IProduct_Hex_Regular_AVX512"), &AVXIProductHex<8>::Create);

// std::string __register_IProduct_Hex_Deformed_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("IProduct_Hex_Deformed_AVX512"), &AVXIProductHex<8,true>::Create);
#endif

} // namespace AVX
} // namespace Nektar