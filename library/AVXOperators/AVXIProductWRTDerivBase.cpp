#include "AVXIProductWRTDerivBase.h"

namespace Nektar
{
namespace AVX
{

#if defined(__AVX2__)
std::string __register_IProductWRTDerivBase_Quad_AVX = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProductWRTDerivBase_Quad_Regular_AVX"), &AVXIProductWRTDerivBaseQuad<4>::Create);

std::string __register_IProductWRTDerivBase_Quad_Deformed_AVX = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProductWRTDerivBase_Quad_Deformed_AVX"), &AVXIProductWRTDerivBaseQuad<4,true>::Create);

// std::string __register_IProductWRTDerivBase_Tri_AVX = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("IProductWRTDerivBase_Tri_Regular_AVX"), &AVXIProductWRTDerivBaseTri<4>::Create);

// std::string __register_IProductWRTDerivBase_Tri_Deformed_AVX = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("IProductWRTDerivBase_Tri_Deformed_AVX"), &AVXIProductWRTDerivBaseTri<4,true>::Create);

// std::string __register_IProductWRTDerivBase_Tet_AVX = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("IProductWRTDerivBase_Tet_Regular_AVX"), &AVXIProductWRTDerivBaseTet<4>::Create);

// std::string __register_IProductWRTDerivBase_Tet_Deformed_AVX = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("IProductWRTDerivBase_Tet_Deformed_AVX"), &AVXIProductWRTDerivBaseTet<4, true>::Create);

// std::string __register_IProductWRTDerivBase_Prism_AVX = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("IProductWRTDerivBase_Prism_Regular_AVX"), &AVXIProductWRTDerivBasePrism<4>::Create);

// std::string __register_IProductWRTDerivBase_Prism_Deformed_AVX = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("IProductWRTDerivBase_Prism_Deformed_AVX"), &AVXIProductWRTDerivBasePrism<4, true>::Create);

// std::string __register_IProductWRTDerivBase_Hex_AVX = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("IProductWRTDerivBase_Hex_Regular_AVX"), &AVXIProductWRTDerivBaseHex<4>::Create);

// std::string __register_IProductWRTDerivBase_Hex_Deformed_AVX = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("IProductWRTDerivBase_Hex_Deformed_AVX"), &AVXIProductWRTDerivBaseHex<4,true>::Create);
#endif

#if defined(__AVX512F__)
std::string __register_IProductWRTDerivBase_Quad_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProductWRTDerivBase_Quad_Regular_AVX512"), &AVXIProductWRTDerivBaseQuad<8>::Create);

std::string __register_IProductWRTDerivBase_Quad_Deformed_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProductWRTDerivBase_Quad_Deformed_AVX512"), &AVXIProductWRTDerivBaseQuad<8,true>::Create);

// std::string __register_IProductWRTDerivBase_Tri_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("IProductWRTDerivBase_Tri_Regular_AVX512"), &AVXIProductWRTDerivBaseTri<8>::Create);

// std::string __register_IProductWRTDerivBase_Tri_Deformed_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("IProductWRTDerivBase_Tri_Deformed_AVX512"), &AVXIProductWRTDerivBaseTri<8,true>::Create);

// std::string __register_IProductWRTDerivBase_Tet_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("IProductWRTDerivBase_Tet_Regular_AVX512"), &AVXIProductWRTDerivBaseTet<8>::Create);

// std::string __register_IProductWRTDerivBase_Tet_Deformed_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("IProductWRTDerivBase_Tet_Deformed_AVX512"), &AVXIProductWRTDerivBaseTet<8, true>::Create);

// std::string __register_IProductWRTDerivBase_Prism_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("IProductWRTDerivBase_Prism_Regular_AVX512"), &AVXIProductWRTDerivBasePrism<8>::Create);

// std::string __register_IProductWRTDerivBase_Prism_Deformed_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("IProductWRTDerivBase_Prism_Deformed_AVX512"), &AVXIProductWRTDerivBasePrism<8, true>::Create);

// std::string __register_IProductWRTDerivBase_Hex_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("IProductWRTDerivBase_Hex_Regular_AVX512"), &AVXIProductWRTDerivBaseHex<8>::Create);

// std::string __register_IProductWRTDerivBase_Hex_Deformed_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("IProductWRTDerivBase_Hex_Deformed_AVX512"), &AVXIProductWRTDerivBaseHex<8,true>::Create);
#endif

} // namespace AVX
} // namespace Nektar