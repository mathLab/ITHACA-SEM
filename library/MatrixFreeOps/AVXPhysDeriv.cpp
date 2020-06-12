#include "AVXPhysDeriv.h"

namespace Nektar
{
namespace AVX
{

std::string __register_PhysDeriv_Quad_AVX = GetOperatorFactory().RegisterCreatorFunction(
    std::string("PhysDeriv_Quad_Regular_AVX"), &AVXPhysDerivQuad<>::Create);

std::string __register_PhysDeriv_Quad_Deformed_AVX = GetOperatorFactory().RegisterCreatorFunction(
    std::string("PhysDeriv_Quad_Deformed_AVX"), &AVXPhysDerivQuad<true>::Create);


#if defined(__AVX2__)

// std::string __register_PhysDeriv_Tri_AVX = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("PhysDeriv_Tri_Regular_AVX"), &AVXPhysDerivTri<4>::Create);

// std::string __register_PhysDeriv_Tri_Deformed_AVX = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("PhysDeriv_Tri_Deformed_AVX"), &AVXPhysDerivTri<4,true>::Create);

// std::string __register_PhysDeriv_Tet_AVX = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("PhysDeriv_Tet_Regular_AVX"), &AVXPhysDerivTet<4>::Create);

// std::string __register_PhysDeriv_Tet_Deformed_AVX = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("PhysDeriv_Tet_Deformed_AVX"), &AVXPhysDerivTet<4, true>::Create);

// std::string __register_PhysDeriv_Hex_AVX = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("PhysDeriv_Hex_Regular_AVX"), &AVXPhysDerivHex<4>::Create);

// std::string __register_PhysDeriv_Hex_Deformed_AVX = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("PhysDeriv_Hex_Deformed_AVX"), &AVXPhysDerivHex<4,true>::Create);

// std::string __register_PhysDeriv_Prism_AVX = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("PhysDeriv_Prism_Regular_AVX"), &AVXPhysDerivPrism<4>::Create);

// std::string __register_PhysDeriv_Prism_Deformed_AVX = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("PhysDeriv_Prism_Deformed_AVX"), &AVXPhysDerivPrism<4, true>::Create);
#endif

#if defined(__AVX512F__)
// std::string __register_PhysDeriv_Quad_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
    // std::string("PhysDeriv_Quad_Regular_AVX512"), &AVXPhysDerivQuad<8>::Create);
//
// std::string __register_PhysDeriv_Quad_Deformed_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
    // std::string("PhysDeriv_Quad_Deformed_AVX512"), &AVXPhysDerivQuad<8, true>::Create);

// std::string __register_PhysDeriv_Tri_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("PhysDeriv_Tri_Regular_AVX512"), &AVXPhysDerivTri<8>::Create);

// std::string __register_PhysDeriv_Tri_Deformed_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("PhysDeriv_Tri_Deformed_AVX512"), &AVXPhysDerivTri<8,true>::Create);

// std::string __register_PhysDeriv_Tet_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("PhysDeriv_Tet_Regular_AVX512"), &AVXPhysDerivTet<8>::Create);

// std::string __register_PhysDeriv_Tet_Deformed_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("PhysDeriv_Tet_Deformed_AVX512"), &AVXPhysDerivTet<8, true>::Create);

// std::string __register_PhysDeriv_Hex_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
    // std::string("PhysDeriv_Hex_Regular_AVX512"), &AVXPhysDerivHex<8>::Create);
//
// std::string __register_PhysDeriv_Hex_Deformed_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
    // std::string("PhysDeriv_Hex_Deformed_AVX512"), &AVXPhysDerivHex<8,true>::Create);

// std::string __register_PhysDeriv_Prism_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("PhysDeriv_Prism_Regular_AVX512"), &AVXPhysDerivPrism<8>::Create);

// std::string __register_PhysDeriv_Prism_Deformed_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("PhysDeriv_Prism_Deformed_AVX512"), &AVXPhysDerivPrism<8, true>::Create);
#endif

} // namespace AVX
} // namespace Nektar