#include "AVXPhysDeriv.h"

namespace Nektar
{
namespace AVX
{

std::string __register_PhysDeriv_Quad_AVX = GetOperatorFactory().RegisterCreatorFunction(
    std::string("PhysDeriv_Quad_Regular_AVX"), &AVXPhysDerivQuad<>::Create);

std::string __register_PhysDeriv_Quad_Deformed_AVX = GetOperatorFactory().RegisterCreatorFunction(
    std::string("PhysDeriv_Quad_Deformed_AVX"), &AVXPhysDerivQuad<true>::Create);


// std::string __register_PhysDeriv_Tri_AVX = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("PhysDeriv_Tri_Regular_AVX"), &AVXPhysDerivTri<4>::Create);

// std::string __register_PhysDeriv_Tri_Deformed_AVX = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("PhysDeriv_Tri_Deformed_AVX"), &AVXPhysDerivTri<4,true>::Create);

// std::string __register_PhysDeriv_Tet_AVX = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("PhysDeriv_Tet_Regular_AVX"), &AVXPhysDerivTet<4>::Create);

// std::string __register_PhysDeriv_Tet_Deformed_AVX = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("PhysDeriv_Tet_Deformed_AVX"), &AVXPhysDerivTet<4, true>::Create);

std::string __register_PhysDeriv_Hex_AVX = GetOperatorFactory().RegisterCreatorFunction(
    std::string("PhysDeriv_Hex_Regular_AVX"), &AVXPhysDerivHex<>::Create);

std::string __register_PhysDeriv_Hex_Deformed_AVX = GetOperatorFactory().RegisterCreatorFunction(
    std::string("PhysDeriv_Hex_Deformed_AVX"), &AVXPhysDerivHex<true>::Create);

// std::string __register_PhysDeriv_Prism_AVX = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("PhysDeriv_Prism_Regular_AVX"), &AVXPhysDerivPrism<4>::Create);

// std::string __register_PhysDeriv_Prism_Deformed_AVX = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("PhysDeriv_Prism_Deformed_AVX"), &AVXPhysDerivPrism<4, true>::Create);

} // namespace AVX
} // namespace Nektar