#include "PhysDeriv.h"

namespace Nektar
{
namespace MatrixFree
{

std::string __register_PhysDeriv_Quad = GetOperatorFactory().RegisterCreatorFunction(
    std::string("PhysDeriv_Quad_Regular"), &PhysDerivQuad<>::Create);

std::string __register_PhysDeriv_Quad_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("PhysDeriv_Quad_Deformed"), &PhysDerivQuad<true>::Create);


std::string __register_PhysDeriv_Tri = GetOperatorFactory().RegisterCreatorFunction(
    std::string("PhysDeriv_Tri_Regular"), &PhysDerivTri<>::Create);

std::string __register_PhysDeriv_Tri_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("PhysDeriv_Tri_Deformed"), &PhysDerivTri<true>::Create);

// std::string __register_PhysDeriv_Tet = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("PhysDeriv_Tet_Regular"), &PhysDerivTet<4>::Create);

// std::string __register_PhysDeriv_Tet_Deformed = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("PhysDeriv_Tet_Deformed"), &PhysDerivTet<4, true>::Create);

std::string __register_PhysDeriv_Hex = GetOperatorFactory().RegisterCreatorFunction(
    std::string("PhysDeriv_Hex_Regular"), &PhysDerivHex<>::Create);

std::string __register_PhysDeriv_Hex_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("PhysDeriv_Hex_Deformed"), &PhysDerivHex<true>::Create);

// std::string __register_PhysDeriv_Prism = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("PhysDeriv_Prism_Regular"), &PhysDerivPrism<4>::Create);

// std::string __register_PhysDeriv_Prism_Deformed = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("PhysDeriv_Prism_Deformed"), &PhysDerivPrism<4, true>::Create);

} // namespace MatrixFree
} // namespace Nektar