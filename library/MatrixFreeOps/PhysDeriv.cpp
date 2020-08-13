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

std::string __register_PhysDeriv_Tet = GetOperatorFactory().RegisterCreatorFunction(
    std::string("PhysDeriv_Tet_Regular"), &PhysDerivTet<>::Create);

std::string __register_PhysDeriv_Tet_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("PhysDeriv_Tet_Deformed"), &PhysDerivTet<true>::Create);

std::string __register_PhysDeriv_Hex = GetOperatorFactory().RegisterCreatorFunction(
    std::string("PhysDeriv_Hex_Regular"), &PhysDerivHex<>::Create);

std::string __register_PhysDeriv_Hex_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("PhysDeriv_Hex_Deformed"), &PhysDerivHex<true>::Create);

std::string __register_PhysDeriv_Prism = GetOperatorFactory().RegisterCreatorFunction(
    std::string("PhysDeriv_Prism_Regular"), &PhysDerivPrism<>::Create);

std::string __register_PhysDeriv_Prism_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("PhysDeriv_Prism_Deformed"), &PhysDerivPrism<true>::Create);

std::string __register_PhysDeriv_Pyr = GetOperatorFactory().RegisterCreatorFunction(
    std::string("PhysDeriv_Pyr_Regular"), &PhysDerivPyr<>::Create);

std::string __register_PhysDeriv_Pyr_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("PhysDeriv_Pyr_Deformed"), &PhysDerivPyr<true>::Create);

} // namespace MatrixFree
} // namespace Nektar
