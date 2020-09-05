#include "Helmholtz.h"

namespace Nektar
{
namespace MatrixFree
{
std::string __register_Helmholtz_Quad = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Quad_Regular"), &HelmholtzQuad<>::Create);

std::string __register_Helmholtz_Quad_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Quad_Deformed"), &HelmholtzQuad<true>::Create);

std::string __register_Helmholtz_Tri = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Tri_Regular"), &HelmholtzTri<>::Create);

std::string __register_Helmholtz_Tri_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Tri_Deformed"), &HelmholtzTri<true>::Create);

std::string __register_Helmholtz_Hex = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Hex_Regular"), &HelmholtzHex<>::Create);

std::string __register_Helmholtz_Hex_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Hex_Deformed"), &HelmholtzHex<true>::Create);

std::string __register_Helmholtz_Prism = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Prism_Regular"), &HelmholtzPrism<>::Create);

std::string __register_Helmholtz_Prism_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Prism_Deformed"), &HelmholtzPrism<true>::Create);

std::string __register_Helmholtz_Pyr = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Pyr_Regular"), &HelmholtzPyr<>::Create);

std::string __register_Helmholtz_Pyr_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Pyr_Deformed"), &HelmholtzPyr<true>::Create);

    std::string __register_Helmholtz_Tet = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Tet_Regular"), &HelmholtzTet<>::Create);

std::string __register_Helmholtz_Tet_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Tet_Deformed"), &HelmholtzTet<true>::Create);
    
} // namespace MatrixFree
} // namespace Nektar
