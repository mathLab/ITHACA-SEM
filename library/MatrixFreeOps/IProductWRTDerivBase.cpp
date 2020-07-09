#include "IProductWRTDerivBase.h"

namespace Nektar
{
namespace MatrixFree
{

std::string __register_IProductWRTDerivBase_Quad = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProductWRTDerivBase_Quad_Regular"), &IProductWRTDerivBaseQuad<>::Create);

std::string __register_IProductWRTDerivBase_Quad_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProductWRTDerivBase_Quad_Deformed"), &IProductWRTDerivBaseQuad<true>::Create);

std::string __register_IProductWRTDerivBase_Tri = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProductWRTDerivBase_Tri_Regular"), &IProductWRTDerivBaseTri<>::Create);

std::string __register_IProductWRTDerivBase_Tri_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProductWRTDerivBase_Tri_Deformed"), &IProductWRTDerivBaseTri<true>::Create);

std::string __register_IProductWRTDerivBase_Tet = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProductWRTDerivBase_Tet_Regular"), &IProductWRTDerivBaseTet<>::Create);

std::string __register_IProductWRTDerivBase_Tet_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProductWRTDerivBase_Tet_Deformed"), &IProductWRTDerivBaseTet<true>::Create);

std::string __register_IProductWRTDerivBase_Prism = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProductWRTDerivBase_Prism_Regular"), &IProductWRTDerivBasePrism<>::Create);

std::string __register_IProductWRTDerivBase_Prism_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProductWRTDerivBase_Prism_Deformed"), &IProductWRTDerivBasePrism<true>::Create);

std::string __register_IProductWRTDerivBase_Hex = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProductWRTDerivBase_Hex_Regular"), &IProductWRTDerivBaseHex<>::Create);

std::string __register_IProductWRTDerivBase_Hex_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProductWRTDerivBase_Hex_Deformed"), &IProductWRTDerivBaseHex<true>::Create);

} // namespace MatrixFree
} // namespace Nektar