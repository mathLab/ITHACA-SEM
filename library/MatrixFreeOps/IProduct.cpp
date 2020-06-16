#include "IProduct.h"

namespace Nektar
{
namespace MatrixFree
{

std::string __register_IProduct_Quad = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProduct_Quad_Regular"), &IProductQuad<>::Create);

std::string __register_IProduct_Quad_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProduct_Quad_Deformed"), &IProductQuad<true>::Create);

std::string __register_IProduct_Tri = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProduct_Tri_Regular"), &IProductTri<>::Create);

std::string __register_IProduct_Tri_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProduct_Tri_Deformed"), &IProductTri<true>::Create);

// std::string __register_IProduct_Tet = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("IProduct_Tet_Regular"), &IProductTet<4>::Create);

// std::string __register_IProduct_Tet_Deformed = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("IProduct_Tet_Deformed"), &IProductTet<4, true>::Create);

// std::string __register_IProduct_Prism = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("IProduct_Prism_Regular"), &IProductPrism<4>::Create);

// std::string __register_IProduct_Prism_Deformed = GetOperatorFactory().RegisterCreatorFunction(
//     std::string("IProduct_Prism_Deformed"), &IProductPrism<4, true>::Create);

std::string __register_IProduct_Hex = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProduct_Hex_Regular"), &IProductHex<>::Create);

std::string __register_IProduct_Hex_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProduct_Hex_Deformed"), &IProductHex<true>::Create);

} // namespace MatrixFree
} // namespace Nektar