#include "IProduct.h"

namespace Nektar
{
namespace MatrixFree
{

std::string __register_IProduct_Quad = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProduct_Quad_Regular"), &IProductQuad<>::Create);

std::string __register_IProduct_Quad_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProduct_Quad_Deformed"), &IProductQuad<true>::Create);

} // namespace MatrixFree
} // namespace Nektar
