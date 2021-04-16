#include "IProduct.h"

namespace Nektar
{
namespace MatrixFree
{

std::string __register_IProduct_Tri = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProduct_Tri_Regular"), &IProductTri<>::Create);

std::string __register_IProduct_Tri_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProduct_Tri_Deformed"), &IProductTri<true>::Create);

} // namespace MatrixFree
} // namespace Nektar
