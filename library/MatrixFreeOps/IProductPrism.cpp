#include "IProduct.h"

namespace Nektar
{
namespace MatrixFree
{

std::string __register_IProduct_Prism = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProduct_Prism_Regular"), &IProductPrism<>::Create);

std::string __register_IProduct_Prism_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProduct_Prism_Deformed"), &IProductPrism<true>::Create);

} // namespace MatrixFree
} // namespace Nektar
