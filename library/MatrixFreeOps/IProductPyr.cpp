#include "IProduct.h"

namespace Nektar
{
namespace MatrixFree
{

std::string __register_IProduct_Pyr = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProduct_Pyr_Regular"), &IProductPyr<>::Create);

std::string __register_IProduct_Pyr_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProduct_Pyr_Deformed"), &IProductPyr<true>::Create);

} // namespace MatrixFree
} // namespace Nektar
