#include "IProduct.h"

namespace Nektar
{
namespace MatrixFree
{

std::string __register_IProduct_Hex = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProduct_Hex_Regular"), &IProductHex<>::Create);

std::string __register_IProduct_Hex_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProduct_Hex_Deformed"), &IProductHex<true>::Create);

} // namespace MatrixFree
} // namespace Nektar
