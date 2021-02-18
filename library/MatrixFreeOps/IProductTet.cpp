#include "IProduct.h"

namespace Nektar
{
namespace MatrixFree
{

std::string __register_IProduct_Tet = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProduct_Tet_Regular"), &IProductTet<>::Create);

std::string __register_IProduct_Tet_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProduct_Tet_Deformed"), &IProductTet<true>::Create);

} // namespace MatrixFree
} // namespace Nektar
