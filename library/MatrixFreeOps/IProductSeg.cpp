#include "IProduct.h"

namespace Nektar
{
namespace MatrixFree
{

std::string __register_IProduct_Seg = GetOperatorFactory().
    RegisterCreatorFunction(
    std::string("IProduct_Seg_Regular"), &IProductSeg<>::Create);

std::string __register_IProduct_Seg_Deformed = GetOperatorFactory().
    RegisterCreatorFunction(
    std::string("IProduct_Seg_Deformed"), &IProductSeg<true>::Create);

} // namespace MatrixFree
} // namespace Nektar
