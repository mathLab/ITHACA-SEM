#include "IProductWRTDerivBase.h"

namespace Nektar
{
namespace MatrixFree
{

std::string __register_IProductWRTDerivBase_Seg = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProductWRTDerivBase_Seg_Regular"), &IProductWRTDerivBaseSeg<>::Create);

std::string __register_IProductWRTDerivBase_Seg_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProductWRTDerivBase_Seg_Deformed"), &IProductWRTDerivBaseSeg<true>::Create);

} // namespace MatrixFree
} // namespace Nektar
