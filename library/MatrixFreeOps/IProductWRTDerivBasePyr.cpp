#include "IProductWRTDerivBase.h"

namespace Nektar
{
namespace MatrixFree
{

std::string __register_IProductWRTDerivBase_Pyr = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProductWRTDerivBase_Pyr_Regular"), &IProductWRTDerivBasePyr<>::Create);

std::string __register_IProductWRTDerivBase_Pyr_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProductWRTDerivBase_Pyr_Deformed"), &IProductWRTDerivBasePyr<true>::Create);

} // namespace MatrixFree
} // namespace Nektar
