#include "IProductWRTDerivBase.h"

namespace Nektar
{
namespace MatrixFree
{

std::string __register_IProductWRTDerivBase_Prism = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProductWRTDerivBase_Prism_Regular"), &IProductWRTDerivBasePrism<>::Create);

std::string __register_IProductWRTDerivBase_Prism_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProductWRTDerivBase_Prism_Deformed"), &IProductWRTDerivBasePrism<true>::Create);

} // namespace MatrixFree
} // namespace Nektar
