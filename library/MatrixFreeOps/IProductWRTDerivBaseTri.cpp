#include "IProductWRTDerivBase.h"

namespace Nektar
{
namespace MatrixFree
{

std::string __register_IProductWRTDerivBase_Tri = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProductWRTDerivBase_Tri_Regular"), &IProductWRTDerivBaseTri<>::Create);

std::string __register_IProductWRTDerivBase_Tri_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProductWRTDerivBase_Tri_Deformed"), &IProductWRTDerivBaseTri<true>::Create);

} // namespace MatrixFree
} // namespace Nektar
