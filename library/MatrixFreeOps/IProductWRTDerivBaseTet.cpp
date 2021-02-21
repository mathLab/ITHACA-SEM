#include "IProductWRTDerivBase.h"

namespace Nektar
{
namespace MatrixFree
{

std::string __register_IProductWRTDerivBase_Tet = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProductWRTDerivBase_Tet_Regular"), &IProductWRTDerivBaseTet<>::Create);

std::string __register_IProductWRTDerivBase_Tet_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProductWRTDerivBase_Tet_Deformed"), &IProductWRTDerivBaseTet<true>::Create);

} // namespace MatrixFree
} // namespace Nektar
