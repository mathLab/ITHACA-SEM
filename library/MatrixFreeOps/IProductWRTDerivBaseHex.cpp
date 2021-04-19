#include "IProductWRTDerivBase.h"

namespace Nektar
{
namespace MatrixFree
{

std::string __register_IProductWRTDerivBase_Hex = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProductWRTDerivBase_Hex_Regular"), &IProductWRTDerivBaseHex<>::Create);

std::string __register_IProductWRTDerivBase_Hex_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProductWRTDerivBase_Hex_Deformed"), &IProductWRTDerivBaseHex<true>::Create);

} // namespace MatrixFree
} // namespace Nektar
