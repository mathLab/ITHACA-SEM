#include "IProductWRTDerivBase.h"

namespace Nektar
{
namespace MatrixFree
{

std::string __register_IProductWRTDerivBase_Quad = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProductWRTDerivBase_Quad_Regular"), &IProductWRTDerivBaseQuad<>::Create);

std::string __register_IProductWRTDerivBase_Quad_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("IProductWRTDerivBase_Quad_Deformed"), &IProductWRTDerivBaseQuad<true>::Create);

} // namespace MatrixFree
} // namespace Nektar
