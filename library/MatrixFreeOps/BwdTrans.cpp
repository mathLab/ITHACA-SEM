#include "BwdTrans.h"

namespace Nektar
{
namespace MatrixFree
{

std::string __register_BwdTrans_Quad = GetOperatorFactory().
    RegisterCreatorFunction( std::string("BwdTrans_Quad_Regular"),
    &BwdTransQuad::Create);

std::string __register_BwdTrans_Tri = GetOperatorFactory().
    RegisterCreatorFunction( std::string("BwdTrans_Tri_Regular"),
    &BwdTransTri::Create);

std::string __register_BwdTrans_Hex = GetOperatorFactory().
    RegisterCreatorFunction( std::string("BwdTrans_Hex_Regular"),
    &BwdTransHex::Create);

std::string __register_BwdTrans_Tet = GetOperatorFactory().
    RegisterCreatorFunction( std::string("BwdTrans_Tet_Regular"),
    &BwdTransTet::Create);

std::string __register_BwdTrans_Prism = GetOperatorFactory().
    RegisterCreatorFunction( std::string("BwdTrans_Prism_Regular"),
    &BwdTransPrism::Create);

} // namespace MatrixFree
} // namespace Nektar