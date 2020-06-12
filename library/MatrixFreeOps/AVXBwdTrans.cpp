#include "AVXBwdTrans.h"

namespace Nektar
{
namespace AVX
{

std::string __register_BwdTrans_Quad_AVX = GetOperatorFactory().
    RegisterCreatorFunction( std::string("BwdTrans_Quad_Regular_AVX"),
    &AVXBwdTransQuad::Create);

// std::string __register_BwdTrans_Tri_AVX = GetOperatorFactory().
//     RegisterCreatorFunction( std::string("BwdTrans_Tri_Regular_AVX"),
//     &AVXBwdTransTri::Create);

std::string __register_BwdTrans_Hex_AVX = GetOperatorFactory().
    RegisterCreatorFunction( std::string("BwdTrans_Hex_Regular_AVX"),
    &AVXBwdTransHex::Create);

// std::string __register_BwdTrans_Tet_AVX = GetOperatorFactory().
//     RegisterCreatorFunction( std::string("BwdTrans_Tet_Regular_AVX"),
//     &AVXBwdTransTet::Create);

// std::string __register_BwdTrans_Prism_AVX = GetOperatorFactory().
//     RegisterCreatorFunction( std::string("BwdTrans_Prism_Regular_AVX"),
//     &AVXBwdTransPrism::Create);

} // namespace AVX
} // namespace Nektar