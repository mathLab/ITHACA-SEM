#include "AVXBwdTrans.h"

namespace Nektar {
namespace AVX {

#if defined(__AVX2__)
std::string __register_BwdTrans_Quad_AVX = GetOperatorFactory().RegisterCreatorFunction(
    std::string("BwdTrans_Quad_Regular_AVX"), &AVXBwdTransQuad<4>::Create);

std::string __register_BwdTrans_Tri_AVX = GetOperatorFactory().RegisterCreatorFunction(
    std::string("BwdTrans_Tri_Regular_AVX"), &AVXBwdTransTri<4>::Create);

std::string __register_BwdTrans_Hex_AVX = GetOperatorFactory().RegisterCreatorFunction(
    std::string("BwdTrans_Hex_Regular_AVX"), &AVXBwdTransHex<4>::Create);

std::string __register_BwdTrans_Tet_AVX = GetOperatorFactory().RegisterCreatorFunction(
    std::string("BwdTrans_Tet_Regular_AVX"), &AVXBwdTransTet<4>::Create);

std::string __register_BwdTrans_Prism_AVX = GetOperatorFactory().RegisterCreatorFunction(
    std::string("BwdTrans_Prism_Regular_AVX"), &AVXBwdTransPrism<4>::Create);
#endif

#if defined(__AVX512F__)
std::string __register_BwdTrans_Quad_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
    std::string("BwdTrans_Quad_Regular_AVX512"), &AVXBwdTransQuad<8>::Create);

std::string __register_BwdTrans_Tri_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
    std::string("BwdTrans_Tri_Regular_AVX512"), &AVXBwdTransTri<8>::Create);

std::string __register_BwdTrans_Hex_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
    std::string("BwdTrans_Hex_Regular_AVX512"), &AVXBwdTransHex<8>::Create);

std::string __register_BwdTrans_Tet_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
    std::string("BwdTrans_Tet_Regular_AVX512"), &AVXBwdTransTet<8>::Create);

std::string __register_BwdTrans_Prism_AVX512 = GetOperatorFactory().RegisterCreatorFunction(
    std::string("BwdTrans_Prism_Regular_AVX512"), &AVXBwdTransPrism<8>::Create);
#endif

} // namespace AVX
} // namespace Nektar