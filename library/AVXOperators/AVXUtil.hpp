#ifndef AVXUTIL_HPP
#define AVXUTIL_HPP

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ShapeType.hpp>
#include <LibUtilities/Foundations/Basis.h>

#include <string>

namespace Nektar {
namespace AVX {

#if defined(__AVX2__)
/// Number of bits in a AVX2 vector
constexpr int SIMD_WIDTH_BITS = 256;
#endif
#if defined(__AVX512F__)
/// Number of bits in a AVX512 vector
constexpr int SIMD_WIDTH_BITS = 512;
#endif

/// Number of bytes in a AVX2 vector
constexpr int SIMD_WIDTH_BYTES = SIMD_WIDTH_BITS / 8;
/// Number of elements in a AVX vector
constexpr int SIMD_WIDTH_SIZE = SIMD_WIDTH_BYTES / sizeof(NekDouble);

/// Get operator string
std::string GetOpstring(LibUtilities::ShapeType shape, bool deformed=false)
{

    std::string op_string = "_";

    if(shape == LibUtilities::eTriangle){
        op_string += "Tri";
    }
    else if(shape == LibUtilities::eQuadrilateral){
        op_string += "Quad";
    }
    else if(shape == LibUtilities::eTetrahedron){
        op_string += "Tet";
    }
    else if(shape == LibUtilities::ePrism){
        op_string += "Prism";
    }
    else if(shape == LibUtilities::eHexahedron){
        op_string += "Hex";
    }

    if (deformed)
    {
        op_string += "_Deformed";
    }
    else
    {
        op_string += "_Regular";
    }

    if (SIMD_WIDTH_SIZE == 4)
    {
        op_string += "_AVX";
    }
    else
    {
        op_string += "_AVX512";
    }

    return op_string;
}


template<typename T>
void TransposeData(const int       nElmt,
                   const int       vecWidth,
                   Array<OneD, T> &inout)
{
    // Make a copy of the data.
    Array<OneD, T> tmp(inout.num_elements(), &inout[0]);

    // Loop over 4 elements at a time and interleave data
    auto *inptr = &tmp[0];
    auto *outptr = &inout[0];

    size_t totSize = inout.num_elements();
    ASSERTL0(totSize % nElmt == 0, "Should be uniform");

    // Size of each element's contributions, e.g. number of coefficients per
    // element
    size_t elmtDataSize = totSize / nElmt;

    // Number of blocks
    size_t numBlocks = nElmt / vecWidth;
    size_t pad = nElmt % vecWidth;

    for (size_t i = 0; i < numBlocks; ++i)
    {
        for (size_t j = 0; j < elmtDataSize; ++j)
        {
            for (size_t k = 0; k < vecWidth; ++k)
            {
                outptr[j*vecWidth + k] = inptr[k*elmtDataSize + j];
            }
        }

        inptr += elmtDataSize * vecWidth;
        outptr += elmtDataSize * vecWidth;
    }
}

template<typename T>
void InvTransposeData(const int       nElmt,
                      const int       vecWidth,
                      Array<OneD, T> &inout)
{
    // Make a copy of the data.
    Array<OneD, T> tmp(inout.num_elements(), &inout[0]);

    // Loop over 4 elements at a time and interleave data
    auto *inptr = &tmp[0];
    auto *outptr = &inout[0];

    size_t totSize = inout.num_elements();
    ASSERTL0(totSize % nElmt == 0, "Should be uniform");

    // Size of each element's contributions, e.g. number of coefficients per
    // element
    size_t elmtDataSize = totSize / nElmt;

    // Number of blocks
    size_t numBlocks = nElmt / vecWidth;
    size_t pad = nElmt % vecWidth;

    for (size_t i = 0; i < numBlocks; ++i)
    {
        for (size_t j = 0; j < elmtDataSize; ++j)
        {
            for (size_t k = 0; k < vecWidth; ++k)
            {
                outptr[k*elmtDataSize + j] = inptr[j*vecWidth + k];
            }
        }

        inptr  += elmtDataSize * vecWidth;
        outptr += elmtDataSize * vecWidth;
    }
}

} // namespace AVX
} // namespace Nektar

#endif // header guard
