#ifndef MF_UTIL_HPP
#define MF_UTIL_HPP

#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include <string>

namespace Nektar {
namespace MatrixFree {


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

} // namespace MatrixFree
} // namespace Nektar

#endif // header guard
