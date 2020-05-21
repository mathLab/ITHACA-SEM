#ifndef AVXUTIL_HPP
#define AVXUTIL_HPP

#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>

#include <string>

namespace Nektar {
namespace AVX {


template<class T, int VW>
void CopyToAlignedVector(Array<OneD, T> &in, AlignedVector<VecData<T, VW>> &out)
{
    size_t nScal = in.num_elements();
    size_t nVec = out.size();

    // check padding and size
    unsigned int pad = nScal % VW;
    ASSERTL1(nVec == nScal + (pad == 0 ? 0 : 1), "incorrect size");

    // copy
    T* inPtr = &in[0];
    for (size_t i = 0; i < size / VW; ++i, inPtr += VW)
    {
        out = inPtr.data();
    }

    // fill padded chuck everything else out
    if (pad > 0)
    {
        std::array<T, VW> tmpArray{};
        for (int i = 0; i < pad; ++i)
        {
            tmpArray[i] = *inPtr;
            ++inPtr;
        }
        out[size] = tmpArray;
    }
}

template<class T, int VW>
void CopyFromAlignedVector(AlignedVector<VecData<T, VW>> &in, Array<OneD, T> &out)
{
    size_t nScal = in.num_elements();
    size_t nVec = out.size();

    // check padding and size
    unsigned int pad = nScal % VW;
    ASSERTL1(nVec == nScal + (pad == 0 ? 0 : 1), "incorrect size");

    // shifted pointers
    T *outPtr0 = out.data();
    T *outPtr1 = out.data() + 1;
    T *outPtr2 = out.data() + 2;
    T *outPtr3 = out.data() + 3;

    size_t i = 0;
    for (; i < nScal/VW; ++i)
    {
        outPtr0[i] = in[i].m_data[0];
        outPtr1[i] = in[i].m_data[1];
        outPtr2[i] = in[i].m_data[2];
        outPtr3[i] = in[i].m_data[3];
    }

    for (size_t j = 0; j < pad; ++j)
    {
        outPtr0[i+j] = in[i].m_data[j];
    }

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
