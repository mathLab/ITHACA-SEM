#ifndef AVXUTIL_HPP
#define AVXUTIL_HPP

#include <SimdOperators/VecData.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include <string>

namespace Nektar {
namespace AVX {

// Copy from nek array to aligned vector with implicit padding
// generic implementation
template<class T, int VW>
void CopyToAlignedVector(const Array<OneD, T> &in, AlignedVector<VecData<T, VW>> &out)
{
    size_t nScal = in.num_elements();
    size_t nVec = out.size();

    // check padding and size
    unsigned short pad = nScal % VW;
    ASSERTL1(nVec == nScal / VW + (pad == 0 ? 0 : 1), "incorrect size");

    // copy
    const T* inPtr = in.data();
    for (size_t i = 0; i < nScal / VW; ++i, inPtr += VW)
    {
        for (size_t j = 0; j < VW; ++j)
        {
            out[i].m_data[j] = inPtr[j];
        }
    }

    // fill padded chunck
    if (pad > 0)
    {
        alignas(VW*sizeof(T)) std::array<T, VW> tmpArray{};
        for (int i = 0; i < pad; ++i)
        {
            tmpArray[i] = inPtr[i];
        }
        out[nVec-1] = tmpArray.data();
    }
}

// Copy from aligned vector to nek array with implicit padding
// generic implementation
template<class T, int VW>
void CopyFromAlignedVector(const AlignedVector<VecData<T, VW>> &in, Array<OneD, T> &out)
{
    size_t nScal = out.num_elements();

    // check padding and size
    unsigned short pad = nScal % VW;
    ASSERTL1(in.size() == nScal / VW + (pad == 0 ? 0 : 1), "incorrect size");

    // copy
    T *outPtr = out.data();
    size_t i = 0;
    for (; i < nScal / VW; ++i)
    {
        for (size_t j = 0; j < VW; ++j)
        {
            outPtr[i*VW+j] = in[i].m_data[j];
        }
    }

    // single step in case of padding
    for (size_t j = 0; j < pad; ++j)
    {
        outPtr[i*VW+j] = in[i].m_data[j];
    }

}

// // Copy from aligned vector to nek array with implicit padding
// // avx2 specialization with loop unroll
// template<>
// void CopyFromAlignedVector<double, 4>(const AlignedVector<VecData<double, 4>> &in,
//     Array<OneD, double> &out)
// {
//     constexpr unsigned short VW = 4;

//     size_t nScal = out.num_elements();
//     size_t nVec = in.size();

//     // check padding and size
//     unsigned short pad = nScal % VW;
//     ASSERTL1(nVec == nScal + (pad == 0 ? 0 : 1), "incorrect size");

//     // copy loop unrolled
//     // shifted pointers
//     double *outPtr0 = out.data();
//     double *outPtr1 = out.data() + 1;
//     double *outPtr2 = out.data() + 2;
//     double *outPtr3 = out.data() + 3;
//     size_t i = 0;
//     for (; i < nScal/VW; ++i)
//     {
//         outPtr0[i*VW] = in[i].m_data[0];
//         outPtr1[i*VW] = in[i].m_data[1];
//         outPtr2[i*VW] = in[i].m_data[2];
//         outPtr3[i*VW] = in[i].m_data[3];
//     }

//     // single step in case of padding
//     for (size_t j = 0; j < pad; ++j)
//     {
//         outPtr0[i*VW+j] = in[i].m_data[j];
//     }

// }

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
