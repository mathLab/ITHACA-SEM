#ifndef NEKTAR_LIBRARY_MF_BWDTRANSKERNELS_HPP
#define NEKTAR_LIBRARY_MF_BWDTRANSKERNELS_HPP

#include <LibUtilities/SimdLib/tinysimd.hpp>

namespace Nektar
{
namespace MatrixFree
{

using namespace tinysimd;
using vec_t = simd<NekDouble>;

template<unsigned short NUMMODE0, unsigned short NUMMODE1,
         unsigned short NUMQUAD0, unsigned short NUMQUAD1>
inline static void BwdTransQuadKernel(
    const std::vector<vec_t, allocator<vec_t>> &in,
    const std::vector<vec_t, allocator<vec_t>> &bdata0,
    const std::vector<vec_t, allocator<vec_t>> &bdata1,
    vec_t *wsp,
    std::vector<vec_t, allocator<vec_t>> &out)
{

    constexpr int nm0 = NUMMODE0, nm1 = NUMMODE1;
    constexpr int nq0 = NUMQUAD0, nq1 = NUMQUAD1;

    for (int i = 0, cnt_iq = 0; i < nq0; ++i)
    {
        for (int q = 0, cnt_pq = 0; q < nm1; ++q, ++cnt_iq)
        {
            vec_t tmp =  in[cnt_pq] * bdata0[i]; //Load 2x
            ++cnt_pq;
            for (int p = 1; p < nm0; ++p, ++cnt_pq)
            {
                tmp.fma(in[cnt_pq], bdata0[p * nq0 + i]); //Load 2x
            }
            wsp[cnt_iq] = tmp; //Store 1x
        }
    }

    for (int j = 0, cnt_ij = 0; j < nq1; ++j)
    {
        for (int i = 0, cnt_iq = 0; i < nq0; ++i, ++cnt_ij)
        {
            vec_t tmp = wsp[cnt_iq] * bdata1[j]; //Load 2x
            ++cnt_iq;
            for (int q = 1; q < nm1; ++q, ++cnt_iq)
            {
                tmp.fma(wsp[cnt_iq], bdata1[q * nq1 + j]); //Load 2x
            }
            out[cnt_ij] = tmp; //Store 1x
        }
    }

}

// template<int NUMMODE0, int NUMMODE1,
//          int NUMQUAD0, int NUMQUAD1,
//          int VW, bool CORRECT>
// inline static void AVXBwdTransTriKernel(
//     const AlignedVector<VecData<NekDouble, VW>> &in,
//     const AlignedVector<VecData<NekDouble, VW>> &bdata0,
//     const AlignedVector<VecData<NekDouble, VW>> &bdata1,
//     VecData<NekDouble, VW> *p_sums,
//     AlignedVector<VecData<NekDouble, VW>> &out)
// {
//     using T = VecData<NekDouble, VW>;

//     constexpr int nm0 = NUMMODE0, nm1 = NUMMODE1;
//     constexpr int nq0 = NUMQUAD0, nq1 = NUMQUAD1;

//     for (int eta1 = 0, eta_idx = 0; eta1 < nq1; ++eta1)
//     {
//         for (int p = 0, mode = 0; p < nm0; ++p)
//         {

//             T p_sum = T(0.0);

//             for (int q = 0; q < (nm1-p); ++q, ++mode)
//             {
//                p_sum.fma(bdata1[mode * nq1 + eta1], in[mode]);
//             }

//             p_sums[p] = p_sum; //Store 1x

//         }

//         //We already have q_sum at each quadrature point in eta 1 for each mode, p.
//         //From this assemble the tensor produce of each quadrature point, eta1
//         for (int eta0 = 0; eta0 < nq0; ++eta0, ++eta_idx)
//         {
//             T p_sum = T(0.0);
//             for (int p = 0; p < nm0; ++p)
//             {
//                 p_sum.fma(p_sums[p], bdata0[p*nq0 + eta0]); //Load 2x
//             }

//             if (CORRECT)
//             {
//                 //p_sum += coef * bdata0 * bdata1
//                 p_sum.fma(in[1] * bdata0[nq0 + eta0], bdata1[nq1 + eta1]);
//             }

//             out[eta_idx] = p_sum;
//         }
//     }
// }

// template<int NUMMODE0, int NUMMODE1, int NUMMODE2,
//          int NUMQUAD0, int NUMQUAD1, int NUMQUAD2,
//          int VW, bool CORRECT, class BasisType>
// inline static void AVXBwdTransPrismKernel(
//     const NekDouble *inptr,
//     const AlignedVector<BasisType> &bdata0,
//     const AlignedVector<BasisType> &bdata1,
//     const AlignedVector<BasisType> &bdata2,
//     VecData<NekDouble, VW> *fpq,
//     VecData<NekDouble, VW> *fp,
//     NekDouble *outptr)
// {
//     using T = VecData<NekDouble, VW>;

//     constexpr int nm0 = NUMMODE0, nm1 = NUMMODE1, nm2 = NUMMODE2;
//     constexpr int nq0 = NUMQUAD0, nq1 = NUMQUAD1, nq2 = NUMQUAD2;

//     int cnt_kji = 0;
//     for(int k = 0; k < nq2; k++)
//     {

//         int mode_pqr = 0, mode_pq = 0, mode_pr = 0;
//         for(int p = 0; p < nm0; p++){
//             for(int q = 0; q < nm1; q++, mode_pq++){
//                 T prod = T(0.0);
//                 for(int r = 0; r < nm2-p; r++, mode_pqr++){
//                     T coef = T(inptr + mode_pqr*VW); //Load 1x
//                     prod.fma(coef, bdata2[(mode_pr + r)*nq2 + k]); //Load 1x
//                 }

//                 fpq[mode_pq] = prod; //Store 1x
//             }

//             mode_pr += nm2 - p;
//         }

//         for(int j = 0; j < nq1; j++){

//             mode_pq = 0;
//             for(int p = 0; p < nm0; p++){
//                 T prod = T(0.0);
//                 for(int q = 0; q < nm1; q++, mode_pq++){
//                     prod.fma(fpq[mode_pq], bdata1[q*nq1 + j]); //Load 2x
//                 }
//                 fp[p] = prod; //Store 1x

//             }

//             for(int i = 0; i < nq0; i++, cnt_kji++){

//                 T val_kji = T(0.0);
//                 for(int p = 0; p < nm0; p++){
//                     val_kji.fma(fp[p], bdata0[p*nq0 + i]); //Load 2x
//                 }

//                 if (CORRECT)
//                 {
//                     T basis_2 = bdata2[nq2 + k]; //Load 1x
//                     T basis_0 = bdata0[nq0 + i]; //Load 1x

//                     for(int q = 0; q < nm1; q++){
//                         T coef_0q1 = T(inptr + (q*nm2 + 1)*VW ); //Load 1x
//                         T basis_1 = bdata1[q*nq1 + j]; //Load 1x
//                         val_kji.fma(basis_2*basis_1, basis_0*coef_0q1);
//                     }
//                 }

//                 val_kji.store(outptr + cnt_kji*VW); //store 1x
//             }

//         }
//     }

// }


template<int NUMMODE0, int NUMMODE1, int NUMMODE2,
         int NUMQUAD0, int NUMQUAD1, int NUMQUAD2>
inline static void BwdTransHexKernel(
    const std::vector<vec_t, allocator<vec_t>> &in,
    const std::vector<vec_t, allocator<vec_t>> &bdata0,
    const std::vector<vec_t, allocator<vec_t>> &bdata1,
    const std::vector<vec_t, allocator<vec_t>> &bdata2,
    vec_t* sum_irq,
    vec_t* sum_jir,
    std::vector<vec_t, allocator<vec_t>> &out)
{
    constexpr auto nm0 = NUMMODE0, nm1 = NUMMODE1, nm2 = NUMMODE2;
    constexpr auto nq0 = NUMQUAD0, nq1 = NUMQUAD1, nq2 = NUMQUAD1;

    for (int i = 0, cnt_irq = 0; i < nq0; ++i)
    {
        for (int r = 0, cnt_rqp = 0; r < nm2; ++r)
        {
            for (int q = 0; q < nm1; ++q, ++cnt_irq)
            {
                vec_t tmp = in[cnt_rqp] * bdata0[i];
                ++cnt_rqp;

                for (int p = 1; p < nm0; ++p, ++cnt_rqp)
                {
                    tmp.fma(in[cnt_rqp], bdata0[p*nq0+i]);
                }

                sum_irq[cnt_irq] = tmp;
            }
        }
    }

    for (int j = 0, cnt_jir = 0; j < nq1; ++j)
    {
        for (int i = 0, cnt_irq = 0; i < nq0; ++i)
        {
            for (int r = 0; r < nm2; ++r, ++cnt_jir)
            {
                vec_t tmp = sum_irq[cnt_irq]* bdata1[j];
                ++cnt_irq;

                for (int q = 1; q < nm1; ++q)
                {
                    tmp.fma(sum_irq[cnt_irq++], bdata1[q*nq1+j]);
                }

                sum_jir[cnt_jir] = tmp;
            }
        }
    }

    for (int k = 0, cnt_kji = 0; k < nq2; ++k)
    {
        for (int j = 0, cnt_jir = 0; j < nq1; ++j)
        {
            for (int i = 0; i < nq0; ++i, ++cnt_kji)
            {
                vec_t tmp = sum_jir[cnt_jir] * bdata2[k];
                ++cnt_jir;

                for (int r = 1; r < nm2; ++r)
                {
                    tmp.fma(sum_jir[cnt_jir++], bdata2[r*nq2+k]);
                }

                out[cnt_kji] = tmp;
            }
        }
    }

}

// template<int NUMMODE0, int NUMMODE1, int NUMMODE2,
//          int NUMQUAD0, int NUMQUAD1, int NUMQUAD2,
//          int VW, bool CORRECT, class BasisType>
// inline static void AVXBwdTransTetKernel(
//     const NekDouble *inptr,
//     const AlignedVector<BasisType> &bdata0,
//     const AlignedVector<BasisType> &bdata1,
//     const AlignedVector<BasisType> &bdata2,
//     VecData<NekDouble, VW> *fpq,
//     VecData<NekDouble, VW> *fp,
//     NekDouble *outptr)
// {
//     using T = VecData<NekDouble, VW>;

//     constexpr int nm0 = NUMMODE0, nm1 = NUMMODE1, nm2 = NUMMODE2;
//     constexpr int nq0 = NUMQUAD0, nq1 = NUMQUAD1, nq2 = NUMQUAD1;

//     int cnt_kji = 0;
//     for (int k = 0; k < nq2; ++k)
//     {
//         int cnt_pqr = 0, cnt_pq = 0, mode = 0;
//         for (int p = 0; p < nm0; ++p)
//         {
//             for (int q = 0; q < nm1-p; ++q, ++cnt_pq)
//             {
//                 T prod = T(inptr + cnt_pqr) * bdata2[k + nq2*mode]; //Load 2x
//                 ++mode;
//                 cnt_pqr += VW;

//                 for (int r = 1; r < nm2-p-q; ++r, ++mode, cnt_pqr += VW)
//                 {
//                     T inxmm = T(inptr + cnt_pqr); //Load 1x
//                     prod.fma(inxmm, bdata2[k + nq2*mode]); //Load 1x
//                 }

//                 fpq[cnt_pq] = prod; //Store 1x
//             }

//             //increment mode in case order1!=order2
//             for(int q = nm1-p; q < nm2-p; ++q)
//             {
//                 mode += nm2-p-q;
//             }
//         }

//         for (int j = 0; j < nq1; ++j)
//         {
//             mode = cnt_pq = 0;
//             for (int p = 0; p < nm0; ++p)
//             {
//                 T prod = fpq[cnt_pq] * bdata1[mode*nq1+j]; //Load 2x
//                 ++cnt_pq;

//                 for (int q = 1; q < nm1 - p; ++q, ++cnt_pq)
//                 {
//                     prod.fma(fpq[cnt_pq], bdata1[(mode+q)*nq1+j]); //Load 2x
//                 }

//                 fp[p] = prod; //Store 1x
//                 mode += nm1 - p;
//             }

//             for (int i = 0; i < nq0; ++i, cnt_kji += VW)
//             {
//                 T tmp = bdata0[i] * fp[0]; //Load 2x

//                 for (int p = 1; p < nm0; ++p)
//                 {
//                     tmp.fma(bdata0[p*nq0+i], fp[p]); //Load 2x
//                 }

//                 if (CORRECT)
//                 {
//                     // top vertex
//                     //
//                     // sum += inarray[1] * base2[nquad2 + k] * (
//                     //     base0[i] * base1[nquad1+j] +
//                     //     base0[nquad0+i] * base1[j] +
//                     //     base0[nquad0+i] * base1[nquad1+j]);

//                     T tmp1 = bdata0[i] * bdata1[nq1+j]; //Load 2x
//                     tmp1.fma(bdata0[nq0+i], bdata1[j]); //Load 2x
//                     tmp1.fma(bdata0[nq0+i], bdata1[nq1+j]); //Load 2x
//                     tmp1 = tmp1 * bdata2[nq2+k]; //Load 1x

//                     T inarray1 = T(inptr + VW); //Load 1x
//                     tmp.fma(tmp1, inarray1);

//                     // bottom vertex
//                     //
//                     // sum += inarray[order2] * base2[k] * (
//                     //     base0[nquad0+i] * base1[nquad1+j]);
//                     tmp1 = bdata0[nq0+i] * bdata1[nq1+j]; //Load 2x
//                     tmp1 = tmp1 * bdata2[k]; //Load 1x
//                     inarray1 = T(inptr + VW*nm2); //Load 1x
//                     tmp.fma(inarray1, tmp1);

//                     // singular edge
//                     for (int r = 1; r < nm2-1; ++r)
//                     {
//                         // sum += inarray[order2+r] * base2[(r+1)*nquad2+k] *
//                         //     base1[nquad1+j] * base0[nquad0+i];
//                         tmp1 = bdata1[nq1+j] * bdata0[nq0+i]; //Load 2x
//                         tmp1 = tmp1 * bdata2[(r+1)*nq2+k]; //Load 1x
//                         inarray1 = T(inptr + (nm2+r)*VW); //Load 1x
//                         tmp.fma(inarray1, tmp1);
//                         // multiply by (1-a)/2
//                     }
//                 }

//                 tmp.store(outptr + cnt_kji); //Store 1x
//             }
//         }
//     }
//  }

} // namespace MatrixFree
} // namespace Nektar

#endif
