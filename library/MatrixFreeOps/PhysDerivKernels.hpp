#ifndef MF_PHYSDERIVKERNELS_HPP
#define MF_PHYSDERIVKERNELS_HPP

#include <LibUtilities/SimdLib/tinysimd.hpp>

namespace Nektar
{
namespace MatrixFree
{
using namespace tinysimd;
using vec_t = simd<NekDouble>;

template<int NUMQUAD0, int NUMQUAD1>
inline static void PhysDerivTensor2DKernel(
    const std::vector<vec_t, allocator<vec_t>> &in,
    const std::vector<vec_t, allocator<vec_t>> &D0,
    const std::vector<vec_t, allocator<vec_t>> &D1,
          std::vector<vec_t, allocator<vec_t>> &outptr_d0,
          std::vector<vec_t, allocator<vec_t>> &outptr_d1)
{
    constexpr auto nq0 = NUMQUAD0, nq1 = NUMQUAD1;

    //All matricies are column major ordered since operators used to be computed via BLAS.

    //D0 * in
    for (int i = 0; i < nq0; ++i)
    { //row index of D0 matrix
        for (int j = 0; j < nq1; ++j)
        { //col index of IN matrix

            vec_t prod_sum = 0.0;
            for (int k = 0; k < nq0; ++k)
            { //Col index of D0, row index of IN
                vec_t v1 = D0[k * nq0 + i]; //Load 1x
                vec_t v2 = in[j * nq0 + k]; //Load 1x

                prod_sum.fma(v1, v2);
            }

            outptr_d0[j * nq0 + i] = prod_sum; //Store 1x
        }
    }

    //in * D1^T
    for (int i = 0; i < nq0; ++i)
    { //row index for grid
        for (int j = 0; j < nq1; ++j)
        { //Column index for D1^T (row idx for D1)

            vec_t prod_sum = 0.0;
            for(int k = 0; k < nq1; ++k)
            {
                vec_t v1 = in[k * nq0 + i]; //Load 1x
                vec_t v2 = D1[k * nq1 + j]; //Load 1x

                prod_sum.fma(v1, v2);
            }

            outptr_d1[j * nq0 + i] = prod_sum; //Store 1x

        }
    }

}

template<int NUMQUAD0, int NUMQUAD1, bool DEFORMED>
static void PhysDerivQuadKernel(
    const std::vector<vec_t, allocator<vec_t>> &in,
    const std::vector<vec_t, allocator<vec_t>> &Z0,
    const std::vector<vec_t, allocator<vec_t>> &Z1,
    const std::vector<vec_t, allocator<vec_t>> &D0,
    const std::vector<vec_t, allocator<vec_t>> &D1,
    const vec_t* df_ptr,
          std::vector<vec_t, allocator<vec_t>> &out_d0,
          std::vector<vec_t, allocator<vec_t>> &out_d1)
{
    boost::ignore_unused(Z0, Z1);

    constexpr auto nq0 = NUMQUAD0, nq1 = NUMQUAD1;
    constexpr auto ndf = 4;

    PhysDerivTensor2DKernel<NUMQUAD0, NUMQUAD1>(
        in,
        D0, D1,
        out_d0, out_d1); //Results written to out_d0, out_d1

    vec_t df0, df1, df2, df3;
    if (!DEFORMED)
    {
        df0 = df_ptr[0];
        df1 = df_ptr[1];
        df2 = df_ptr[2];
        df3 = df_ptr[3];
    }

    int cnt_ji = 0;
    for (int j = 0; j < nq0; ++j)
    {
        for (int i = 0; i < nq1; ++i, ++cnt_ji)
        {

            vec_t d0 = out_d0[cnt_ji];// Load 1x
            vec_t d1 = out_d1[cnt_ji]; //Load 1x

            if (DEFORMED)
            {
                df0 = df_ptr[cnt_ji * ndf];
                df1 = df_ptr[cnt_ji * ndf + 1];
                df2 = df_ptr[cnt_ji * ndf + 2];
                df3 = df_ptr[cnt_ji * ndf + 3];
            }

            //Multiply by derivative factors
            vec_t out0 = d0 * df0; //d0 * df00 + d1 * df10
            out0.fma(d1, df1);
            out_d0[cnt_ji] = out0; //Store 1x

            vec_t out1 = d0 * df2; //d0 * df2 + d1 * df30
            out1.fma(d1, df3);
            out_d1[cnt_ji] = out1; //Store 1x
        }
    }

}

template<int NUMQUAD0, int NUMQUAD1, bool DEFORMED>
static void PhysDerivTriKernel(
    const std::vector<vec_t, allocator<vec_t>> &in,
    const std::vector<vec_t, allocator<vec_t>> &Z0,
    const std::vector<vec_t, allocator<vec_t>> &Z1,
    const std::vector<vec_t, allocator<vec_t>> &D0,
    const std::vector<vec_t, allocator<vec_t>> &D1,
    const vec_t* df_ptr,
    std::vector<vec_t, allocator<vec_t>> &out_d0,
    std::vector<vec_t, allocator<vec_t>> &out_d1)
{

    constexpr auto nq0 = NUMQUAD0, nq1 = NUMQUAD1;
    constexpr auto ndf = 4;

    PhysDerivTensor2DKernel<NUMQUAD0, NUMQUAD1>(
        in,
        D0, D1,
        out_d0, out_d1); //Results written to out_d0, out_d1

    vec_t df0, df1, df2, df3;
    if (!DEFORMED) //Optimized out by compiler
    {
        df0 = df_ptr[0];
        df1 = df_ptr[1];
        df2 = df_ptr[2];
        df3 = df_ptr[3];
    }

    int cnt_ji = 0;
    for (int j = 0; j < nq1; ++j)
    {

        vec_t xfrm0 = 2.0 / (1.0 - Z1[j]); //Load 1x

        for (int i = 0; i < nq0; ++i, ++cnt_ji)
        {

            if (DEFORMED)
            {
                df0 = df_ptr[cnt_ji * ndf];
                df1 = df_ptr[cnt_ji * ndf + 1];
                df2 = df_ptr[cnt_ji * ndf + 2];
                df3 = df_ptr[cnt_ji * ndf + 3];
            }

            //moving from standard to collapsed coordinates
            vec_t d0 = xfrm0 * out_d0[cnt_ji]; //Load 1x
            vec_t d1 = out_d1[cnt_ji]; //Load 1x
            vec_t xfrm1 = 0.5 * (1.0 + Z0[i]); //Load 1x
            d1.fma(d0, xfrm1);

            //Multiply by derivative factors
            vec_t out0 = d0 * df0; //d0 * df0 + d1 * df10
            out0.fma(d1, df1);
            out_d0[cnt_ji] = out0; // store 1x

            vec_t out1 = d0 * df2; //d0 * df2 + d1 * df3
            out1.fma(d1, df3);
            out_d1[cnt_ji] = out1; //Store 1x

        }
    }

}

template<int NUMQUAD0, int NUMQUAD1, int NUMQUAD2>
inline static void PhysDerivTensor3DKernel(
    const std::vector<vec_t, allocator<vec_t>> &in,
    const std::vector<vec_t, allocator<vec_t>> &D0,
    const std::vector<vec_t, allocator<vec_t>> &D1,
    const std::vector<vec_t, allocator<vec_t>> &D2,
          std::vector<vec_t, allocator<vec_t>> &out_d0,
          std::vector<vec_t, allocator<vec_t>> &out_d1,
          std::vector<vec_t, allocator<vec_t>> &out_d2)
{
    constexpr auto nq0 = NUMQUAD0, nq1 = NUMQUAD1, nq2 = NUMQUAD2;

    //Direction 1
    for (int i = 0; i < nq0; ++i)
    {
        for (int j = 0; j < nq1*nq2; ++j)
        {
            vec_t prod_sum = 0.0;
            for (int k = 0; k < nq0; ++k)
            {
                vec_t v1 = D0[k*nq0+i];
                vec_t v2 = in[j*nq0+k];

                prod_sum.fma(v1, v2);
            }

            out_d0[j*nq0+i] = prod_sum;

        }
    }

    // Direction 2
    for (int block = 0; block < nq2; ++block)
    {
        int start = block * nq0 * nq1;

        for (int i = 0; i < nq0; ++i)
        {
            for (int j = 0; j < nq1; ++j)
            {

                vec_t prod_sum = 0.0;
                for (int k = 0; k < nq1; ++k)
                {
                    vec_t v1 = in[start + k*nq0+i];
                    vec_t v2 = D1[k*nq1 + j];

                    prod_sum.fma(v1,v2);
                }

                out_d1[ start + j*nq0+i] = prod_sum;
            }
        }

    }

    //Direction 3
    for (int i = 0; i < nq0 * nq1; ++i)
    {
        for (int j = 0; j < nq2; ++j)
        {

            vec_t prod_sum = 0.0;
            for (int k = 0; k < nq2; ++k)
            {
                vec_t v1 = in[k*nq0*nq1 + i];
                vec_t v2 = D2[k*nq2 + j];

                prod_sum.fma(v1,v2);
            }

            out_d2[j*nq0*nq1 + i] = prod_sum;
        }
    }

}

template<int NUMQUAD0, int NUMQUAD1, int NUMQUAD2, bool DEFORMED>
static void PhysDerivHexKernel(
    const std::vector<vec_t, allocator<vec_t>> &inptr,
    const std::vector<vec_t, allocator<vec_t>> &Z0,
    const std::vector<vec_t, allocator<vec_t>> &Z1,
    const std::vector<vec_t, allocator<vec_t>> &Z2,
    const std::vector<vec_t, allocator<vec_t>> &D0,
    const std::vector<vec_t, allocator<vec_t>> &D1,
    const std::vector<vec_t, allocator<vec_t>> &D2,
    const vec_t* df_ptr,
    std::vector<vec_t, allocator<vec_t>> &outptr_d0,
    std::vector<vec_t, allocator<vec_t>> &outptr_d1,
    std::vector<vec_t, allocator<vec_t>> &outptr_d2)
{
    boost::ignore_unused(Z0, Z1, Z2);
    constexpr auto nq0 = NUMQUAD0, nq1 = NUMQUAD1, nq2 = NUMQUAD2;
    constexpr auto ndf = 9;

    PhysDerivTensor3DKernel<NUMQUAD0, NUMQUAD1, NUMQUAD2>(
        inptr,
        D0, D1, D2,
        outptr_d0, outptr_d1, outptr_d2);

    vec_t df0, df1, df2, df3, df4, df5, df6, df7, df8;
    if(!DEFORMED){
        df0 = df_ptr[0];
        df1 = df_ptr[1];
        df2 = df_ptr[2];
        df3 = df_ptr[3];
        df4 = df_ptr[4];
        df5 = df_ptr[5];
        df6 = df_ptr[6];
        df7 = df_ptr[7];
        df8 = df_ptr[8];
    }

    int cnt_ijk = 0;
    for (int k = 0; k < nq2; ++k)
    {
        for (int j = 0; j < nq1; ++j)
        {
            for (int i = 0; i < nq0; ++i, ++cnt_ijk)
            {
                vec_t d0 = outptr_d0[cnt_ijk]; //Load 1x
                vec_t d1 = outptr_d1[cnt_ijk]; //Load 1x
                vec_t d2 = outptr_d2[cnt_ijk]; //Load 1x

                if(DEFORMED){
                    df0 = df_ptr[cnt_ijk*ndf + 0];
                    df1 = df_ptr[cnt_ijk*ndf + 1];
                    df2 = df_ptr[cnt_ijk*ndf + 2];
                    df3 = df_ptr[cnt_ijk*ndf + 3];
                    df4 = df_ptr[cnt_ijk*ndf + 4];
                    df5 = df_ptr[cnt_ijk*ndf + 5];
                    df6 = df_ptr[cnt_ijk*ndf + 6];
                    df7 = df_ptr[cnt_ijk*ndf + 7];
                    df8 = df_ptr[cnt_ijk*ndf + 8];
                }

                vec_t out0 = d0 * df0;
                out0.fma(d1, df1);
                out0.fma(d2, df2);
                outptr_d0[cnt_ijk] = out0; //Store 1x

                vec_t out1 = d0 * df3;
                out1.fma(d1, df4);
                out1.fma(d2, df5);
                outptr_d1[cnt_ijk] = out1; //Store 1x

                vec_t out2 = d0 * df6;
                out2.fma(d1, df7);
                out2.fma(d2, df8);
                outptr_d2[cnt_ijk] = out2; //Store 1x

            }
        }
    }

}

// template<int NUMQUAD0, int NUMQUAD1, int NUMQUAD2,
//          int VW, bool DEFORMED, class BasisType>
// static void PhysDerivPrismKernel(
//     const NekDouble *inptr,
//     const AlignedVector<BasisType> &Z0,
//     const AlignedVector<BasisType> &Z1,
//     const AlignedVector<BasisType> &Z2,
//     const AlignedVector<BasisType> &D0,
//     const AlignedVector<BasisType> &D1,
//     const AlignedVector<BasisType> &D2,
//     const VecData<double, VW> *df_ptr,
//     NekDouble *outptr_d0,
//     NekDouble *outptr_d1,
//     NekDouble *outptr_d2)
// {
//     constexpr int nq0 = NUMQUAD0, nq1 = NUMQUAD1, nq2 = NUMQUAD2;
//     constexpr int ndf = 9;
//     using T = VecData<double, VW>;

//     PhysDerivTensor3DKernel<NUMQUAD0, NUMQUAD1, NUMQUAD2, VW>(
//         inptr,
//         D0, D1, D2,
//         outptr_d0, outptr_d1, outptr_d2);

//     T df0, df1, df2, df3, df4, df5, df6, df7, df8;
//     if(!DEFORMED){
//         df0 = df_ptr[0];
//         df1 = df_ptr[1];
//         df2 = df_ptr[2];
//         df3 = df_ptr[3];
//         df4 = df_ptr[4];
//         df5 = df_ptr[5];
//         df6 = df_ptr[6];
//         df7 = df_ptr[7];
//         df8 = df_ptr[8];
//     }

//     int cnt_ijk = 0;
//     for(int k = 0; k < nq2; k++){
//         T xfrm_eta2 = T(2.0) / (T(1.0) - Z2[k]); //Load 1x
//         for(int j = 0; j < nq1; j++){
//             for(int i = 0; i < nq0; i++, cnt_ijk++){
//                 //Chain-rule for eta_0 and eta_2
//                 T d0 = T(outptr_d0 + cnt_ijk*VW) * xfrm_eta2; //Load 1x

//                 T d2 = T(outptr_d2 + cnt_ijk*VW);//Load 1x
//                 T xfrm_eta0 = T(0.5) * (T(1.0) + Z0[i]); //Load 1x
//                 d2.fma(xfrm_eta0, d0);

//                 T d1 = T(outptr_d1 + cnt_ijk*VW); //Load 1x

//                 if(DEFORMED){
//                     df0 = df_ptr[cnt_ijk*ndf + 0];
//                     df1 = df_ptr[cnt_ijk*ndf + 1];
//                     df2 = df_ptr[cnt_ijk*ndf + 2];
//                     df3 = df_ptr[cnt_ijk*ndf + 3];
//                     df4 = df_ptr[cnt_ijk*ndf + 4];
//                     df5 = df_ptr[cnt_ijk*ndf + 5];
//                     df6 = df_ptr[cnt_ijk*ndf + 6];
//                     df7 = df_ptr[cnt_ijk*ndf + 7];
//                     df8 = df_ptr[cnt_ijk*ndf + 8];
//                 }

//                 //Metric for eta_0, xi_1, eta_2
//                 T out0 = d0 * df0;
//                 out0.fma(d1, df1);
//                 out0.fma(d2, df2);
//                 out0.store(outptr_d0 + cnt_ijk*VW); //Store 1x

//                 T out1 = d0 * df3;
//                 out1.fma(d1, df4);
//                 out1.fma(d2, df5);
//                 out1.store(outptr_d1 + cnt_ijk*VW); //Store 1x

//                 T out2 = d0 * df6;
//                 out2.fma(d1, df7);
//                 out2.fma(d2, df8);
//                 out2.store(outptr_d2 + cnt_ijk*VW); //Store 1x

//             }
//         }
//     }
// }

// template<int NUMQUAD0, int NUMQUAD1, int NUMQUAD2,
//          int VW, bool DEFORMED, class BasisType>
// static void PhysDerivTetKernel(
//     const NekDouble *inptr,
//     const AlignedVector<BasisType> &Z0,
//     const AlignedVector<BasisType> &Z1,
//     const AlignedVector<BasisType> &Z2,
//     const AlignedVector<BasisType> &D0,
//     const AlignedVector<BasisType> &D1,
//     const AlignedVector<BasisType> &D2,
//     const VecData<double, VW> *df_ptr,
//     NekDouble *diff0,
//     NekDouble *diff1,
//     NekDouble *diff2,
//     NekDouble *outptr_d0,
//     NekDouble *outptr_d1,
//     NekDouble *outptr_d2)
// {

//     constexpr int nq0 = NUMQUAD0, nq1 = NUMQUAD1, nq2 = NUMQUAD2;
//     constexpr int ndf = 9;
//     using T = VecData<double, VW>;

//     PhysDerivTensor3DKernel<NUMQUAD0, NUMQUAD1, NUMQUAD2, VW>(
//         inptr,
//         D0, D1, D2,
//         diff0, diff1, diff2);

//     T df0, df1, df2, df3, df4, df5, df6, df7, df8;
//     if(!DEFORMED){
//         df0 = df_ptr[0];
//         df1 = df_ptr[1];
//         df2 = df_ptr[2];
//         df3 = df_ptr[3];
//         df4 = df_ptr[4];
//         df5 = df_ptr[5];
//         df6 = df_ptr[6];
//         df7 = df_ptr[7];
//         df8 = df_ptr[8];
//     }

//     int eta0 = 0;
//     for(int k=0; k < nq2; k++){
//         T xfrm_eta2 = T(2.0) / (T(1.0) - Z2[k]); //Load 1x

//         for(int j = 0; j < nq1; j++){

//             T xfrm_eta1 = T(2.0) / (T(1.0) -  Z1[j]); //Load 1x
//             T xfrm = xfrm_eta1 * xfrm_eta2;

//             for(int i = 0; i < nq0; i++, eta0++)
//             {
//                 T d0 = xfrm * T(diff0 + eta0*VW); //Load 1x

//                 d0.store(outptr_d0 + eta0*VW); //Store 2x
//                 d0.store(diff0 + eta0*VW); //Store partial form for reuse

//             }
//         }
//     }

//     eta0 = 0;
//     for(int k = 0; k < nq2; k++){
//         T xfrm_eta2 = T(2.0) / (T(1.0) - Z2[k]); //Load 1x

//         for(int j = 0; j < nq1; j++){
//             for(int i = 0; i < nq0; i++, eta0++){
//                 T xfrm_eta0 = T(0.5) * (T(1.0) + Z0[i]); //Load 1x

//                 T out0 = xfrm_eta0 * T(diff0 + eta0*VW); //Load 1x
//                 out0.store(diff0 + eta0*VW); //2 * (1 + eta_0) / (1 - eta_1)(1-eta2) | store 1x

//                 T d1 = T(diff1 + eta0*VW); //Load 1x
//                 d1 = xfrm_eta2 * d1;

//                 (out0 + d1).store(outptr_d1 + eta0*VW); //store 1x
//                 d1.store(diff1 + eta0*VW); //store 1x
//             }
//         }
//     }

//     eta0 = 0;
//     for(int k = 0; k < nq2; k++){
//         for(int j = 0; j < nq1; j++){
//             T xfrm_eta1 = T(0.5) * (T(1.0) + Z1[j]); //Load 1x

//             for(int i = 0; i < nq0; i++, eta0++){
//                 T out = T(diff0 + eta0*VW);//Load 1x

//                 T d1 = T(diff1 + eta0*VW); //Load1x
//                 out.fma(d1, xfrm_eta1);

//                 out = out + T(diff2 + eta0*VW); //Load 1x

//                 out.store(outptr_d2 + eta0*VW); //Store 1x
//             }
//         }
//     }

//     int cnt_kji = 0;
//     for(int k = 0; k < nq2; k++){
//         for(int j = 0; j < nq1; j++){
//             for(int i = 0; i < nq0; i++, cnt_kji++){

//                 T d0 = T(outptr_d0 + cnt_kji*VW); //Load 1x
//                 T d1 = T(outptr_d1 + cnt_kji*VW); //Load 1x
//                 T d2 = T(outptr_d2 + cnt_kji*VW); //Load 1x

//                 if(DEFORMED){
//                     df0 = df_ptr[cnt_kji*ndf + 0];
//                     df1 = df_ptr[cnt_kji*ndf + 1];
//                     df2 = df_ptr[cnt_kji*ndf + 2];
//                     df3 = df_ptr[cnt_kji*ndf + 3];
//                     df4 = df_ptr[cnt_kji*ndf + 4];
//                     df5 = df_ptr[cnt_kji*ndf + 5];
//                     df6 = df_ptr[cnt_kji*ndf + 6];
//                     df7 = df_ptr[cnt_kji*ndf + 7];
//                     df8 = df_ptr[cnt_kji*ndf + 8];
//                 }

//                 T out0 = d0 * df0;
//                 out0.fma(d1, df1);
//                 out0.fma(d2, df2);
//                 out0.store(outptr_d0 + cnt_kji*VW);//Store 1x

//                 T out1 = d0 * df3;
//                 out1.fma(d1, df4);
//                 out1.fma(d2, df5);
//                 out1.store(outptr_d1 + cnt_kji*VW);//Store 1x

//                 T out2 = d0 * df6;
//                 out2.fma(d1, df7);
//                 out2.fma(d2, df8);
//                 out2.store(outptr_d2 + cnt_kji*VW);//Store 1x
//             }
//         }
//     }

// }

} // namespace MatrixFree
} // namespace Nektar

#endif
