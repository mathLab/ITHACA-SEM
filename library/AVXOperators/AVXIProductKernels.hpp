#ifndef NEKTAR_LIBRARY_AVXIPRODUCTKERNELS_HPP
#define NEKTAR_LIBRARY_AVXIPRODUCTKERNELS_HPP

#include "VecData.hpp"

namespace Nektar
{
namespace AVX
{

template<int VW, bool SCALE, bool APPEND>
inline static void ScaleAppend(
    double *storeptr,
    VecData<double, VW> &pos,
    double scale)
{
    using T = VecData<double, VW>;
    if (SCALE && APPEND)
    {
        T tmp(storeptr);
        tmp.fma(pos, T(scale));
        tmp.store(storeptr);
    }
    else if (APPEND)
    {
        T tmp = T(storeptr) + pos;
        tmp.store(storeptr);
    }
    else if (SCALE)
    {
        pos = pos * T(scale);
        pos.store(storeptr);
    }
    else
    {
        pos.store(storeptr);
    }
}

template<int NUMMODE0, int NUMMODE1,
         int NUMQUAD0, int NUMQUAD1,
         int VW, bool SCALE, bool APPEND, bool DEFORMED,
         class BasisType>
inline static void AVXIProductQuadKernel(
    const double *inptr,
    const AlignedVector<BasisType> &basis0,
    const AlignedVector<BasisType> &basis1,
    const AlignedVector<BasisType> &w0,
    const AlignedVector<BasisType> &w1,
    const VecData<double, VW> *jac,
    VecData<double, VW> *sums_j,
    double *outptr,
    double scale = 1.0)
{
    using T = VecData<double, VW>;

    constexpr int nm0 = NUMMODE0, nm1 = NUMMODE1;
    constexpr int nq0 = NUMQUAD0, nq1 = NUMQUAD1;

    for(int p = 0; p < nm0; p++){
        int cnt_ji = 0;
        for(int j = 0; j < nq1; j++){

            T sum_j = T(0.0);
            for(int i = 0; i < nq0; i++, cnt_ji++){
                T jac_val;
                if(DEFORMED){
                    jac_val = jac[j*nq0 + i]; //J for each quadrature point.
                }
                else{
                    jac_val = jac[0];
                }

                T prod = T(inptr + cnt_ji*VW) * basis0[p*nq0 + i] * jac_val; //Load 2x
                sum_j.fma(prod, w0[i]); //Load 1x
            }

            sums_j[j] = sum_j; //Store 1
        }

        for(int q = 0; q < nm1; q++){
            T sum = T(0.0);
            for(int j = 0; j < nq1; j++){
                T prod = sums_j[j] * basis1[q*nq1 + j]; //Load 2x
                sum.fma(prod, w1[j]); // Load 1x
            }

            //Modes are reversed from what they normally are for tris, tets etc.
            ScaleAppend<VW, SCALE, APPEND>(outptr + (q*nm1 + p)*VW, sum, scale); // Store x1
        }
    }

}

template<int NUMMODE0, int NUMMODE1,
         int NUMQUAD0, int NUMQUAD1,
         int VW, bool CORRECT, bool SCALE, bool APPEND,
         bool DEFORMED, class BasisType>
inline static void AVXIProductTriKernel(
    const double *inptr,
    const AlignedVector<BasisType> &basis0,
    const AlignedVector<BasisType> &basis1,
    const AlignedVector<BasisType> &w0,
    const AlignedVector<BasisType> &w1,
    const VecData<double, VW> *jac,
    VecData<double, VW> *eta0_sums,
    double *outptr,
    double scale = 1.0)
{
    using T = VecData<double, VW>;

    constexpr int nm0 = NUMMODE0, nm1 = NUMMODE1;
    constexpr int nq0 = NUMQUAD0, nq1 = NUMQUAD1;

    int mode = 0;
    for (int p = 0; p < nm0; p++){

        int eta_idx = 0;
        //Our inner loop is phi_p not phi_pq since we want to put as much work
        //as we can in the p-only loop instead of the full pq loop.
        for(int eta1 = 0; eta1 < nq1; eta1++){

            T eta0_sum = T(0.0);

            for(int eta0 = 0; eta0 < nq0; eta0++, eta_idx += VW){
                //eta0_sum += phi_p(eta0) * fn(eta0, eta1) * J * w0(eta0)
                T jac_val;
                if(DEFORMED){
                    jac_val = jac[eta1*nq0 + eta0];
                }
                else{
                    jac_val = jac[0];
                }

                T prod = T(inptr + eta_idx) * basis0[p*nq0 + eta0] * jac_val; //Load 2x
                eta0_sum.fma(prod, w0[eta0]); //Load 1x
            }

            eta0_sums[eta1] = eta0_sum;

        }

        for(int q = 0; q < nm1 - p; q++, mode++)
        {
            T sum_eta1 = 0.0;
            for(int eta1 = 0; eta1 < nq1; eta1++)
            {
                T prod = eta0_sums[eta1] * basis1[mode*nq1 + eta1]; //Load 2x
                sum_eta1.fma(prod, w1[eta1]); //Load 1x
            }

            ScaleAppend<VW, SCALE, APPEND>(outptr + mode*VW, sum_eta1, scale); //Store x1
        }
    }

    //Correction for singular vertex in collpased coordinates.
    //Basically we add phi_1 * phi_01 * (weighting, etc) to mode 00
    //With contributions from every quadrature point
    if (CORRECT)
    {
        int eta_idx = 0;
        T iprod_01 = 0.0;//T(outptr + VW); //Load 1x
        for(int eta1 = 0; eta1 < nq1; eta1++){

            T preweight_eta1;
            if(DEFORMED){
                preweight_eta1 = w1[eta1] * basis1[nq1 + eta1];
            }
            else{
                preweight_eta1 = w1[eta1] * jac[0] * basis1[nq1 + eta1];
            }

            for(int eta0 = 0; eta0 < nq0; eta0++, eta_idx += VW){
                T prod = T(inptr + eta_idx) * preweight_eta1 * w0[eta0];

                if(DEFORMED){
                    prod = prod * jac[eta1*nq0 + eta0];
                }

                T basis_val1 = basis0[nq0 + eta0];
                iprod_01.fma(prod, basis_val1);
            }
        }

        ScaleAppend<VW, SCALE, true>(outptr + VW, iprod_01, scale);
    }
}

template<int NUMMODE0, int NUMMODE1, int NUMMODE2,
         int NUMQUAD0, int NUMQUAD1, int NUMQUAD2,
         int VW, bool SCALE, bool APPEND, bool DEFORMED, class BasisType>
inline static void AVXIProductHexKernel(
    const double *inptr,
    const AlignedVector<BasisType> &bdata0,
    const AlignedVector<BasisType> &bdata1,
    const AlignedVector<BasisType> &bdata2,
    const AlignedVector<BasisType> &w0,
    const AlignedVector<BasisType> &w1,
    const AlignedVector<BasisType> &w2,
    const VecData<double, VW> *jac,
    VecData<double, VW> *sums_kj,
    VecData<double, VW> *sums_k,
    double *outptr,
    double scale = 1.0)
{
    using T = VecData<double, VW>;

    constexpr int nm0 = NUMMODE0, nm1 = NUMMODE1, nm2 = NUMMODE2;
    constexpr int nq0 = NUMQUAD0, nq1 = NUMQUAD1, nq2 = NUMQUAD2;

    for (int p = 0; p < nm0; p++)
    {
        int cnt_kji = 0, cnt_kj = 0;
        for(int k = 0; k < nq2; k++){
            for(int j = 0; j < nq1; j++, cnt_kj++){
                T sum_kj = T(0.0);

                for(int i = 0; i < nq0; i++, cnt_kji++){

                    T jac_val;
                    if(DEFORMED){
                        jac_val = jac[nq0*nq1*k + nq0*j + i];
                    }
                    else{
                        jac_val = jac[0];
                    }

                    T prod = T(inptr + cnt_kji*VW); // Load 1x
                    prod = prod * bdata0[i + nq0*p] * jac_val; // load 1x
                    sum_kj.fma(prod, w0[i]); //Load 1x
                }

                sums_kj[cnt_kj] = sum_kj;
            }
        }

        for(int q = 0; q < nm1; q++){
            cnt_kj = 0;
            for(int k = 0; k < nq2; k++){
                T sum_k = T(0.0);
                for(int j = 0; j < nq1; j++, cnt_kj++){
                    T prod = sums_kj[cnt_kj] * bdata1[q*nq1 + j]; //Load 2x
                    sum_k.fma(prod, w1[j]); //Load 1x
                }

                sums_k[k] = sum_k;
            }

            for(int r = 0; r < nm2; r++){
                T sum = T(0.0);

                for(int k = 0; k < nq2; k++){
                    T prod = sums_k[k] * bdata2[r*nq2 + k]; //Load 2x
                    sum.fma(prod,w2[k]); //Load 1x
                }

                ScaleAppend<VW, SCALE, APPEND>(outptr + (r*nm0*nm1 + q*nm0 + p) * VW, sum, scale); // Store x1
            }
        }
    }
}

template<int NUMMODE0, int NUMMODE1, int NUMMODE2,
         int NUMQUAD0, int NUMQUAD1, int NUMQUAD2,
         int VW, bool CORRECT, bool SCALE, bool APPEND, bool DEFORMED,
         class BasisType>
inline static void AVXIProductPrismKernel(
    const double *inptr,
    const AlignedVector<BasisType> &bdata0,
    const AlignedVector<BasisType> &bdata1,
    const AlignedVector<BasisType> &bdata2,
    const AlignedVector<BasisType> &w0,
    const AlignedVector<BasisType> &w1,
    const AlignedVector<BasisType> &w2,
    const VecData<double, VW> *jac,
    VecData<double, VW> *sums_kj,
    VecData<double, VW> *sums_k,
    VecData<double, VW> *corr_q,
    double *outptr,
    double scale = 1.0)
{
    using T = VecData<double, VW>;

    constexpr int nm0 = NUMMODE0, nm1 = NUMMODE1, nm2 = NUMMODE2;
    constexpr int nq0 = NUMQUAD0, nq1 = NUMQUAD1, nq2 = NUMQUAD2;

    int mode_pr = 0, mode_pqr = 0;
    for(int p = 0; p < nm0; p++){

        int cnt_kji = 0, cnt_kj = 0;
        for(int k = 0; k < nq2; k++){
            for(int j = 0; j < nq1; j++, cnt_kj++){
                T sum_kj = T(0.0);

                for(int i = 0; i < nq0; i++, cnt_kji++){
                    T jac_val;
                    if(DEFORMED){
                        jac_val = jac[nq0*nq1*k + nq0*j + i];
                    }
                    else{
                        jac_val = jac[0];
                    }

                    T prod = bdata0[nq0*p + i] * jac_val * w0[i]; // load 2x
                    T fn = T(inptr + cnt_kji*VW); //load 1x
                    sum_kj.fma(prod, fn);
                }

                sums_kj[cnt_kj] = sum_kj; //store 1x
            }
        }

        for(int q = 0; q < nm1; q++){
            cnt_kj = 0;
            for(int k = 0; k < nq2; k++){

                T sum_k = T(0.0);
                for(int j = 0; j < nq1; j++, cnt_kj++)
                {
                    sum_k.fma(bdata1[q*nq1 + j] * w1[j], sums_kj[cnt_kj]); //Load 3x
                }

                sums_k[k] = sum_k; //Store 1x
            }

            //Start with nesting. Should be able to move out of q loop and sotre identical copies...
            for(int r = 0; r < nm2-p; r++, mode_pqr++){
                T k_sum = T(0.0);
                for(int k = 0; k < nq2; k++){
                    k_sum.fma(bdata2[(mode_pr + r)*nq2 + k] * w2[k], sums_k[k]); //Load 3x
                }

                ScaleAppend<VW, SCALE, APPEND>(outptr + mode_pqr*VW, k_sum, scale); //Store 1x
            }
        }

        mode_pr += nm2-p;

    }

    if (CORRECT)
    {
        //Corrections for singular edge
        for(int q = 0; q < nm1; q++){
            corr_q[q] = 0.0;//T(outptr + (nm2*q + 1)*VW);
        }

        int cnt_kji = 0;
        for(int k = 0; k < nq2; k++){

            T k_weight = w2[k];
            if(!DEFORMED){
                k_weight = k_weight * jac[0];
            }

            for(int j = 0; j < nq1; j++){
                T kj_weight = k_weight * w1[j];
                for(int i = 0; i < nq0; i++, cnt_kji++){

                    T kji_weight = kj_weight * w0[i];
                    T prod = kji_weight * T(inptr + cnt_kji*VW);

                    if(DEFORMED){
                        prod = prod * jac[k*nq1*nq0 + j*nq0 + i];
                    }

                    T basis_2 = bdata2[nq2 + k];
                    T basis_0 = bdata0[nq0 + i];
                    //Add phi_1q1 to phi_0q1
                    for(int q = 0; q < nm1; q++)
                    {
                        T basis_1 = bdata1[q*nq1 + j];

                        corr_q[q].fma(basis_2*basis_1, basis_0*prod);
                    }
                }
            }
        }

        for(int q = 0; q < nm1; q++){
            ScaleAppend<VW, SCALE, true>(outptr + (nm2*q + 1)*VW, corr_q[q], scale);
        }
    }
}

template<int NUMMODE0, int NUMMODE1, int NUMMODE2,
         int NUMQUAD0, int NUMQUAD1, int NUMQUAD2,
         int VW, bool CORRECT, bool SCALE, bool APPEND,
         bool DEFORMED, class BasisType>
inline static void AVXIProductTetKernel(
    const double *in,
    const AlignedVector<BasisType> &bdata0,
    const AlignedVector<BasisType> &bdata1,
    const AlignedVector<BasisType> &bdata2,
    const AlignedVector<BasisType> &w0,
    const AlignedVector<BasisType> &w1,
    const AlignedVector<BasisType> &w2,
    const VecData<double, VW> *jac,
    VecData<double, VW> *wsp,
    double *out,
    double scale = 1.0)
{
    using T = VecData<double, VW>;

    constexpr int nm0 = NUMMODE0, nm1 = NUMMODE1, nm2 = NUMMODE2;
    constexpr int nq0 = NUMQUAD0, nq1 = NUMQUAD1, nq2 = NUMQUAD2;

    T *f = wsp, *fb = wsp + nq1 * nq2;
    int mode = 0, mode2 = 0, cnt_pqr = 0;

    for (int p = 0; p < nm0; ++p)
    {
        int cnt_kji = 0, cnt_kj = 0;
        for (int k = 0; k < nq2; ++k)
        {
            for (int j = 0; j < nq1; ++j, ++cnt_kj)
            {
                // Unroll first entry of for loop below and multiply by
                // quadrature weights in dir0 & jacobian.
                T jac_val;
                if(DEFORMED){
                    jac_val = jac[nq0*nq1*k + nq0*j];
                }
                else{
                    jac_val = jac[0];
                }

                T f_kj = T(in + cnt_kji) * bdata0[nq0*p] * jac_val * w0[0]; //Load 3x
                cnt_kji += VW;

                for (int i = 1; i < nq0; ++i, cnt_kji += VW)
                {
                    if(DEFORMED){
                        jac_val = jac[nq0*nq1*k + nq0*j + i];
                    }
                    else{
                        jac_val = jac[0];
                    }

                    T inxmm = T(in + cnt_kji) * bdata0[i + nq0*p] * jac_val; //Load 2x
                    f_kj.fma(inxmm, w0[i]); //Load 1x
                }

                f[cnt_kj] =  f_kj; //Store 1x
            }
        }

        for (int q = 0; q < nm1-p; ++q, ++mode)
        {
            cnt_kj = 0;
            for (int k = 0; k < nq2; ++k)
            {
                T f_k = bdata1[mode*nq1] * f[cnt_kj] * w1[0]; //Load 3x
                ++cnt_kj;

                for (int j = 1; j < nq1; ++j, ++cnt_kj)
                {
                    T tmp2 = bdata1[mode*nq1 + j] * f[cnt_kj]; //Load 2x
                    f_k.fma(tmp2, w1[j]); //Load 1x
                }

                fb[k] = f_k; //Store 1x
            }

            // TODO increase mode somehow here

            for (int r = 0; r < nm2-p-q; ++r, ++mode2, cnt_pqr += VW)
            {
                T tmp = fb[0] * bdata2[mode2*nq2] * w2[0]; //Load 3x

                for (int k = 1; k < nq2; ++k)
                {
                    T tmp2 = fb[k] * bdata2[mode2*nq2 + k]; //Load 2x
                    tmp.fma(tmp2, w2[k]); //Load 1x
                }

                ScaleAppend<VW, SCALE, APPEND>(out + cnt_pqr, tmp, scale); //Store 1x

                // TODO increase mode2 somehow here
            }
        }
    }

    if (CORRECT)
    {
        for (int k = 0, cnt = 0; k < nq2; ++k)
        {

            T tmpQ2 = w2[k]; //Load 1x
            if(!DEFORMED){
                tmpQ2 = tmpQ2 * jac[0];
            }

            for (int j = 0; j < nq1; ++j)
            {
                T tmpQ1 = tmpQ2 * w1[j]; //Load 1x

                for (int i = 0; i < nq0; ++i, cnt += VW)
                {
                    // Store jac * quadrature weight
                    T tmpQ = tmpQ1 * w0[i]; //Load 1x
                    T tmpIn = T(in + cnt); //Load 1x

                    if(DEFORMED){
                        tmpQ = tmpQ * jac[k*nq0*nq1 + j*nq0 + i];
                    }

                    // top vertex
                    //
                    // outarray[1] += inarray[cnt] * bdata2[nq2 + k] * (
                    //     bdata0[i]*bdata1[nq1+j] + bdata0[nq0+i]*bdata1[j] +
                    //     bdata0[nq0+i]*bdata1[nq1+j]);

                    T tmp = bdata0[i] * bdata1[nq1+j]; //Load 2x
                    tmp.fma(bdata0[nq0+i], bdata1[j]); //Load 2x
                    tmp.fma(bdata0[nq0+i], bdata1[nq1+j]); //Load 2x
                    tmp = tmp * bdata2[nq2+k]; //Load 1x
                    tmp = tmp * tmpIn;

                    // add to existing entry
                    T tmpOut = tmp * tmpQ;
                    ScaleAppend<VW, SCALE, true>(out + VW, tmpOut, scale); //Store 1x

                    // bottom vertex
                    //
                    // outarray[nm2] += inarray[cnt] * bdata2[k] * (
                    //    bdata0[nq0+i] * bdata1[nq1+j]);

                    tmp = bdata0[nq0+i] * bdata1[nq1+j] * bdata2[k] * tmpIn; //Load 3x
                    tmpOut = tmp * tmpQ;
                    ScaleAppend<VW, SCALE, true>(out + nm2 * VW, tmpOut, scale); //Store 1x

                    // singular edge
                    for (int r = 1; r < nm2-1; ++r)
                    {
                        // outarray[nm2+r] += inarray[cnt] *
                        //     bdata2[(r+1)*nq2+k] * bdata1[nq1+j] * bdata0[nq0+i];
                        tmp = bdata2[(r+1)*nq2+k] * bdata1[nq1+j] * bdata0[nq0+i] * tmpIn; //Load 3x
                        tmpOut = tmp * tmpQ;
                        ScaleAppend<VW, SCALE, true>(out + (nm2+r) * VW, tmpOut, scale); //Store 1x
                    }
                }
            }
        }
    }
}

} // namespace AVX
} // namespace Nektar

#endif
