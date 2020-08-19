#ifndef NEKTAR_LIBRARY_MF_IPRODUCTKERNELS_HPP
#define NEKTAR_LIBRARY_MF_IPRODUCTKERNELS_HPP

#include <LibUtilities/SimdLib/tinysimd.hpp>

namespace Nektar
{
namespace MatrixFree
{

using namespace tinysimd;
using vec_t = simd<NekDouble>;

template<bool SCALE, bool APPEND>
inline static void ScaleAppend(
    vec_t &store,
    vec_t &pos,
    NekDouble scale)
{

    if (SCALE && APPEND)
    {
        store.fma(pos, vec_t(scale));
    }
    else if (APPEND)
    {
        store = store + pos;
    }
    else if (SCALE)
    {
        store = pos * scale;
    }
    else
    {
        store = pos;
    }
}

template<int NUMMODE0, int NUMQUAD0,
         bool SCALE, bool APPEND, bool DEFORMED>
inline static void IProductSegKernel(
    const std::vector<vec_t, allocator<vec_t>> &in,
    const std::vector<vec_t, allocator<vec_t>> &basis0,
    const std::vector<vec_t, allocator<vec_t>> &w0,
    const vec_t *jac,
    std::vector<vec_t, allocator<vec_t>> &out,
    NekDouble scale = 1.0)
{
    constexpr auto nm0 = NUMMODE0;
    constexpr auto nq0 = NUMQUAD0;

    for (int p = 0; p < nm0; ++p)
    {
        vec_t sum = 0.0;
        for (int i = 0; i < nq0; ++i)
        {
            vec_t jac_val;
            if(DEFORMED)
            {
                jac_val = jac[i]; //J for each quadrature point.
            }
            else
            {
                jac_val = jac[0];
            }

            vec_t prod = in[i] * basis0[p*nq0 + i] * jac_val; //Load 2x
            sum.fma(prod, w0[i]); //Load 1x
        }

        //Modes are reversed from what they normally are for tris, tets etc.
        ScaleAppend<SCALE, APPEND>(out[p], sum, scale); // Store x1
    }
}

template<int NUMMODE0, int NUMMODE1,
         int NUMQUAD0, int NUMQUAD1,
         bool SCALE, bool APPEND, bool DEFORMED>
inline static void IProductQuadKernel(
    const std::vector<vec_t, allocator<vec_t>> &in,
    const std::vector<vec_t, allocator<vec_t>> &basis0,
    const std::vector<vec_t, allocator<vec_t>> &basis1,
    const std::vector<vec_t, allocator<vec_t>> &w0,
    const std::vector<vec_t, allocator<vec_t>> &w1,
    const vec_t *jac,
    vec_t *sums_j,
    std::vector<vec_t, allocator<vec_t>> &out,
    NekDouble scale = 1.0)
{
    constexpr auto nm0 = NUMMODE0, nm1 = NUMMODE1;
    constexpr auto nq0 = NUMQUAD0, nq1 = NUMQUAD1;

    for (int p = 0; p < nm0; ++p)
    {
        int cnt_ji = 0;
        for (int j = 0; j < nq1; ++j)
        {

            vec_t sum_j = 0.0;
            for (int i = 0; i < nq0; ++i, ++cnt_ji)
            {
                vec_t jac_val;
                if(DEFORMED)
                {
                    jac_val = jac[j*nq0 + i]; //J for each quadrature point.
                }
                else
                {
                    jac_val = jac[0];
                }

                vec_t prod = in[cnt_ji] * basis0[p*nq0 + i] * jac_val; //Load 2x
                sum_j.fma(prod, w0[i]); //Load 1x
            }

            sums_j[j] = sum_j; //Store 1
        }

        for (int q = 0; q < nm1; ++q)
        {
            vec_t sum = 0.0;
            for (int j = 0; j < nq1; ++j)
            {
                vec_t prod = sums_j[j] * basis1[q*nq1 + j]; //Load 2x
                sum.fma(prod, w1[j]); // Load 1x
            }

            //Modes are reversed from what they normally are for tris, tets etc.
            ScaleAppend<SCALE, APPEND>(out[q*nm1 + p], sum, scale); // Store x1
        }
    }

}
    
template<int NUMMODE0, int NUMMODE1,
         int NUMQUAD0, int NUMQUAD1,
         bool CORRECT, bool SCALE, bool APPEND,
         bool DEFORMED>
inline static void IProductTriKernel(
    const std::vector<vec_t, allocator<vec_t>> &in,
    const std::vector<vec_t, allocator<vec_t>> &basis0,
    const std::vector<vec_t, allocator<vec_t>> &basis1,
    const std::vector<vec_t, allocator<vec_t>> &w0,
    const std::vector<vec_t, allocator<vec_t>> &w1,
    const vec_t* jac,
    vec_t* eta0_sums,
    std::vector<vec_t, allocator<vec_t>> &out,
    NekDouble scale = 1.0)
{
    constexpr auto nm0 = NUMMODE0, nm1 = NUMMODE1;
    constexpr auto nq0 = NUMQUAD0, nq1 = NUMQUAD1;

    int mode = 0;
    for (int p = 0; p < nm0; ++p)
    {

        int eta_idx = 0;
        //Our inner loop is phi_p not phi_pq since we want to put as much work
        //as we can in the p-only loop instead of the full pq loop.
        for (int eta1 = 0; eta1 < nq1; ++eta1)
        {

            vec_t eta0_sum = 0.0;

            for (int eta0 = 0; eta0 < nq0; ++eta0, ++eta_idx)
            {
                //eta0_sum += phi_p(eta0) * fn(eta0, eta1) * J * w0(eta0)
                vec_t jac_val;
                if (DEFORMED)
                {
                    jac_val = jac[eta1*nq0 + eta0];
                }
                else
                {
                    jac_val = jac[0];
                }

                vec_t prod = in[eta_idx] * basis0[p*nq0 + eta0] * jac_val; //Load 2x
                eta0_sum.fma(prod, w0[eta0]); //Load 1x
            }

            eta0_sums[eta1] = eta0_sum;

        }

        for (int q = 0; q < nm1 - p; ++q, ++mode)
        {
            vec_t sum_eta1 = 0.0;
            for (int eta1 = 0; eta1 < nq1; ++eta1)
            {
                vec_t prod = eta0_sums[eta1] * basis1[mode*nq1 + eta1]; //Load 2x
                sum_eta1.fma(prod, w1[eta1]); //Load 1x
            }

            ScaleAppend<SCALE, APPEND>(out[mode], sum_eta1, scale); //Store x1
        }
    }

    //Correction for singular vertex in collpased coordinates.
    //Basically we add phi_1 * phi_01 * (weighting, etc) to mode 00
    //With contributions from every quadrature point
    if (CORRECT)
    {
        int eta_idx = 0;
        vec_t iprod_01 = 0.0;//T(outptr + VW); //Load 1x
        for (int eta1 = 0; eta1 < nq1; ++eta1)
        {
            vec_t preweight_eta1;
            if (DEFORMED)
            {
                preweight_eta1 = w1[eta1] * basis1[nq1 + eta1];
            }
            else
            {
                preweight_eta1 = w1[eta1] * jac[0] * basis1[nq1 + eta1];
            }

            for (int eta0 = 0; eta0 < nq0; ++eta0, ++eta_idx)
            {
                vec_t prod = in[eta_idx] * preweight_eta1 * w0[eta0];

                if (DEFORMED)
                {
                    prod = prod * jac[eta1*nq0 + eta0];
                }

                vec_t basis_val1 = basis0[nq0 + eta0];
                iprod_01.fma(prod, basis_val1);
            }
        }

        ScaleAppend<SCALE, true>(out[1], iprod_01, scale);
    }
}

template<int NUMMODE0, int NUMMODE1, int NUMMODE2,
         int NUMQUAD0, int NUMQUAD1, int NUMQUAD2,
         bool SCALE, bool APPEND, bool DEFORMED>
inline static void IProductHexKernel(
    const std::vector<vec_t, allocator<vec_t>> &in,
    const std::vector<vec_t, allocator<vec_t>> &bdata0,
    const std::vector<vec_t, allocator<vec_t>> &bdata1,
    const std::vector<vec_t, allocator<vec_t>> &bdata2,
    const std::vector<vec_t, allocator<vec_t>> &w0,
    const std::vector<vec_t, allocator<vec_t>> &w1,
    const std::vector<vec_t, allocator<vec_t>> &w2,
    const vec_t *jac,
    vec_t *sums_kj,
    vec_t *sums_k,
    std::vector<vec_t, allocator<vec_t>> &out,
    NekDouble scale = 1.0)
{
    constexpr auto nm0 = NUMMODE0, nm1 = NUMMODE1, nm2 = NUMMODE2;
    constexpr auto nq0 = NUMQUAD0, nq1 = NUMQUAD1, nq2 = NUMQUAD2;

    for (int p = 0; p < nm0; ++p)
    {
        int cnt_kji = 0, cnt_kj = 0;
        for (int k = 0; k < nq2; ++k)
        {
            for (int j = 0; j < nq1; ++j, ++cnt_kj)
            {
                vec_t sum_kj = 0.0;
                for (int i = 0; i < nq0; ++i, ++cnt_kji)
                {

                    vec_t jac_val;
                    if(DEFORMED)
                    {
                        jac_val = jac[nq0*nq1*k + nq0*j + i];
                    }
                    else
                    {
                        jac_val = jac[0];
                    }

                    vec_t prod = in[cnt_kji] * bdata0[i + nq0*p] * jac_val; // load 2x
                    sum_kj.fma(prod, w0[i]); //Load 1x
                }

                sums_kj[cnt_kj] = sum_kj;
            }
        }

        for (int q = 0; q < nm1; ++q)
        {
            cnt_kj = 0;
            for (int k = 0; k < nq2; ++k)
            {
                vec_t sum_k = 0.0;
                for (int j = 0; j < nq1; ++j, ++cnt_kj)
                {
                    vec_t prod = sums_kj[cnt_kj] * bdata1[q*nq1 + j]; //Load 2x
                    sum_k.fma(prod, w1[j]); //Load 1x
                }

                sums_k[k] = sum_k;
            }

            for (int r = 0; r < nm2; ++r)
            {
                vec_t sum = 0.0;

                for(int k = 0; k < nq2; ++k)
                {
                    vec_t prod = sums_k[k] * bdata2[r*nq2 + k]; //Load 2x
                    sum.fma(prod,w2[k]); //Load 1x
                }

                ScaleAppend<SCALE, APPEND>(out[r*nm0*nm1 + q*nm0 + p], sum, scale); // Store x1
            }
        }
    }
}

template<int NUMMODE0, int NUMMODE1, int NUMMODE2,
         int NUMQUAD0, int NUMQUAD1, int NUMQUAD2,
         bool CORRECT, bool SCALE, bool APPEND, bool DEFORMED>
inline static void IProductPrismKernel(
    const std::vector<vec_t, allocator<vec_t>>& in,
    const std::vector<vec_t, allocator<vec_t>>& bdata0,
    const std::vector<vec_t, allocator<vec_t>>& bdata1,
    const std::vector<vec_t, allocator<vec_t>>& bdata2,
    const std::vector<vec_t, allocator<vec_t>>& w0,
    const std::vector<vec_t, allocator<vec_t>>& w1,
    const std::vector<vec_t, allocator<vec_t>>& w2,
    const vec_t* jac,
    vec_t* sums_kj,
    vec_t* sums_k,
    vec_t* corr_q,
    std::vector<vec_t, allocator<vec_t>>& out,
    NekDouble scale = 1.0)
{
    constexpr auto nm0 = NUMMODE0, nm1 = NUMMODE1, nm2 = NUMMODE2;
    constexpr auto nq0 = NUMQUAD0, nq1 = NUMQUAD1, nq2 = NUMQUAD2;

    int mode_pr = 0, mode_pqr = 0;
    for (int p = 0; p < nm0; ++p)
    {

        int cnt_kji = 0, cnt_kj = 0;
        for (int k = 0; k < nq2; ++k)
        {
            for (int j = 0; j < nq1; ++j, ++cnt_kj)
            {
                vec_t sum_kj = 0.0;

                for (int i = 0; i < nq0; ++i, ++cnt_kji)
                {
                    vec_t jac_val;
                    if (DEFORMED)
                    {
                        jac_val = jac[nq0*nq1*k + nq0*j + i];
                    }
                    else{
                        jac_val = jac[0];
                    }

                    vec_t prod = bdata0[nq0*p + i] * jac_val * w0[i]; // load 2x
                    vec_t fn = in[cnt_kji]; //load 1x
                    sum_kj.fma(prod, fn);
                }

                sums_kj[cnt_kj] = sum_kj; //store 1x
            }
        }

        for (int q = 0; q < nm1; ++q)
        {
            cnt_kj = 0;
            for (int k = 0; k < nq2; ++k)
            {

                vec_t sum_k = 0.0;
                for (int j = 0; j < nq1; ++j, ++cnt_kj)
                {
                    sum_k.fma(bdata1[q*nq1 + j] * w1[j], sums_kj[cnt_kj]); //Load 3x
                }

                sums_k[k] = sum_k; //Store 1x
            }

            //Start with nesting. Should be able to move out of q loop and sotre identical copies...
            for (int r = 0; r < nm2-p; ++r, ++mode_pqr)
            {
                vec_t k_sum = 0.0;
                for (int k = 0; k < nq2; ++k)
                {
                    k_sum.fma(bdata2[(mode_pr + r)*nq2 + k] * w2[k], sums_k[k]); //Load 3x
                }

                ScaleAppend<SCALE, APPEND>(out[mode_pqr], k_sum, scale); //Store 1x
            }
        }

        mode_pr += nm2 - p;

    }

    if (CORRECT)
    {
        // Corrections for singular edge
        for (int q = 0; q < nm1; ++q)
        {
            corr_q[q] = 0.0;//T(outptr + (nm2*q + 1)*VW);
        }

        int cnt_kji = 0;
        for (int k = 0; k < nq2; ++k)
        {

            vec_t k_weight = w2[k];
            if (!DEFORMED)
            {
                k_weight = k_weight * jac[0];
            }

            for (int j = 0; j < nq1; ++j)
            {
                vec_t kj_weight = k_weight * w1[j];
                for (int i = 0; i < nq0; ++i, ++cnt_kji)
                {

                    vec_t kji_weight = kj_weight * w0[i];
                    vec_t prod = kji_weight * in[cnt_kji];

                    if (DEFORMED)
                    {
                        prod *= jac[k*nq1*nq0 + j*nq0 + i];
                    }

                    vec_t basis_2 = bdata2[nq2 + k];
                    vec_t basis_0 = bdata0[nq0 + i];
                    //Add phi_1q1 to phi_0q1
                    for (int q = 0; q < nm1; ++q)
                    {
                        vec_t basis_1 = bdata1[q*nq1 + j];

                        corr_q[q].fma(basis_2*basis_1, basis_0*prod);
                    }
                }
            }
        }

        for (int q = 0; q < nm1; ++q)
        {
            ScaleAppend<SCALE, true>(out[nm2*q + 1], corr_q[q], scale);
        }
    }
}

template<int NUMMODE0, int NUMMODE1, int NUMMODE2,
         int NUMQUAD0, int NUMQUAD1, int NUMQUAD2,
         bool CORRECT, bool SCALE, bool APPEND, bool DEFORMED>
inline static void IProductPyrKernel(
    const std::vector<vec_t, allocator<vec_t>>& in,
    const std::vector<vec_t, allocator<vec_t>>& bdata0,
    const std::vector<vec_t, allocator<vec_t>>& bdata1,
    const std::vector<vec_t, allocator<vec_t>>& bdata2,
    const std::vector<vec_t, allocator<vec_t>>& w0,
    const std::vector<vec_t, allocator<vec_t>>& w1,
    const std::vector<vec_t, allocator<vec_t>>& w2,
    const vec_t* jac,
    vec_t* sums_kj,
    vec_t* sums_k,
    std::vector<vec_t, allocator<vec_t>>& out,
    NekDouble scale = 1.0)
{
    constexpr auto nm0 = NUMMODE0, nm1 = NUMMODE1, nm2 = NUMMODE2;
    constexpr auto nq0 = NUMQUAD0, nq1 = NUMQUAD1, nq2 = NUMQUAD2;

    int mode_pqr = 0;
    for (int p = 0; p < nm0; ++p)
    {
        int cnt_kji = 0, cnt_kj = 0;
        for (int k = 0; k < nq2; ++k)
        {
            for (int j = 0; j < nq1; ++j, ++cnt_kj)
            {
                vec_t sum_kj = 0.0;

                for (int i = 0; i < nq0; ++i, ++cnt_kji)
                {
                    vec_t jac_val;
                    if (DEFORMED)
                    {
                        jac_val = jac[nq0*nq1*k + nq0*j + i];
                    }
                    else{
                        jac_val = jac[0];
                    }
                    vec_t prod = bdata0[nq0*p + i] * jac_val * w0[i]; // load 2x
                    vec_t fn = in[cnt_kji]; //load 1x
                    sum_kj.fma(prod, fn);
                }

                sums_kj[cnt_kj] = sum_kj; //store 1x
            }
        }

        for (int q = 0; q < p; ++q)
        {
            cnt_kj = 0;
            for (int k = 0; k < nq2; ++k)
            {
                vec_t sum_k = 0.0;
                for (int j = 0; j < nq1; ++j, ++cnt_kj)
                {
                    sum_k.fma(bdata1[q*nq1+j]*w1[j], sums_kj[cnt_kj]); //Load 3x
                }

                sums_k[k] = sum_k; //Store 1x
            }

            for (int r = 0; r < nm2-p; ++r, ++mode_pqr)
            {
                vec_t k_sum = 0.0;
                for (int k = 0; k < nq2; ++k)
                {
                    k_sum.fma(bdata2[mode_pqr*nq2+k]*w2[k], sums_k[k]); //Load 3x
                }
                
                ScaleAppend<SCALE, APPEND>(out[mode_pqr], k_sum, scale); //Store 1x
            }
        }

        for (int q = p; q < nm1; ++q)
        {
            cnt_kj = 0;
            for (int k = 0; k < nq2; ++k)
            {
                vec_t sum_k = 0.0;
                for (int j = 0; j < nq1; ++j, ++cnt_kj)
                {
                    sum_k.fma(bdata1[q*nq1+j]*w1[j], sums_kj[cnt_kj]); //Load 3x
                }

                sums_k[k] = sum_k; //Store 1x
            }

            for (int r = 0; r < nm2-q; ++r, ++mode_pqr)
            {
                vec_t k_sum = 0.0;
                for (int k = 0; k < nq2; ++k)
                {
                    k_sum.fma(bdata2[mode_pqr*nq2+k]*w2[k], sums_k[k]); //Load 3x
                }

                ScaleAppend<SCALE, APPEND>(out[mode_pqr], k_sum, scale); //Store 1x
            }
        }
    }

    if (CORRECT)
    {
        for (int k = 0, cnt = 0; k < nq2; ++k)
        {
            vec_t tmpQ2 = w2[k]; //Load 1x
            if (!DEFORMED)
            {
                tmpQ2 = tmpQ2 * jac[0];
            }
            
            for (int j = 0; j < nq1; ++j)
            { 
                vec_t tmpQ1 = tmpQ2 * w1[j]; //Load 1x

                for (int i = 0; i < nq0; ++i, ++cnt)
                {
                    // Store jac * quadrature weight
                    vec_t tmpQ = tmpQ1 * w0[i]; //Load 1x
                    vec_t tmpIn = in[cnt]; //Load 1x

                    if (DEFORMED)
                    {
                        tmpQ = tmpQ * jac[k*nq0*nq1 + j*nq0 + i];
                    }

                    // top vertex
                    //
                    // outarray[1] += inarray[cnt] * bdata2[nq2 + k] * (
                    //     bdata0[i]*bdata1[nq1+j] + bdata0[nq0+i]*bdata1[j] +
                    //     bdata0[nq0+i]*bdata1[nq1+j]);

                    vec_t tmp = bdata0[i] * bdata1[nq1+j]; //Load 2x
                    tmp.fma(bdata0[nq0+i], bdata1[j]); //Load 2x
                    tmp.fma(bdata0[nq0+i], bdata1[nq1+j]); //Load 2x
                    tmp = tmp * bdata2[nq2+k]; //Load 1x
                    tmp = tmp * tmpIn;

                    // add to existing entry
                    vec_t tmpOut = tmp * tmpQ;
                    ScaleAppend<SCALE, true>(out[1], tmpOut, scale); //Store 1x
                }
            }
        }
    }
}

template<int NUMMODE0, int NUMMODE1, int NUMMODE2,
         int NUMQUAD0, int NUMQUAD1, int NUMQUAD2,
         bool CORRECT, bool SCALE, bool APPEND,
         bool DEFORMED>
inline static void IProductTetKernel(
    const std::vector<vec_t, allocator<vec_t>>& in,
    const std::vector<vec_t, allocator<vec_t>>& bdata0,
    const std::vector<vec_t, allocator<vec_t>>& bdata1,
    const std::vector<vec_t, allocator<vec_t>>& bdata2,
    const std::vector<vec_t, allocator<vec_t>>& w0,
    const std::vector<vec_t, allocator<vec_t>>& w1,
    const std::vector<vec_t, allocator<vec_t>>& w2,
    const vec_t* jac,
    vec_t* wsp,
    std::vector<vec_t, allocator<vec_t>>& out,
    NekDouble scale = 1.0)
{
    constexpr auto nm0 = NUMMODE0, nm1 = NUMMODE1, nm2 = NUMMODE2;
    constexpr auto nq0 = NUMQUAD0, nq1 = NUMQUAD1, nq2 = NUMQUAD2;

    vec_t* f = wsp;
    vec_t* fb = wsp + nq1 * nq2;

    for (int p = 0, mode = 0, mode2 = 0, cnt_pqr = 0; p < nm0; ++p)
    {
        int cnt_kji = 0;
        for (int k = 0, cnt_kj = 0; k < nq2; ++k)
        {
            for (int j = 0; j < nq1; ++j, ++cnt_kj)
            {
                // Unroll first entry of for loop below and multiply by
                // quadrature weights in dir0 & jacobian.
                vec_t jac_val;
                if (DEFORMED)
                {
                    jac_val = jac[nq0*nq1*k + nq0*j];
                }
                else
                {
                    jac_val = jac[0];
                }

                vec_t f_kj = in[cnt_kji] * bdata0[nq0*p] * jac_val * w0[0]; //Load 3x
                ++cnt_kji;

                for (int i = 1; i < nq0; ++i, ++cnt_kji)
                {
                    if (DEFORMED)
                    {
                        jac_val = jac[nq0*nq1*k + nq0*j + i];
                    }
                    else
                    {
                        jac_val = jac[0];
                    }

                    vec_t inxmm = in[cnt_kji] * bdata0[i + nq0*p] * jac_val; //Load 2x
                    f_kj.fma(inxmm, w0[i]); //Load 1x
                }

                f[cnt_kj] =  f_kj; //Store 1x
            }
        }

        for (int q = 0; q < nm1-p; ++q, ++mode)
        {
            int cnt_kj = 0;
            for (int k = 0; k < nq2; ++k)
            {
                vec_t f_k = bdata1[mode*nq1] * f[cnt_kj] * w1[0]; //Load 3x
                ++cnt_kj;

                for (int j = 1; j < nq1; ++j, ++cnt_kj)
                {
                    vec_t tmp2 = bdata1[mode*nq1 + j] * f[cnt_kj]; //Load 2x
                    f_k.fma(tmp2, w1[j]); //Load 1x
                }

                fb[k] = f_k; //Store 1x
            }

            for (int r = 0; r < nm2-p-q; ++r, ++mode2, ++cnt_pqr)
            {
                vec_t tmp = fb[0] * bdata2[mode2*nq2] * w2[0]; //Load 3x

                for (int k = 1; k < nq2; ++k)
                {
                    vec_t tmp2 = fb[k] * bdata2[mode2*nq2 + k]; //Load 2x
                    tmp.fma(tmp2, w2[k]); //Load 1x
                }

                ScaleAppend<SCALE, APPEND>(out[cnt_pqr], tmp, scale); //Store 1x

            }
        }
    }

    if (CORRECT)
    {
        for (int k = 0, cnt = 0; k < nq2; ++k)
        {

            vec_t tmpQ2 = w2[k]; //Load 1x
            if (!DEFORMED)
            {
                tmpQ2 = tmpQ2 * jac[0];
            }

            for (int j = 0; j < nq1; ++j)
            {
                vec_t tmpQ1 = tmpQ2 * w1[j]; //Load 1x

                for (int i = 0; i < nq0; ++i, ++cnt)
                {
                    // Store jac * quadrature weight
                    vec_t tmpQ = tmpQ1 * w0[i]; //Load 1x
                    vec_t tmpIn = in[cnt]; //Load 1x

                    if (DEFORMED)
                    {
                        tmpQ = tmpQ * jac[k*nq0*nq1 + j*nq0 + i];
                    }

                    // top vertex
                    //
                    // outarray[1] += inarray[cnt] * bdata2[nq2 + k] * (
                    //     bdata0[i]*bdata1[nq1+j] + bdata0[nq0+i]*bdata1[j] +
                    //     bdata0[nq0+i]*bdata1[nq1+j]);

                    vec_t tmp = bdata0[i] * bdata1[nq1+j]; //Load 2x
                    tmp.fma(bdata0[nq0+i], bdata1[j]); //Load 2x
                    tmp.fma(bdata0[nq0+i], bdata1[nq1+j]); //Load 2x
                    tmp = tmp * bdata2[nq2+k]; //Load 1x
                    tmp = tmp * tmpIn;

                    // add to existing entry
                    vec_t tmpOut = tmp * tmpQ;
                    ScaleAppend<SCALE, true>(out[1], tmpOut, scale); //Store 1x

                    // bottom vertex
                    //
                    // outarray[nm2] += inarray[cnt] * bdata2[k] * (
                    //    bdata0[nq0+i] * bdata1[nq1+j]);

                    tmp = bdata0[nq0+i] * bdata1[nq1+j] * bdata2[k] * tmpIn; //Load 3x
                    tmpOut = tmp * tmpQ;
                    ScaleAppend<SCALE, true>(out[nm2], tmpOut, scale); //Store 1x

                    // singular edge
                    for (int r = 1; r < nm2-1; ++r)
                    {
                        // outarray[nm2+r] += inarray[cnt] *
                        //     bdata2[(r+1)*nq2+k] * bdata1[nq1+j] * bdata0[nq0+i];
                        tmp = bdata2[(r+1)*nq2+k] * bdata1[nq1+j] * bdata0[nq0+i] * tmpIn; //Load 3x
                        tmpOut = tmp * tmpQ;
                        ScaleAppend<SCALE, true>(out[nm2+r], tmpOut, scale); //Store 1x
                    }
                }
            }
        }
    }
}

} // namespace MatrixFree
} // namespace Nektar

#endif
