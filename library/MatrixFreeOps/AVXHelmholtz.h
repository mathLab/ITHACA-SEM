#ifndef AVXHELMHOLTZ_H
#define AVXHELMHOLTZ_H

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ShapeType.hpp>
#include <LibUtilities/Foundations/Basis.h>

#include "VecData.hpp"
#include "Operator.hpp"

#include "AVXIProduct.h"
#include "AVXBwdTrans.h"
#include "AVXPhysDeriv.h"
#include "AVXIProductKernels.hpp"
#include "AVXBwdTransKernels.hpp"
#include "AVXPhysDerivKernels.hpp"

template<int VW, bool DEFORMED=false>
struct AVXHelmholtzQuad : public Helmholtz, public AVXHelper<VW, 2, DEFORMED>
{
    AVXHelmholtzQuad(std::vector<LibUtilities::BasisSharedPtr> basis,
                     int nElmt)
    : Helmholtz(basis, nElmt),
      AVXHelper<VW, 2, DEFORMED>(basis, nElmt),
      m_nmTot(LibUtilities::StdQuadData::getNumberOfCoefficients(
                  this->m_nm[0], this->m_nm[1])),
      m_lambda(1.0)
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<AVXHelmholtzQuad<VW, DEFORMED>>(basis, nElmt);
    }

    static NekDouble FlopsPerElement(
        const int nm,
        const int nq0,
        const int nq1)
    {

        const double bwdTrans = AVXBwdTransQuad<VW>::FlopsPerElement(nm,nq0,nq1);
        const double iprod1 = AVXIProductQuad<VW>::FlopsPerElement(nm,nq0,nq1);
        const double physDeriv = AVXPhysDerivQuad<VW>::FlopsPerElement(nq0,nq1);

        double metrics;
        if(DEFORMED){
            metrics = nq1 * nq0 * (9 + 6);
        }
        else{
            const double metric_calc = 9;
            const double metrics_x0 = nq0 * nq1 * 3;
            const double metrics_x1 = nq0 * nq1 * 3;
            metrics = metric_calc + metrics_x0 + metrics_x1;
        }

        const double iprod2 = AVXIProductQuad<VW>::FlopsPerElement(nm,nq0,nq1);
        const double iprod3 = AVXIProductQuad<VW>::FlopsPerElement(nm, nq0, nq1);

        return bwdTrans + iprod1 + physDeriv +
            metrics + iprod2 + iprod3;
    }

    virtual double GFlops() override
    {
        const int nm = m_basis[0]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();

        const double flops = m_nElmt * AVXHelmholtzQuad::FlopsPerElement(nm, nq0, nq1);
        return flops * 1e-9;
    }

    virtual NekDouble Ndof() override
    {
        return m_nmTot * this->m_nElmt;
    }

    virtual void operator()(const Array<OneD, const NekDouble> &in,
                                  Array<OneD,       NekDouble> &out) override
    {
        switch(m_basis[0]->GetNumModes())
        {
            case 2:  AVXHelmholtzQuadImpl<2 ,2 ,3 ,3 >(in, out); break;
            case 3:  AVXHelmholtzQuadImpl<3 ,3 ,4 ,4 >(in, out); break;
            case 4:  AVXHelmholtzQuadImpl<4 ,4 ,5 ,5 >(in, out); break;
            case 5:  AVXHelmholtzQuadImpl<5 ,5 ,6 ,6 >(in, out); break;
            case 6:  AVXHelmholtzQuadImpl<6 ,6 ,7 ,7 >(in, out); break;
            case 7:  AVXHelmholtzQuadImpl<7 ,7 ,8 ,8 >(in, out); break;
            case 8:  AVXHelmholtzQuadImpl<8 ,8 ,9 ,9 >(in, out); break;
            case 9:  AVXHelmholtzQuadImpl<9 ,9 ,10,10>(in, out); break;
            case 10: AVXHelmholtzQuadImpl<10,10,11,11>(in, out); break;
            case 11: AVXHelmholtzQuadImpl<11,11,12,12>(in, out); break;
        }
    }

    template<int NM0, int NM1, int NQ0, int NQ1>
    void AVXHelmholtzQuadImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &out)
    {
        using T = VecData<double, VW>;

        // Pointers for input/output across blocks
        const NekDouble *inptr = &input[0];
        NekDouble *outptr = &out[0];

        constexpr int nqBlocks = NQ0 * NQ1 * VW;
        constexpr int nmBlocks = NM0 * NM1 * VW;
        constexpr int ndf = 4;
        constexpr int nq = NQ0 * NQ1;

        // Allocate sufficient workspace for backwards transform and inner
        // product kernels.
        constexpr int wspInnerProd = NQ1;
        constexpr int wspBwdTrans = NM0 * NQ0;
        constexpr int wspSize = wspInnerProd > wspBwdTrans ?
            wspInnerProd : wspBwdTrans;

        // Workspace for kernels
        T wsp[wspSize];

        // Aligned stack storage
        std::aligned_storage<sizeof(double), SIMD_WIDTH_BYTES>::type bwd_storage[nqBlocks];
        std::aligned_storage<sizeof(double), SIMD_WIDTH_BYTES>::type deriv0_storage[nqBlocks];
        std::aligned_storage<sizeof(double), SIMD_WIDTH_BYTES>::type deriv1_storage[nqBlocks];
        double *bwd = reinterpret_cast<double *>(&bwd_storage[0]);
        double *deriv0 = reinterpret_cast<double *>(&deriv0_storage[0]);
        double *deriv1 = reinterpret_cast<double *>(&deriv1_storage[0]);

        for (int e = 0; e < this->m_nBlocks; e++)
        {
            // Precompute Laplacian metrics
            T metric00, metric01, metric11;
            const VecData<double, VW> *jac_ptr;
            if(!DEFORMED){
                metric00 = this->m_df[e*ndf] * this->m_df[e*ndf];
                metric00.fma(this->m_df[e*ndf + 2], this->m_df[e*ndf + 2]);
                metric01 = this->m_df[e*ndf] * this->m_df[e*ndf + 1];
                metric01.fma(this->m_df[e*ndf + 2], this->m_df[e*ndf + 3]);
                metric11 = this->m_df[e*ndf + 1] * this->m_df[e*ndf + 1];
                metric11.fma(this->m_df[e*ndf + 3], this->m_df[e*ndf + 3]);
                jac_ptr = &(this->m_jac[e]);
            }
            else{
                jac_ptr = &(this->m_jac[e*nq]);
            }


            // Step 1: BwdTrans
            AVXBwdTransQuadKernel<NM0, NM1, NQ0, NQ1, VW>(
                inptr, this->m_bdata[0], this->m_bdata[1], wsp, bwd);

            // Step 2: inner product for mass matrix operation
            AVXIProductQuadKernel<NM0, NM1, NQ0, NQ1, VW, true, false, DEFORMED>(
                bwd, this->m_bdata[0], this->m_bdata[1], this->m_w[0],
                this->m_w[1], jac_ptr, wsp, outptr, m_lambda);

            // Step 3: take derivatives in standard space
            AVXPhysDerivTensor2DKernel<NQ0, NQ1, VW>(
                bwd, this->m_D[0], this->m_D[1], deriv0, deriv1);

            // Step 4: Apply Laplacian metrics & inner product
            if(DEFORMED){
                //Since we need a unique metric for each i,j, we can't take advantage of the factorizaiton used in the uniform case
                // Direction 1
                int cnt = 0;
                for(int j = 0; j < NQ1; j++)
                {
                    for (int i = 0; i < NQ0; i++, cnt += VW)
                    {
                        VecData<double, VW> *df_start = &(this->m_df[e*ndf*nq]);
                        int ji = j*NQ0 + i;
                        metric00 = df_start[ji*ndf] * df_start[ji*ndf];
                        metric00.fma(df_start[ji*ndf + 2], df_start[ji*ndf + 2]);
                        metric01 = df_start[ji*ndf] * df_start[ji*ndf + 1];
                        metric01.fma(df_start[ji*ndf + 2], df_start[ji*ndf + 3]);
                        metric11 = df_start[ji*ndf + 1] * df_start[ji*ndf + 1];
                        metric11.fma(df_start[ji*ndf + 3], df_start[ji*ndf + 3]);

                        T d0(deriv0 + cnt), d1(deriv1 + cnt);
                        T tmp = metric00 * d0;

                        tmp.fma(metric01, d1);
                        tmp.store(bwd + cnt);

                        tmp = metric01 * d0;
                        tmp.fma(metric11, d1);
                        tmp.store(deriv0 + cnt);
                    }
                }

                AVXIProductQuadKernel<NM0, NM1, NQ0, NQ1, VW, false, true, true>(
                    bwd, this->m_dbdata[0], this->m_bdata[1], this->m_w[0],
                    this->m_w[1], jac_ptr, wsp, outptr);

                AVXIProductQuadKernel<NM0, NM1, NQ0, NQ1, VW, false, true, true>(
                    deriv0, this->m_bdata[0], this->m_dbdata[1], this->m_w[0],
                    this->m_w[1], jac_ptr, wsp, outptr);
            }
            else{
                // Direction 1
                for (int i = 0; i < nqBlocks; i += VW)
                {
                    T tmp = metric00 * T(deriv0 + i);
                    tmp.fma(metric01, T(deriv1 + i));
                    tmp.store(bwd + i);
                }

                AVXIProductQuadKernel<NM0, NM1, NQ0, NQ1, VW, false, true, false>(
                    bwd, this->m_dbdata[0], this->m_bdata[1], this->m_w[0],
                    this->m_w[1], jac_ptr, wsp, outptr);

                // Direction 2
                for (int i = 0; i < nqBlocks; i += VW)
                {
                    T tmp = metric01 * T(deriv0 + i);
                    tmp.fma(metric11, T(deriv1 + i));
                    tmp.store(bwd + i);
                }

                AVXIProductQuadKernel<NM0, NM1, NQ0, NQ1, VW, false, true, false>(
                    bwd, this->m_bdata[0], this->m_dbdata[1], this->m_w[0],
                    this->m_w[1], jac_ptr, wsp, outptr);
            }

            inptr += nmBlocks;
            outptr += nmBlocks;
        }

    }

private:
    int m_nmTot;
    double m_lambda;
};

template<int VW, bool DEFORMED = false>
struct AVXHelmholtzTri : public Helmholtz, public AVXHelper<VW, 2, DEFORMED>
{
    AVXHelmholtzTri(std::vector<LibUtilities::BasisSharedPtr> basis,
                     int nElmt)
    : Helmholtz(basis, nElmt),
      AVXHelper<VW, 2, DEFORMED>(basis, nElmt),
      m_nmTot(LibUtilities::StdTriData::getNumberOfCoefficients(
                  this->m_nm[0], this->m_nm[1])),
      m_lambda(1.0),
      m_h0(basis[0]->GetNumPoints()),
      m_h1(basis[1]->GetNumPoints())
    {
        const int nq0 = basis[0]->GetNumPoints();
        const int nq1 = basis[1]->GetNumPoints();

        const Array<OneD, const NekDouble> &z0 = basis[0]->GetZ();
        const Array<OneD, const NekDouble> &z1 = basis[1]->GetZ();

        for (int i = 0; i < nq0; ++i)
        {
            m_h0[i] = 0.5 * (1 + z0[i]);
        }

        for (int j = 0; j < nq1; ++j)
        {
            m_h1[j] = 2.0 / (1 - z1[j]);
        }
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<AVXHelmholtzTri<VW,DEFORMED>>(basis, nElmt);
    }

    static NekDouble FlopsPerElement(
        const int nm,
        const int nq0,
        const int nq1)
    {

        const double bwdTrans = AVXBwdTransTri<VW>::FlopsPerElement(nm,nq0,nq1);
        const double iprod1 = AVXIProductTri<VW>::FlopsPerElement(nm, nq0, nq1);
        const double physDeriv = AVXPhysDerivTri<VW>::FlopsPerElement(nq0,nq1);
        const double metric = nq1 * nq0 * (5 + 7 + 3 + 3 + 3);
        const double iprod2 = AVXIProductTri<VW>::FlopsPerElement(nm, nq0, nq1);
        const double iprod3 = AVXIProductTri<VW>::FlopsPerElement(nm, nq0, nq1);

        return bwdTrans + iprod1 + physDeriv + metric + iprod2 + iprod3;
    }

    virtual double GFlops() override
    {
        const int nm = m_basis[0]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();

        const double flops = m_nElmt * AVXHelmholtzTri::FlopsPerElement(nm, nq0, nq1);
        return flops * 1e-9;
    }

    virtual NekDouble Ndof() override
    {
        return m_nmTot * this->m_nElmt;
    }

    virtual void operator()(const Array<OneD, const NekDouble> &in,
                                  Array<OneD,       NekDouble> &out) override
    {
        if (m_basis[0]->GetBasisType() == LibUtilities::eModified_A)
        {
            switch(m_basis[0]->GetNumModes())
            {
                case 2:  AVXHelmholtzTriImpl<2 ,2 ,3 ,2 ,true>(in, out); break;
                case 3:  AVXHelmholtzTriImpl<3 ,3 ,4 ,3 ,true>(in, out); break;
                case 4:  AVXHelmholtzTriImpl<4 ,4 ,5 ,4 ,true>(in, out); break;
                case 5:  AVXHelmholtzTriImpl<5 ,5 ,6 ,5 ,true>(in, out); break;
                case 6:  AVXHelmholtzTriImpl<6 ,6 ,7 ,6 ,true>(in, out); break;
                case 7:  AVXHelmholtzTriImpl<7 ,7 ,8 ,7 ,true>(in, out); break;
                case 8:  AVXHelmholtzTriImpl<8 ,8 ,9 ,8 ,true>(in, out); break;
                case 9:  AVXHelmholtzTriImpl<9 ,9 ,10,9 ,true>(in, out); break;
                case 10: AVXHelmholtzTriImpl<10,10,11,10,true>(in, out); break;
                case 11: AVXHelmholtzTriImpl<11,11,12,11,true>(in, out); break;
            }
        }
        else
        {
            switch(m_basis[0]->GetNumModes())
            {
                case 2:  AVXHelmholtzTriImpl<2 ,2 ,3 ,2 ,false>(in, out); break;
                case 3:  AVXHelmholtzTriImpl<3 ,3 ,4 ,3 ,false>(in, out); break;
                case 4:  AVXHelmholtzTriImpl<4 ,4 ,5 ,4 ,false>(in, out); break;
                case 5:  AVXHelmholtzTriImpl<5 ,5 ,6 ,5 ,false>(in, out); break;
                case 6:  AVXHelmholtzTriImpl<6 ,6 ,7 ,6 ,false>(in, out); break;
                case 7:  AVXHelmholtzTriImpl<7 ,7 ,8 ,7 ,false>(in, out); break;
                case 8:  AVXHelmholtzTriImpl<8 ,8 ,9 ,8 ,false>(in, out); break;
                case 9:  AVXHelmholtzTriImpl<9 ,9 ,10,9 ,false>(in, out); break;
                case 10: AVXHelmholtzTriImpl<10,10,11,10,false>(in, out); break;
                case 11: AVXHelmholtzTriImpl<11,11,12,11,false>(in, out); break;
            }
        }
    }

    template<int NM0, int NM1, int NQ0, int NQ1, bool CORRECT>
    void AVXHelmholtzTriImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &out)
    {
        using T = VecData<double, VW>;

        // Pointers for input/output across blocks
        const NekDouble *inptr = &input[0];
        NekDouble *outptr = &out[0];

        constexpr int nqBlocks = NQ0 * NQ1 * VW;
        const int nmBlocks = m_nmTot * VW;
        constexpr int ndf = 4;
        constexpr int nq = NQ0*NQ1;

        // Allocate sufficient workspace for backwards transform and inner
        // product kernels.
        constexpr int wspInnerProd = NQ1;
        constexpr int wspBwdTrans = NM0;
        constexpr int wspSize = wspInnerProd > wspBwdTrans ?
            wspInnerProd : wspBwdTrans;

        // Workspace for kernels
        T wsp[wspSize];

        // Aligned stack storage
        std::aligned_storage<sizeof(double), 64>::type bwd_storage[nqBlocks];
        std::aligned_storage<sizeof(double), 64>::type deriv0_storage[nqBlocks];
        std::aligned_storage<sizeof(double), 64>::type deriv1_storage[nqBlocks];
        std::aligned_storage<sizeof(double), 64>::type metric_storage[nqBlocks];
        double *bwd = reinterpret_cast<double *>(&bwd_storage[0]);
        double *deriv0 = reinterpret_cast<double *>(&deriv0_storage[0]);
        double *deriv1 = reinterpret_cast<double *>(&deriv1_storage[0]);
        double *metric = reinterpret_cast<double *>(&metric_storage[0]);

        T df0, df1, df2, df3;

        for (int e = 0; e < this->m_nBlocks; e++)
        {
            const VecData<double, VW> *jac_ptr;
            if(DEFORMED){
                jac_ptr = &(this->m_jac[e*nq]);
            }
            else{
                jac_ptr = &(this->m_jac[e]);
            }

            // Step 1: BwdTrans
            AVXBwdTransTriKernel<NM0, NM1, NQ0, NQ1, VW, CORRECT>(
                inptr, this->m_bdata[0], this->m_bdata[1], wsp, bwd);

            // Step 2: inner product for mass matrix operation
            AVXIProductTriKernel<NM0, NM1, NQ0, NQ1, VW, CORRECT, true, false, DEFORMED>(
                bwd, this->m_bdata[0], this->m_bdata[1], this->m_w[0],
                this->m_w[1], jac_ptr, wsp, outptr, m_lambda);

            // Step 3: take derivatives in collapsed coordinate space
            AVXPhysDerivTensor2DKernel<NQ0, NQ1, VW>(
                bwd, this->m_D[0], this->m_D[1], deriv0, deriv1);

            if(!DEFORMED){
                df0 = this->m_df[e*ndf];
                df1 = this->m_df[e*ndf + 1];
                df2 = this->m_df[e*ndf + 2];
                df3 = this->m_df[e*ndf + 3];
            }

            // Step 4a: Construct Laplacian metrics
            for (int j = 0, cnt = 0; j < NQ1; ++j)
            {
                T h1j = m_h1[j];
                for (int i = 0; i < NQ0; ++i, cnt += VW)
                {
                    if(DEFORMED){
                        int ji = j*NQ0 + i;
                        VecData<double, VW> *df_start = &(this->m_df[e*ndf*nq]);
                        df0 = df_start[ji*ndf];
                        df1 = df_start[ji*ndf + 1];
                        df2 = df_start[ji*ndf + 2];
                        df3 = df_start[ji*ndf + 3];
                    }

                    T h0i = m_h0[i];
                    T metric00 = h1j * (df0 + h0i * df1);
                    T metric01 = metric00 * df1;
                    metric00 = metric00 * metric00;

                    T tmp = h1j * (df2 + h0i * df3);
                    metric01.fma(tmp, df3);
                    metric00.fma(tmp, tmp);

                    T metric11 = df1 * df1;
                    metric11.fma(df3, df3);

                    T d0(deriv0 + cnt), d1(deriv1 + cnt);
                    tmp = metric00 * d0;
                    tmp.fma(metric01, d1);
                    tmp.store(bwd + cnt);

                    tmp = metric01 * d0;
                    tmp.fma(metric11, d1);
                    tmp.store(deriv0 + cnt);
                }
            }

            // Step 4b: Take inner products
            AVXIProductTriKernel<NM0, NM1, NQ0, NQ1, VW, CORRECT, false, true, DEFORMED>(
                bwd, this->m_dbdata[0], this->m_bdata[1], this->m_w[0],
                this->m_w[1], jac_ptr, wsp, outptr);

            AVXIProductTriKernel<NM0, NM1, NQ0, NQ1, VW, CORRECT, false, true, DEFORMED>(
                deriv0, this->m_bdata[0], this->m_dbdata[1], this->m_w[0],
                this->m_w[1], jac_ptr, wsp, outptr);

            inptr += nmBlocks;
            outptr += nmBlocks;
        }

    }

private:
    int m_nmTot;
    double m_lambda;
    AlignedVector<VecData<double, VW>> m_h0;
    AlignedVector<VecData<double, VW>> m_h1;
};

template<int VW, bool DEFORMED = false>
struct AVXHelmholtzPrism : public Helmholtz, public AVXHelper<VW, 3, DEFORMED>
{
    AVXHelmholtzPrism(std::vector<LibUtilities::BasisSharedPtr> basis,
                      int nElmt)
    : Helmholtz(basis, nElmt),
      AVXHelper<VW, 3, DEFORMED>(basis, nElmt),
      m_nmTot(LibUtilities::StdPrismData::getNumberOfCoefficients(
                  this->m_nm[0], this->m_nm[1], this->m_nm[2])),
      m_lambda(1.0),
      m_h0(basis[0]->GetNumPoints()),
      m_h1(basis[2]->GetNumPoints())
    {
        const int nq0 = basis[0]->GetNumPoints();
        const int nq2 = basis[2]->GetNumPoints();

        const Array<OneD, const NekDouble> &z0 = basis[0]->GetZ();
        const Array<OneD, const NekDouble> &z2 = basis[2]->GetZ();

        for (int i = 0; i < nq0; ++i)
        {
            m_h0[i] = 0.5 * (1 + z0[i]);
        }

        for (int k = 0; k < nq2; ++k)
        {
            m_h1[k] = 2.0 / (1 - z2[k]);
        }
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<AVXHelmholtzPrism<VW, DEFORMED>>(basis, nElmt);
    }

    static NekDouble FlopsPerElement(
        const int nm,
        const int nq0,
        const int nq1,
        const int nq2)
    {

        const double bwdTrans = AVXBwdTransPrism<VW>::FlopsPerElement(nm,nq0,nq1, nq2);
        const double iprod1 = AVXIProductPrism<VW>::FlopsPerElement(nm, nq0, nq1, nq2);
        const double physDeriv = AVXPhysDerivPrism<VW>::FlopsPerElement(nq0,nq1, nq2);
        const double metric = nq2 * nq1 * nq0 * (9 + 5*9);
        const double iprod2 = AVXIProductPrism<VW>::FlopsPerElement(nm, nq0, nq1, nq2);
        const double iprod3 = AVXIProductPrism<VW>::FlopsPerElement(nm, nq0, nq1, nq2);
        const double iprod4 = AVXIProductPrism<VW>::FlopsPerElement(nm, nq0, nq1, nq2);

        return bwdTrans + iprod1 + physDeriv + metric + iprod2 + iprod3 + iprod4;
    }

    virtual double GFlops() override
    {
        const int nm = m_basis[0]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        const double flops = m_nElmt * AVXHelmholtzPrism::FlopsPerElement(nm, nq0, nq1, nq2);
        return flops * 1e-9;
    }

    virtual NekDouble Ndof() override
    {
        return m_nmTot * this->m_nElmt;
    }

    virtual void operator()(const Array<OneD, const NekDouble> &in,
                                  Array<OneD,       NekDouble> &out) override
    {
        if (m_basis[0]->GetBasisType() == LibUtilities::eModified_A)
        {
            switch(m_basis[0]->GetNumModes())
            {
                case 2:  AVXHelmholtzPrismImpl<2 ,2 ,2 ,3 ,3 ,2 ,true>(in, out); break;
                case 3:  AVXHelmholtzPrismImpl<3 ,3 ,3 ,4 ,4 ,3 ,true>(in, out); break;
                case 4:  AVXHelmholtzPrismImpl<4 ,4 ,4 ,5 ,5 ,4 ,true>(in, out); break;
                case 5:  AVXHelmholtzPrismImpl<5 ,5 ,5 ,6 ,6 ,5 ,true>(in, out); break;
                case 6:  AVXHelmholtzPrismImpl<6 ,6 ,6 ,7 ,7 ,6 ,true>(in, out); break;
                case 7:  AVXHelmholtzPrismImpl<7 ,7 ,7 ,8 ,8 ,7 ,true>(in, out); break;
                case 8:  AVXHelmholtzPrismImpl<8 ,8 ,8 ,9 ,9 ,8 ,true>(in, out); break;
                case 9:  AVXHelmholtzPrismImpl<9 ,9 ,9 ,10,10,9 ,true>(in, out); break;
                case 10: AVXHelmholtzPrismImpl<10,10,10,11,11,10,true>(in, out); break;
                case 11: AVXHelmholtzPrismImpl<11,11,11,12,12,11,true>(in, out); break;
            }
        }
        else
        {
            switch(m_basis[0]->GetNumModes())
            {
                case 2:  AVXHelmholtzPrismImpl<2 ,2 ,2 ,3 ,3 ,2 ,false>(in, out); break;
                case 3:  AVXHelmholtzPrismImpl<3 ,3 ,3 ,4 ,4 ,3 ,false>(in, out); break;
                case 4:  AVXHelmholtzPrismImpl<4 ,4 ,4 ,5 ,5 ,4 ,false>(in, out); break;
                case 5:  AVXHelmholtzPrismImpl<5 ,5 ,5 ,6 ,6 ,5 ,false>(in, out); break;
                case 6:  AVXHelmholtzPrismImpl<6 ,6 ,6 ,7 ,7 ,6 ,false>(in, out); break;
                case 7:  AVXHelmholtzPrismImpl<7 ,7 ,7 ,8 ,8 ,7 ,false>(in, out); break;
                case 8:  AVXHelmholtzPrismImpl<8 ,8 ,8 ,9 ,9 ,8 ,false>(in, out); break;
                case 9:  AVXHelmholtzPrismImpl<9 ,9 ,9 ,10,10,9 ,false>(in, out); break;
                case 10: AVXHelmholtzPrismImpl<10,10,10,11,11,10,false>(in, out); break;
                case 11: AVXHelmholtzPrismImpl<11,11,11,12,12,11,false>(in, out); break;
            }
        }
    }

    template<int NM0, int NM1, int NM2, int NQ0, int NQ1, int NQ2, bool CORRECT>
    void AVXHelmholtzPrismImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &out)
    {
        using T = VecData<double, VW>;

        // Pointers for input/output across blocks
        const NekDouble *inptr = &input[0];
        NekDouble *outptr = &out[0];

        constexpr int nqBlocks = NQ0 * NQ1 * NQ2 * VW;
        const int nmBlocks = m_nmTot * VW;
        constexpr int ndf = 9;
        constexpr int nq = NQ0 * NQ1 * NQ2;

        // Workspace for kernels
        T wsp1[NQ1 * NQ2], wsp2[NQ2], wsp3[NM1];

        // Aligned stack storage
        std::aligned_storage<sizeof(double), 64>::type bwd_storage[nqBlocks];
        std::aligned_storage<sizeof(double), 64>::type deriv0_storage[nqBlocks];
        std::aligned_storage<sizeof(double), 64>::type deriv1_storage[nqBlocks];
        std::aligned_storage<sizeof(double), 64>::type deriv2_storage[nqBlocks];
        double *bwd    = reinterpret_cast<double *>(&bwd_storage[0]);
        double *deriv0 = reinterpret_cast<double *>(&deriv0_storage[0]);
        double *deriv1 = reinterpret_cast<double *>(&deriv1_storage[0]);
        double *deriv2 = reinterpret_cast<double *>(&deriv2_storage[0]);

        for (int e = 0; e < this->m_nBlocks; e++)
        {
            const VecData<double, VW> *jac_ptr;
            if(DEFORMED){
                jac_ptr = &(this->m_jac[e*nq]);
            }
            else{
                jac_ptr = &(this->m_jac[e]);
            }

            T df0, df1, df2, df3, df4, df5, df6, df7, df8;
            if(!DEFORMED){
                df0 = this->m_df[e*ndf];
                df1 = this->m_df[e*ndf + 1];
                df2 = this->m_df[e*ndf + 2];
                df3 = this->m_df[e*ndf + 3];
                df4 = this->m_df[e*ndf + 4];
                df5 = this->m_df[e*ndf + 5];
                df6 = this->m_df[e*ndf + 6];
                df7 = this->m_df[e*ndf + 7];
                df8 = this->m_df[e*ndf + 8];
            }

            // Step 1: BwdTrans
            AVXBwdTransPrismKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, VW, CORRECT>(
                inptr, this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
                wsp1, wsp2, bwd);

            // Step 2: inner product for mass matrix operation
            AVXIProductPrismKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, VW, CORRECT, true, false, DEFORMED>(
                bwd, this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
                this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                wsp1, wsp2, wsp3, outptr, m_lambda);

            // Step 3: take derivatives in standard space
            AVXPhysDerivTensor3DKernel<NQ0, NQ1, NQ2, VW>(
                bwd, this->m_D[0], this->m_D[1], this->m_D[2], deriv0, deriv1,
                deriv2);

            // Step 4: Apply Laplacian metrics & inner product

            // Step 4a: Construct Laplacian metrics
            for (int k = 0, cnt = 0; k < NQ2; ++k)
            {
                T h1 = m_h1[k];
                for (int j = 0; j < NQ1; ++j)
                {
                    for (int i = 0; i < NQ0; ++i, cnt += VW)
                    {
                        T h0 = m_h0[i];

                        if(DEFORMED){
                            int kji = k*NQ0*NQ1 + j*NQ0 + i;
                            VecData<double, VW> *df_start = &(this->m_df[e*ndf*nq]);
                            df0 = df_start[kji*ndf];
                            df1 = df_start[kji*ndf + 1];
                            df2 = df_start[kji*ndf + 2];
                            df3 = df_start[kji*ndf + 3];
                            df4 = df_start[kji*ndf + 4];
                            df5 = df_start[kji*ndf + 5];
                            df6 = df_start[kji*ndf + 6];
                            df7 = df_start[kji*ndf + 7];
                            df8 = df_start[kji*ndf + 8];
                        }

                        T tmp1 = h1 * (h0 * df2 + df0);
                        T tmp2 = h1 * (h0 * df5 + df3);
                        T tmp3 = h1 * (h0 * df8 + df6);

                        T g0 = tmp1 * tmp1;
                        g0.fma(tmp2, tmp2);
                        g0.fma(tmp3, tmp3);

                        T g3 = df1 * tmp1;
                        g3.fma(df4, tmp2);
                        g3.fma(df7, tmp3);

                        T g4 = df2 * tmp1;
                        g4.fma(df5, tmp2);
                        g4.fma(df8, tmp3);

                        T g1 = df1 * df1;
                        g1.fma(df4, df4);
                        g1.fma(df7, df7);

                        T g2 = df2 * df2;
                        g2.fma(df5, df5);
                        g2.fma(df8, df8);

                        T g5 = df1 * df2;
                        g5.fma(df4, df5);
                        g5.fma(df7, df8);

                        // metric00 = g0
                        // metric11 = g1
                        // metric22 = g2
                        // metric01 = g3
                        // metric02 = g4
                        // metric12 = g5

                        T d0(deriv0 + cnt), d1(deriv1 + cnt), d2(deriv2 + cnt);
                        tmp1 = g0 * d0;
                        tmp1.fma(g3, d1);
                        tmp1.fma(g4, d2);

                        tmp2 = g3 * d0;
                        tmp2.fma(g1, d1);
                        tmp2.fma(g5, d2);

                        tmp3 = g4 * d0;
                        tmp3.fma(g5, d1);
                        tmp3.fma(g2, d2);

                        tmp1.store(deriv0 + cnt);
                        tmp2.store(deriv1 + cnt);
                        tmp3.store(deriv2 + cnt);
                    }
                }
            }

            AVXIProductPrismKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, VW, CORRECT, false, true, DEFORMED>(
                deriv0, this->m_dbdata[0], this->m_bdata[1], this->m_bdata[2],
                this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                wsp1, wsp2, wsp3, outptr);

            AVXIProductPrismKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, VW, CORRECT, false, true, DEFORMED>(
                deriv1, this->m_bdata[0], this->m_dbdata[1], this->m_bdata[2],
                this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                wsp1, wsp2, wsp3, outptr);

            AVXIProductPrismKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, VW, CORRECT, false, true, DEFORMED>(
                deriv2, this->m_bdata[0], this->m_bdata[1], this->m_dbdata[2],
                this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                wsp1, wsp2, wsp3, outptr);

            inptr += nmBlocks;
            outptr += nmBlocks;
        }
    }

private:
    int m_nmTot;
    double m_lambda;
    AlignedVector<VecData<double, VW>> m_h0;
    AlignedVector<VecData<double, VW>> m_h1;
};

template<int VW, bool DEFORMED = false>
struct AVXHelmholtzTet : public Helmholtz, public AVXHelper<VW, 3, DEFORMED>
{
    AVXHelmholtzTet(std::vector<LibUtilities::BasisSharedPtr> basis,
                    int nElmt)
    : Helmholtz(basis, nElmt),
      AVXHelper<VW, 3, DEFORMED>(basis, nElmt),
      m_nmTot(LibUtilities::StdTetData::getNumberOfCoefficients(
                  this->m_nm[0], this->m_nm[1], this->m_nm[2])),
      m_lambda(1.0),
      m_h0(basis[0]->GetNumPoints()),
      m_h1(basis[1]->GetNumPoints()),
      m_h2(basis[1]->GetNumPoints()),
      m_h3(basis[2]->GetNumPoints())
    {
        const int nq0 = basis[0]->GetNumPoints();
        const int nq1 = basis[1]->GetNumPoints();
        const int nq2 = basis[2]->GetNumPoints();

        const Array<OneD, const NekDouble> &z0 = basis[0]->GetZ();
        const Array<OneD, const NekDouble> &z1 = basis[1]->GetZ();
        const Array<OneD, const NekDouble> &z2 = basis[2]->GetZ();

        for (int i = 0; i < nq0; ++i)
        {
            m_h0[i] = 0.5 * (1 + z0[i]);
        }

        for (int j = 0; j < nq1; ++j)
        {
            m_h1[j] = 0.5 * (1 + z1[j]);
            m_h2[j] = 2.0 / (1 - z1[j]);
        }

        for (int k = 0; k < nq2; ++k)
        {
            m_h3[k] = 2.0 / (1 - z2[k]);
        }
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<AVXHelmholtzTet<VW, DEFORMED>>(basis, nElmt);
    }

    static NekDouble FlopsPerElement(
        const int nm,
        const int nq0,
        const int nq1,
        const int nq2)
    {

        const double bwdTrans = AVXBwdTransTet<VW>::FlopsPerElement(nm,nq0,nq1, nq2);
        const double iprod1 = AVXIProductTet<VW>::FlopsPerElement(nm, nq0, nq1, nq2);
        const double physDeriv = AVXPhysDerivTet<VW>::FlopsPerElement(nq0,nq1, nq2);
        const double metric = nq2 * nq1 * (2 + nq0 * (1 + 12 + 9 + 5*9));
        const double iprod2 = AVXIProductTet<VW>::FlopsPerElement(nm, nq0, nq1, nq2);
        const double iprod3 = AVXIProductTet<VW>::FlopsPerElement(nm, nq0, nq1, nq2);
        const double iprod4 = AVXIProductTet<VW>::FlopsPerElement(nm, nq0, nq1, nq2);

        return bwdTrans + iprod1 + physDeriv + metric + iprod2 + iprod3 + iprod4;
    }

    virtual double GFlops() override
    {
        const int nm = m_basis[0]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        const double flops = m_nElmt * AVXHelmholtzTet::FlopsPerElement(nm, nq0, nq1, nq2);
        return flops * 1e-9;
    }

    virtual NekDouble Ndof() override
    {
        return m_nmTot * this->m_nElmt;
    }

    virtual void operator()(const Array<OneD, const NekDouble> &in,
                                  Array<OneD,       NekDouble> &out) override
    {
        if (m_basis[0]->GetBasisType() == LibUtilities::eModified_A)
        {
            switch(m_basis[0]->GetNumModes())
            {
                case 2:  AVXHelmholtzTetImpl<2 ,2 ,2 ,3 ,2 ,2 ,true>(in, out); break;
                case 3:  AVXHelmholtzTetImpl<3 ,3 ,3 ,4 ,3 ,3 ,true>(in, out); break;
                case 4:  AVXHelmholtzTetImpl<4 ,4 ,4 ,5 ,4 ,4 ,true>(in, out); break;
                case 5:  AVXHelmholtzTetImpl<5 ,5 ,5 ,6 ,5 ,5 ,true>(in, out); break;
                case 6:  AVXHelmholtzTetImpl<6 ,6 ,6 ,7 ,6 ,6 ,true>(in, out); break;
                case 7:  AVXHelmholtzTetImpl<7 ,7 ,7 ,8 ,7 ,7 ,true>(in, out); break;
                case 8:  AVXHelmholtzTetImpl<8 ,8 ,8 ,9 ,8 ,8 ,true>(in, out); break;
                case 9:  AVXHelmholtzTetImpl<9 ,9 ,9 ,10,9 ,9 ,true>(in, out); break;
                case 10: AVXHelmholtzTetImpl<10,10,10,11,10,10,true>(in, out); break;
                case 11: AVXHelmholtzTetImpl<11,11,11,12,11,11,true>(in, out); break;
            }
        }
        else
        {
            switch(m_basis[0]->GetNumModes())
            {
                case 2:  AVXHelmholtzTetImpl<2 ,2 ,2 ,3 ,2 ,2 ,false>(in, out); break;
                case 3:  AVXHelmholtzTetImpl<3 ,3 ,3 ,4 ,3 ,3 ,false>(in, out); break;
                case 4:  AVXHelmholtzTetImpl<4 ,4 ,4 ,5 ,4 ,4 ,false>(in, out); break;
                case 5:  AVXHelmholtzTetImpl<5 ,5 ,5 ,6 ,5 ,5 ,false>(in, out); break;
                case 6:  AVXHelmholtzTetImpl<6 ,6 ,6 ,7 ,6 ,6 ,false>(in, out); break;
                case 7:  AVXHelmholtzTetImpl<7 ,7 ,7 ,8 ,7 ,7 ,false>(in, out); break;
                case 8:  AVXHelmholtzTetImpl<8 ,8 ,8 ,9 ,8 ,8 ,false>(in, out); break;
                case 9:  AVXHelmholtzTetImpl<9 ,9 ,9 ,10,9 ,9 ,false>(in, out); break;
                case 10: AVXHelmholtzTetImpl<10,10,10,11,10,10,false>(in, out); break;
                case 11: AVXHelmholtzTetImpl<11,11,11,12,11,11,false>(in, out); break;
            }
        }
    }

    template<int NM0, int NM1, int NM2, int NQ0, int NQ1, int NQ2, bool CORRECT>
    void AVXHelmholtzTetImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &out)
    {
        using T = VecData<double, VW>;

        // Pointers for input/output across blocks
        const NekDouble *inptr = &input[0];
        NekDouble *outptr = &out[0];

        constexpr int nqBlocks = NQ0 * NQ1 * NQ2 * VW;
        constexpr int ndf = 9;
        constexpr int nq = NQ0 * NQ1 * NQ2;
        const int nmBlocks = m_nmTot * VW;

        // Workspace for kernels
        T wsp1[NQ1 * NQ2 + NQ2], wsp2[NM0];

        // Aligned stack storage
        std::aligned_storage<sizeof(double), 64>::type bwd_storage[nqBlocks];
        std::aligned_storage<sizeof(double), 64>::type deriv0_storage[nqBlocks];
        std::aligned_storage<sizeof(double), 64>::type deriv1_storage[nqBlocks];
        std::aligned_storage<sizeof(double), 64>::type deriv2_storage[nqBlocks];
        double *bwd    = reinterpret_cast<double *>(&bwd_storage[0]);
        double *deriv0 = reinterpret_cast<double *>(&deriv0_storage[0]);
        double *deriv1 = reinterpret_cast<double *>(&deriv1_storage[0]);
        double *deriv2 = reinterpret_cast<double *>(&deriv2_storage[0]);

        for (int e = 0; e < this->m_nBlocks; e++)
        {
            const VecData<double, VW> *jac_ptr;
            if(DEFORMED){
                jac_ptr = &(this->m_jac[e*nq]);
            }
            else{
                jac_ptr = &(this->m_jac[e]);
            }

            T df0, df1, df2, df3, df4, df5, df6, df7, df8;
            if(!DEFORMED){
                df0 = this->m_df[e*ndf];
                df1 = this->m_df[e*ndf + 1];
                df2 = this->m_df[e*ndf + 2];
                df3 = this->m_df[e*ndf + 3];
                df4 = this->m_df[e*ndf + 4];
                df5 = this->m_df[e*ndf + 5];
                df6 = this->m_df[e*ndf + 6];
                df7 = this->m_df[e*ndf + 7];
                df8 = this->m_df[e*ndf + 8];
            }

            // Step 1: BwdTrans
            AVXBwdTransTetKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, VW, CORRECT>(
                inptr, this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
                wsp1, wsp2, bwd);

            // Step 2: inner product for mass matrix operation
            AVXIProductTetKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, VW, CORRECT, true, false, DEFORMED>(
                bwd, this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
                this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                wsp1, outptr, m_lambda);

            // Step 3: take derivatives in standard space
            AVXPhysDerivTensor3DKernel<NQ0, NQ1, NQ2, VW>(
                bwd, this->m_D[0], this->m_D[1], this->m_D[2], deriv0, deriv1,
                deriv2);

            // Step 4: Apply Laplacian metrics & inner product

            // Step 4a: Construct Laplacian metrics
            for (int k = 0, cnt = 0; k < NQ2; ++k)
            {
                T h3 = m_h3[k];
                for (int j = 0; j < NQ1; ++j)
                {
                    T h1 = m_h1[j];
                    T h2 = m_h2[j];
                    T h2h3 = h2 * h3;
                    T h1h3 = h1 * h3;

                    for (int i = 0; i < NQ0; ++i, cnt += VW)
                    {
                        if(DEFORMED){
                            int kji = k*NQ0*NQ1 + j*NQ0 + i;
                            VecData<double, VW> *df_start = &(this->m_df[e*ndf*nq]);
                            df0 = df_start[kji*ndf];
                            df1 = df_start[kji*ndf + 1];
                            df2 = df_start[kji*ndf + 2];
                            df3 = df_start[kji*ndf + 3];
                            df4 = df_start[kji*ndf + 4];
                            df5 = df_start[kji*ndf + 5];
                            df6 = df_start[kji*ndf + 6];
                            df7 = df_start[kji*ndf + 7];
                            df8 = df_start[kji*ndf + 8];
                        }

                        T h0h2h3 = m_h0[i] * h2h3;

                        T tmp1 = h0h2h3 * (df1 + df2);
                        tmp1.fma(df0, h2h3);
                        T tmp2 = h0h2h3 * (df4 + df5);
                        tmp2.fma(df3, h2h3);
                        T tmp3 = h0h2h3 * (df7 + df8);
                        tmp3.fma(df6, h2h3);

                        T g0 = tmp1 * tmp1;
                        g0.fma(tmp2, tmp2);
                        g0.fma(tmp3, tmp3);

                        T g4 = df2 * tmp1;
                        g4.fma(df5, tmp2);
                        g4.fma(df8, tmp3);

                        T tmp4 = df1 * h3;
                        tmp4.fma(df2, h1h3);
                        T tmp5 = df4 * h3;
                        tmp5.fma(df5, h1h3);
                        T tmp6 = df7 * h3;
                        tmp6.fma(df8, h1h3);

                        T g3 = tmp1 * tmp4;
                        g3.fma(tmp2, tmp5);
                        g3.fma(tmp3, tmp6);

                        T g1 = tmp4 * tmp4;
                        g1.fma(tmp5, tmp5);
                        g1.fma(tmp6, tmp6);

                        T g5 = df2 * tmp4;
                        g5.fma(df5, tmp5);
                        g5.fma(df8, tmp6);

                        T g2 = df2 * df2;
                        g2.fma(df5, df5);
                        g2.fma(df8, df8);

                        // metric00 = g0
                        // metric11 = g1
                        // metric22 = g2
                        // metric01 = g3
                        // metric02 = g4
                        // metric12 = g5

                        T d0(deriv0 + cnt), d1(deriv1 + cnt), d2(deriv2 + cnt);

                        tmp1 = g0 * d0;
                        tmp1.fma(g3, d1);
                        tmp1.fma(g4, d2);

                        tmp2 = g3 * d0;
                        tmp2.fma(g1, d1);
                        tmp2.fma(g5, d2);

                        tmp3 = g4 * d0;
                        tmp3.fma(g5, d1);
                        tmp3.fma(g2, d2);

                        tmp1.store(deriv0 + cnt);
                        tmp2.store(deriv1 + cnt);
                        tmp3.store(deriv2 + cnt);
                    }
                }
            }

            AVXIProductTetKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, VW, CORRECT, false, true, DEFORMED>(
                deriv0, this->m_dbdata[0], this->m_bdata[1], this->m_bdata[2],
                this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                wsp1, outptr);

            AVXIProductTetKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, VW, CORRECT, false, true, DEFORMED>(
                deriv1, this->m_bdata[0], this->m_dbdata[1], this->m_bdata[2],
                this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                wsp1, outptr);

            AVXIProductTetKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, VW, CORRECT, false, true, DEFORMED>(
                deriv2, this->m_bdata[0], this->m_bdata[1], this->m_dbdata[2],
                this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                wsp1, outptr);

            inptr += nmBlocks;
            outptr += nmBlocks;
        }
    }

private:
    int m_nmTot;
    double m_lambda;
    AlignedVector<VecData<double, VW>> m_h0, m_h1, m_h2, m_h3;
};

template<int VW, bool DEFORMED = false>
struct AVXHelmholtzHex : public Helmholtz, public AVXHelper<VW, 3, DEFORMED>
{
    AVXHelmholtzHex(std::vector<LibUtilities::BasisSharedPtr> basis,
                    int nElmt)
    : Helmholtz(basis, nElmt),
      AVXHelper<VW, 3, DEFORMED>(basis, nElmt),
      m_nmTot(LibUtilities::StdHexData::getNumberOfCoefficients(
                  this->m_nm[0], this->m_nm[1], this->m_nm[2])),
      m_lambda(1.0)
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<AVXHelmholtzHex<VW, DEFORMED>>(basis, nElmt);
    }

    static NekDouble FlopsPerElement(
        const int nm,
        const int nq0,
        const int nq1,
        const int nq2)
    {

        const double bwdTrans = AVXBwdTransHex<VW>::FlopsPerElement(nm,nq0,nq1, nq2);
        const double iprod1 = AVXIProductHex<VW>::FlopsPerElement(nm, nq0, nq1, nq2);
        const double physDeriv = AVXPhysDerivHex<VW>::FlopsPerElement(nq0,nq1, nq2);

        const double metric_calc = 30;
        double metric;
        if(DEFORMED){
            metric = nq2 *nq1 * nq0 * (metric_calc + 15);
        }
        else{
            metric = (nq2 * nq1 * nq0 * 15) + metric_calc;
        }

        const double iprod2 = AVXIProductHex<VW>::FlopsPerElement(nm, nq0, nq1, nq2);
        const double iprod3 = AVXIProductHex<VW>::FlopsPerElement(nm, nq0, nq1, nq2);
        const double iprod4 = AVXIProductHex<VW>::FlopsPerElement(nm, nq0, nq1, nq2);

        return bwdTrans + iprod1 + physDeriv + metric + iprod2 + iprod3 + iprod4;
    }

    virtual double GFlops() override
    {
        const int nm = m_basis[0]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        const double flops = m_nElmt * AVXHelmholtzHex::FlopsPerElement(nm, nq0, nq1, nq2);
        return flops * 1e-9;
    }

    virtual NekDouble Ndof() override
    {
        return m_nmTot * this->m_nElmt;
    }

    virtual void operator()(const Array<OneD, const NekDouble> &in,
                                  Array<OneD,       NekDouble> &out) override
    {
        switch(m_basis[0]->GetNumModes())
        {
            case 2:  AVXHelmholtzHexImpl<2 ,2 ,2 ,3 ,3 ,3 >(in, out); break;
            case 3:  AVXHelmholtzHexImpl<3 ,3 ,3 ,4 ,4 ,4 >(in, out); break;
            case 4:  AVXHelmholtzHexImpl<4 ,4 ,4 ,5 ,5 ,5 >(in, out); break;
            case 5:  AVXHelmholtzHexImpl<5 ,5 ,5 ,6 ,6 ,6 >(in, out); break;
            case 6:  AVXHelmholtzHexImpl<6 ,6 ,6 ,7 ,7 ,7 >(in, out); break;
            case 7:  AVXHelmholtzHexImpl<7 ,7 ,7 ,8 ,8 ,8 >(in, out); break;
            case 8:  AVXHelmholtzHexImpl<8 ,8 ,8 ,9 ,9 ,9 >(in, out); break;
            case 9:  AVXHelmholtzHexImpl<9 ,9 ,9 ,10,10,10>(in, out); break;
            case 10: AVXHelmholtzHexImpl<10,10,10,11,11,11>(in, out); break;
            case 11: AVXHelmholtzHexImpl<11,11,11,12,12,12>(in, out); break;
        }
    }

#if 0
    template<int NM0, int NM1, int NM2, int NQ0, int NQ1, int NQ2>
    void AVXHelmholtzHexImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &out)
    {
        using T = VecData<double, VW>;

        // Pointers for input/output across blocks
        const NekDouble *inptr = &input[0];
        NekDouble *outptr = &out[0];

        constexpr int nqBlocks = NQ0 * NQ1 * NQ2 * VW;
        constexpr int nmBlocks = NM0 * NM1 * NM2 * VW;
        constexpr int ndf = 9;
        constexpr int nq = NQ0 * NQ1 * NQ2;

        // Workspace for kernels
        //T wsp1_0[NM2 * NQ0 * NM1];//, wsp1_1[NM2 * NQ0 * NM1];
        //T wsp2_0[NM2 * NQ1 * NQ0];//, wsp2_1[NM2 * NQ1 * NQ0], wsp2_2[NM2 * NQ1 * NQ0];
        T wsp1[NQ0 * NQ1 * NQ2], wsp2[NQ0 * NQ1 * NQ2];

        // Aligned stack storage
        std::aligned_storage<sizeof(double), 64>::type bwd_storage[nqBlocks];
        std::aligned_storage<sizeof(double), 64>::type out_storage[nmBlocks];
        std::aligned_storage<sizeof(double), 64>::type deriv0_storage[nqBlocks];
        std::aligned_storage<sizeof(double), 64>::type deriv1_storage[nqBlocks];
        std::aligned_storage<sizeof(double), 64>::type deriv2_storage[nqBlocks];
        double *bwd    = reinterpret_cast<double *>(&bwd_storage[0]);
        double *outtmp = reinterpret_cast<double *>(&out_storage[0]);
        double *deriv0 = reinterpret_cast<double *>(&deriv0_storage[0]);
        double *deriv1 = reinterpret_cast<double *>(&deriv1_storage[0]);
        double *deriv2 = reinterpret_cast<double *>(&deriv2_storage[0]);

        VecData<double, VW> *jac_ptr = &this->m_jac[0];
        VecData<double, VW> *df_ptr = &this->m_df[0];

        // Precompute Laplacian metrics
        T metric00, metric01, metric02, metric11, metric12, metric22;
        metric00 = df_ptr[0] * df_ptr[0];
        metric00.fma(df_ptr[3], df_ptr[3]);
        metric00.fma(df_ptr[6], df_ptr[6]);

        metric01 = df_ptr[0] * df_ptr[1];
        metric01.fma(df_ptr[3], df_ptr[4]);
        metric01.fma(df_ptr[6], df_ptr[7]);

        metric02 = df_ptr[0] * df_ptr[2];
        metric02.fma(df_ptr[3], df_ptr[5]);
        metric02.fma(df_ptr[6], df_ptr[8]);

        metric11 = df_ptr[1] * df_ptr[1];
        metric11.fma(df_ptr[4], df_ptr[4]);
        metric11.fma(df_ptr[7], df_ptr[7]);

        metric12 = df_ptr[1] * df_ptr[2];
        metric12.fma(df_ptr[4], df_ptr[5]);
        metric12.fma(df_ptr[7], df_ptr[8]);

        metric22 = df_ptr[2] * df_ptr[2];
        metric22.fma(df_ptr[5], df_ptr[5]);
        metric22.fma(df_ptr[8], df_ptr[8]);

        for (int e = 0; e < this->m_nBlocks; e++)
        {
            // x-gradient
            TensorContractDealII<VW, 3, 0, NM0, NQ0, true, false>(this->m_bdata[0], (T *)inptr, wsp1);
            TensorContractDealII<VW, 3, 1, NM0, NQ0, true, false>(this->m_dbdata[0], wsp1, wsp2);
            TensorContractDealII<VW, 3, 2, NM0, NQ0, true, false>(this->m_dbdata[0], wsp2, (T *)deriv0);

            // y-gradient
            TensorContractDealII<VW, 3, 0, NM0, NQ0, true, false>(this->m_bdata[0], (T *)inptr, wsp1);
            TensorContractDealII<VW, 3, 1, NM0, NQ0, true, false>(this->m_dbdata[0], wsp1, wsp2);
            TensorContractDealII<VW, 3, 2, NM0, NQ0, true, false>(this->m_bdata[0], wsp2, (T *)deriv1);

            // z-gradient
            TensorContractDealII<VW, 3, 1, NM0, NQ0, true, false>(this->m_bdata[0], wsp1, wsp2);
            TensorContractDealII<VW, 3, 2, NM0, NQ0, true, false>(this->m_dbdata[0], wsp2, (T *)deriv2);

            // values
            TensorContractDealII<VW, 3, 2, NM0, NQ0, true, false>(this->m_bdata[0], wsp2, (T *)bwd);

            // Step 4: Apply Laplacian metrics & inner product
            if (DEFORMED)
            {
                int cnt = 0;
                for (int k = 0; k < NQ2; k++)
                {
                    for (int j = 0; j < NQ1; j++)
                    {
                        for (int i = 0; i < NQ0; i++, cnt += VW, df_ptr += 9)
                        {
                            metric00 = df_ptr[0] * df_ptr[0];
                            metric00.fma(df_ptr[3], df_ptr[3]);
                            metric00.fma(df_ptr[6], df_ptr[6]);

                            metric01 = df_ptr[0] * df_ptr[1];
                            metric01.fma(df_ptr[3], df_ptr[4]);
                            metric01.fma(df_ptr[6], df_ptr[7]);

                            metric02 = df_ptr[0] * df_ptr[2];
                            metric02.fma(df_ptr[3], df_ptr[5]);
                            metric02.fma(df_ptr[6], df_ptr[8]);

                            metric11 = df_ptr[1] * df_ptr[1];
                            metric11.fma(df_ptr[4], df_ptr[4]);
                            metric11.fma(df_ptr[7], df_ptr[7]);

                            metric12 = df_ptr[1] * df_ptr[2];
                            metric12.fma(df_ptr[4], df_ptr[5]);
                            metric12.fma(df_ptr[7], df_ptr[8]);

                            metric22 = df_ptr[2] * df_ptr[2];
                            metric22.fma(df_ptr[5], df_ptr[5]);
                            metric22.fma(df_ptr[8], df_ptr[8]);

                            T d0 = T(deriv0 + cnt), d1 = T(deriv1 + cnt), d2 = T(deriv2 + cnt);

                            T tmp = metric00 * d0;
                            tmp.fma(metric01, d1);
                            tmp.fma(metric02, d2);
                            tmp.store(deriv0 + cnt);

                            tmp = metric01 * d0;
                            tmp.fma(metric11, d1);
                            tmp.fma(metric12, d2);
                            tmp.store(deriv1 + cnt);

                            tmp = metric02 * d0;
                            tmp.fma(metric12, d1);
                            tmp.fma(metric22, d2);
                            tmp.store(deriv2 + cnt);
                        }
                    }
                }

                AVXIProductHexKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, VW, false, true, true>(
                    deriv0, this->m_dbdata[0], this->m_bdata[1], this->m_bdata[2],
                    this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                    wsp1, wsp2, outtmp);

                AVXIProductHexKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, VW, false, true, true>(
                    deriv1, this->m_bdata[0], this->m_dbdata[1], this->m_bdata[2],
                    this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                    wsp1, wsp2, outtmp);

                AVXIProductHexKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, VW, false, true, true>(
                    deriv2, this->m_bdata[0], this->m_bdata[1], this->m_dbdata[2],
                    this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                    wsp1, wsp2, outtmp);

                // Increment Jacobian pointer
                jac_ptr += nq;
            }
            else
            {
                // Direction 1
                for (int i = 0; i < nqBlocks; i += VW)
                {
                    T deriv0i(deriv0 + i);
                    T deriv1i(deriv1 + i);
                    T deriv2i(deriv2 + i);

                    T *d0new = (T *)(deriv0 + i);
                    T *d1new = (T *)(deriv1 + i);
                    T *d2new = (T *)(deriv2 + i);

                    *d0new = metric00 * deriv0i;
                    d0new->fma(metric01, deriv1i);
                    d0new->fma(metric02, deriv2i);

                    *d1new = metric01 * deriv0i;
                    d1new->fma(metric11, deriv1i);
                    d1new->fma(metric12, deriv2i);

                    *d2new = metric02 * deriv0i;
                    d2new->fma(metric12, deriv1i);
                    d2new->fma(metric22, deriv2i);
                }

                TensorContractDealII<VW, 3, 2, NM0, NQ0, false, false>(this->m_dbdata[0], (T *)deriv2, wsp1);
                TensorContractDealII<VW, 3, 2, NM0, NQ0, false, true>(this->m_bdata[0], (T *)bwd, wsp1);
                TensorContractDealII<VW, 3, 1, NM0, NQ0, false, false>(this->m_bdata[0], wsp1, wsp2);
                TensorContractDealII<VW, 3, 2, NM0, NQ0, false, false>(this->m_bdata[0], (T *)deriv1, wsp1);

                TensorContractDealII<VW, 3, 1, NM0, NQ0, false, true>(this->m_dbdata[0], wsp1, wsp2);
                TensorContractDealII<VW, 3, 0, NM0, NQ0, false, true>(this->m_bdata[0], wsp2, (T *)outptr);

                TensorContractDealII<VW, 3, 2, NM0, NQ0, false, false>(this->m_bdata[0], (T *)deriv0, wsp1);
                TensorContractDealII<VW, 3, 1, NM0, NQ0, false, false>(this->m_bdata[0], wsp1, wsp2);
                TensorContractDealII<VW, 3, 0, NM0, NQ0, false, true>(this->m_dbdata[0], wsp1, (T *)outptr);

                /*
                AVXIProductHexKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, VW, false, true, false>(
                    bwd, this->m_dbdata[0], this->m_bdata[0], this->m_bdata[0],
                    this->m_w[0], this->m_w[0], this->m_w[0], jac_ptr,
                    wsp1, wsp2, outtmp);

                AVXIProductHexKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, VW, false, true, false>(
                    bwd, this->m_bdata[0], this->m_dbdata[0], this->m_bdata[0],
                    this->m_w[0], this->m_w[0], this->m_w[0], jac_ptr,
                    wsp1, wsp2, outtmp);

                AVXIProductHexKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, VW, false, true, false>(
                    bwd, this->m_bdata[0], this->m_bdata[0], this->m_dbdata[0],
                    this->m_w[0], this->m_w[0], this->m_w[0], jac_ptr,
                    wsp1, wsp2, outtmp);

                  /*
                df_ptr += 9;
                jac_ptr++;
                */
            }

            for (int i = 0; i < NM0 * NM1 * NM2; ++i)
            {
                T outxmm(outtmp + i * VW);
                outxmm.store(outptr + i * VW);
            }

            inptr += nmBlocks;
            outptr += nmBlocks;
        }
    }
#else
    template<int NM0, int NM1, int NM2, int NQ0, int NQ1, int NQ2>
    void AVXHelmholtzHexImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &out)
    {
        using T = VecData<double, VW>;

        // Pointers for input/output across blocks
        const NekDouble *inptr = &input[0];
        NekDouble *outptr = &out[0];

        constexpr int nqBlocks = NQ0 * NQ1 * NQ2 * VW;
        constexpr int nmBlocks = NM0 * NM1 * NM2 * VW;
        constexpr int ndf = 9;
        constexpr int nq = NQ0 * NQ1 * NQ2;

        // Workspace for kernels
        T wsp1[NQ0 * NQ1 * NQ2], wsp2[NQ0 * NQ1 * NQ2];

        // Aligned stack storage
        std::aligned_storage<sizeof(double), 64>::type bwd_storage[nqBlocks];
        std::aligned_storage<sizeof(double), 64>::type out_storage[nmBlocks];
        std::aligned_storage<sizeof(double), 64>::type deriv0_storage[nqBlocks];
        std::aligned_storage<sizeof(double), 64>::type deriv1_storage[nqBlocks];
        std::aligned_storage<sizeof(double), 64>::type deriv2_storage[nqBlocks];
        double *bwd    = reinterpret_cast<double *>(&bwd_storage[0]);
        double *outtmp = reinterpret_cast<double *>(&out_storage[0]);
        double *deriv0 = reinterpret_cast<double *>(&deriv0_storage[0]);
        double *deriv1 = reinterpret_cast<double *>(&deriv1_storage[0]);
        double *deriv2 = reinterpret_cast<double *>(&deriv2_storage[0]);

        VecData<double, VW> *jac_ptr = &this->m_jac[0];
        VecData<double, VW> *df_ptr = &this->m_df[0];

        T metric00, metric01, metric02, metric11, metric12, metric22;

        for (int e = 0; e < this->m_nBlocks; e++)
        {
            // Precompute Laplacian metrics

            // Step 1: BwdTrans
            AVXBwdTransHexKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, VW>(
               inptr, this->m_bdata[0], this->m_bdata[0], this->m_bdata[0],
               wsp1, wsp2, bwd);

            // Step 2: inner product for mass matrix operation
            AVXIProductHexKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, VW, true, false, DEFORMED>(
                bwd, this->m_bdata[0], this->m_bdata[0], this->m_bdata[0],
                this->m_w[0], this->m_w[0], this->m_w[0], jac_ptr,
                wsp1, wsp2, outtmp, m_lambda);

            // Step 3: take derivatives in standard space
            AVXPhysDerivTensor3DKernel<NQ0, NQ1, NQ2, VW>(
                bwd, this->m_D[0], this->m_D[0], this->m_D[0], deriv0, deriv1,
                deriv2);

            // Step 4: Apply Laplacian metrics & inner product
            if (DEFORMED)
            {
                int cnt = 0;
                for (int k = 0; k < NQ2; k++)
                {
                    for (int j = 0; j < NQ1; j++)
                    {
                        for (int i = 0; i < NQ0; i++, cnt += VW, df_ptr += 9)
                        {
                            metric00 = df_ptr[0] * df_ptr[0];
                            metric00.fma(df_ptr[3], df_ptr[3]);
                            metric00.fma(df_ptr[6], df_ptr[6]);

                            metric01 = df_ptr[0] * df_ptr[1];
                            metric01.fma(df_ptr[3], df_ptr[4]);
                            metric01.fma(df_ptr[6], df_ptr[7]);

                            metric02 = df_ptr[0] * df_ptr[2];
                            metric02.fma(df_ptr[3], df_ptr[5]);
                            metric02.fma(df_ptr[6], df_ptr[8]);

                            metric11 = df_ptr[1] * df_ptr[1];
                            metric11.fma(df_ptr[4], df_ptr[4]);
                            metric11.fma(df_ptr[7], df_ptr[7]);

                            metric12 = df_ptr[1] * df_ptr[2];
                            metric12.fma(df_ptr[4], df_ptr[5]);
                            metric12.fma(df_ptr[7], df_ptr[8]);

                            metric22 = df_ptr[2] * df_ptr[2];
                            metric22.fma(df_ptr[5], df_ptr[5]);
                            metric22.fma(df_ptr[8], df_ptr[8]);

                            T d0 = T(deriv0 + cnt), d1 = T(deriv1 + cnt), d2 = T(deriv2 + cnt);

                            T tmp = metric00 * d0;
                            tmp.fma(metric01, d1);
                            tmp.fma(metric02, d2);
                            tmp.store(deriv0 + cnt);

                            tmp = metric01 * d0;
                            tmp.fma(metric11, d1);
                            tmp.fma(metric12, d2);
                            tmp.store(deriv1 + cnt);

                            tmp = metric02 * d0;
                            tmp.fma(metric12, d1);
                            tmp.fma(metric22, d2);
                            tmp.store(deriv2 + cnt);
                        }
                    }
                }

                AVXIProductHexKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, VW, false, true, true>(
                    deriv0, this->m_dbdata[0], this->m_bdata[1], this->m_bdata[2],
                    this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                    wsp1, wsp2, outtmp);

                AVXIProductHexKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, VW, false, true, true>(
                    deriv1, this->m_bdata[0], this->m_dbdata[1], this->m_bdata[2],
                    this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                    wsp1, wsp2, outtmp);

                AVXIProductHexKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, VW, false, true, true>(
                    deriv2, this->m_bdata[0], this->m_bdata[1], this->m_dbdata[2],
                    this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                    wsp1, wsp2, outtmp);

                // Increment Jacobian pointer
                jac_ptr += nq;
            }
            else
            {
                metric00 = df_ptr[0] * df_ptr[0];
                metric00.fma(df_ptr[3], df_ptr[3]);
                metric00.fma(df_ptr[6], df_ptr[6]);

                metric01 = df_ptr[0] * df_ptr[1];
                metric01.fma(df_ptr[3], df_ptr[4]);
                metric01.fma(df_ptr[6], df_ptr[7]);

                metric02 = df_ptr[0] * df_ptr[2];
                metric02.fma(df_ptr[3], df_ptr[5]);
                metric02.fma(df_ptr[6], df_ptr[8]);

                metric11 = df_ptr[1] * df_ptr[1];
                metric11.fma(df_ptr[4], df_ptr[4]);
                metric11.fma(df_ptr[7], df_ptr[7]);

                metric12 = df_ptr[1] * df_ptr[2];
                metric12.fma(df_ptr[4], df_ptr[5]);
                metric12.fma(df_ptr[7], df_ptr[8]);

                metric22 = df_ptr[2] * df_ptr[2];
                metric22.fma(df_ptr[5], df_ptr[5]);
                metric22.fma(df_ptr[8], df_ptr[8]);

                for (int i = 0; i < nqBlocks; i += VW)
                {
                    T deriv0i(deriv0 + i);
                    T deriv1i(deriv1 + i);
                    T deriv2i(deriv2 + i);

                    T *d0new = (T *)(deriv0 + i);
                    T *d1new = (T *)(deriv1 + i);
                    T *d2new = (T *)(deriv2 + i);

                    *d0new = metric00 * deriv0i;
                    d0new->fma(metric01, deriv1i);
                    d0new->fma(metric02, deriv2i);

                    *d1new = metric01 * deriv0i;
                    d1new->fma(metric11, deriv1i);
                    d1new->fma(metric12, deriv2i);

                    *d2new = metric02 * deriv0i;
                    d2new->fma(metric12, deriv1i);
                    d2new->fma(metric22, deriv2i);
                }

                AVXIProductHexKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, VW, false, true, false>(
                    deriv0, this->m_dbdata[0], this->m_bdata[0], this->m_bdata[0],
                    this->m_w[0], this->m_w[0], this->m_w[0], jac_ptr,
                    wsp1, wsp2, outtmp);

                AVXIProductHexKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, VW, false, true, false>(
                    deriv1, this->m_bdata[0], this->m_dbdata[0], this->m_bdata[0],
                    this->m_w[0], this->m_w[0], this->m_w[0], jac_ptr,
                    wsp1, wsp2, outtmp);

                AVXIProductHexKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, VW, false, true, false>(
                    deriv2, this->m_bdata[0], this->m_bdata[0], this->m_dbdata[0],
                    this->m_w[0], this->m_w[0], this->m_w[0], jac_ptr,
                    wsp1, wsp2, outtmp);

                df_ptr += 9;
                jac_ptr++;
            }

            for (int i = 0; i < NM0 * NM1 * NM2; ++i)
            {
                T outxmm(outtmp + i * VW);
                outxmm.store(outptr + i * VW);
            }

            inptr += nmBlocks;
            outptr += nmBlocks;
        }
    }
#endif

private:
    int m_nmTot;
    double m_lambda;
};

#endif
