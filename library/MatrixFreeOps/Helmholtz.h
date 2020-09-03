#ifndef NEKTAR_LIBRARY_MF_HELMHOLTZ_H
#define NEKTAR_LIBRARY_MF_HELMHOLTZ_H

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ShapeType.hpp>
#include <LibUtilities/Foundations/Basis.h>

#include "IProduct.h"
#include "IProductKernels.hpp"
#include "BwdTrans.h"
#include "PhysDeriv.h"

namespace Nektar
{
namespace MatrixFree
{

template<bool DEFORMED = false>
struct HelmholtzQuad : public Helmholtz, public Helper<2, DEFORMED>
{
    HelmholtzQuad(std::vector<LibUtilities::BasisSharedPtr> basis,
                     int nElmt)
    : Helmholtz(basis, nElmt),
      Helper<2, DEFORMED>(basis, nElmt),
      m_nmTot(LibUtilities::StdQuadData::getNumberOfCoefficients(
                                 this->m_nm[0], this->m_nm[1]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<HelmholtzQuad<DEFORMED>>(basis, nElmt);
    }

    virtual void operator()(const Array<OneD, const NekDouble> &in,
                                  Array<OneD,       NekDouble> &out) override
    {
        // Check preconditions
        ASSERTL0(m_basis[0]->GetNumModes() == m_basis[1]->GetNumModes() &&
            m_basis[0]->GetNumPoints() == m_basis[1]->GetNumPoints(),
            "MatrixFree requires homogenous modes/points");

        switch(m_basis[0]->GetNumModes())
        {
        case 2:  
            switch(m_basis[0]->GetNumPoints())
            {
            case 2: HelmholtzQuadImpl<2 ,2 ,2 ,2 >(in, out); break;
            case 3: HelmholtzQuadImpl<2 ,2 ,3 ,3 >(in, out); break;
            case 4: HelmholtzQuadImpl<2 ,2 ,4 ,4 >(in, out); break;
            default: NEKERROR(ErrorUtil::efatal,
                              "HelmholtzQuad: # of modes / points combo not implemented.");
            } break;
        case 3:  HelmholtzQuadImpl<3 ,3 ,4 ,4 >(in, out); break;
        case 4:  HelmholtzQuadImpl<4 ,4 ,5 ,5 >(in, out); break;
        case 5:  HelmholtzQuadImpl<5 ,5 ,6 ,6 >(in, out); break;
        case 6:  HelmholtzQuadImpl<6 ,6 ,7 ,7 >(in, out); break;
        case 7:  HelmholtzQuadImpl<7 ,7 ,8 ,8 >(in, out); break;
        case 8:  HelmholtzQuadImpl<8 ,8 ,9 ,9 >(in, out); break;
        case 9:  HelmholtzQuadImpl<9 ,9 ,10,10>(in, out); break;
        case 10: HelmholtzQuadImpl<10,10,11,11>(in, out); break;
        case 11: HelmholtzQuadImpl<11,11,12,12>(in, out); break;
        }
    }
    
    template<int NM0, int NM1, int NQ0, int NQ1>
    void HelmholtzQuadImpl(
          const Array<OneD, const NekDouble> &input, Array<OneD, NekDouble> &out)
    {
        auto *inptr  = &input[0];
        auto *outptr = &out[0];

        constexpr auto ndf = 4;
        constexpr auto nqTot = NQ0 * NQ1;
        //constexpr auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        // Allocate sufficient workspace for backwards transform and inner
        // product kernels.
        constexpr auto wspInnerProd = NQ1;
        const auto wspSize = wspInnerProd > m_nmTot ? wspInnerProd : m_nmTot;

        vec_t wsp[wspSize]; // workspace for kernels

        std::vector<vec_t, allocator<vec_t>> tmpIn(m_nmTot), tmpOut(m_nmTot);
        std::vector<vec_t, allocator<vec_t>> bwd(nqTot),  deriv0(nqTot), deriv1(nqTot);

        const vec_t* jac_ptr;
        const vec_t* df_ptr;

        vec_t df0,df1,df2,df3;
        vec_t metric00,metric01,metric11; 
        
        // Get size of derivative factor block
        auto dfSize = ndf;
        if (DEFORMED)
        {
            dfSize *= nqTot;
        }

        for (int e = 0; e < this->m_nBlocks; e++)
        {
            df_ptr = &(this->m_df[e*dfSize]);

            // Load and transpose data
            load_interleave(inptr, m_nmTot, tmpIn);

            // Precompute Laplacian metricsp
            if(!DEFORMED)
            {
                df0 = df_ptr[0];  df1 = df_ptr[1];
                df2 = df_ptr[2];  df3 = df_ptr[3];

                metric00 = df0*df0;
                metric00.fma(df2,df2);
                metric01 = df0 * df1; 
                metric01.fma(df2,df3);
                metric11 = df1*df1;
                metric11.fma(df3,df3);

                jac_ptr = &(this->m_jac[e]);
            }
            else
            {
                jac_ptr = &(this->m_jac[e*nqTot]);
            }

            // Step 1: BwdTrans
            BwdTransQuadKernel<NM0, NM1, NQ0, NQ1>
                (tmpIn, this->m_bdata[0], this->m_bdata[1], wsp, bwd);

            // Step 2: inner product for mass matrix operation
            IProductQuadKernel<NM0, NM1, NQ0, NQ1, true, false, DEFORMED>
                (bwd, this->m_bdata[0], this->m_bdata[1], this->m_w[0],
                 this->m_w[1], jac_ptr, wsp, tmpOut, m_lambda);

            // Step 3: take derivatives in quadrature space
            PhysDerivTensor2DKernel<NQ0, NQ1>
                (bwd, this->m_D[0], this->m_D[1], deriv0, deriv1);

            // Step 4: Apply Laplacian metrics & inner product
            if(DEFORMED)
            {
                for(size_t j = 0, cnt = 0; j < NQ1; ++j)
                {
                    for (size_t i = 0; i < NQ0; ++i, ++cnt)
                    {
                        df0 = df_ptr[cnt * ndf];
                        df1 = df_ptr[cnt * ndf + 1];
                        df2 = df_ptr[cnt * ndf + 2];
                        df3 = df_ptr[cnt * ndf + 3];


                        metric00 = df0*df0;
                        metric00.fma(df2,df2);
                        metric01 = df0 * df1; 
                        metric01.fma(df2,df3);
                        metric11 = df1*df1;
                        metric11.fma(df3,df3);

                        vec_t d0 = deriv0[cnt];
                        vec_t d1 = deriv1[cnt];

                        vec_t tmp = metric00 * d0;
                        tmp.fma(metric01, d1);
                        bwd[cnt]  = tmp;
                        
                        tmp = metric01 * d0;
                        tmp.fma(metric11, d1);
                        deriv0[cnt] = tmp;
                    }
                }
            }
            else
            {
                for (int i = 0; i < nqTot; ++i)
                {
                    vec_t d0 = deriv0[i];
                    vec_t d1 = deriv1[i];

                    vec_t tmp = metric00 * d0;
                    tmp.fma(metric01, d1);
                    bwd[i] = tmp;
                    
                    tmp = metric01 * d0;
                    tmp.fma(metric11, d1);
                    deriv0[i] = tmp; 
                }
            }
            
            IProductQuadKernel<NM0, NM1, NQ0, NQ1, false, true, DEFORMED>
                (bwd, this->m_dbdata[0], this->m_bdata[1], this->m_w[0],
                 this->m_w[1], jac_ptr, wsp, tmpOut);
            
            IProductQuadKernel<NM0, NM1, NQ0, NQ1, false, true, DEFORMED>
                (deriv0, this->m_bdata[0], this->m_dbdata[1], this->m_w[0],
                 this->m_w[1], jac_ptr, wsp, tmpOut);
           
            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);
            
            inptr  += nmBlocks;
            outptr += nmBlocks;
        }
    }

public:
    static NekDouble FlopsPerElement(const int nm, const int nq0, const int nq1)
    {
        const NekDouble bwdTrans  = BwdTransQuad::FlopsPerElement(nm,nq0,nq1);
        const NekDouble iprod1    = IProductQuad<DEFORMED>::FlopsPerElement(nm,nq0,nq1);
        const NekDouble physDeriv = PhysDerivQuad<DEFORMED>::FlopsPerElement(nq0,nq1);

        NekDouble metrics;
        if(DEFORMED)
        {
            metrics = nq1 * nq0 * (9 + 6);
        }
        else
        {
            const NekDouble metric_calc = 9;
            const NekDouble metrics_x0 = nq0 * nq1 * 3;
            const NekDouble metrics_x1 = nq0 * nq1 * 3;
            metrics = metric_calc + metrics_x0 + metrics_x1;
        }

        const NekDouble iprod2 = IProductQuad<DEFORMED>::FlopsPerElement(nm,nq0,nq1);
        const NekDouble iprod3 = IProductQuad<DEFORMED>::FlopsPerElement(nm,nq0,nq1);

        return bwdTrans + iprod1 + physDeriv + metrics + iprod2 + iprod3;
    }

    virtual NekDouble GFlops() override
    {
        const int nm  = m_basis[0]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();

        const NekDouble flops = m_nElmt * HelmholtzQuad::FlopsPerElement(nm, nq0, nq1);
        return flops * 1e-9;
    }

    virtual NekDouble Ndof() override
    {
        return m_nmTot * this->m_nElmt;
    }

private:
    int m_nmTot;
};


} // namespace MatrixFree
} // namespace Nektar

#endif
