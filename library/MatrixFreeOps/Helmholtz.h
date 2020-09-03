#ifndef NEKTAR_LIBRARY_MF_HELMHOLTZ_H
#define NEKTAR_LIBRARY_MF_HELMHOLTZ_H

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ShapeType.hpp>
#include <LibUtilities/Foundations/Basis.h>

#include "IProduct.h"
#include "IProductKernels.hpp"
#include "BwdTransKernels.hpp"
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
                              "HelmholtzQuad: # of modes / points combo not initialised.");
            } break;
        case 3:
            switch(m_basis[0]->GetNumPoints())
            {
            case 3: HelmholtzQuadImpl<3 ,3 ,3 ,3 >(in, out); break;
            case 4: HelmholtzQuadImpl<3 ,3 ,4 ,4 >(in, out); break;
            case 5: HelmholtzQuadImpl<3 ,3 ,5 ,5 >(in, out); break;
            case 6: HelmholtzQuadImpl<3 ,3 ,6 ,6 >(in, out); break;
            default: NEKERROR(ErrorUtil::efatal,
                              "HelmholtzQuad: # of modes / points combo not initialised.");
            } break;
        case 4:
            switch(m_basis[0]->GetNumPoints())
            {
            case 4: HelmholtzQuadImpl<4 ,4 ,4 ,4 >(in, out); break;
            case 5: HelmholtzQuadImpl<4 ,4 ,5 ,5 >(in, out); break;
            case 6: HelmholtzQuadImpl<4 ,4 ,6 ,6 >(in, out); break;
            case 7: HelmholtzQuadImpl<4 ,4 ,7 ,7 >(in, out); break;
            case 8: HelmholtzQuadImpl<4 ,4 ,8 ,8 >(in, out); break;
            default: NEKERROR(ErrorUtil::efatal,
                              "HelmholtzQuad: # of modes / points combo not initialised.");
            } break;
        case 5:
            switch(m_basis[0]->GetNumPoints())
            {
            case 5: HelmholtzQuadImpl<5 ,5 ,5 ,5 >(in, out); break;
            case 6: HelmholtzQuadImpl<5 ,5 ,6 ,6 >(in, out); break;
            case 7: HelmholtzQuadImpl<5 ,5 ,7 ,7 >(in, out); break;
            case 8: HelmholtzQuadImpl<5 ,5 ,8 ,8 >(in, out); break;
            case 9: HelmholtzQuadImpl<5 ,5 ,9 ,9 >(in, out); break;
            case 10: HelmholtzQuadImpl<5 ,5 ,10 ,10 >(in, out); break;
            default: NEKERROR(ErrorUtil::efatal,
                              "HelmholtzQuad: # of modes / points combo not initialised.");
            } break;
        case 6:
            switch(m_basis[0]->GetNumPoints())
            {
            case 6: HelmholtzQuadImpl<6 ,6 ,6 ,6 >(in, out); break;
            case 7: HelmholtzQuadImpl<6 ,6 ,7 ,7 >(in, out); break;
            case 8: HelmholtzQuadImpl<6 ,6 ,8 ,8 >(in, out); break;
            case 9: HelmholtzQuadImpl<6 ,6 ,9 ,9 >(in, out); break;
            case 10: HelmholtzQuadImpl<6 ,6 ,10 ,10 >(in, out); break;
            case 11: HelmholtzQuadImpl<6 ,6 ,11 ,11 >(in, out); break;
            case 12: HelmholtzQuadImpl<6 ,6 ,12 ,12 >(in, out); break;
            default: NEKERROR(ErrorUtil::efatal,
                              "HelmholtzQuad: # of modes / points combo not initialised.");
            } break;
        case 7:
            switch(m_basis[0]->GetNumPoints())
            {
            case 7: HelmholtzQuadImpl<7 ,7 ,7 ,7 >(in, out); break;
            case 8: HelmholtzQuadImpl<7 ,7 ,8 ,8 >(in, out); break;
            case 9: HelmholtzQuadImpl<7 ,7 ,9 ,9 >(in, out); break;
            case 10: HelmholtzQuadImpl<7 ,7 ,10 ,10 >(in, out); break;
            case 11: HelmholtzQuadImpl<7 ,7 ,11 ,11 >(in, out); break;
            case 12: HelmholtzQuadImpl<7 ,7 ,12 ,12 >(in, out); break;
            case 13: HelmholtzQuadImpl<7 ,7 ,13 ,13 >(in, out); break;
            case 14: HelmholtzQuadImpl<7 ,7 ,14 ,14 >(in, out); break;
            default: NEKERROR(ErrorUtil::efatal,
                              "HelmholtzQuad: # of modes / points combo not initialised.");
            } break;
        case 8:
            switch(m_basis[0]->GetNumPoints())
            {
            case 8: HelmholtzQuadImpl<8 ,8 ,8 ,8 >(in, out); break;
            case 9: HelmholtzQuadImpl<8 ,8 ,9 ,9 >(in, out); break;
            case 10: HelmholtzQuadImpl<8 ,8 ,10 ,10 >(in, out); break;
            case 11: HelmholtzQuadImpl<8 ,8 ,11 ,11 >(in, out); break;
            case 12: HelmholtzQuadImpl<8 ,8 ,12 ,12 >(in, out); break;
            case 13: HelmholtzQuadImpl<8 ,8 ,13 ,13 >(in, out); break;
            case 14: HelmholtzQuadImpl<8 ,8 ,14 ,14 >(in, out); break;
            case 15: HelmholtzQuadImpl<8 ,8 ,15 ,15 >(in, out); break;
            case 16: HelmholtzQuadImpl<8 ,8 ,16 ,16 >(in, out); break;
            default: NEKERROR(ErrorUtil::efatal,
                              "HelmholtzQuad: # of modes / points combo not initialised.");
            } break;;
        default: NEKERROR(ErrorUtil::efatal,
                          "HelmholtzQuad: # of modes / points combo not yet initialised.");
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
        const auto nmBlocks = m_nmTot * vec_t::width;

        // Allocate sufficient workspace for backwards transform and inner
        // product kernels.
        constexpr auto wspInnerProd = NQ1;
        constexpr auto wspBwdTrans = NQ0 * NM0;
        constexpr auto wspSize = wspInnerProd > wspBwdTrans ? wspInnerProd : wspBwdTrans;

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
    const int m_nmTot;
};


template<bool DEFORMED = false>
struct HelmholtzTri : public Helmholtz, public Helper<2, DEFORMED>
{
    HelmholtzTri(std::vector<LibUtilities::BasisSharedPtr> basis,
                     int nElmt)
    : Helmholtz(basis, nElmt),
      Helper<2, DEFORMED>(basis, nElmt),
      m_nmTot(LibUtilities::StdTriData::getNumberOfCoefficients(
                  this->m_nm[0], this->m_nm[1])),
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
        return std::make_shared<HelmholtzTri<DEFORMED>>(basis, nElmt);
    }

    virtual void operator()(const Array<OneD, const NekDouble> &in,
                                  Array<OneD,       NekDouble> &out) override
    {
        // Check preconditions
        ASSERTL0(m_basis[0]->GetNumModes() == m_basis[1]->GetNumModes() &&
            m_basis[0]->GetNumPoints() == m_basis[1]->GetNumPoints()+1,
            "MatrixFree requires homogenous modes/points");

        if (m_basis[0]->GetBasisType() == LibUtilities::eModified_A)
        {
            switch(m_basis[0]->GetNumModes())
            {
                case 2: switch(m_basis[0]->GetNumPoints())
                {
                    case 3: HelmholtzTriImpl<2, 2, 3, 2, true>(in, out); break;
                    case 4: HelmholtzTriImpl<2, 2, 4, 3, true>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "HelmholtzTri: # of modes / points combo not initialised.");
                } break;
            case 3:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 4: HelmholtzTriImpl<3, 3, 4, 3, true>(in, out); break;
                    case 5: HelmholtzTriImpl<3, 3, 5, 4, true>(in, out); break;
                    case 6: HelmholtzTriImpl<3, 3, 6, 5, true>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "HelmholtzTri: # of modes / points combo not initialised.");
                } break;
            case 4:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 5: HelmholtzTriImpl<4, 4, 5, 4, true>(in, out); break;
                    case 6: HelmholtzTriImpl<4, 4, 6, 5, true>(in, out); break;
                    case 7: HelmholtzTriImpl<4, 4, 7, 6, true>(in, out); break;
                    case 8: HelmholtzTriImpl<4, 4, 8, 7, true>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "HelmholtzTri: # of modes / points combo not initialised.");
                } break;
            case 5:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 6: HelmholtzTriImpl<5, 5, 6, 5, true>(in, out); break;
                    case 7: HelmholtzTriImpl<5, 5, 7, 6, true>(in, out); break;
                    case 8: HelmholtzTriImpl<5, 5, 8, 7, true>(in, out); break;
                    case 9: HelmholtzTriImpl<5, 5, 9, 8, true>(in, out); break;
                    case 10: HelmholtzTriImpl<5, 5, 10, 9, true>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "HelmholtzTri: # of modes / points combo not initialised.");
                } break;
            case 6:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 7: HelmholtzTriImpl<6, 6, 7, 6, true>(in, out); break;
                    case 8: HelmholtzTriImpl<6, 6, 8, 7, true>(in, out); break;
                    case 9: HelmholtzTriImpl<6, 6, 9, 8, true>(in, out); break;
                    case 10: HelmholtzTriImpl<6, 6, 10, 9, true>(in, out); break;
                    case 11: HelmholtzTriImpl<6, 6, 11, 10, true>(in, out); break;
                    case 12: HelmholtzTriImpl<6, 6, 12, 11, true>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "HelmholtzTri: # of modes / points combo not initialised.");
                } break;
            case 7:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 8: HelmholtzTriImpl<7, 7, 8, 7, true>(in, out); break;
                    case 9: HelmholtzTriImpl<7, 7, 9, 8, true>(in, out); break;
                    case 10: HelmholtzTriImpl<7, 7, 10, 9, true>(in, out); break;
                    case 11: HelmholtzTriImpl<7, 7, 11, 10, true>(in, out); break;
                    case 12: HelmholtzTriImpl<7, 7, 12, 11, true>(in, out); break;
                    case 13: HelmholtzTriImpl<7, 7, 13, 12, true>(in, out); break;
                    case 14: HelmholtzTriImpl<7, 7, 14, 13, true>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "HelmholtzTri: # of modes / points combo not initialised.");
                } break;
            case 8:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 9: HelmholtzTriImpl<8, 8, 9, 8, true>(in, out); break;
                    case 10: HelmholtzTriImpl<8, 8, 10, 9, true>(in, out); break;
                    case 11: HelmholtzTriImpl<8, 8, 11, 10, true>(in, out); break;
                    case 12: HelmholtzTriImpl<8, 8, 12, 11, true>(in, out); break;
                    case 13: HelmholtzTriImpl<8, 8, 13, 12, true>(in, out); break;
                    case 14: HelmholtzTriImpl<8, 8, 14, 13, true>(in, out); break;
                    case 15: HelmholtzTriImpl<8, 8, 15, 14, true>(in, out); break;
                    case 16: HelmholtzTriImpl<8, 8, 16, 15, true>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "HelmholtzTri: # of modes / points combo not initialised.");
                } break;
            default: NEKERROR(ErrorUtil::efatal,
                "HelmholtzTri: # of modes / points combo not initialised.");

            }
        }
        else
        {
            switch(m_basis[0]->GetNumModes())
            {
                case 2: switch(m_basis[0]->GetNumPoints())
                {
                    case 3: HelmholtzTriImpl<2 ,2 ,3 ,2, false >(in, out); break;
                    case 4: HelmholtzTriImpl<2 ,2 ,4 ,3, false >(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "HelmholtzTri: # of modes / points combo not initialised.");
                } break;
            case 3:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 4: HelmholtzTriImpl<3 ,3 ,4 ,3, false >(in, out); break;
                    case 5: HelmholtzTriImpl<3 ,3 ,5 ,4, false >(in, out); break;
                    case 6: HelmholtzTriImpl<3 ,3 ,6 ,5, false >(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "HelmholtzTri: # of modes / points combo not initialised.");
                } break;
            case 4:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 5: HelmholtzTriImpl<4, 4, 5, 4, false>(in, out); break;
                    case 6: HelmholtzTriImpl<4, 4, 6, 5, false>(in, out); break;
                    case 7: HelmholtzTriImpl<4, 4, 7, 6, false>(in, out); break;
                    case 8: HelmholtzTriImpl<4, 4, 8, 7, false>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "HelmholtzTri: # of modes / points combo not initialised.");
                } break;
            case 5:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 6: HelmholtzTriImpl<5, 5, 6, 5, false>(in, out); break;
                    case 7: HelmholtzTriImpl<5, 5, 7, 6, false>(in, out); break;
                    case 8: HelmholtzTriImpl<5, 5, 8, 7, false>(in, out); break;
                    case 9: HelmholtzTriImpl<5, 5, 9, 8, false>(in, out); break;
                    case 10: HelmholtzTriImpl<5, 5, 10, 9, false>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "HelmholtzTri: # of modes / points combo not initialised.");
                } break;
            case 6:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 7: HelmholtzTriImpl<6, 6, 7, 6, false>(in, out); break;
                    case 8: HelmholtzTriImpl<6, 6, 8, 7, false>(in, out); break;
                    case 9: HelmholtzTriImpl<6, 6, 9, 8, false>(in, out); break;
                    case 10: HelmholtzTriImpl<6, 6, 10, 9, false>(in, out); break;
                    case 11: HelmholtzTriImpl<6, 6, 11, 10, false>(in, out); break;
                    case 12: HelmholtzTriImpl<6, 6, 12, 11, false>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "HelmholtzTri: # of modes / points combo not initialised.");
                } break;
            case 7:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 8: HelmholtzTriImpl<7, 7, 8, 7, false>(in, out); break;
                    case 9: HelmholtzTriImpl<7, 7, 9, 8, false>(in, out); break;
                    case 10: HelmholtzTriImpl<7, 7, 10, 9, false>(in, out); break;
                    case 11: HelmholtzTriImpl<7, 7, 11, 10, false>(in, out); break;
                    case 12: HelmholtzTriImpl<7, 7, 12, 11, false>(in, out); break;
                    case 13: HelmholtzTriImpl<7, 7, 13, 12, false>(in, out); break;
                    case 14: HelmholtzTriImpl<7, 7, 14, 13, false>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "HelmholtzTri: # of modes / points combo not initialised.");
                } break;
            case 8:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 9: HelmholtzTriImpl<8, 8, 9, 8, false>(in, out); break;
                    case 10: HelmholtzTriImpl<8, 8, 10, 9, false>(in, out); break;
                    case 11: HelmholtzTriImpl<8, 8, 11, 10, false>(in, out); break;
                    case 12: HelmholtzTriImpl<8, 8, 12, 11, false>(in, out); break;
                    case 13: HelmholtzTriImpl<8, 8, 13, 12, false>(in, out); break;
                    case 14: HelmholtzTriImpl<8, 8, 14, 13, false>(in, out); break;
                    case 15: HelmholtzTriImpl<8, 8, 15, 14, false>(in, out); break;
                    case 16: HelmholtzTriImpl<8, 8, 16, 15, false>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "HelmholtzTri: # of modes / points combo not initialised.");
                } break;
            default: NEKERROR(ErrorUtil::efatal,
                "HelmholtzTri: # of modes / points combo not initialised.");

            }
        }
    }

    template<int NM0, int NM1, int NQ0, int NQ1, bool CORRECT>
    void HelmholtzTriImpl(
        const Array<OneD, const NekDouble> &input, Array<OneD,NekDouble> &out)
    {
        auto *inptr  = &input[0];
        auto *outptr = &out[0];

        constexpr auto ndf = 4;
        constexpr auto nqTot = NQ0 * NQ1;
        const auto nmBlocks = m_nmTot * vec_t::width;


        // Allocate sufficient workspace for backwards transform and inner
        // product kernels.
        constexpr auto wspInnerProd = NQ1;
        constexpr auto wspBwdTrans = NM0;
        constexpr auto wspSize = wspInnerProd > wspBwdTrans ? wspInnerProd : wspBwdTrans;

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
            BwdTransTriKernel<NM0, NM1, NQ0, NQ1, CORRECT>
                (tmpIn, this->m_bdata[0], this->m_bdata[1], wsp, bwd);

            // Step 2: inner product for mass matrix operation
            IProductTriKernel<NM0, NM1, NQ0, NQ1, CORRECT, true, false, DEFORMED>
                (bwd, this->m_bdata[0], this->m_bdata[1], this->m_w[0],
                this->m_w[1], jac_ptr, wsp, tmpOut, m_lambda);

            // Step 3: take derivatives in collapsed coordinate space
            PhysDerivTensor2DKernel<NQ0, NQ1>(
                bwd, this->m_D[0], this->m_D[1], deriv0, deriv1);

            // Step 4a: Construct Laplacian metrics
            for (size_t j = 0, cnt = 0; j < NQ1; ++j)
            {
                vec_t h1j = m_h1[j];
                for (size_t i = 0; i < NQ0; ++i, ++cnt)
                {
                    if(DEFORMED)
                    {
                        df0 = df_ptr[cnt * ndf];
                        df1 = df_ptr[cnt * ndf + 1];
                        df2 = df_ptr[cnt * ndf + 2];
                        df3 = df_ptr[cnt * ndf + 3];
                    }

                    vec_t h0i = m_h0[i];
                    metric00 = h1j * (df0 + h0i * df1);
                    metric01 = metric00 * df1;
                    metric00 = metric00 * metric00;

                    vec_t tmp = h1j * (df2 + h0i * df3);
                    metric01.fma(tmp, df3);
                    metric00.fma(tmp, tmp);

                    metric11 = df1 * df1;
                    metric11.fma(df3, df3);

                    vec_t d0 = deriv0[cnt];
                    vec_t d1 = deriv1[cnt];

                    tmp = metric00 * d0;
                    tmp.fma(metric01, d1);
                    bwd[cnt]  = tmp;
                    
                    tmp = metric01 * d0;
                    tmp.fma(metric11, d1);
                    deriv0[cnt] = tmp;
                }
            }

            // Step 4b: Take inner products
            IProductTriKernel<NM0, NM1, NQ0, NQ1, CORRECT, false, true, DEFORMED>(
                bwd, this->m_dbdata[0], this->m_bdata[1], this->m_w[0],
                this->m_w[1], jac_ptr, wsp, tmpOut);

            IProductTriKernel<NM0, NM1, NQ0, NQ1, CORRECT, false, true, DEFORMED>(
                deriv0, this->m_bdata[0], this->m_dbdata[1], this->m_w[0],
                this->m_w[1], jac_ptr, wsp, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr  += nmBlocks;
            outptr += nmBlocks;
        }
    }
public:

        static NekDouble FlopsPerElement(
        const int nm,
        const int nq0,
        const int nq1)
    {

        const double bwdTrans = BwdTransTri::FlopsPerElement(nm,nq0,nq1);
        const double iprod1 =   IProductTri<DEFORMED>::FlopsPerElement(nm, nq0, nq1);
        const double physDeriv = PhysDerivTri<DEFORMED>::FlopsPerElement(nq0,nq1);
        const double metric = nq1 * nq0 * (5 + 7 + 3 + 3 + 3);
        const double iprod2 = IProductTri<DEFORMED>::FlopsPerElement(nm, nq0, nq1);
        const double iprod3 = IProductTri<DEFORMED>::FlopsPerElement(nm, nq0, nq1);

        return bwdTrans + iprod1 + physDeriv + metric + iprod2 + iprod3;
    }

    virtual double GFlops() override
    {
        const int nm  = m_basis[0]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();

        const double flops = m_nElmt * HelmholtzTri::FlopsPerElement(nm, nq0, nq1);
        return flops * 1e-9;
    }

    virtual NekDouble Ndof() override
    {
        return m_nmTot * this->m_nElmt;
    }
    
private:
    const int m_nmTot;
    std::vector<vec_t, allocator<vec_t>> m_h0;
    std::vector<vec_t, allocator<vec_t>> m_h1;
};

} // namespace MatrixFree
} // namespace Nektar

#endif
