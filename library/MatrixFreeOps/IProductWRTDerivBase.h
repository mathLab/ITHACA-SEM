#ifndef NEKTAR_LIBRARY_MF_IPRODUCTWRTDERIVBASE_H
#define NEKTAR_LIBRARY_MF_IPRODUCTWRTDERIVBASE_H

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ShapeType.hpp>
#include <LibUtilities/Foundations/Basis.h>

#include "Operator.hpp"

#include "IProductKernels.hpp"

namespace Nektar
{
namespace MatrixFree
{

template<bool DEFORMED = false>
struct IProductWRTDerivBaseQuad : public IProductWRTDerivBase, public Helper<2, DEFORMED>
{
public:

    IProductWRTDerivBaseQuad(std::vector<LibUtilities::BasisSharedPtr> basis,
                   int nElmt)
        : IProductWRTDerivBase(basis, nElmt),
          Helper<2, DEFORMED>(basis, nElmt),
          m_nmTot(LibUtilities::StdQuadData::getNumberOfCoefficients(
                      this->m_nm[0], this->m_nm[1]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<IProductWRTDerivBaseQuad<DEFORMED>>(basis, nElmt);
    }

    void operator()( const Array<OneD, Array<OneD, NekDouble>> &in,
        Array<OneD, NekDouble> &out) final
    {
        switch(m_basis[0]->GetNumModes())
        {
            case 2:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 2: IProductWRTDerivBaseQuadImpl<2 ,2 ,2 ,2 >
                    (in[0], in[1], out); break;
                    case 3: IProductWRTDerivBaseQuadImpl<2 ,2 ,3 ,3 >
                    (in[0], in[1], out); break;
                    case 4: IProductWRTDerivBaseQuadImpl<2 ,2 ,4 ,4 >
                    (in[0], in[1], out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductWRTDerivBaseQuad: # of modes / points combo not"
                " implemented.");
                } break;
            case 3:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 3: IProductWRTDerivBaseQuadImpl<3 ,3 ,3 ,3 >
                    (in[0], in[1], out); break;
                    case 4: IProductWRTDerivBaseQuadImpl<3 ,3 ,4 ,4 >
                    (in[0], in[1], out); break;
                    case 5: IProductWRTDerivBaseQuadImpl<3 ,3 ,5 ,5 >
                    (in[0], in[1], out); break;
                    case 6: IProductWRTDerivBaseQuadImpl<3 ,3 ,6 ,6 >
                    (in[0], in[1], out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductWRTDerivBaseQuad: # of modes / points combo not"
                " implemented.");
                } break;
            case 4:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 4: IProductWRTDerivBaseQuadImpl<4 ,4 ,4 ,4 >
                    (in[0], in[1], out); break;
                    case 5: IProductWRTDerivBaseQuadImpl<4 ,4 ,5 ,5 >
                    (in[0], in[1], out); break;
                    case 6: IProductWRTDerivBaseQuadImpl<4 ,4 ,6 ,6 >
                    (in[0], in[1], out); break;
                    case 7: IProductWRTDerivBaseQuadImpl<4 ,4 ,7 ,7 >
                    (in[0], in[1], out); break;
                    case 8: IProductWRTDerivBaseQuadImpl<4 ,4 ,8 ,8 >
                    (in[0], in[1], out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductWRTDerivBaseQuad: # of modes / points combo not"
                " implemented.");
                } break;
            case 5:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 5: IProductWRTDerivBaseQuadImpl<5 ,5 ,5 ,5 >
                    (in[0], in[1], out); break;
                    case 6: IProductWRTDerivBaseQuadImpl<5 ,5 ,6 ,6 >
                    (in[0], in[1], out); break;
                    case 7: IProductWRTDerivBaseQuadImpl<5 ,5 ,7 ,7 >
                    (in[0], in[1], out); break;
                    case 8: IProductWRTDerivBaseQuadImpl<5 ,5 ,8 ,8 >
                    (in[0], in[1], out); break;
                    case 9: IProductWRTDerivBaseQuadImpl<5 ,5 ,9 ,9 >
                    (in[0], in[1], out); break;
                    case 10: IProductWRTDerivBaseQuadImpl<5 ,5 ,10 ,10 >
                    (in[0], in[1], out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductWRTDerivBaseQuad: # of modes / points combo not"
                " implemented.");
                } break;
            case 6:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 6: IProductWRTDerivBaseQuadImpl<6 ,6 ,6 ,6 >
                    (in[0], in[1], out); break;
                    case 7: IProductWRTDerivBaseQuadImpl<6 ,6 ,7 ,7 >
                    (in[0], in[1], out); break;
                    case 8: IProductWRTDerivBaseQuadImpl<6 ,6 ,8 ,8 >
                    (in[0], in[1], out); break;
                    case 9: IProductWRTDerivBaseQuadImpl<6 ,6 ,9 ,9 >
                    (in[0], in[1], out); break;
                    case 10: IProductWRTDerivBaseQuadImpl<6 ,6 ,10 ,10 >
                    (in[0], in[1], out); break;
                    case 11: IProductWRTDerivBaseQuadImpl<6 ,6 ,11 ,11 >
                    (in[0], in[1], out); break;
                    case 12: IProductWRTDerivBaseQuadImpl<6 ,6 ,12 ,12 >
                    (in[0], in[1], out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductWRTDerivBaseQuad: # of modes / points combo not"
                " implemented.");
                } break;
            case 7:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 7: IProductWRTDerivBaseQuadImpl<7 ,7 ,7 ,7 >
                    (in[0], in[1], out); break;
                    case 8: IProductWRTDerivBaseQuadImpl<7 ,7 ,8 ,8 >
                    (in[0], in[1], out); break;
                    case 9: IProductWRTDerivBaseQuadImpl<7 ,7 ,9 ,9 >
                    (in[0], in[1], out); break;
                    case 10: IProductWRTDerivBaseQuadImpl<7 ,7 ,10 ,10 >
                    (in[0], in[1], out); break;
                    case 11: IProductWRTDerivBaseQuadImpl<7 ,7 ,11 ,11 >
                    (in[0], in[1], out); break;
                    case 12: IProductWRTDerivBaseQuadImpl<7 ,7 ,12 ,12 >
                    (in[0], in[1], out); break;
                    case 13: IProductWRTDerivBaseQuadImpl<7 ,7 ,13 ,13 >
                    (in[0], in[1], out); break;
                    case 14: IProductWRTDerivBaseQuadImpl<7 ,7 ,14 ,14 >
                    (in[0], in[1], out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductWRTDerivBaseQuad: # of modes / points combo not"
                " implemented.");
                } break;
            case 8:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 8: IProductWRTDerivBaseQuadImpl<8 ,8 ,8 ,8 >
                    (in[0], in[1], out); break;
                    case 9: IProductWRTDerivBaseQuadImpl<8 ,8 ,9 ,9 >
                    (in[0], in[1], out); break;
                    case 10: IProductWRTDerivBaseQuadImpl<8 ,8 ,10 ,10 >
                    (in[0], in[1], out); break;
                    case 11: IProductWRTDerivBaseQuadImpl<8 ,8 ,11 ,11 >
                    (in[0], in[1], out); break;
                    case 12: IProductWRTDerivBaseQuadImpl<8 ,8 ,12 ,12 >
                    (in[0], in[1], out); break;
                    case 13: IProductWRTDerivBaseQuadImpl<8 ,8 ,13 ,13 >
                    (in[0], in[1], out); break;
                    case 14: IProductWRTDerivBaseQuadImpl<8 ,8 ,14 ,14 >
                    (in[0], in[1], out); break;
                    case 15: IProductWRTDerivBaseQuadImpl<8 ,8 ,15 ,15 >
                    (in[0], in[1], out); break;
                    case 16: IProductWRTDerivBaseQuadImpl<8 ,8 ,16 ,16 >
                    (in[0], in[1], out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductWRTDerivBaseQuad: # of modes / points combo not"
                " implemented.");
                } break;;
            default: NEKERROR(ErrorUtil::efatal,
                "IProductWRTDerivBaseQuad: # of modes / points combo not"
                " implemented.");
        }
    }

    template<int NM0, int NM1, int NQ0, int NQ1>
    void IProductWRTDerivBaseQuadImpl(
        const Array<OneD, const NekDouble> &input0,
        const Array<OneD, const NekDouble> &input1,
              Array<OneD,       NekDouble> &output)
    {
        auto* inptr0 = input0.data();
        auto* inptr1 = input1.data();
        auto* outptr = output.data();

        constexpr auto ndf = 4u;
        constexpr auto nqTot = NQ0 * NQ1;
        constexpr auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        // Get size of jacobian factor block
        auto dJSize = 1u;
        auto dfSize = ndf;
        if(DEFORMED)
        {
            dJSize = nqTot;
            dfSize = ndf*nqTot;
        }

        vec_t sums_j[NQ1]; //Sums over eta0 for each value of eta1;
        std::vector<vec_t, allocator<vec_t>> tmpIn0(nqTot), tmpIn1(nqTot),
            tmp0(nqTot), tmp1(nqTot), tmpOut(m_nmTot);

        const vec_t* df_ptr;
        const vec_t* jac_ptr;
        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            // Jacobian
            jac_ptr = &(this->m_jac[dJSize*e]);

            // Derivative factor
            df_ptr = &(this->m_df[e*dfSize]);

            // Load and transpose data
            load_interleave(inptr0, nqTot, tmpIn0);
            load_interleave(inptr1, nqTot, tmpIn1);

            // Calculate dx/dxi in[0] + dy/dxi in[1]
            vec_t df0, df1, df2, df3;
            if(!DEFORMED)
            {
                df0 = df_ptr[0];
                df1 = df_ptr[1];
                df2 = df_ptr[2];
                df3 = df_ptr[3];
            }
            for (int i = 0; i < nqTot; ++i)
            {
                if(DEFORMED)
                {
                    df0 = df_ptr[i * ndf];
                    df1 = df_ptr[i * ndf + 1];
                    df2 = df_ptr[i * ndf + 2];
                    df3 = df_ptr[i * ndf + 3];
                }
                tmp0[i] = df0 * tmpIn0[i] + df2 * tmpIn1[i];
                tmp1[i] = df1 * tmpIn0[i] + df3 * tmpIn1[i];
            }

            // IP DB0 B1
            IProductQuadKernel<NM0, NM1, NQ0, NQ1, false, false, DEFORMED>(
                tmp0, this->m_dbdata[0], this->m_bdata[1],
                this->m_w[0], this->m_w[1], jac_ptr,
                sums_j, tmpOut);

            // IP DB1 B0
            IProductQuadKernel<NM0, NM1, NQ0, NQ1, false, true, DEFORMED>(
                tmp1, this->m_bdata[0], this->m_dbdata[1],
                this->m_w[0], this->m_w[1], jac_ptr,
                sums_j, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr0 += nqBlocks;
            inptr1 += nqBlocks;
            outptr += nmBlocks;
        }
    }


private:
    int m_nmTot;
};

template<bool DEFORMED = false>
struct IProductWRTDerivBaseTri : public IProductWRTDerivBase, public Helper<2, DEFORMED>
{
    IProductWRTDerivBaseTri(std::vector<LibUtilities::BasisSharedPtr> basis,
                   int nElmt)
        : IProductWRTDerivBase(basis, nElmt),
          Helper<2, DEFORMED>(basis, nElmt),
          m_nmTot(LibUtilities::StdTriData::getNumberOfCoefficients(
                      this->m_nm[0], this->m_nm[1]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<IProductWRTDerivBaseTri<DEFORMED>>(basis, nElmt);
    }

   void operator()( const Array<OneD, Array<OneD, NekDouble>> &in,
        Array<OneD, NekDouble> &out) final
    {
        // Check preconditions
        ASSERTL0(m_basis[0]->GetNumModes() == m_basis[1]->GetNumModes() &&
            m_basis[0]->GetNumPoints() == m_basis[1]->GetNumPoints()+1,
            "MatrixFree requires homogenous modes/points");

        if (m_basis[0]->GetBasisType() == LibUtilities::eModified_A)
        {
            switch(m_basis[0]->GetNumModes())
            {
                case 2:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 3: IProductWRTDerivBaseTriImpl<2 ,2 ,3 ,2 ,true>(in[0], in[1], out); break;
                        case 4: IProductWRTDerivBaseTriImpl<2 ,2 ,4 ,3 ,true>(in[0], in[1], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTri: # of modes / points combo not implemented.");
                    } break;
                case 3:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 4: IProductWRTDerivBaseTriImpl<3 ,3 ,4 ,3 ,true>(in[0], in[1], out); break;
                        case 5: IProductWRTDerivBaseTriImpl<3 ,3 ,5 ,4 ,true>(in[0], in[1], out); break;
                        case 6: IProductWRTDerivBaseTriImpl<3 ,3 ,6 ,5 ,true>(in[0], in[1], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTri: # of modes / points combo not implemented.");
                    } break;
                case 4:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 5: IProductWRTDerivBaseTriImpl<4 ,4 ,5 ,4 ,true>(in[0], in[1], out); break;
                        case 6: IProductWRTDerivBaseTriImpl<4 ,4 ,6 ,5 ,true>(in[0], in[1], out); break;
                        case 7: IProductWRTDerivBaseTriImpl<4 ,4 ,7 ,6 ,true>(in[0], in[1], out); break;
                        case 8: IProductWRTDerivBaseTriImpl<4 ,4 ,8 ,7 ,true>(in[0], in[1], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTri: # of modes / points combo not implemented.");
                    } break;
                case 5:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 6: IProductWRTDerivBaseTriImpl<5 ,5 ,6 ,5 ,true>(in[0], in[1], out); break;
                        case 7: IProductWRTDerivBaseTriImpl<5 ,5 ,7 ,6 ,true>(in[0], in[1], out); break;
                        case 8: IProductWRTDerivBaseTriImpl<5 ,5 ,8 ,7 ,true>(in[0], in[1], out); break;
                        case 9: IProductWRTDerivBaseTriImpl<5 ,5 ,9 ,8 ,true>(in[0], in[1], out); break;
                        case 10: IProductWRTDerivBaseTriImpl<5 ,5 ,10 ,9 ,true>(in[0], in[1], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTri: # of modes / points combo not implemented.");
                    } break;
                case 6:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 7: IProductWRTDerivBaseTriImpl<6 ,6 ,7 ,6 ,true>(in[0], in[1], out); break;
                        case 8: IProductWRTDerivBaseTriImpl<6 ,6 ,8 ,7 ,true>(in[0], in[1], out); break;
                        case 9: IProductWRTDerivBaseTriImpl<6 ,6 ,9 ,8 ,true>(in[0], in[1], out); break;
                        case 10: IProductWRTDerivBaseTriImpl<6 ,6 ,10 ,9 ,true>(in[0], in[1], out); break;
                        case 11: IProductWRTDerivBaseTriImpl<6 ,6 ,11 ,10 ,true>(in[0], in[1], out); break;
                        case 12: IProductWRTDerivBaseTriImpl<6 ,6 ,12 ,11 ,true>(in[0], in[1], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTri: # of modes / points combo not implemented.");
                    } break;
                case 7:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 8: IProductWRTDerivBaseTriImpl<7 ,7 ,8 ,7 ,true>(in[0], in[1], out); break;
                        case 9: IProductWRTDerivBaseTriImpl<7 ,7 ,9 ,8 ,true>(in[0], in[1], out); break;
                        case 10: IProductWRTDerivBaseTriImpl<7 ,7 ,10 ,9 ,true>(in[0], in[1], out); break;
                        case 11: IProductWRTDerivBaseTriImpl<7 ,7 ,11 ,10 ,true>(in[0], in[1], out); break;
                        case 12: IProductWRTDerivBaseTriImpl<7 ,7 ,12 ,11 ,true>(in[0], in[1], out); break;
                        case 13: IProductWRTDerivBaseTriImpl<7 ,7 ,13 ,12 ,true>(in[0], in[1], out); break;
                        case 14: IProductWRTDerivBaseTriImpl<7 ,7 ,14 ,13 ,true>(in[0], in[1], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTri: # of modes / points combo not implemented.");
                    } break;
                case 8:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 9: IProductWRTDerivBaseTriImpl<8 ,8 ,9 ,8 ,true>(in[0], in[1], out); break;
                        case 10: IProductWRTDerivBaseTriImpl<8 ,8 ,10 ,9 ,true>(in[0], in[1], out); break;
                        case 11: IProductWRTDerivBaseTriImpl<8 ,8 ,11 ,10 ,true>(in[0], in[1], out); break;
                        case 12: IProductWRTDerivBaseTriImpl<8 ,8 ,12 ,11 ,true>(in[0], in[1], out); break;
                        case 13: IProductWRTDerivBaseTriImpl<8 ,8 ,13 ,12 ,true>(in[0], in[1], out); break;
                        case 14: IProductWRTDerivBaseTriImpl<8 ,8 ,14 ,13 ,true>(in[0], in[1], out); break;
                        case 15: IProductWRTDerivBaseTriImpl<8 ,8 ,15 ,14 ,true>(in[0], in[1], out); break;
                        case 16: IProductWRTDerivBaseTriImpl<8 ,8 ,16 ,15 ,true>(in[0], in[1], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTri: # of modes / points combo not implemented.");
                    } break;
                default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTri: # of modes / points combo not implemented.");
            }
        }
        else
        {
            switch(m_basis[0]->GetNumModes())
            {
                case 2:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 3: IProductWRTDerivBaseTriImpl<2 ,2 ,3 ,2 ,false>(in[0], in[1], out); break;
                        case 4: IProductWRTDerivBaseTriImpl<2 ,2 ,4 ,3 ,false>(in[0], in[1], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTri: # of modes / points combo not implemented.");
                    } break;
                case 3:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 4: IProductWRTDerivBaseTriImpl<3 ,3 ,4 ,3 ,false>(in[0], in[1], out); break;
                        case 5: IProductWRTDerivBaseTriImpl<3 ,3 ,5 ,4 ,false>(in[0], in[1], out); break;
                        case 6: IProductWRTDerivBaseTriImpl<3 ,3 ,6 ,5 ,false>(in[0], in[1], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTri: # of modes / points combo not implemented.");
                    } break;
                case 4:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 5: IProductWRTDerivBaseTriImpl<4 ,4 ,5 ,4 ,false>(in[0], in[1], out); break;
                        case 6: IProductWRTDerivBaseTriImpl<4 ,4 ,6 ,5 ,false>(in[0], in[1], out); break;
                        case 7: IProductWRTDerivBaseTriImpl<4 ,4 ,7 ,6 ,false>(in[0], in[1], out); break;
                        case 8: IProductWRTDerivBaseTriImpl<4 ,4 ,8 ,7 ,false>(in[0], in[1], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTri: # of modes / points combo not implemented.");
                    } break;
                case 5:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 6: IProductWRTDerivBaseTriImpl<5 ,5 ,6 ,5 ,false>(in[0], in[1], out); break;
                        case 7: IProductWRTDerivBaseTriImpl<5 ,5 ,7 ,6 ,false>(in[0], in[1], out); break;
                        case 8: IProductWRTDerivBaseTriImpl<5 ,5 ,8 ,7 ,false>(in[0], in[1], out); break;
                        case 9: IProductWRTDerivBaseTriImpl<5 ,5 ,9 ,8 ,false>(in[0], in[1], out); break;
                        case 10: IProductWRTDerivBaseTriImpl<5 ,5 ,10 ,9 ,false>(in[0], in[1], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTri: # of modes / points combo not implemented.");
                    } break;
                case 6:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 7: IProductWRTDerivBaseTriImpl<6 ,6 ,7 ,6 ,false>(in[0], in[1], out); break;
                        case 8: IProductWRTDerivBaseTriImpl<6 ,6 ,8 ,7 ,false>(in[0], in[1], out); break;
                        case 9: IProductWRTDerivBaseTriImpl<6 ,6 ,9 ,8 ,false>(in[0], in[1], out); break;
                        case 10: IProductWRTDerivBaseTriImpl<6 ,6 ,10 ,9 ,false>(in[0], in[1], out); break;
                        case 11: IProductWRTDerivBaseTriImpl<6 ,6 ,11 ,10 ,false>(in[0], in[1], out); break;
                        case 12: IProductWRTDerivBaseTriImpl<6 ,6 ,12 ,11 ,false>(in[0], in[1], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTri: # of modes / points combo not implemented.");
                    } break;
                case 7:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 8: IProductWRTDerivBaseTriImpl<7 ,7 ,8 ,7 ,false>(in[0], in[1], out); break;
                        case 9: IProductWRTDerivBaseTriImpl<7 ,7 ,9 ,8 ,false>(in[0], in[1], out); break;
                        case 10: IProductWRTDerivBaseTriImpl<7 ,7 ,10 ,9 ,false>(in[0], in[1], out); break;
                        case 11: IProductWRTDerivBaseTriImpl<7 ,7 ,11 ,10 ,false>(in[0], in[1], out); break;
                        case 12: IProductWRTDerivBaseTriImpl<7 ,7 ,12 ,11 ,false>(in[0], in[1], out); break;
                        case 13: IProductWRTDerivBaseTriImpl<7 ,7 ,13 ,12 ,false>(in[0], in[1], out); break;
                        case 14: IProductWRTDerivBaseTriImpl<7 ,7 ,14 ,13 ,false>(in[0], in[1], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTri: # of modes / points combo not implemented.");
                    } break;
                case 8:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 9: IProductWRTDerivBaseTriImpl<8 ,8 ,9 ,8 ,false>(in[0], in[1], out); break;
                        case 10: IProductWRTDerivBaseTriImpl<8 ,8 ,10 ,9 ,false>(in[0], in[1], out); break;
                        case 11: IProductWRTDerivBaseTriImpl<8 ,8 ,11 ,10 ,false>(in[0], in[1], out); break;
                        case 12: IProductWRTDerivBaseTriImpl<8 ,8 ,12 ,11 ,false>(in[0], in[1], out); break;
                        case 13: IProductWRTDerivBaseTriImpl<8 ,8 ,13 ,12 ,false>(in[0], in[1], out); break;
                        case 14: IProductWRTDerivBaseTriImpl<8 ,8 ,14 ,13 ,false>(in[0], in[1], out); break;
                        case 15: IProductWRTDerivBaseTriImpl<8 ,8 ,15 ,14 ,false>(in[0], in[1], out); break;
                        case 16: IProductWRTDerivBaseTriImpl<8 ,8 ,16 ,15 ,false>(in[0], in[1], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTri: # of modes / points combo not implemented.");
                    } break;
                default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTri: # of modes / points combo not implemented.");
            }

            switch(m_basis[0]->GetNumModes())
            {
                case 2:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 3: IProductWRTDerivBaseTriImpl<2 ,2 ,3 ,2 ,true>(in[0], in[1], out); break;
                        case 4: IProductWRTDerivBaseTriImpl<2 ,2 ,4 ,3 ,true>(in[0], in[1], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTri: # of modes / points combo not implemented.");
                    } break;
                case 3:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 4: IProductWRTDerivBaseTriImpl<3 ,3 ,4 ,3 ,true>(in[0], in[1], out); break;
                        case 5: IProductWRTDerivBaseTriImpl<3 ,3 ,5 ,4 ,true>(in[0], in[1], out); break;
                        case 6: IProductWRTDerivBaseTriImpl<3 ,3 ,6 ,5 ,true>(in[0], in[1], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTri: # of modes / points combo not implemented.");
                    } break;
                case 4:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 5: IProductWRTDerivBaseTriImpl<4 ,4 ,5 ,4 ,true>(in[0], in[1], out); break;
                        case 6: IProductWRTDerivBaseTriImpl<4 ,4 ,6 ,5 ,true>(in[0], in[1], out); break;
                        case 7: IProductWRTDerivBaseTriImpl<4 ,4 ,7 ,6 ,true>(in[0], in[1], out); break;
                        case 8: IProductWRTDerivBaseTriImpl<4 ,4 ,8 ,7 ,true>(in[0], in[1], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTri: # of modes / points combo not implemented.");
                    } break;
                case 5:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 6: IProductWRTDerivBaseTriImpl<5 ,5 ,6 ,5 ,true>(in[0], in[1], out); break;
                        case 7: IProductWRTDerivBaseTriImpl<5 ,5 ,7 ,6 ,true>(in[0], in[1], out); break;
                        case 8: IProductWRTDerivBaseTriImpl<5 ,5 ,8 ,7 ,true>(in[0], in[1], out); break;
                        case 9: IProductWRTDerivBaseTriImpl<5 ,5 ,9 ,8 ,true>(in[0], in[1], out); break;
                        case 10: IProductWRTDerivBaseTriImpl<5 ,5 ,10 ,9 ,true>(in[0], in[1], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTri: # of modes / points combo not implemented.");
                    } break;
                case 6:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 7: IProductWRTDerivBaseTriImpl<6 ,6 ,7 ,6 ,true>(in[0], in[1], out); break;
                        case 8: IProductWRTDerivBaseTriImpl<6 ,6 ,8 ,7 ,true>(in[0], in[1], out); break;
                        case 9: IProductWRTDerivBaseTriImpl<6 ,6 ,9 ,8 ,true>(in[0], in[1], out); break;
                        case 10: IProductWRTDerivBaseTriImpl<6 ,6 ,10 ,9 ,true>(in[0], in[1], out); break;
                        case 11: IProductWRTDerivBaseTriImpl<6 ,6 ,11 ,10 ,true>(in[0], in[1], out); break;
                        case 12: IProductWRTDerivBaseTriImpl<6 ,6 ,12 ,11 ,true>(in[0], in[1], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTri: # of modes / points combo not implemented.");
                    } break;
                case 7:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 8: IProductWRTDerivBaseTriImpl<7 ,7 ,8 ,7 ,true>(in[0], in[1], out); break;
                        case 9: IProductWRTDerivBaseTriImpl<7 ,7 ,9 ,8 ,true>(in[0], in[1], out); break;
                        case 10: IProductWRTDerivBaseTriImpl<7 ,7 ,10 ,9 ,true>(in[0], in[1], out); break;
                        case 11: IProductWRTDerivBaseTriImpl<7 ,7 ,11 ,10 ,true>(in[0], in[1], out); break;
                        case 12: IProductWRTDerivBaseTriImpl<7 ,7 ,12 ,11 ,true>(in[0], in[1], out); break;
                        case 13: IProductWRTDerivBaseTriImpl<7 ,7 ,13 ,12 ,true>(in[0], in[1], out); break;
                        case 14: IProductWRTDerivBaseTriImpl<7 ,7 ,14 ,13 ,true>(in[0], in[1], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTri: # of modes / points combo not implemented.");
                    } break;
                case 8:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 9: IProductWRTDerivBaseTriImpl<8 ,8 ,9 ,8 ,true>(in[0], in[1], out); break;
                        case 10: IProductWRTDerivBaseTriImpl<8 ,8 ,10 ,9 ,true>(in[0], in[1], out); break;
                        case 11: IProductWRTDerivBaseTriImpl<8 ,8 ,11 ,10 ,true>(in[0], in[1], out); break;
                        case 12: IProductWRTDerivBaseTriImpl<8 ,8 ,12 ,11 ,true>(in[0], in[1], out); break;
                        case 13: IProductWRTDerivBaseTriImpl<8 ,8 ,13 ,12 ,true>(in[0], in[1], out); break;
                        case 14: IProductWRTDerivBaseTriImpl<8 ,8 ,14 ,13 ,true>(in[0], in[1], out); break;
                        case 15: IProductWRTDerivBaseTriImpl<8 ,8 ,15 ,14 ,true>(in[0], in[1], out); break;
                        case 16: IProductWRTDerivBaseTriImpl<8 ,8 ,16 ,15 ,true>(in[0], in[1], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTri: # of modes / points combo not implemented.");
                    } break;
                default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTri: # of modes / points combo not implemented.");
            }
        }
    }


    template<int NM0, int NM1, int NQ0, int NQ1, bool CORRECT>
    void IProductWRTDerivBaseTriImpl(
        const Array<OneD, const NekDouble> &input0,
        const Array<OneD, const NekDouble> &input1,
              Array<OneD,       NekDouble> &output)
    {
        auto* inptr0 = input0.data();
        auto* inptr1 = input1.data();
        auto* outptr = output.data();

        constexpr auto ndf = 4u;
        constexpr auto nqTot = NQ0 * NQ1;
        constexpr auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;


        // Get size of jacobian factor block
        auto dJSize = 1u;
        auto dfSize = ndf;
        if(DEFORMED)
        {
            dJSize = nqTot;
            dfSize = ndf*nqTot;
        }

        vec_t sums_j[NQ1]; //Sums over eta0 for each value of eta1;
        std::vector<vec_t, allocator<vec_t>> tmpIn0(nqTot), tmpIn1(nqTot),
            tmp0(nqTot), tmp1(nqTot), tmpOut(m_nmTot);

        const vec_t* df_ptr;
        const vec_t* jac_ptr;

        std::vector<vec_t, allocator<vec_t>> Z0 = this->m_Z[0];
        std::vector<vec_t, allocator<vec_t>> Z1 = this->m_Z[1];

        for (int e =0; e < this->m_nBlocks; ++e)
        {
            // Jacobian
            jac_ptr = &(this->m_jac[dJSize*e]);

            // Derivative factor
            df_ptr = &(this->m_df[dfSize*e]);

            // Load and transpose data
            load_interleave(inptr0, nqTot, tmpIn0);
            load_interleave(inptr1, nqTot, tmpIn1);

            vec_t df0, df1, df2, df3;
            if (!DEFORMED)
            {
                df0 = df_ptr[0];
                df1 = df_ptr[1];
                df2 = df_ptr[2];
                df3 = df_ptr[3];
            }

            size_t cnt_ji = 0;
            for (size_t j = 0; j < NQ1; ++j)
            {
                vec_t f0 = 2.0 / (1.0 - Z1[j]); //Load 1x

                for (size_t i = 0; i < NQ0; ++i, ++cnt_ji)
                {
                    if (DEFORMED)
                    {
                        df0 = df_ptr[cnt_ji * ndf];
                        df1 = df_ptr[cnt_ji * ndf + 1];
                        df2 = df_ptr[cnt_ji * ndf + 2];
                        df3 = df_ptr[cnt_ji * ndf + 3];
                    }

                    // Calculate dx/dxi in[0] + dy/dxi in[1]
                    vec_t t0 = df0 * tmpIn0[cnt_ji] + df2 * tmpIn1[cnt_ji]; // load 2x
                    vec_t t1 = df1 * tmpIn0[cnt_ji] + df3 * tmpIn1[cnt_ji]; // load 2x

                    // Multiply by geometric factors
                    vec_t hf1 = 0.5 * (1.0 + Z0[i]); //Load 1x
                    // Scale by geometric factor 2/(1-z1)
                    t0 *= f0;
                    // Scale by geometric factor (1+z0)/(1-z1)
                    vec_t c1 = hf1 * t1;
                    t0.fma(c1, f0);

                    // Store
                    tmp0[cnt_ji] = t0; // store 1x
                    tmp1[cnt_ji] = t1; // store 1x
                }
            }



            // IP DB0 B1
            IProductTriKernel<NM0, NM1, NQ0, NQ1, CORRECT, false, false, DEFORMED>(
                tmp0, this->m_dbdata[0], this->m_bdata[1],
                this->m_w[0], this->m_w[1], jac_ptr,
                sums_j, tmpOut);

            // IP DB1 B0
            IProductTriKernel<NM0, NM1, NQ0, NQ1, CORRECT, false, true, DEFORMED>(
                tmp1, this->m_bdata[0], this->m_dbdata[1],
                this->m_w[0], this->m_w[1], jac_ptr,
                sums_j, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr0 += nqBlocks;
            inptr1 += nqBlocks;
            outptr += nmBlocks;
        }
    }

private:
    int m_nmTot;
};

template<bool DEFORMED = false>
struct IProductWRTDerivBaseHex : public IProductWRTDerivBase, public Helper<3, DEFORMED>
{
public:
    IProductWRTDerivBaseHex(std::vector<LibUtilities::BasisSharedPtr> basis,
                   int nElmt)
        : IProductWRTDerivBase(basis, nElmt),
          Helper<3, DEFORMED>(basis, nElmt),
          m_nmTot(LibUtilities::StdHexData::getNumberOfCoefficients(
                      this->m_nm[0], this->m_nm[1], this->m_nm[2]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<IProductWRTDerivBaseHex<DEFORMED>>(basis, nElmt);
    }

    void operator()(const Array<OneD, Array<OneD, NekDouble>> &in,
        Array<OneD,       NekDouble> &out) final
    {
        switch(m_basis[0]->GetNumModes())
        {
            case 2:  IProductWRTDerivBaseHexImpl<2 ,2 ,2 ,3 ,3 ,3 >
                (in[0], in[1], in[2], out); break;
            case 3:  IProductWRTDerivBaseHexImpl<3 ,3 ,3 ,4 ,4 ,4 >
                (in[0], in[1], in[2], out); break;
            case 4:  IProductWRTDerivBaseHexImpl<4 ,4 ,4 ,5 ,5 ,5 >
                (in[0], in[1], in[2], out); break;
            case 5:  IProductWRTDerivBaseHexImpl<5 ,5 ,5 ,6 ,6 ,6 >
                (in[0], in[1], in[2], out); break;
            case 6:  IProductWRTDerivBaseHexImpl<6 ,6 ,6 ,7 ,7 ,7 >
                (in[0], in[1], in[2], out); break;
            case 7:  IProductWRTDerivBaseHexImpl<7 ,7 ,7 ,8 ,8 ,8 >
                (in[0], in[1], in[2], out); break;
            case 8:  IProductWRTDerivBaseHexImpl<8 ,8 ,8 ,9 ,9 ,9 >
                (in[0], in[1], in[2], out); break;
            case 9:  IProductWRTDerivBaseHexImpl<9 ,9 ,9 ,10,10,10>
                (in[0], in[1], in[2], out); break;
            case 10: IProductWRTDerivBaseHexImpl<10,10,10,11,11,11>
                (in[0], in[1], in[2], out); break;
            case 11: IProductWRTDerivBaseHexImpl<11,11,11,12,12,12>
                (in[0], in[1], in[2], out); break;
        }
    }

    template<int NM0, int NM1, int NM2, int NQ0, int NQ1, int NQ2>
    void IProductWRTDerivBaseHexImpl(
        const Array<OneD, NekDouble> &input0,
        const Array<OneD, NekDouble> &input1,
        const Array<OneD, NekDouble> &input2,
              Array<OneD, NekDouble> &output)
    {
        using namespace tinysimd;
        using vec_t = simd<NekDouble>;

        auto* inptr0 = input0.data();
        auto* inptr1 = input1.data();
        auto* inptr2 = input2.data();
        auto* outptr = output.data();

        constexpr auto ndf = 9;
        constexpr auto nqTot = NQ0 * NQ1 * NQ2;
        constexpr auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        // Get size of jacobian factor block
        int dJSize{}, dfSize{};
        if(DEFORMED)
        {
            dJSize = nqTot;
            dfSize = ndf*nqTot;
        }
        else
        {
            dJSize = 1;
            dfSize = ndf;
        }

        vec_t sums_kj[NQ1 * NQ2];
        vec_t sums_k[NQ2];

        std::vector<vec_t, allocator<vec_t>> tmpIn0(nqTot), tmpIn1(nqTot),
            tmpIn2(nqTot),  tmp0(nqTot), tmp1(nqTot), tmp2(nqTot),
            tmpOut(m_nmTot);

        const vec_t* df_ptr;
        const vec_t* jac_ptr;
        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            // Jacobian
            jac_ptr = &(this->m_jac[dJSize*e]);

            // Derivative factor
            df_ptr = &(this->m_df[e*dfSize]);

            // Load and transpose data
            load_interleave(inptr0, nqTot, tmpIn0);
            load_interleave(inptr1, nqTot, tmpIn1);
            load_interleave(inptr2, nqTot, tmpIn2);

            // Calculate dx/dxi in[0] + dy/dxi in[1] + dz/dxi in[2]
            vec_t df0, df1, df2, df3, df4, df5, df6, df7, df8;
            if(!DEFORMED)
            {
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
            for (int i = 0; i < nqTot; ++i)
            {
                if(DEFORMED)
                {
                    df0 = df_ptr[i * ndf];
                    df1 = df_ptr[i * ndf + 1];
                    df2 = df_ptr[i * ndf + 2];
                    df3 = df_ptr[i * ndf + 3];
                    df4 = df_ptr[i * ndf + 4];
                    df5 = df_ptr[i * ndf + 5];
                    df6 = df_ptr[i * ndf + 6];
                    df7 = df_ptr[i * ndf + 7];
                    df8 = df_ptr[i * ndf + 8];
                }
                tmp0[i] = df0 * tmpIn0[i] + df3 * tmpIn1[i] + df6 *tmpIn2[i];
                tmp1[i] = df1 * tmpIn0[i] + df4 * tmpIn1[i] + df7 *tmpIn2[i];
                tmp2[i] = df2 * tmpIn0[i] + df5 * tmpIn1[i] + df8 *tmpIn2[i];
            }

            // IP DB0 B1 B2
            IProductHexKernel
            <NM0, NM1, NM2, NQ0, NQ1, NQ2, false, false, DEFORMED>(
                tmp0, this->m_dbdata[0], this->m_bdata[1], this->m_bdata[2],
                this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                sums_kj, sums_k,
                tmpOut);

            // IP B0 DB1 B2
            IProductHexKernel
            <NM0, NM1, NM2, NQ0, NQ1, NQ2, false, true, DEFORMED>(
                tmp1, this->m_bdata[0], this->m_dbdata[1], this->m_bdata[2],
                this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                sums_kj, sums_k,
                tmpOut);

            // IP B0 B1 DB2
            IProductHexKernel
            <NM0, NM1, NM2, NQ0, NQ1, NQ2, false, true, DEFORMED>(
                tmp2, this->m_bdata[0], this->m_bdata[1], this->m_dbdata[2],
                this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                sums_kj, sums_k,
                tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr0 += nqBlocks;
            inptr1 += nqBlocks;
            inptr2 += nqBlocks;
            outptr += nmBlocks;
        }
    }



private:
    /// Padded basis
    int m_nmTot;
};

// template<int VW, bool DEFORMED = false>
// struct ProductWRTDerivBasePrism : public IProduct, public Helper<VW, 3, DEFORMED>
// {
//     AVXIProductPrism(std::vector<LibUtilities::BasisSharedPtr> basis,
//                      int nElmt)
//         : IProduct(basis, nElmt),
//           Helper<VW, 3, DEFORMED>(basis, nElmt),
//           m_nmTot(LibUtilities::StdPrismData::getNumberOfCoefficients(
//                       this->m_nm[0], this->m_nm[1], this->m_nm[2]))
//     {
//     }

//     static std::shared_ptr<Operator> Create(
//         std::vector<LibUtilities::BasisSharedPtr> basis,
//         int nElmt)
//     {
//         return std::make_shared<AVXIProductPrism<VW, DEFORMED>>(basis, nElmt);
//     }

//     virtual void operator()(const Array<OneD, const NekDouble> &in,
//                                   Array<OneD,       NekDouble> &out)
//     {
//         if (m_basis[0]->GetBasisType() == LibUtilities::eModified_A)
//         {
//             switch(m_basis[0]->GetNumModes())
//             {
//                 case 2:  AVXIProductPrismImpl<2 ,2 ,2 ,3 ,3 ,2 ,true>(in, out); break;
//                 case 3:  AVXIProductPrismImpl<3 ,3 ,3 ,4 ,4 ,3 ,true>(in, out); break;
//                 case 4:  AVXIProductPrismImpl<4 ,4 ,4 ,5 ,5 ,4 ,true>(in, out); break;
//                 case 5:  AVXIProductPrismImpl<5 ,5 ,5 ,6 ,6 ,5 ,true>(in, out); break;
//                 case 6:  AVXIProductPrismImpl<6 ,6 ,6 ,7 ,7 ,6 ,true>(in, out); break;
//                 case 7:  AVXIProductPrismImpl<7 ,7 ,7 ,8 ,8 ,7 ,true>(in, out); break;
//                 case 8:  AVXIProductPrismImpl<8 ,8 ,8 ,9 ,9 ,8 ,true>(in, out); break;
//                 case 9:  AVXIProductPrismImpl<9 ,9 ,9 ,10,10,9 ,true>(in, out); break;
//                 case 10: AVXIProductPrismImpl<10,10,10,11,11,10,true>(in, out); break;
//                 case 11: AVXIProductPrismImpl<11,11,11,12,12,11,true>(in, out); break;
//             }
//         }
//         else
//         {
//             switch(m_basis[0]->GetNumModes())
//             {
//                 case 2:  AVXIProductPrismImpl<2 ,2 ,2 ,3 ,3 ,2 ,false>(in, out); break;
//                 case 3:  AVXIProductPrismImpl<3 ,3 ,3 ,4 ,4 ,3 ,false>(in, out); break;
//                 case 4:  AVXIProductPrismImpl<4 ,4 ,4 ,5 ,5 ,4 ,false>(in, out); break;
//                 case 5:  AVXIProductPrismImpl<5 ,5 ,5 ,6 ,6 ,5 ,false>(in, out); break;
//                 case 6:  AVXIProductPrismImpl<6 ,6 ,6 ,7 ,7 ,6 ,false>(in, out); break;
//                 case 7:  AVXIProductPrismImpl<7 ,7 ,7 ,8 ,8 ,7 ,false>(in, out); break;
//                 case 8:  AVXIProductPrismImpl<8 ,8 ,8 ,9 ,9 ,8 ,false>(in, out); break;
//                 case 9:  AVXIProductPrismImpl<9 ,9 ,9 ,10,10,9 ,false>(in, out); break;
//                 case 10: AVXIProductPrismImpl<10,10,10,11,11,10,false>(in, out); break;
//                 case 11: AVXIProductPrismImpl<11,11,11,12,12,11,false>(in, out); break;
//             }
//         }
//     }

//     template<int NM0, int NM1, int NM2, int NQ0, int NQ1, int NQ2, bool CORRECT>
//     void AVXIProductPrismImpl(
//         const Array<OneD, const NekDouble> &input,
//               Array<OneD,       NekDouble> &output)
//     {
//         using T = VecData<double, VW>;
//         auto *inptr = &input[0];
//         auto *outptr = &output[0];

//         constexpr int nqBlocks = NQ0 * NQ1 * NQ2 * VW;
//         const int nmBlocks = m_nmTot * VW;

//         T sums_kj[NQ1 * NQ2];
//         T sums_k[NQ2];

//         T corr_q[NM1];

//         for (int e = 0; e < this->m_nBlocks; ++e)
//         {
//             VecData<NekDouble, VW> *jac_ptr;
//             if(DEFORMED){
//                 jac_ptr = &(this->m_jac[NQ0*NQ1*NQ2*e]);
//             }
//             else{
//                 jac_ptr = &(this->m_jac[e]);
//             }

//             AVXIProductPrismKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, VW, CORRECT,
//                                    false, false, DEFORMED>(
//                 inptr, this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
//                 this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
//                 sums_kj, sums_k,
//                 corr_q,
//                 outptr);

//             inptr += nqBlocks;
//             outptr += nmBlocks;
//         }

//     }
// public:

//     static NekDouble FlopsPerElement(
//         const int nm,
//         const int nq0,
//         const int nq1,
//         const int nq2)
//     {
//         int loop_ijk = nm * nq0*nq1*nq2*4;
//         int loop_kj = nm*nm*nq1*nq2*3;
//         int loop_k = nm*nq2* 3 * nm*(nm+1) / 2.0;

//         int corr;
//         if(DEFORMED){
//             int corr_inner = 3 + nm*4;
//             corr = nq2*nq1*(1 + nq0*corr_inner);
//         }
//         else{
//             int corr_inner = 2 + nm*4;
//             corr = nq2*(1 + nq1*(1 + nq0*corr_inner));
//         }

//         return  (loop_ijk + loop_kj + loop_k + corr);
//     }

//     virtual NekDouble GFlops() override
//     {
//         const int nm = m_basis[0]->GetNumModes();
//         const int nq0 = m_basis[0]->GetNumPoints();
//         const int nq1 = m_basis[1]->GetNumPoints();
//         const int nq2 = m_basis[2]->GetNumPoints();

//         int flops = m_nElmt * AVXIProductPrism::FlopsPerElement(nm, nq0, nq1, nq2);
//         return flops * 1e-9;
//     }

//     virtual NekDouble NLoads() override
//     {
//         const int nm0 = m_basis[0]->GetNumModes();
//         const int nm1 = m_basis[1]->GetNumModes();
//         const int nm2 = m_basis[2]->GetNumModes();
//         const int nq0 = m_basis[0]->GetNumPoints();
//         const int nq1 = m_basis[1]->GetNumPoints();
//         const int nq2 = m_basis[2]->GetNumPoints();

//         const int load_pkji = nm0 * nq2 * nq1 * nq0 * 3;
//         const int load_pqkj = nm0 * nm1 * nq2 * nq1 * 3;
//         const int load_pqrk = (nm0 *(nm0+1)/2)* nm1 * nq2 * 3;
//         const int load_corr_setup = nm1;
//         const int load_corr = nq2*(1 + nq1*(1 + nq0*(4 + nm1)));
//         const int load_corr_finish = nm1;
//         const int load_expected = load_pkji + load_pqkj + load_pqrk + load_corr_setup + load_corr + load_corr_finish;

//         return load_expected * this->m_nElmt;
//     }

//     virtual NekDouble NStores() override
//     {
//         const int nm0 = m_basis[0]->GetNumModes();
//         const int nm1 = m_basis[1]->GetNumModes();
//         const int nm2 = m_basis[2]->GetNumModes();
//         const int nq0 = m_basis[0]->GetNumPoints();
//         const int nq1 = m_basis[1]->GetNumPoints();
//         const int nq2 = m_basis[2]->GetNumPoints();

//         const int store_pkj = nm0 * nq2 * nq1;
//         const int store_pqk = nm0 * nm1 * nq2;
//         const int store_pqr = (nm0*(nm0+1) / 2) * nm1;
//         const int store_corr_start_finish = nm1*2;
//         const int store_expected = store_pkj + store_pqk + store_pqr + store_corr_start_finish;

//         return store_expected * this->m_nElmt;
//     }

//     virtual NekDouble Ndof() override
//     {
//         return m_nmTot * this->m_nElmt;
//     }

// private:
//     /// Padded basis
//     int m_nmTot;
// };

// template<int VW, bool DEFORMED = false>
// struct AVXIProductTet : public IProduct, public Helper<VW, 3, DEFORMED>
// {
//     AVXIProductTet(std::vector<LibUtilities::BasisSharedPtr> basis,
//                    int nElmt)
//         : IProduct(basis, nElmt),
//           Helper<VW, 3, DEFORMED>(basis, nElmt),
//           m_nmTot(LibUtilities::StdTetData::getNumberOfCoefficients(
//                       this->m_nm[0], this->m_nm[1], this->m_nm[2]))
//     {
//     }

//     static std::shared_ptr<Operator> Create(
//         std::vector<LibUtilities::BasisSharedPtr> basis,
//         int nElmt)
//     {
//         return std::make_shared<AVXIProductTet<VW, DEFORMED>>(basis, nElmt);
//     }

//     virtual void operator()(const Array<OneD, const NekDouble> &in,
//                                   Array<OneD,       NekDouble> &out)
//     {
//         if (m_basis[0]->GetBasisType() == LibUtilities::eModified_A)
//         {
//             switch(m_basis[0]->GetNumModes())
//             {
//                 case 2:  AVXIProductTetImpl<2 ,2 ,2 ,3 ,2 ,2 ,true>(in, out); break;
//                 case 3:  AVXIProductTetImpl<3 ,3 ,3 ,4 ,3 ,3 ,true>(in, out); break;
//                 case 4:  AVXIProductTetImpl<4 ,4 ,4 ,5 ,4 ,4 ,true>(in, out); break;
//                 case 5:  AVXIProductTetImpl<5 ,5 ,5 ,6 ,5 ,5 ,true>(in, out); break;
//                 case 6:  AVXIProductTetImpl<6 ,6 ,6 ,7 ,6 ,6 ,true>(in, out); break;
//                 case 7:  AVXIProductTetImpl<7 ,7 ,7 ,8 ,7 ,7 ,true>(in, out); break;
//                 case 8:  AVXIProductTetImpl<8 ,8 ,8 ,9 ,8 ,8 ,true>(in, out); break;
//                 case 9:  AVXIProductTetImpl<9 ,9 ,9 ,10,9 ,9 ,true>(in, out); break;
//                 case 10: AVXIProductTetImpl<10,10,10,11,10,10,true>(in, out); break;
//                 case 11: AVXIProductTetImpl<11,11,11,12,11,11,true>(in, out); break;
//             }
//         }
//         else
//         {
//             switch(m_basis[0]->GetNumModes())
//             {
//                 case 2:  AVXIProductTetImpl<2 ,2 ,2 ,3 ,2 ,2 ,false>(in, out); break;
//                 case 3:  AVXIProductTetImpl<3 ,3 ,3 ,4 ,3 ,3 ,false>(in, out); break;
//                 case 4:  AVXIProductTetImpl<4 ,4 ,4 ,5 ,4 ,4 ,false>(in, out); break;
//                 case 5:  AVXIProductTetImpl<5 ,5 ,5 ,6 ,5 ,5 ,false>(in, out); break;
//                 case 6:  AVXIProductTetImpl<6 ,6 ,6 ,7 ,6 ,6 ,false>(in, out); break;
//                 case 7:  AVXIProductTetImpl<7 ,7 ,7 ,8 ,7 ,7 ,false>(in, out); break;
//                 case 8:  AVXIProductTetImpl<8 ,8 ,8 ,9 ,8 ,8 ,false>(in, out); break;
//                 case 9:  AVXIProductTetImpl<9 ,9 ,9 ,10,9 ,9 ,false>(in, out); break;
//                 case 10: AVXIProductTetImpl<10,10,10,11,10,10,false>(in, out); break;
//                 case 11: AVXIProductTetImpl<11,11,11,12,11,11,false>(in, out); break;
//             }
//         }
//     }

//     template<int NM0, int NM1, int NM2, int NQ0, int NQ1, int NQ2, bool CORRECT>
//     void AVXIProductTetImpl(
//         const Array<OneD, const NekDouble> &input,
//               Array<OneD,       NekDouble> &output)
//     {
//         using T = VecData<double, VW>;
//         auto *inptr = &input[0];
//         auto *outptr = &output[0];

//         constexpr int nqBlocks = NQ0 * NQ1 * NQ2 * VW;
//         const int nmBlocks = m_nmTot * VW;

//         T wsp[NQ1 * NQ2 + NQ2];

//         for (int e = 0; e < this->m_nBlocks; ++e)
//         {
//             VecData<NekDouble, VW> *jac_ptr;
//             if(DEFORMED){
//                 jac_ptr = &(this->m_jac[NQ0*NQ1*NQ2*e]);
//             }
//             else{
//                 jac_ptr = &(this->m_jac[e]);
//             }
//             AVXIProductTetKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, VW, CORRECT,
//                                  false, false, DEFORMED>(
//                 inptr, this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
//                 this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
//                 wsp, outptr);

//             inptr += nqBlocks;
//             outptr += nmBlocks;
//         }
//     }
// public:

//     static NekDouble FlopsPerElement(
//         const int nm,
//         const int nq0,
//         const int nq1,
//         const int nq2)
//     {
//         double loop_p = nm * nq2 * nq1 * nq0 * 4;
//         double loop_pq = nm * (nm + 1) * nq2*nq1 * 3 / 2;
//         double loop_pqr = nm * (nm + 1) * (nm + 2) * nq2 * 3 / 6;

//         double loop_ijk;
//         if(DEFORMED){
//             double inner = 11 + 5 * (nm - 1);
//             loop_ijk = nq2 * nq1 * (1 + nq0 * inner);
//         }
//         else{
//             double inner = 10 + 5 * (nm - 1);
//             loop_ijk = nq2 * (1 + nq1 * (1 + nq0 * inner));
//         }

//         return (loop_p + loop_pq + loop_pqr + loop_ijk);
//     }

//     virtual NekDouble GFlops() override
//     {
//         const int nm = m_basis[0]->GetNumModes();
//         const int nq0 = m_basis[0]->GetNumPoints();
//         const int nq1 = m_basis[1]->GetNumPoints();
//         const int nq2 = m_basis[2]->GetNumPoints();

//         double flops = m_nElmt * AVXIProductTet::FlopsPerElement(nm, nq0, nq1, nq2);
//         return flops * 1e-9;
//     }

//     virtual NekDouble NLoads() override
//     {
//         const int nm0 = m_basis[0]->GetNumModes();
//         const int nm1 = m_basis[1]->GetNumModes();
//         const int nm2 = m_basis[2]->GetNumModes();
//         const int nq0 = m_basis[0]->GetNumPoints();
//         const int nq1 = m_basis[1]->GetNumPoints();
//         const int nq2 = m_basis[2]->GetNumPoints();

//         int load_pkji = nm0 * nq2 * nq1 * nq0 * 3;
//         int load_pqkj = (nm0*(nm0+1)/2) * nq2 * nq1 * 3;
//         int load_pqrk = (nm0*(nm0+1)*(nm0+2)/6)*nq2 * 3;
//         int load_corr_inner = 9 + (nm2-1)*3;
//         int load_corr = nq2 * (1 + nq1 * (1 + nq0*load_corr_inner));
//         int load_expected = load_pkji + load_pqkj + load_pqrk + load_corr;

//         return load_expected * this->m_nElmt;
//     }

//     virtual NekDouble NStores() override
//     {
//         const int nm0 = m_basis[0]->GetNumModes();
//         const int nm1 = m_basis[1]->GetNumModes();
//         const int nm2 = m_basis[2]->GetNumModes();
//         const int nq0 = m_basis[0]->GetNumPoints();
//         const int nq1 = m_basis[1]->GetNumPoints();
//         const int nq2 = m_basis[2]->GetNumPoints();

//         int store_pkj = nm0 * nq2 * nq1;
//         int store_pqk = (nm0*(nm0+1)/2) * nq2;
//         int store_pqr = (nm0*(nm0+1)*(nm0+2)/6);
//         int store_corr_inner = 1 + (nm2-1);
//         int store_corr = nq2 * nq1 * nq0 * store_corr_inner;
//         int store_expected =  store_pkj + store_pqk + store_pqr + store_corr;

//         return store_expected * this->m_nElmt;
//     }

//     virtual NekDouble Ndof() override
//     {
//         return m_nmTot * this->m_nElmt;
//     }

// private:
//     /// Padded basis
//     int m_nmTot;
// };

} // namespace MatrixFree
} // namespace Nektar

#endif
