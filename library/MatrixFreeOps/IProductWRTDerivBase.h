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
        // Check preconditions
        ASSERTL0(m_basis[0]->GetNumModes() == m_basis[1]->GetNumModes() &&
            m_basis[0]->GetNumPoints() == m_basis[1]->GetNumPoints(),
            "MatrixFree requires homogenous modes/points");

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

        std::vector<vec_t, allocator<vec_t>>& Z0 = this->m_Z[0];
        std::vector<vec_t, allocator<vec_t>>& Z1 = this->m_Z[1];

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
                    vec_t tI0 = tmpIn0[cnt_ji]; // load 1x
                    vec_t tI1 = tmpIn1[cnt_ji]; // load 1x
                    vec_t t0 = df0 * tI0 + df2 * tI1;
                    vec_t t1 = df1 * tI0 + df3 * tI1;

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
        auto nm0 = m_basis[0]->GetNumModes();
        auto nm1 = m_basis[1]->GetNumModes();
        auto nm2 = m_basis[2]->GetNumModes();
        ASSERTL0( nm0 == nm1 && nm0 == nm2,
            "IProductWRTDBHex: anisotropy not implemented.");

        auto np0 = m_basis[0]->GetNumPoints();

        switch(nm0)
        {
            case 2:
                switch(np0)
                {
                    case 2: IProductWRTDerivBaseHexImpl<2, 2, 2, 2, 2, 2>
                        (in[0], in[1], in[2], out); break;
                    case 3: IProductWRTDerivBaseHexImpl<2, 2, 2, 3, 3, 3>
                        (in[0], in[1], in[2], out); break;
                    case 4: IProductWRTDerivBaseHexImpl<2, 2, 2, 4, 4, 4>
                        (in[0], in[1], in[2], out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductWRTDBHex: # of modes / points combo not implemented.");
                } break;
            case 3:
                switch(np0)
                {
                    case 3: IProductWRTDerivBaseHexImpl<3, 3, 3, 3, 3, 3>
                        (in[0], in[1], in[2], out); break;
                    case 4: IProductWRTDerivBaseHexImpl<3, 3, 3, 4, 4, 4>
                        (in[0], in[1], in[2], out); break;
                    case 5: IProductWRTDerivBaseHexImpl<3, 3, 3, 5, 5, 5>
                        (in[0], in[1], in[2], out); break;
                    case 6: IProductWRTDerivBaseHexImpl<3, 3, 3, 6, 6, 6>
                        (in[0], in[1], in[2], out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductWRTDBHex: # of modes / points combo not implemented.");
                } break;
            case 4:
                switch(np0)
                {
                    case 4: IProductWRTDerivBaseHexImpl<4, 4, 4, 4, 4 ,4>
                        (in[0], in[1], in[2], out); break;
                    case 5: IProductWRTDerivBaseHexImpl<4, 4, 4, 5, 5 ,5>
                        (in[0], in[1], in[2], out); break;
                    case 6: IProductWRTDerivBaseHexImpl<4, 4, 4, 6, 6 ,6>
                        (in[0], in[1], in[2], out); break;
                    case 7: IProductWRTDerivBaseHexImpl<4, 4, 4, 7, 7 ,7>
                        (in[0], in[1], in[2], out); break;
                    case 8: IProductWRTDerivBaseHexImpl<4, 4, 4, 8, 8 ,8>
                        (in[0], in[1], in[2], out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductWRTDBHex: # of modes / points combo not implemented.");
                } break;
            case 5:
                switch(np0)
                {
                    case 5: IProductWRTDerivBaseHexImpl<5, 5, 5, 5, 5, 5>
                        (in[0], in[1], in[2], out); break;
                    case 6: IProductWRTDerivBaseHexImpl<5, 5, 5, 6, 6, 6>
                        (in[0], in[1], in[2], out); break;
                    case 7: IProductWRTDerivBaseHexImpl<5, 5, 5, 7, 7, 7>
                        (in[0], in[1], in[2], out); break;
                    case 8: IProductWRTDerivBaseHexImpl<5, 5, 5, 8, 8, 8>
                        (in[0], in[1], in[2], out); break;
                    case 9: IProductWRTDerivBaseHexImpl<5, 5, 5, 9, 9, 9>
                        (in[0], in[1], in[2], out); break;
                    case 10: IProductWRTDerivBaseHexImpl<5, 5, 5, 10, 10, 10>
                        (in[0], in[1], in[2], out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductWRTDBHex: # of modes / points combo not implemented.");
                } break;
            case 6:
                switch(np0)
                {
                    case 6: IProductWRTDerivBaseHexImpl<6, 6, 6, 6, 6, 6>
                        (in[0], in[1], in[2], out); break;
                    case 7: IProductWRTDerivBaseHexImpl<6, 6, 6, 7, 7, 7>
                        (in[0], in[1], in[2], out); break;
                    case 8: IProductWRTDerivBaseHexImpl<6, 6, 6, 8, 8, 8>
                        (in[0], in[1], in[2], out); break;
                    case 9: IProductWRTDerivBaseHexImpl<6, 6, 6, 9, 9, 9>
                        (in[0], in[1], in[2], out); break;
                    case 10: IProductWRTDerivBaseHexImpl<6, 6, 6, 10, 10, 10>
                        (in[0], in[1], in[2], out); break;
                    case 11: IProductWRTDerivBaseHexImpl<6, 6, 6, 11, 11, 11>
                        (in[0], in[1], in[2], out); break;
                    case 12: IProductWRTDerivBaseHexImpl<6, 6, 6, 12, 12, 12>
                        (in[0], in[1], in[2], out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductWRTDBHex: # of modes / points combo not implemented.");
                } break;
            case 7:
                switch(np0)
                {
                    case 7: IProductWRTDerivBaseHexImpl<7, 7, 7, 7, 7, 7>
                        (in[0], in[1], in[2], out); break;
                    case 8: IProductWRTDerivBaseHexImpl<7, 7, 7, 8, 8, 8>
                        (in[0], in[1], in[2], out); break;
                    case 9: IProductWRTDerivBaseHexImpl<7, 7, 7, 9, 9, 9>
                        (in[0], in[1], in[2], out); break;
                    case 10: IProductWRTDerivBaseHexImpl<7, 7, 7, 10, 10, 10>
                        (in[0], in[1], in[2], out); break;
                    case 11: IProductWRTDerivBaseHexImpl<7, 7, 7, 11, 11, 11>
                        (in[0], in[1], in[2], out); break;
                    case 12: IProductWRTDerivBaseHexImpl<7, 7, 7, 12, 12, 12>
                        (in[0], in[1], in[2], out); break;
                    case 13: IProductWRTDerivBaseHexImpl<7, 7, 7, 13, 13, 13>
                        (in[0], in[1], in[2], out); break;
                    case 14: IProductWRTDerivBaseHexImpl<7, 7, 7, 14, 14, 14>
                        (in[0], in[1], in[2], out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductWRTDBHex: # of modes / points combo not implemented.");
                } break;
            case 8:
                switch(np0)
                {
                    case 8: IProductWRTDerivBaseHexImpl<8, 8, 8, 8, 8, 8>
                        (in[0], in[1], in[2], out); break;
                    case 9: IProductWRTDerivBaseHexImpl<8, 8, 8, 9, 9, 9>
                        (in[0], in[1], in[2], out); break;
                    case 10: IProductWRTDerivBaseHexImpl<8, 8, 8, 10, 10, 10>
                        (in[0], in[1], in[2], out); break;
                    case 11: IProductWRTDerivBaseHexImpl<8, 8, 8, 11, 11, 11>
                        (in[0], in[1], in[2], out); break;
                    case 12: IProductWRTDerivBaseHexImpl<8, 8, 8, 12, 12, 12>
                        (in[0], in[1], in[2], out); break;
                    case 13: IProductWRTDerivBaseHexImpl<8, 8, 8, 13, 13, 13>
                        (in[0], in[1], in[2], out); break;
                    case 14: IProductWRTDerivBaseHexImpl<8, 8, 8, 14, 14, 14>
                        (in[0], in[1], in[2], out); break;
                    case 15: IProductWRTDerivBaseHexImpl<8, 8, 8, 15, 15, 15>
                        (in[0], in[1], in[2], out); break;
                    case 16: IProductWRTDerivBaseHexImpl<8, 8, 8, 16, 16, 16>
                        (in[0], in[1], in[2], out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductWRTDBHex: # of modes / points combo not implemented.");
                } break;;
            default: NEKERROR(ErrorUtil::efatal,
                "IProductWRTDBHex: # of modes / points combo not implemented.");
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

        constexpr auto ndf = 9u;
        constexpr auto nqTot = NQ0 * NQ1 * NQ2;
        constexpr auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        // Get size of jacobian factor block
        auto dJSize = 1u;
        auto dfSize = ndf;
        if(DEFORMED)
        {
            dJSize *= nqTot;
            dfSize *= nqTot;
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
            if (!DEFORMED)
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
                if (DEFORMED)
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
                vec_t tI0 = tmpIn0[i]; // load 1x
                vec_t tI1 = tmpIn1[i]; // load 1x
                vec_t tI2 = tmpIn2[i]; // load 1x
                tmp0[i] = df0 * tI0 + df3 * tI1 + df6 * tI2;
                tmp1[i] = df1 * tI0 + df4 * tI1 + df7 * tI2;
                tmp2[i] = df2 * tI0 + df5 * tI1 + df8 * tI2;
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

template<bool DEFORMED = false>
struct IProductWRTDerivBasePrism : public IProductWRTDerivBase, public Helper<3, DEFORMED>
{
    IProductWRTDerivBasePrism(std::vector<LibUtilities::BasisSharedPtr> basis,
                   int nElmt)
        : IProductWRTDerivBase(basis, nElmt),
          Helper<3, DEFORMED>(basis, nElmt),
          m_nmTot(LibUtilities::StdPrismData::getNumberOfCoefficients(
                      this->m_nm[0], this->m_nm[1], this->m_nm[2]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<IProductWRTDerivBasePrism<DEFORMED>>(basis, nElmt);
    }

   void operator()( const Array<OneD, Array<OneD, NekDouble>> &in,
        Array<OneD, NekDouble> &out) final
    {
        // Check preconditions
        ASSERTL0(m_basis[0]->GetNumModes() == m_basis[1]->GetNumModes() &&
            m_basis[0]->GetNumModes() == m_basis[2]->GetNumModes() &&
            m_basis[0]->GetNumPoints() == m_basis[1]->GetNumPoints() &&
            m_basis[0]->GetNumPoints() == m_basis[2]->GetNumPoints()+1,
            "MatrixFree requires homogenous modes/points");

        if (m_basis[0]->GetBasisType() == LibUtilities::eModified_A)
        {
            switch(m_basis[0]->GetNumModes())
            {
                case 2:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 3: IProductWRTDerivBasePrismImpl<2, 2, 2, 3, 3, 2, true>
                            (in[0], in[1], in[2], out); break;
                        case 4: IProductWRTDerivBasePrismImpl<2, 2, 2, 4, 4, 3, true>
                            (in[0], in[1], in[2], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBPrism: # of modes / points combo not implemented.");
                    } break;
                case 3:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 4: IProductWRTDerivBasePrismImpl<3, 3, 3, 4, 4, 3, true>
                            (in[0], in[1], in[2], out); break;
                        case 5: IProductWRTDerivBasePrismImpl<3, 3, 3, 5, 5, 4, true>
                            (in[0], in[1], in[2], out); break;
                        case 6: IProductWRTDerivBasePrismImpl<3, 3, 3, 6, 6, 5, true>
                            (in[0], in[1], in[2], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBPrism: # of modes / points combo not implemented.");
                    } break;
                case 4:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 5: IProductWRTDerivBasePrismImpl<4, 4, 4, 5, 5, 4, true>
                            (in[0], in[1], in[2], out); break;
                        case 6: IProductWRTDerivBasePrismImpl<4, 4, 4, 6, 6, 5, true>
                            (in[0], in[1], in[2], out); break;
                        case 7: IProductWRTDerivBasePrismImpl<4, 4, 4, 7, 7, 6, true>
                            (in[0], in[1], in[2], out); break;
                        case 8: IProductWRTDerivBasePrismImpl<4, 4, 4, 8, 8, 7, true>
                            (in[0], in[1], in[2], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBPrism: # of modes / points combo not implemented.");
                    } break;
                case 5:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 6: IProductWRTDerivBasePrismImpl<5 ,5, 5, 6, 6, 5, true>
                            (in[0], in[1], in[2], out); break;
                        case 7: IProductWRTDerivBasePrismImpl<5 ,5, 5, 7, 7, 6, true>
                            (in[0], in[1], in[2], out); break;
                        case 8: IProductWRTDerivBasePrismImpl<5 ,5, 5, 8, 8, 7, true>
                            (in[0], in[1], in[2], out); break;
                        case 9: IProductWRTDerivBasePrismImpl<5 ,5, 5, 9, 9, 8, true>
                            (in[0], in[1], in[2], out); break;
                        case 10: IProductWRTDerivBasePrismImpl<5 ,5, 5, 10, 10, 9, true>
                            (in[0], in[1], in[2], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBPrism: # of modes / points combo not implemented.");
                    } break;
                case 6:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 7: IProductWRTDerivBasePrismImpl<6, 6, 6, 7, 7, 6, true>
                            (in[0], in[1], in[2], out); break;
                        case 8: IProductWRTDerivBasePrismImpl<6, 6, 6, 8, 8, 7, true>
                            (in[0], in[1], in[2], out); break;
                        case 9: IProductWRTDerivBasePrismImpl<6, 6, 6, 9, 9, 8, true>
                            (in[0], in[1], in[2], out); break;
                        case 10: IProductWRTDerivBasePrismImpl<6, 6, 6, 10, 10, 9, true>
                            (in[0], in[1], in[2], out); break;
                        case 11: IProductWRTDerivBasePrismImpl<6, 6, 6, 11, 11, 10, true>
                            (in[0], in[1], in[2], out); break;
                        case 12: IProductWRTDerivBasePrismImpl<6, 6, 6, 12, 12, 11, true>
                            (in[0], in[1], in[2], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBPrism: # of modes / points combo not implemented.");
                    } break;
                case 7:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 8: IProductWRTDerivBasePrismImpl<7, 7, 7, 8, 8, 7, true>
                            (in[0], in[1], in[2], out); break;
                        case 9: IProductWRTDerivBasePrismImpl<7, 7, 7, 9, 9, 8, true>
                            (in[0], in[1], in[2], out); break;
                        case 10: IProductWRTDerivBasePrismImpl<7, 7, 7, 10, 10, 9, true>
                            (in[0], in[1], in[2], out); break;
                        case 11: IProductWRTDerivBasePrismImpl<7, 7, 7, 11, 11, 10, true>
                            (in[0], in[1], in[2], out); break;
                        case 12: IProductWRTDerivBasePrismImpl<7, 7, 7, 12, 12, 11, true>
                            (in[0], in[1], in[2], out); break;
                        case 13: IProductWRTDerivBasePrismImpl<7, 7, 7, 13, 13, 12, true>
                            (in[0], in[1], in[2], out); break;
                        case 14: IProductWRTDerivBasePrismImpl<7, 7, 7, 14, 14, 13, true>
                            (in[0], in[1], in[2], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBPrism: # of modes / points combo not implemented.");
                    } break;
                case 8:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 9: IProductWRTDerivBasePrismImpl<8, 8, 8, 9, 9, 8, true>
                            (in[0], in[1], in[2], out); break;
                        case 10: IProductWRTDerivBasePrismImpl<8, 8, 8, 10, 10, 9, true>
                            (in[0], in[1], in[2], out); break;
                        case 11: IProductWRTDerivBasePrismImpl<8, 8, 8, 11, 11, 10, true>
                            (in[0], in[1], in[2], out); break;
                        case 12: IProductWRTDerivBasePrismImpl<8, 8, 8, 12, 12, 11, true>
                            (in[0], in[1], in[2], out); break;
                        case 13: IProductWRTDerivBasePrismImpl<8, 8, 8, 13, 13, 12, true>
                            (in[0], in[1], in[2], out); break;
                        case 14: IProductWRTDerivBasePrismImpl<8, 8, 8, 14, 14, 13, true>
                            (in[0], in[1], in[2], out); break;
                        case 15: IProductWRTDerivBasePrismImpl<8, 8, 8, 15, 15, 14, true>
                            (in[0], in[1], in[2], out); break;
                        case 16: IProductWRTDerivBasePrismImpl<8, 8, 8, 16, 16, 15, true>
                            (in[0], in[1], in[2], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBPrism: # of modes / points combo not implemented.");
                    } break;
                default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBPrism: # of modes / points combo not implemented.");
            }
        }
        else
        {
            switch(m_basis[0]->GetNumModes())
            {
                case 2:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 3: IProductWRTDerivBasePrismImpl<2, 2, 2, 3, 3, 2, false>
                            (in[0], in[1], in[2], out); break;
                        case 4: IProductWRTDerivBasePrismImpl<2, 2, 2, 4, 4, 3, false>
                            (in[0], in[1], in[2], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBPrism: # of modes / points combo not implemented.");
                    } break;
                case 3:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 4: IProductWRTDerivBasePrismImpl<3, 3, 3, 4, 4, 3, false>
                            (in[0], in[1], in[2], out); break;
                        case 5: IProductWRTDerivBasePrismImpl<3, 3, 3, 5, 5, 4, false>
                            (in[0], in[1], in[2], out); break;
                        case 6: IProductWRTDerivBasePrismImpl<3, 3, 3, 6, 6, 5, false>
                            (in[0], in[1], in[2], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBPrism: # of modes / points combo not implemented.");
                    } break;
                case 4:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 5: IProductWRTDerivBasePrismImpl<4, 4, 4, 5, 5, 4, false>
                            (in[0], in[1], in[2], out); break;
                        case 6: IProductWRTDerivBasePrismImpl<4, 4, 4, 6, 6, 5, false>
                            (in[0], in[1], in[2], out); break;
                        case 7: IProductWRTDerivBasePrismImpl<4, 4, 4, 7, 7, 6, false>
                            (in[0], in[1], in[2], out); break;
                        case 8: IProductWRTDerivBasePrismImpl<4, 4, 4, 8, 8, 7, false>
                            (in[0], in[1], in[2], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBPrism: # of modes / points combo not implemented.");
                    } break;
                case 5:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 6: IProductWRTDerivBasePrismImpl<5 ,5, 5, 6, 6, 5, false>
                            (in[0], in[1], in[2], out); break;
                        case 7: IProductWRTDerivBasePrismImpl<5 ,5, 5, 7, 7, 6, false>
                            (in[0], in[1], in[2], out); break;
                        case 8: IProductWRTDerivBasePrismImpl<5 ,5, 5, 8, 8, 7, false>
                            (in[0], in[1], in[2], out); break;
                        case 9: IProductWRTDerivBasePrismImpl<5 ,5, 5, 9, 9, 8, false>
                            (in[0], in[1], in[2], out); break;
                        case 10: IProductWRTDerivBasePrismImpl<5 ,5, 5, 10, 10, 9, false>
                            (in[0], in[1], in[2], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBPrism: # of modes / points combo not implemented.");
                    } break;
                case 6:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 7: IProductWRTDerivBasePrismImpl<6, 6, 6, 7, 7, 6, false>
                            (in[0], in[1], in[2], out); break;
                        case 8: IProductWRTDerivBasePrismImpl<6, 6, 6, 8, 8, 7, false>
                            (in[0], in[1], in[2], out); break;
                        case 9: IProductWRTDerivBasePrismImpl<6, 6, 6, 9, 9, 8, false>
                            (in[0], in[1], in[2], out); break;
                        case 10: IProductWRTDerivBasePrismImpl<6, 6, 6, 10, 10, 9, false>
                            (in[0], in[1], in[2], out); break;
                        case 11: IProductWRTDerivBasePrismImpl<6, 6, 6, 11, 11, 10, false>
                            (in[0], in[1], in[2], out); break;
                        case 12: IProductWRTDerivBasePrismImpl<6, 6, 6, 12, 12, 11, false>
                            (in[0], in[1], in[2], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBPrism: # of modes / points combo not implemented.");
                    } break;
                case 7:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 8: IProductWRTDerivBasePrismImpl<7, 7, 7, 8, 8, 7, false>
                            (in[0], in[1], in[2], out); break;
                        case 9: IProductWRTDerivBasePrismImpl<7, 7, 7, 9, 9, 8, false>
                            (in[0], in[1], in[2], out); break;
                        case 10: IProductWRTDerivBasePrismImpl<7, 7, 7, 10, 10, 9, false>
                            (in[0], in[1], in[2], out); break;
                        case 11: IProductWRTDerivBasePrismImpl<7, 7, 7, 11, 11, 10, false>
                            (in[0], in[1], in[2], out); break;
                        case 12: IProductWRTDerivBasePrismImpl<7, 7, 7, 12, 12, 11, false>
                            (in[0], in[1], in[2], out); break;
                        case 13: IProductWRTDerivBasePrismImpl<7, 7, 7, 13, 13, 12, false>
                            (in[0], in[1], in[2], out); break;
                        case 14: IProductWRTDerivBasePrismImpl<7, 7, 7, 14, 14, 13, false>
                            (in[0], in[1], in[2], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBPrism: # of modes / points combo not implemented.");
                    } break;
                case 8:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 9: IProductWRTDerivBasePrismImpl<8, 8, 8, 9, 9, 8, false>
                            (in[0], in[1], in[2], out); break;
                        case 10: IProductWRTDerivBasePrismImpl<8, 8, 8, 10, 10, 9, false>
                            (in[0], in[1], in[2], out); break;
                        case 11: IProductWRTDerivBasePrismImpl<8, 8, 8, 11, 11, 10, false>
                            (in[0], in[1], in[2], out); break;
                        case 12: IProductWRTDerivBasePrismImpl<8, 8, 8, 12, 12, 11, false>
                            (in[0], in[1], in[2], out); break;
                        case 13: IProductWRTDerivBasePrismImpl<8, 8, 8, 13, 13, 12, false>
                            (in[0], in[1], in[2], out); break;
                        case 14: IProductWRTDerivBasePrismImpl<8, 8, 8, 14, 14, 13, false>
                            (in[0], in[1], in[2], out); break;
                        case 15: IProductWRTDerivBasePrismImpl<8, 8, 8, 15, 15, 14, false>
                            (in[0], in[1], in[2], out); break;
                        case 16: IProductWRTDerivBasePrismImpl<8, 8, 8, 16, 16, 15, false>
                            (in[0], in[1], in[2], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBPrism: # of modes / points combo not implemented.");
                    } break;
                default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBPrism: # of modes / points combo not implemented.");
            }
        }
    }


    template<int NM0, int NM1, int NM2, int NQ0, int NQ1, int NQ2, bool CORRECT>
    void IProductWRTDerivBasePrismImpl(
        const Array<OneD, const NekDouble> &input0,
        const Array<OneD, const NekDouble> &input1,
        const Array<OneD, const NekDouble> &input2,
              Array<OneD,       NekDouble> &output)
    {
        auto* inptr0 = input0.data();
        auto* inptr1 = input1.data();
        auto* inptr2 = input2.data();
        auto* outptr = output.data();

        constexpr auto ndf = 9u;
        constexpr auto nqTot = NQ0 * NQ1 * NQ2;
        constexpr auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;


        // Get size of jacobian factor block
        auto dJSize = 1u;
        auto dfSize = ndf;
        if (DEFORMED)
        {
            dJSize = nqTot;
            dfSize = ndf*nqTot;
        }

        vec_t sums_kj[NQ1 * NQ2];
        vec_t sums_k[NQ2];
        vec_t corr_q[NM1];

        std::vector<vec_t, allocator<vec_t>>
            tmpIn0(nqTot), tmpIn1(nqTot), tmpIn2(nqTot),
            tmp0(nqTot), tmp1(nqTot), tmp2(nqTot), tmpOut(m_nmTot);

        const vec_t* df_ptr;
        const vec_t* jac_ptr;

        std::vector<vec_t, allocator<vec_t>>& Z0 = this->m_Z[0];
        // std::vector<vec_t, allocator<vec_t>>& Z1 = this->m_Z[1];
        std::vector<vec_t, allocator<vec_t>>& Z2 = this->m_Z[2];

        for (int e =0; e < this->m_nBlocks; ++e)
        {
            // Jacobian
            jac_ptr = &(this->m_jac[dJSize*e]);

            // Derivative factor
            df_ptr = &(this->m_df[dfSize*e]);

            // Load and transpose data
            load_interleave(inptr0, nqTot, tmpIn0);
            load_interleave(inptr1, nqTot, tmpIn1);
            load_interleave(inptr2, nqTot, tmpIn2);

            vec_t df0, df1, df2, df3, df4, df5, df6, df7, df8;
            if (!DEFORMED)
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

            for (size_t k = 0, cnt_kji = 0; k < NQ2; ++k)
            {
                // div in most external loop
                vec_t f0 = 2.0 / (1.0 - Z2[k]); // Load 1x

                for (size_t j = 0; j < NQ1; ++j)
                {
                    for (size_t i = 0; i < NQ0; ++i, ++cnt_kji)
                    {
                        if (DEFORMED)
                        {
                            df0 = df_ptr[cnt_kji * ndf];
                            df1 = df_ptr[cnt_kji * ndf + 1];
                            df2 = df_ptr[cnt_kji * ndf + 2];
                            df3 = df_ptr[cnt_kji * ndf + 3];
                            df4 = df_ptr[cnt_kji * ndf + 4];
                            df5 = df_ptr[cnt_kji * ndf + 5];
                            df6 = df_ptr[cnt_kji * ndf + 6];
                            df7 = df_ptr[cnt_kji * ndf + 7];
                            df8 = df_ptr[cnt_kji * ndf + 8];
                        }

                        // Calculate dx/dxi in[0] + dy/dxi in[1] + dz/dxi in[2]
                        vec_t tI0 = tmpIn0[cnt_kji]; // load 1x
                        vec_t tI1 = tmpIn1[cnt_kji]; // load 1x
                        vec_t tI2 = tmpIn2[cnt_kji]; // load 1x
                        vec_t t0 = df0 * tI0 + df3 * tI1 + df6 * tI2;
                        vec_t t1 = df1 * tI0 + df4 * tI1 + df7 * tI2;
                        vec_t t2 = df2 * tI0 + df5 * tI1 + df8 * tI2;

                        // Multiply by geometric factors
                        vec_t hf1 = 0.5 * (1.0 + Z0[i]); //Load 1x
                        // Scale by geometric factor 2/(1-z2)
                        t0 *= f0;
                        // Scale by geometric factor (1+z0)/(1-z2)
                        vec_t f1t2 = hf1 * t2;
                        t0.fma(f1t2, f0);

                        // Store
                        tmp0[cnt_kji] = t0; // store 1x
                        tmp1[cnt_kji] = t1; // store 1x
                        tmp2[cnt_kji] = t2; // store 1x
                    }
                }
            }



            // IP DB0 B1 B2
            IProductPrismKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, CORRECT, false,
                false, DEFORMED>(
                tmp0, this->m_dbdata[0], this->m_bdata[1], this->m_bdata[2],
                this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                sums_kj, sums_k,
                corr_q,
                tmpOut);

            // IP B0 DB1 B2
            IProductPrismKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, CORRECT, false,
                true, DEFORMED>(
                tmp1, this->m_bdata[0], this->m_dbdata[1], this->m_bdata[2],
                this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                sums_kj, sums_k,
                corr_q,
                tmpOut);

            // IP B0 B1 DB2
            IProductPrismKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, CORRECT, false,
                true, DEFORMED>(
                tmp2, this->m_bdata[0], this->m_bdata[1], this->m_dbdata[2],
                this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                sums_kj, sums_k,
                corr_q,
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
    int m_nmTot;
};

template<bool DEFORMED = false>
struct IProductWRTDerivBaseTet : public IProductWRTDerivBase, public Helper<3, DEFORMED>
{
    IProductWRTDerivBaseTet(std::vector<LibUtilities::BasisSharedPtr> basis,
                   int nElmt)
        : IProductWRTDerivBase(basis, nElmt),
          Helper<3, DEFORMED>(basis, nElmt),
          m_nmTot(LibUtilities::StdTetData::getNumberOfCoefficients(
                      this->m_nm[0], this->m_nm[1], this->m_nm[2]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<IProductWRTDerivBaseTet<DEFORMED>>(basis, nElmt);
    }

   void operator()( const Array<OneD, Array<OneD, NekDouble>> &in,
        Array<OneD, NekDouble> &out) final
    {
        // Check preconditions
        ASSERTL0(m_basis[0]->GetNumModes() == m_basis[1]->GetNumModes() &&
            m_basis[0]->GetNumModes() == m_basis[2]->GetNumModes() &&
            m_basis[0]->GetNumPoints() == m_basis[1]->GetNumPoints()+1 &&
            m_basis[0]->GetNumPoints() == m_basis[2]->GetNumPoints()+1,
            "MatrixFree requires homogenous modes/points");

        if (m_basis[0]->GetBasisType() == LibUtilities::eModified_A)
        {
            switch(m_basis[0]->GetNumModes())
            {
                case 2:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 3: IProductWRTDerivBaseTetImpl<2, 2, 2, 3, 2, 2, true>
                            (in[0], in[1], in[2], out); break;
                        case 4: IProductWRTDerivBaseTetImpl<2, 2, 2, 4, 3, 3, true>
                            (in[0], in[1], in[2], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTet: # of modes / points combo not implemented.");
                    } break;
                case 3:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 4: IProductWRTDerivBaseTetImpl<3, 3, 3, 4, 3, 3, true>
                            (in[0], in[1], in[2], out); break;
                        case 5: IProductWRTDerivBaseTetImpl<3, 3, 3, 5, 4, 4, true>
                            (in[0], in[1], in[2], out); break;
                        case 6: IProductWRTDerivBaseTetImpl<3, 3, 3, 6, 5, 5, true>
                            (in[0], in[1], in[2], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTet: # of modes / points combo not implemented.");
                    } break;
                case 4:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 5: IProductWRTDerivBaseTetImpl<4, 4, 4, 5, 4, 4, true>
                            (in[0], in[1], in[2], out); break;
                        case 6: IProductWRTDerivBaseTetImpl<4, 4, 4, 6, 5, 5, true>
                            (in[0], in[1], in[2], out); break;
                        case 7: IProductWRTDerivBaseTetImpl<4, 4, 4, 7, 6, 6, true>
                            (in[0], in[1], in[2], out); break;
                        case 8: IProductWRTDerivBaseTetImpl<4, 4, 4, 8, 7, 7, true>
                            (in[0], in[1], in[2], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTet: # of modes / points combo not implemented.");
                    } break;
                case 5:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 6: IProductWRTDerivBaseTetImpl<5 ,5, 5, 6, 5, 5, true>
                            (in[0], in[1], in[2], out); break;
                        case 7: IProductWRTDerivBaseTetImpl<5 ,5, 5, 7, 6, 6, true>
                            (in[0], in[1], in[2], out); break;
                        case 8: IProductWRTDerivBaseTetImpl<5 ,5, 5, 8, 7, 7, true>
                            (in[0], in[1], in[2], out); break;
                        case 9: IProductWRTDerivBaseTetImpl<5 ,5, 5, 9, 8, 8, true>
                            (in[0], in[1], in[2], out); break;
                        case 10: IProductWRTDerivBaseTetImpl<5 ,5, 5, 10, 9, 9, true>
                            (in[0], in[1], in[2], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTet: # of modes / points combo not implemented.");
                    } break;
                case 6:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 7: IProductWRTDerivBaseTetImpl<6, 6, 6, 7, 6, 6, true>
                            (in[0], in[1], in[2], out); break;
                        case 8: IProductWRTDerivBaseTetImpl<6, 6, 6, 8, 7, 7, true>
                            (in[0], in[1], in[2], out); break;
                        case 9: IProductWRTDerivBaseTetImpl<6, 6, 6, 9, 8, 8, true>
                            (in[0], in[1], in[2], out); break;
                        case 10: IProductWRTDerivBaseTetImpl<6, 6, 6, 10, 9, 9, true>
                            (in[0], in[1], in[2], out); break;
                        case 11: IProductWRTDerivBaseTetImpl<6, 6, 6, 11, 10, 10, true>
                            (in[0], in[1], in[2], out); break;
                        case 12: IProductWRTDerivBaseTetImpl<6, 6, 6, 12, 11, 11, true>
                            (in[0], in[1], in[2], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTet: # of modes / points combo not implemented.");
                    } break;
                case 7:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 8: IProductWRTDerivBaseTetImpl<7, 7, 7, 8, 7, 7, true>
                            (in[0], in[1], in[2], out); break;
                        case 9: IProductWRTDerivBaseTetImpl<7, 7, 7, 9, 8, 8, true>
                            (in[0], in[1], in[2], out); break;
                        case 10: IProductWRTDerivBaseTetImpl<7, 7, 7, 10, 9, 9, true>
                            (in[0], in[1], in[2], out); break;
                        case 11: IProductWRTDerivBaseTetImpl<7, 7, 7, 11, 10, 10, true>
                            (in[0], in[1], in[2], out); break;
                        case 12: IProductWRTDerivBaseTetImpl<7, 7, 7, 12, 11, 11, true>
                            (in[0], in[1], in[2], out); break;
                        case 13: IProductWRTDerivBaseTetImpl<7, 7, 7, 13, 12, 12, true>
                            (in[0], in[1], in[2], out); break;
                        case 14: IProductWRTDerivBaseTetImpl<7, 7, 7, 14, 13, 13, true>
                            (in[0], in[1], in[2], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTet: # of modes / points combo not implemented.");
                    } break;
                case 8:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 9: IProductWRTDerivBaseTetImpl<8, 8, 8, 9, 8, 8, true>
                            (in[0], in[1], in[2], out); break;
                        case 10: IProductWRTDerivBaseTetImpl<8, 8, 8, 10, 9, 9, true>
                            (in[0], in[1], in[2], out); break;
                        case 11: IProductWRTDerivBaseTetImpl<8, 8, 8, 11, 10, 10, true>
                            (in[0], in[1], in[2], out); break;
                        case 12: IProductWRTDerivBaseTetImpl<8, 8, 8, 12, 11, 11, true>
                            (in[0], in[1], in[2], out); break;
                        case 13: IProductWRTDerivBaseTetImpl<8, 8, 8, 13, 12, 12, true>
                            (in[0], in[1], in[2], out); break;
                        case 14: IProductWRTDerivBaseTetImpl<8, 8, 8, 14, 13, 13, true>
                            (in[0], in[1], in[2], out); break;
                        case 15: IProductWRTDerivBaseTetImpl<8, 8, 8, 15, 14, 14, true>
                            (in[0], in[1], in[2], out); break;
                        case 16: IProductWRTDerivBaseTetImpl<8, 8, 8, 16, 15, 15, true>
                            (in[0], in[1], in[2], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTet: # of modes / points combo not implemented.");
                    } break;
                default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTet: # of modes / points combo not implemented.");
            }
        }
        else
        {
            switch(m_basis[0]->GetNumModes())
            {
                case 2:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 3: IProductWRTDerivBaseTetImpl<2, 2, 2, 3, 2, 2, false>
                            (in[0], in[1], in[2], out); break;
                        case 4: IProductWRTDerivBaseTetImpl<2, 2, 2, 4, 3, 3, false>
                            (in[0], in[1], in[2], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTet: # of modes / points combo not implemented.");
                    } break;
                case 3:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 4: IProductWRTDerivBaseTetImpl<3, 3, 3, 4, 3, 3, false>
                            (in[0], in[1], in[2], out); break;
                        case 5: IProductWRTDerivBaseTetImpl<3, 3, 3, 5, 4, 4, false>
                            (in[0], in[1], in[2], out); break;
                        case 6: IProductWRTDerivBaseTetImpl<3, 3, 3, 6, 5, 5, false>
                            (in[0], in[1], in[2], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTet: # of modes / points combo not implemented.");
                    } break;
                case 4:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 5: IProductWRTDerivBaseTetImpl<4, 4, 4, 5, 4, 4, false>
                            (in[0], in[1], in[2], out); break;
                        case 6: IProductWRTDerivBaseTetImpl<4, 4, 4, 6, 5, 5, false>
                            (in[0], in[1], in[2], out); break;
                        case 7: IProductWRTDerivBaseTetImpl<4, 4, 4, 7, 6, 6, false>
                            (in[0], in[1], in[2], out); break;
                        case 8: IProductWRTDerivBaseTetImpl<4, 4, 4, 8, 7, 7, false>
                            (in[0], in[1], in[2], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTet: # of modes / points combo not implemented.");
                    } break;
                case 5:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 6: IProductWRTDerivBaseTetImpl<5 ,5, 5, 6, 5, 5, false>
                            (in[0], in[1], in[2], out); break;
                        case 7: IProductWRTDerivBaseTetImpl<5 ,5, 5, 7, 6, 6, false>
                            (in[0], in[1], in[2], out); break;
                        case 8: IProductWRTDerivBaseTetImpl<5 ,5, 5, 8, 7, 7, false>
                            (in[0], in[1], in[2], out); break;
                        case 9: IProductWRTDerivBaseTetImpl<5 ,5, 5, 9, 8, 8, false>
                            (in[0], in[1], in[2], out); break;
                        case 10: IProductWRTDerivBaseTetImpl<5 ,5, 5, 10, 9, 9, false>
                            (in[0], in[1], in[2], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTet: # of modes / points combo not implemented.");
                    } break;
                case 6:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 7: IProductWRTDerivBaseTetImpl<6, 6, 6, 7, 6, 6, false>
                            (in[0], in[1], in[2], out); break;
                        case 8: IProductWRTDerivBaseTetImpl<6, 6, 6, 8, 7, 7, false>
                            (in[0], in[1], in[2], out); break;
                        case 9: IProductWRTDerivBaseTetImpl<6, 6, 6, 9, 8, 8, false>
                            (in[0], in[1], in[2], out); break;
                        case 10: IProductWRTDerivBaseTetImpl<6, 6, 6, 10, 9, 9, false>
                            (in[0], in[1], in[2], out); break;
                        case 11: IProductWRTDerivBaseTetImpl<6, 6, 6, 11, 10, 10, false>
                            (in[0], in[1], in[2], out); break;
                        case 12: IProductWRTDerivBaseTetImpl<6, 6, 6, 12, 11, 11, false>
                            (in[0], in[1], in[2], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTet: # of modes / points combo not implemented.");
                    } break;
                case 7:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 8: IProductWRTDerivBaseTetImpl<7, 7, 7, 8, 7, 7, false>
                            (in[0], in[1], in[2], out); break;
                        case 9: IProductWRTDerivBaseTetImpl<7, 7, 7, 9, 8, 8, false>
                            (in[0], in[1], in[2], out); break;
                        case 10: IProductWRTDerivBaseTetImpl<7, 7, 7, 10, 9, 9, false>
                            (in[0], in[1], in[2], out); break;
                        case 11: IProductWRTDerivBaseTetImpl<7, 7, 7, 11, 10, 10, false>
                            (in[0], in[1], in[2], out); break;
                        case 12: IProductWRTDerivBaseTetImpl<7, 7, 7, 12, 11, 11, false>
                            (in[0], in[1], in[2], out); break;
                        case 13: IProductWRTDerivBaseTetImpl<7, 7, 7, 13, 12, 12, false>
                            (in[0], in[1], in[2], out); break;
                        case 14: IProductWRTDerivBaseTetImpl<7, 7, 7, 14, 13, 13, false>
                            (in[0], in[1], in[2], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTet: # of modes / points combo not implemented.");
                    } break;
                case 8:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 9: IProductWRTDerivBaseTetImpl<8, 8, 8, 9, 8, 8, false>
                            (in[0], in[1], in[2], out); break;
                        case 10: IProductWRTDerivBaseTetImpl<8, 8, 8, 10, 9, 9, false>
                            (in[0], in[1], in[2], out); break;
                        case 11: IProductWRTDerivBaseTetImpl<8, 8, 8, 11, 10, 10, false>
                            (in[0], in[1], in[2], out); break;
                        case 12: IProductWRTDerivBaseTetImpl<8, 8, 8, 12, 11, 11, false>
                            (in[0], in[1], in[2], out); break;
                        case 13: IProductWRTDerivBaseTetImpl<8, 8, 8, 13, 12, 12, false>
                            (in[0], in[1], in[2], out); break;
                        case 14: IProductWRTDerivBaseTetImpl<8, 8, 8, 14, 13, 13, false>
                            (in[0], in[1], in[2], out); break;
                        case 15: IProductWRTDerivBaseTetImpl<8, 8, 8, 15, 14, 14, false>
                            (in[0], in[1], in[2], out); break;
                        case 16: IProductWRTDerivBaseTetImpl<8, 8, 8, 16, 15, 15, false>
                            (in[0], in[1], in[2], out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTet: # of modes / points combo not implemented.");
                    } break;
                default: NEKERROR(ErrorUtil::efatal,
                    "IProductWRTDBTet: # of modes / points combo not implemented.");
            }
        }
    }


    template<int NM0, int NM1, int NM2, int NQ0, int NQ1, int NQ2, bool CORRECT>
    void IProductWRTDerivBaseTetImpl(
        const Array<OneD, const NekDouble> &input0,
        const Array<OneD, const NekDouble> &input1,
        const Array<OneD, const NekDouble> &input2,
              Array<OneD,       NekDouble> &output)
    {
        auto* inptr0 = input0.data();
        auto* inptr1 = input1.data();
        auto* inptr2 = input2.data();
        auto* outptr = output.data();

        constexpr auto ndf = 9u;
        constexpr auto nqTot = NQ0 * NQ1 * NQ2;
        constexpr auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;


        // Get size of jacobian factor block
        auto dJSize = 1u;
        auto dfSize = ndf;
        if (DEFORMED)
        {
            dJSize = nqTot;
            dfSize = ndf*nqTot;
        }

        vec_t wsp[NQ1 * NQ2 + NQ2];

        std::vector<vec_t, allocator<vec_t>>
            tmpIn0(nqTot), tmpIn1(nqTot), tmpIn2(nqTot),
            tmp0(nqTot), tmp1(nqTot), tmp2(nqTot), tmpOut(m_nmTot);

        const vec_t* df_ptr;
        const vec_t* jac_ptr;

        std::vector<vec_t, allocator<vec_t>>& Z0 = this->m_Z[0];
        std::vector<vec_t, allocator<vec_t>>& Z1 = this->m_Z[1];
        std::vector<vec_t, allocator<vec_t>>& Z2 = this->m_Z[2];

        for (int e =0; e < this->m_nBlocks; ++e)
        {
            // Jacobian
            jac_ptr = &(this->m_jac[dJSize*e]);

            // Derivative factor
            df_ptr = &(this->m_df[dfSize*e]);

            // Load and transpose data
            load_interleave(inptr0, nqTot, tmpIn0);
            load_interleave(inptr1, nqTot, tmpIn1);
            load_interleave(inptr2, nqTot, tmpIn2);

            vec_t df0, df1, df2, df3, df4, df5, df6, df7, df8;
            if (!DEFORMED)
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

            for (size_t k = 0, cnt_kji = 0; k < NQ2; ++k)
            {
                vec_t f2 = 2.0 / (1.0 - Z2[k]); // Load 1x

                for (size_t j = 0; j < NQ1; ++j)
                {
                    vec_t Z1Load = Z1[j]; // Load1x
                    vec_t f3 = 0.5 * (1.0 + Z1Load);
                    vec_t f0 = 2.0 * f2 / (1.0 - Z1Load);

                    for (size_t i = 0; i < NQ0; ++i, ++cnt_kji)
                    {
                        if (DEFORMED)
                        {
                            df0 = df_ptr[cnt_kji * ndf];
                            df1 = df_ptr[cnt_kji * ndf + 1];
                            df2 = df_ptr[cnt_kji * ndf + 2];
                            df3 = df_ptr[cnt_kji * ndf + 3];
                            df4 = df_ptr[cnt_kji * ndf + 4];
                            df5 = df_ptr[cnt_kji * ndf + 5];
                            df6 = df_ptr[cnt_kji * ndf + 6];
                            df7 = df_ptr[cnt_kji * ndf + 7];
                            df8 = df_ptr[cnt_kji * ndf + 8];
                        }

                        // Calculate dx/dxi in[0] + dy/dxi in[1] + dz/dxi in[2]
                        vec_t tI0 = tmpIn0[cnt_kji]; // load 1x
                        vec_t tI1 = tmpIn1[cnt_kji]; // load 1x
                        vec_t tI2 = tmpIn2[cnt_kji]; // load 1x
                        vec_t t0 = df0 * tI0 + df3 * tI1 + df6 * tI2;
                        vec_t t1 = df1 * tI0 + df4 * tI1 + df7 * tI2;
                        vec_t t2 = df2 * tI0 + df5 * tI1 + df8 * tI2;

                        // Multiply by geometric factors
                        vec_t f1 = 0.5 * (1.0 + Z0[i]); // Load 1x

                        // Scale by geometric factor 1 and add to t0
                        t0.fma(t1+t2, f1);
                        // Scale by geometric factor 0
                        t0 *= f0;

                        // Scale by geometric factor 3 and add to t1
                        t1.fma(t2, f3);
                        // Scale by geometric factor 2
                        t1 *= f2;

                        // Store
                        tmp0[cnt_kji] = t0; // store 1x
                        tmp1[cnt_kji] = t1; // store 1x
                        tmp2[cnt_kji] = t2; // store 1x
                    }
                }
            }

            // IP DB0 B1 B2
            IProductTetKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, CORRECT, false,
                false, DEFORMED>(
                tmp0, this->m_dbdata[0], this->m_bdata[1], this->m_bdata[2],
                this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                wsp, tmpOut);

            // IP B0 DB1 B2
            IProductTetKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, CORRECT, false,
                true, DEFORMED>(
                tmp1, this->m_bdata[0], this->m_dbdata[1], this->m_bdata[2],
                this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                wsp, tmpOut);

            // IP B0 B1 DB2
            IProductTetKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, CORRECT, false,
                true, DEFORMED>(
                tmp2, this->m_bdata[0], this->m_bdata[1], this->m_dbdata[2],
                this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                wsp, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr0 += nqBlocks;
            inptr1 += nqBlocks;
            inptr2 += nqBlocks;
            outptr += nmBlocks;
        }
    }

private:
    int m_nmTot;
};

} // namespace MatrixFree
} // namespace Nektar

#endif
