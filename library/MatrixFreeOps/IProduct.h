#ifndef NEKTAR_LIBRARY_MF_IPRODUCT_H
#define NEKTAR_LIBRARY_MF_IPRODUCT_H

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
struct IProductQuad : public IProduct, public Helper<2, DEFORMED>
{
    IProductQuad(std::vector<LibUtilities::BasisSharedPtr> basis,
                   int nElmt)
        : IProduct(basis, nElmt),
          Helper<2, DEFORMED>(basis, nElmt),
          m_nmTot(LibUtilities::StdQuadData::getNumberOfCoefficients(
                      this->m_nm[0], this->m_nm[1]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<IProductQuad<DEFORMED>>(basis, nElmt);
    }

    void operator()(const Array<OneD, const NekDouble> &in,
                                  Array<OneD,       NekDouble> &out) final
    {
        // Check preconditions
        ASSERTL0(m_basis[0]->GetNumModes() == m_basis[1]->GetNumModes() &&
            m_basis[0]->GetNumPoints() == m_basis[1]->GetNumPoints(),
            " requires homogenous modes/points");

        switch(m_basis[0]->GetNumModes())
        {
            case 2:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 2: IProductQuadImpl<2 ,2 ,2 ,2 >(in, out); break;
                    case 3: IProductQuadImpl<2 ,2 ,3 ,3 >(in, out); break;
                    case 4: IProductQuadImpl<2 ,2 ,4 ,4 >(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductQuad: # of modes / points combo not implemented.");
                } break;
            case 3:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 3: IProductQuadImpl<3 ,3 ,3 ,3 >(in, out); break;
                    case 4: IProductQuadImpl<3 ,3 ,4 ,4 >(in, out); break;
                    case 5: IProductQuadImpl<3 ,3 ,5 ,5 >(in, out); break;
                    case 6: IProductQuadImpl<3 ,3 ,6 ,6 >(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductQuad: # of modes / points combo not implemented.");
                } break;
            case 4:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 4: IProductQuadImpl<4 ,4 ,4 ,4 >(in, out); break;
                    case 5: IProductQuadImpl<4 ,4 ,5 ,5 >(in, out); break;
                    case 6: IProductQuadImpl<4 ,4 ,6 ,6 >(in, out); break;
                    case 7: IProductQuadImpl<4 ,4 ,7 ,7 >(in, out); break;
                    case 8: IProductQuadImpl<4 ,4 ,8 ,8 >(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductQuad: # of modes / points combo not implemented.");
                } break;
            case 5:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 5: IProductQuadImpl<5 ,5 ,5 ,5 >(in, out); break;
                    case 6: IProductQuadImpl<5 ,5 ,6 ,6 >(in, out); break;
                    case 7: IProductQuadImpl<5 ,5 ,7 ,7 >(in, out); break;
                    case 8: IProductQuadImpl<5 ,5 ,8 ,8 >(in, out); break;
                    case 9: IProductQuadImpl<5 ,5 ,9 ,9 >(in, out); break;
                    case 10: IProductQuadImpl<5 ,5 ,10 ,10 >(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductQuad: # of modes / points combo not implemented.");
                } break;
            case 6:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 6: IProductQuadImpl<6 ,6 ,6 ,6 >(in, out); break;
                    case 7: IProductQuadImpl<6 ,6 ,7 ,7 >(in, out); break;
                    case 8: IProductQuadImpl<6 ,6 ,8 ,8 >(in, out); break;
                    case 9: IProductQuadImpl<6 ,6 ,9 ,9 >(in, out); break;
                    case 10: IProductQuadImpl<6 ,6 ,10 ,10 >(in, out); break;
                    case 11: IProductQuadImpl<6 ,6 ,11 ,11 >(in, out); break;
                    case 12: IProductQuadImpl<6 ,6 ,12 ,12 >(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductQuad: # of modes / points combo not implemented.");
                } break;
            case 7:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 7: IProductQuadImpl<7 ,7 ,7 ,7 >(in, out); break;
                    case 8: IProductQuadImpl<7 ,7 ,8 ,8 >(in, out); break;
                    case 9: IProductQuadImpl<7 ,7 ,9 ,9 >(in, out); break;
                    case 10: IProductQuadImpl<7 ,7 ,10 ,10 >(in, out); break;
                    case 11: IProductQuadImpl<7 ,7 ,11 ,11 >(in, out); break;
                    case 12: IProductQuadImpl<7 ,7 ,12 ,12 >(in, out); break;
                    case 13: IProductQuadImpl<7 ,7 ,13 ,13 >(in, out); break;
                    case 14: IProductQuadImpl<7 ,7 ,14 ,14 >(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductQuad: # of modes / points combo not implemented.");
                } break;
            case 8:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 8: IProductQuadImpl<8 ,8 ,8 ,8 >(in, out); break;
                    case 9: IProductQuadImpl<8 ,8 ,9 ,9 >(in, out); break;
                    case 10: IProductQuadImpl<8 ,8 ,10 ,10 >(in, out); break;
                    case 11: IProductQuadImpl<8 ,8 ,11 ,11 >(in, out); break;
                    case 12: IProductQuadImpl<8 ,8 ,12 ,12 >(in, out); break;
                    case 13: IProductQuadImpl<8 ,8 ,13 ,13 >(in, out); break;
                    case 14: IProductQuadImpl<8 ,8 ,14 ,14 >(in, out); break;
                    case 15: IProductQuadImpl<8 ,8 ,15 ,15 >(in, out); break;
                    case 16: IProductQuadImpl<8 ,8 ,16 ,16 >(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductQuad: # of modes / points combo not implemented.");
                } break;;
            default: NEKERROR(ErrorUtil::efatal,
                "IProductQuad: # of modes / points combo not implemented.");
        }
    }

    template<int NM0, int NM1, int NQ0, int NQ1>
    void IProductQuadImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &output)
    {
        auto* inptr = input.data();
        auto* outptr = output.data();

        constexpr auto nqTot = NQ0 * NQ1;
        constexpr auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        vec_t sums_j[NQ1]; //Sums over eta0 for each value of eta1;

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpOut(m_nmTot);
        vec_t* jac_ptr;
        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            if (DEFORMED)
            {
                jac_ptr = &((*this->m_jac)[nqTot*e]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e]);
            }

            // Load and transpose data
            load_interleave(inptr, nqTot, tmpIn);

            IProductQuadKernel<NM0, NM1, NQ0, NQ1, false, false, DEFORMED>(
                tmpIn, this->m_bdata[0], this->m_bdata[1],
                this->m_w[0], this->m_w[1], jac_ptr,
                sums_j, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr += nqBlocks;
            outptr += nmBlocks;
        }
    }
public:
    NekDouble Ndof() final
    {
        return m_nmTot * this->m_nElmt;
    }

private:
    int m_nmTot;
};

template<bool DEFORMED = false>
struct IProductTri : public IProduct, public Helper<2, DEFORMED>
{
    IProductTri(std::vector<LibUtilities::BasisSharedPtr> basis,
                   int nElmt)
        : IProduct(basis, nElmt),
          Helper<2, DEFORMED>(basis, nElmt),
          m_nmTot(LibUtilities::StdTriData::getNumberOfCoefficients(
                      this->m_nm[0], this->m_nm[1]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<IProductTri<DEFORMED>>(basis, nElmt);
    }

    void operator()(const Array<OneD, const NekDouble> &in,
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
                        case 3: IProductTriImpl<2 ,2 ,3 ,2 ,true>(in, out); break;
                        case 4: IProductTriImpl<2 ,2 ,4 ,3 ,true>(in, out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductTri: # of modes / points combo not implemented.");
                    } break;
                case 3:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 4: IProductTriImpl<3 ,3 ,4 ,3 ,true>(in, out); break;
                        case 5: IProductTriImpl<3 ,3 ,5 ,4 ,true>(in, out); break;
                        case 6: IProductTriImpl<3 ,3 ,6 ,5 ,true>(in, out); break;
                        case 7: IProductTriImpl<3 ,3 ,7 ,6 ,true>(in, out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductTri: # of modes / points combo not implemented.");
                    } break;
                case 4:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 5: IProductTriImpl<4 ,4 ,5 ,4 ,true>(in, out); break;
                        case 6: IProductTriImpl<4 ,4 ,6 ,5 ,true>(in, out); break;
                        case 7: IProductTriImpl<4 ,4 ,7 ,6 ,true>(in, out); break;
                        case 8: IProductTriImpl<4 ,4 ,8 ,7 ,true>(in, out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductTri: # of modes / points combo not implemented.");
                    } break;
                case 5:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 6: IProductTriImpl<5 ,5 ,6 ,5 ,true>(in, out); break;
                        case 7: IProductTriImpl<5 ,5 ,7 ,6 ,true>(in, out); break;
                        case 8: IProductTriImpl<5 ,5 ,8 ,7 ,true>(in, out); break;
                        case 9: IProductTriImpl<5 ,5 ,9 ,8 ,true>(in, out); break;
                        case 10: IProductTriImpl<5 ,5 ,10 ,9 ,true>(in, out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductTri: # of modes / points combo not implemented.");
                    } break;
                case 6:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 7: IProductTriImpl<6 ,6 ,7 ,6 ,true>(in, out); break;
                        case 8: IProductTriImpl<6 ,6 ,8 ,7 ,true>(in, out); break;
                        case 9: IProductTriImpl<6 ,6 ,9 ,8 ,true>(in, out); break;
                        case 10: IProductTriImpl<6 ,6 ,10 ,9 ,true>(in, out); break;
                        case 11: IProductTriImpl<6 ,6 ,11 ,10 ,true>(in, out); break;
                        case 12: IProductTriImpl<6 ,6 ,12 ,11 ,true>(in, out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductTri: # of modes / points combo not implemented.");
                    } break;
                case 7:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 8: IProductTriImpl<7 ,7 ,8 ,7 ,true>(in, out); break;
                        case 9: IProductTriImpl<7 ,7 ,9 ,8 ,true>(in, out); break;
                        case 10: IProductTriImpl<7 ,7 ,10 ,9 ,true>(in, out); break;
                        case 11: IProductTriImpl<7 ,7 ,11 ,10 ,true>(in, out); break;
                        case 12: IProductTriImpl<7 ,7 ,12 ,11 ,true>(in, out); break;
                        case 13: IProductTriImpl<7 ,7 ,13 ,12 ,true>(in, out); break;
                        case 14: IProductTriImpl<7 ,7 ,14 ,13 ,true>(in, out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductTri: # of modes / points combo not implemented.");
                    } break;
                case 8:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 9: IProductTriImpl<8 ,8 ,9 ,8 ,true>(in, out); break;
                        case 10: IProductTriImpl<8 ,8 ,10 ,9 ,true>(in, out); break;
                        case 11: IProductTriImpl<8 ,8 ,11 ,10 ,true>(in, out); break;
                        case 12: IProductTriImpl<8 ,8 ,12 ,11 ,true>(in, out); break;
                        case 13: IProductTriImpl<8 ,8 ,13 ,12 ,true>(in, out); break;
                        case 14: IProductTriImpl<8 ,8 ,14 ,13 ,true>(in, out); break;
                        case 15: IProductTriImpl<8 ,8 ,15 ,14 ,true>(in, out); break;
                        case 16: IProductTriImpl<8 ,8 ,16 ,15 ,true>(in, out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductTri: # of modes / points combo not implemented.");
                    } break;
                default: NEKERROR(ErrorUtil::efatal,
                    "IProductTri: # of modes / points combo not implemented.");
            }
        }
        else
        {
            switch(m_basis[0]->GetNumModes())
            {
                case 2:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 3: IProductTriImpl<2 ,2 ,3 ,2 ,false>(in, out); break;
                        case 4: IProductTriImpl<2 ,2 ,4 ,3 ,false>(in, out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductTri: # of modes / points combo not implemented.");
                    } break;
                case 3:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 4: IProductTriImpl<3 ,3 ,4 ,3 ,false>(in, out); break;
                        case 5: IProductTriImpl<3 ,3 ,5 ,4 ,false>(in, out); break;
                        case 6: IProductTriImpl<3 ,3 ,6 ,5 ,false>(in, out); break;
                        case 7: IProductTriImpl<3 ,3 ,7 ,6 ,false>(in, out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductTri: # of modes / points combo not implemented.");
                    } break;
                case 4:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 5: IProductTriImpl<4 ,4 ,5 ,4 ,false>(in, out); break;
                        case 6: IProductTriImpl<4 ,4 ,6 ,5 ,false>(in, out); break;
                        case 7: IProductTriImpl<4 ,4 ,7 ,6 ,false>(in, out); break;
                        case 8: IProductTriImpl<4 ,4 ,8 ,7 ,false>(in, out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductTri: # of modes / points combo not implemented.");
                    } break;
                case 5:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 6: IProductTriImpl<5 ,5 ,6 ,5 ,false>(in, out); break;
                        case 7: IProductTriImpl<5 ,5 ,7 ,6 ,false>(in, out); break;
                        case 8: IProductTriImpl<5 ,5 ,8 ,7 ,false>(in, out); break;
                        case 9: IProductTriImpl<5 ,5 ,9 ,8 ,false>(in, out); break;
                        case 10: IProductTriImpl<5 ,5 ,10 ,9 ,false>(in, out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductTri: # of modes / points combo not implemented.");
                    } break;
                case 6:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 7: IProductTriImpl<6 ,6 ,7 ,6 ,false>(in, out); break;
                        case 8: IProductTriImpl<6 ,6 ,8 ,7 ,false>(in, out); break;
                        case 9: IProductTriImpl<6 ,6 ,9 ,8 ,false>(in, out); break;
                        case 10: IProductTriImpl<6 ,6 ,10 ,9 ,false>(in, out); break;
                        case 11: IProductTriImpl<6 ,6 ,11 ,10 ,false>(in, out); break;
                        case 12: IProductTriImpl<6 ,6 ,12 ,11 ,false>(in, out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductTri: # of modes / points combo not implemented.");
                    } break;
                case 7:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 8: IProductTriImpl<7 ,7 ,8 ,7 ,false>(in, out); break;
                        case 9: IProductTriImpl<7 ,7 ,9 ,8 ,false>(in, out); break;
                        case 10: IProductTriImpl<7 ,7 ,10 ,9 ,false>(in, out); break;
                        case 11: IProductTriImpl<7 ,7 ,11 ,10 ,false>(in, out); break;
                        case 12: IProductTriImpl<7 ,7 ,12 ,11 ,false>(in, out); break;
                        case 13: IProductTriImpl<7 ,7 ,13 ,12 ,false>(in, out); break;
                        case 14: IProductTriImpl<7 ,7 ,14 ,13 ,false>(in, out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductTri: # of modes / points combo not implemented.");
                    } break;
                case 8:
                    switch(m_basis[0]->GetNumPoints())
                    {
                        case 9: IProductTriImpl<8 ,8 ,9 ,8 ,false>(in, out); break;
                        case 10: IProductTriImpl<8 ,8 ,10 ,9 ,false>(in, out); break;
                        case 11: IProductTriImpl<8 ,8 ,11 ,10 ,false>(in, out); break;
                        case 12: IProductTriImpl<8 ,8 ,12 ,11 ,false>(in, out); break;
                        case 13: IProductTriImpl<8 ,8 ,13 ,12 ,false>(in, out); break;
                        case 14: IProductTriImpl<8 ,8 ,14 ,13 ,false>(in, out); break;
                        case 15: IProductTriImpl<8 ,8 ,15 ,14 ,false>(in, out); break;
                        case 16: IProductTriImpl<8 ,8 ,16 ,15 ,false>(in, out); break;
                        default: NEKERROR(ErrorUtil::efatal,
                    "IProductTri: # of modes / points combo not implemented.");
                    } break;
                default: NEKERROR(ErrorUtil::efatal,
                    "IProductTri: # of modes / points combo not implemented.");
            }
        }
    }

    template<int NM0, int NM1, int NQ0, int NQ1, bool CORRECT>
    void IProductTriImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &output)
    {
        auto *inptr = input.data();
        auto *outptr = output.data();

        constexpr auto nqTot = NQ0 * NQ1;
        constexpr auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        vec_t eta0_sums[NQ1]; //Sums over eta0 for each value of eta1;

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpOut(m_nmTot);
        vec_t* jac_ptr;
        for (int e =0; e < this->m_nBlocks; ++e)
        {

            if (DEFORMED)
            {
                jac_ptr = &((*this->m_jac)[nqTot*e]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e]);
            }

            // Load and transpose data
            load_interleave(inptr, nqTot, tmpIn);

            IProductTriKernel<NM0, NM1, NQ0, NQ1, CORRECT, false,
                false, DEFORMED>(
                tmpIn, this->m_bdata[0], this->m_bdata[1],
                this->m_w[0], this->m_w[1], jac_ptr,
                eta0_sums, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr += nqBlocks;
            outptr += nmBlocks;
        }
    }
public:

    NekDouble Ndof() final
    {
        return m_nmTot * this->m_nElmt;
    }

private:
    int m_nmTot;
};

template<bool DEFORMED = false>
struct IProductHex : public IProduct, public Helper<3, DEFORMED>
{
    IProductHex(std::vector<LibUtilities::BasisSharedPtr> basis,
                   int nElmt)
        : IProduct(basis, nElmt),
          Helper<3, DEFORMED>(basis, nElmt),
          m_nmTot(LibUtilities::StdHexData::getNumberOfCoefficients(
                      this->m_nm[0], this->m_nm[1], this->m_nm[2]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<IProductHex<DEFORMED>>(basis, nElmt);
    }

    void operator()(const Array<OneD, const NekDouble> &in,
                                  Array<OneD,       NekDouble> &out)
    {
        auto nm0 = m_basis[0]->GetNumModes();
        auto nm1 = m_basis[1]->GetNumModes();
        auto nm2 = m_basis[2]->GetNumModes();
        ASSERTL0( nm0 == nm1 && nm0 == nm2,
            "IProductHex: anisotropy not implemented.");

        auto np0 = m_basis[0]->GetNumPoints();


        switch(nm0)
        {
            case 2:
                switch(np0)
                {
                    case 2: IProductHexImpl<2, 2, 2, 2, 2, 2>(in, out); break;
                    case 3: IProductHexImpl<2, 2, 2, 3, 3, 3>(in, out); break;
                    case 4: IProductHexImpl<2, 2, 2, 4, 4, 4>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductHex: # of modes / points combo not implemented.");
                } break;
            case 3:
                switch(np0)
                {
                    case 3: IProductHexImpl<3, 3, 3, 3, 3, 3>(in, out); break;
                    case 4: IProductHexImpl<3, 3, 3, 4, 4, 4>(in, out); break;
                    case 5: IProductHexImpl<3, 3, 3, 5, 5, 5>(in, out); break;
                    case 6: IProductHexImpl<3, 3, 3, 6, 6, 6>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductHex: # of modes / points combo not implemented.");
                } break;
            case 4:
                switch(np0)
                {
                    case 4: IProductHexImpl<4, 4, 4, 4, 4 ,4>(in, out); break;
                    case 5: IProductHexImpl<4, 4, 4, 5, 5 ,5>(in, out); break;
                    case 6: IProductHexImpl<4, 4, 4, 6, 6 ,6>(in, out); break;
                    case 7: IProductHexImpl<4, 4, 4, 7, 7 ,7>(in, out); break;
                    case 8: IProductHexImpl<4, 4, 4, 8, 8 ,8>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductHex: # of modes / points combo not implemented.");
                } break;
            case 5:
                switch(np0)
                {
                    case 5: IProductHexImpl<5, 5, 5, 5, 5, 5>(in, out); break;
                    case 6: IProductHexImpl<5, 5, 5, 6, 6, 6>(in, out); break;
                    case 7: IProductHexImpl<5, 5, 5, 7, 7, 7>(in, out); break;
                    case 8: IProductHexImpl<5, 5, 5, 8, 8, 8>(in, out); break;
                    case 9: IProductHexImpl<5, 5, 5, 9, 9, 9>(in, out); break;
                    case 10: IProductHexImpl<5, 5, 5, 10, 10, 10>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductHex: # of modes / points combo not implemented.");
                } break;
            case 6:
                switch(np0)
                {
                    case 6: IProductHexImpl<6, 6, 6, 6, 6, 6>(in, out); break;
                    case 7: IProductHexImpl<6, 6, 6, 7, 7, 7>(in, out); break;
                    case 8: IProductHexImpl<6, 6, 6, 8, 8, 8>(in, out); break;
                    case 9: IProductHexImpl<6, 6, 6, 9, 9, 9>(in, out); break;
                    case 10: IProductHexImpl<6, 6, 6, 10, 10, 10>(in, out); break;
                    case 11: IProductHexImpl<6, 6, 6, 11, 11, 11>(in, out); break;
                    case 12: IProductHexImpl<6, 6, 6, 12, 12, 12>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductHex: # of modes / points combo not implemented.");
                } break;
            case 7:
                switch(np0)
                {
                    case 7: IProductHexImpl<7, 7, 7, 7, 7, 7>(in, out); break;
                    case 8: IProductHexImpl<7, 7, 7, 8, 8, 8>(in, out); break;
                    case 9: IProductHexImpl<7, 7, 7, 9, 9, 9>(in, out); break;
                    case 10: IProductHexImpl<7, 7, 7, 10, 10, 10>(in, out); break;
                    case 11: IProductHexImpl<7, 7, 7, 11, 11, 11>(in, out); break;
                    case 12: IProductHexImpl<7, 7, 7, 12, 12, 12>(in, out); break;
                    case 13: IProductHexImpl<7, 7, 7, 13, 13, 13>(in, out); break;
                    case 14: IProductHexImpl<7, 7, 7, 14, 14, 14>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductHex: # of modes / points combo not implemented.");
                } break;
            case 8:
                switch(np0)
                {
                    case 8: IProductHexImpl<8, 8, 8, 8, 8, 8>(in, out); break;
                    case 9: IProductHexImpl<8, 8, 8, 9, 9, 9>(in, out); break;
                    case 10: IProductHexImpl<8, 8, 8, 10, 10, 10>(in, out); break;
                    case 11: IProductHexImpl<8, 8, 8, 11, 11, 11>(in, out); break;
                    case 12: IProductHexImpl<8, 8, 8, 12, 12, 12>(in, out); break;
                    case 13: IProductHexImpl<8, 8, 8, 13, 13, 13>(in, out); break;
                    case 14: IProductHexImpl<8, 8, 8, 14, 14, 14>(in, out); break;
                    case 15: IProductHexImpl<8, 8, 8, 15, 15, 15>(in, out); break;
                    case 16: IProductHexImpl<8, 8, 8, 16, 16, 16>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductHex: # of modes / points combo not implemented.");
                } break;;
            default: NEKERROR(ErrorUtil::efatal,
                "IProductHex: # of modes / points combo not implemented.");
        }
    }

    template<int NM0, int NM1, int NM2, int NQ0, int NQ1, int NQ2>
    void IProductHexImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &output)
    {
        using namespace tinysimd;
        using vec_t = simd<NekDouble>;

        auto* inptr = input.data();
        auto* outptr = output.data();

        constexpr auto nqTot = NQ0 * NQ1 * NQ2;
        constexpr auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        vec_t sums_kj[NQ1 * NQ2];
        vec_t sums_k[NQ2];

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpOut(m_nmTot);
        vec_t* jac_ptr;
        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            if (DEFORMED)
            {
                jac_ptr = &((*this->m_jac)[nqTot*e]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e]);
            }

            // Load and transpose data
            load_interleave(inptr, nqTot, tmpIn);

            IProductHexKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, false, false, DEFORMED>(
                tmpIn, this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
                this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                sums_kj, sums_k,
                tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr += nqBlocks;
            outptr += nmBlocks;
        }
    }

public:

    NekDouble Ndof() final
    {
        return m_nmTot * this->m_nElmt;
    }


private:
    /// Padded basis
    int m_nmTot;
};

template<bool DEFORMED = false>
struct IProductPrism : public IProduct, public Helper<3, DEFORMED>
{
    IProductPrism(std::vector<LibUtilities::BasisSharedPtr> basis,
                     int nElmt)
        : IProduct(basis, nElmt),
          Helper<3, DEFORMED>(basis, nElmt),
          m_nmTot(LibUtilities::StdPrismData::getNumberOfCoefficients(
                      this->m_nm[0], this->m_nm[1], this->m_nm[2]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<IProductPrism<DEFORMED>>(basis, nElmt);
    }

    void operator()(const Array<OneD, const NekDouble> &in,
        Array<OneD,       NekDouble> &out) final
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
                case 2: switch(m_basis[0]->GetNumPoints())
                {
                    case 3: IProductPrismImpl<2, 2, 2, 3, 3, 2, true>
                        (in, out); break;
                    case 4: IProductPrismImpl<2, 2, 2, 4, 4, 3, true>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductPrism: # of modes / points combo not implemented.");
                } break;
            case 3:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 4: IProductPrismImpl<3, 3, 3, 4, 4, 3, true>
                        (in, out); break;
                    case 5: IProductPrismImpl<3, 3, 3, 5, 5, 4, true>
                        (in, out); break;
                    case 6: IProductPrismImpl<3, 3, 3, 6, 6, 5, true>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductPrism: # of modes / points combo not implemented.");
                } break;
            case 4:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 5: IProductPrismImpl<4, 4, 4, 5, 5, 4, true>
                        (in, out); break;
                    case 6: IProductPrismImpl<4, 4, 4, 6, 6, 5, true>
                        (in, out); break;
                    case 7: IProductPrismImpl<4, 4, 4, 7, 7, 6, true>
                        (in, out); break;
                    case 8: IProductPrismImpl<4, 4, 4, 8, 8, 7, true>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductPrism: # of modes / points combo not implemented.");
                } break;
            case 5:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 6: IProductPrismImpl<5, 5, 5, 6, 6, 5, true>
                        (in, out); break;
                    case 7: IProductPrismImpl<5, 5, 5, 7, 7, 6, true>
                        (in, out); break;
                    case 8: IProductPrismImpl<5, 5, 5, 8, 8, 7, true>
                        (in, out); break;
                    case 9: IProductPrismImpl<5, 5, 5, 9, 9, 8, true>
                        (in, out); break;
                    case 10: IProductPrismImpl<5, 5, 5, 10, 10, 9, true>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductPrism: # of modes / points combo not implemented.");
                } break;
            case 6:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 7: IProductPrismImpl<6, 6, 6, 7, 7, 6, true>
                        (in, out); break;
                    case 8: IProductPrismImpl<6, 6, 6, 8, 8, 7, true>
                        (in, out); break;
                    case 9: IProductPrismImpl<6, 6, 6, 9, 9, 8, true>
                        (in, out); break;
                    case 10: IProductPrismImpl<6, 6, 6, 10, 10, 9, true>
                        (in, out); break;
                    case 11: IProductPrismImpl<6, 6, 6, 11, 11, 10, true>
                        (in, out); break;
                    case 12: IProductPrismImpl<6, 6, 6, 12, 12, 11, true>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductPrism: # of modes / points combo not implemented.");
                } break;
            case 7:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 8: IProductPrismImpl<7, 7, 7, 8, 8, 7, true>
                        (in, out); break;
                    case 9: IProductPrismImpl<7, 7, 7, 9, 9, 8, true>
                        (in, out); break;
                    case 10: IProductPrismImpl<7, 7, 7, 10, 10, 9, true>
                        (in, out); break;
                    case 11: IProductPrismImpl<7, 7, 7, 11, 11, 10, true>
                        (in, out); break;
                    case 12: IProductPrismImpl<7, 7, 7, 12, 12, 11, true>
                        (in, out); break;
                    case 13: IProductPrismImpl<7, 7, 7, 13, 13, 12, true>
                        (in, out); break;
                    case 14: IProductPrismImpl<7, 7, 7, 14, 14, 13, true>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductPrism: # of modes / points combo not implemented.");
                } break;
            case 8:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 9: IProductPrismImpl<8, 8, 8, 9, 9, 8, true>
                        (in, out); break;
                    case 10: IProductPrismImpl<8, 8, 8, 10, 10, 9, true>
                        (in, out); break;
                    case 11: IProductPrismImpl<8, 8, 8, 11, 11, 10, true>
                        (in, out); break;
                    case 12: IProductPrismImpl<8, 8, 8, 12, 12, 11, true>
                        (in, out); break;
                    case 13: IProductPrismImpl<8, 8, 8, 13, 13, 12, true>
                        (in, out); break;
                    case 14: IProductPrismImpl<8, 8, 8, 14, 14, 13, true>
                        (in, out); break;
                    case 15: IProductPrismImpl<8, 8, 8, 15, 15, 14, true>
                        (in, out); break;
                    case 16: IProductPrismImpl<8, 8, 8, 16, 16, 15, true>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductPrism: # of modes / points combo not implemented.");
                } break;
            default: NEKERROR(ErrorUtil::efatal,
                "IProductPrism: # of modes / points combo not implemented.");

            }
        }
        else
        {
            switch(m_basis[0]->GetNumModes())
            {
                case 2: switch(m_basis[0]->GetNumPoints())
                {
                    case 3: IProductPrismImpl<2, 2, 2, 3, 3, 2, false>
                        (in, out); break;
                    case 4: IProductPrismImpl<2, 2, 2, 4, 4, 3, false>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductPrism: # of modes / points combo not implemented.");
                } break;
            case 3:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 4: IProductPrismImpl<3, 3, 3, 4, 4, 3, false>
                        (in, out); break;
                    case 5: IProductPrismImpl<3, 3, 3, 5, 5, 4, false>
                        (in, out); break;
                    case 6: IProductPrismImpl<3, 3, 3, 6, 6, 5, false>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductPrism: # of modes / points combo not implemented.");
                } break;
            case 4:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 5: IProductPrismImpl<4, 4, 4, 5, 5, 4, false>
                        (in, out); break;
                    case 6: IProductPrismImpl<4, 4, 4, 6, 6, 5, false>
                        (in, out); break;
                    case 7: IProductPrismImpl<4, 4, 4, 7, 7, 6, false>
                        (in, out); break;
                    case 8: IProductPrismImpl<4, 4, 4, 8, 8, 7, false>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductPrism: # of modes / points combo not implemented.");
                } break;
            case 5:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 6: IProductPrismImpl<5, 5, 5, 6, 6, 5, false>
                        (in, out); break;
                    case 7: IProductPrismImpl<5, 5, 5, 7, 7, 6, false>
                        (in, out); break;
                    case 8: IProductPrismImpl<5, 5, 5, 8, 8, 7, false>
                        (in, out); break;
                    case 9: IProductPrismImpl<5, 5, 5, 9, 9, 8, false>
                        (in, out); break;
                    case 10: IProductPrismImpl<5, 5, 5, 10, 10, 9, false>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductPrism: # of modes / points combo not implemented.");
                } break;
            case 6:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 7: IProductPrismImpl<6, 6, 6, 7, 7, 6, false>
                        (in, out); break;
                    case 8: IProductPrismImpl<6, 6, 6, 8, 8, 7, false>
                        (in, out); break;
                    case 9: IProductPrismImpl<6, 6, 6, 9, 9, 8, false>
                        (in, out); break;
                    case 10: IProductPrismImpl<6, 6, 6, 10, 10, 9, false>
                        (in, out); break;
                    case 11: IProductPrismImpl<6, 6, 6, 11, 11, 10, false>
                        (in, out); break;
                    case 12: IProductPrismImpl<6, 6, 6, 12, 12, 11, false>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductPrism: # of modes / points combo not implemented.");
                } break;
            case 7:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 8: IProductPrismImpl<7, 7, 7, 8, 8, 7, false>
                        (in, out); break;
                    case 9: IProductPrismImpl<7, 7, 7, 9, 9, 8, false>
                        (in, out); break;
                    case 10: IProductPrismImpl<7, 7, 7, 10, 10, 9, false>
                        (in, out); break;
                    case 11: IProductPrismImpl<7, 7, 7, 11, 11, 10, false>
                        (in, out); break;
                    case 12: IProductPrismImpl<7, 7, 7, 12, 12, 11, false>
                        (in, out); break;
                    case 13: IProductPrismImpl<7, 7, 7, 13, 13, 12, false>
                        (in, out); break;
                    case 14: IProductPrismImpl<7, 7, 7, 14, 14, 13, false>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductPrism: # of modes / points combo not implemented.");
                } break;
            case 8:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 9: IProductPrismImpl<8, 8, 8, 9, 9, 8, false>
                        (in, out); break;
                    case 10: IProductPrismImpl<8, 8, 8, 10, 10, 9, false>
                        (in, out); break;
                    case 11: IProductPrismImpl<8, 8, 8, 11, 11, 10, false>
                        (in, out); break;
                    case 12: IProductPrismImpl<8, 8, 8, 12, 12, 11, false>
                        (in, out); break;
                    case 13: IProductPrismImpl<8, 8, 8, 13, 13, 12, false>
                        (in, out); break;
                    case 14: IProductPrismImpl<8, 8, 8, 14, 14, 13, false>
                        (in, out); break;
                    case 15: IProductPrismImpl<8, 8, 8, 15, 15, 14, false>
                        (in, out); break;
                    case 16: IProductPrismImpl<8, 8, 8, 16, 16, 15, false>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductPrism: # of modes / points combo not implemented.");
                } break;
            default: NEKERROR(ErrorUtil::efatal,
                "IProductPrism: # of modes / points combo not implemented.");

            }

        }
    }

    template<int NM0, int NM1, int NM2, int NQ0, int NQ1, int NQ2, bool CORRECT>
    void IProductPrismImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &output)
    {
        auto* inptr = input.data();
        auto* outptr = output.data();

        constexpr auto nqTot = NQ0 * NQ1 * NQ2;
        constexpr auto nqBlocks = NQ0 * NQ1 * NQ2 * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        vec_t sums_kj[NQ1 * NQ2];
        vec_t sums_k[NQ2];
        vec_t corr_q[NM1];

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpOut(m_nmTot);
        vec_t* jac_ptr;
        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            if (DEFORMED)
            {
                jac_ptr = &((*this->m_jac)[nqTot*e]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e]);
            }

            // Load and transpose data
            load_interleave(inptr, nqTot, tmpIn);

            IProductPrismKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, CORRECT, false,
                false, DEFORMED>(
                tmpIn, this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
                this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                sums_kj, sums_k,
                corr_q,
                tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr += nqBlocks;
            outptr += nmBlocks;
        }

    }
public:

    NekDouble Ndof() final
    {
        return m_nmTot * this->m_nElmt;
    }

private:
    /// Padded basis
    int m_nmTot;
};

template<bool DEFORMED = false>
struct IProductTet : public IProduct, public Helper<3, DEFORMED>
{
    IProductTet(std::vector<LibUtilities::BasisSharedPtr> basis,
                   int nElmt)
        : IProduct(basis, nElmt),
          Helper<3, DEFORMED>(basis, nElmt),
          m_nmTot(LibUtilities::StdTetData::getNumberOfCoefficients(
                      this->m_nm[0], this->m_nm[1], this->m_nm[2]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<IProductTet<DEFORMED>>(basis, nElmt);
    }

    void operator()(const Array<OneD, const NekDouble> &in,
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
                case 2: switch(m_basis[0]->GetNumPoints())
                {
                    case 3: IProductTetImpl<2, 2, 2, 3, 2, 2, true>
                        (in, out); break;
                    case 4: IProductTetImpl<2, 2, 2, 4, 3, 3, true>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductTet: # of modes / points combo not implemented.");
                } break;
            case 3:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 4: IProductTetImpl<3, 3, 3, 4, 3, 3, true>
                        (in, out); break;
                    case 5: IProductTetImpl<3, 3, 3, 5, 4, 4, true>
                        (in, out); break;
                    case 6: IProductTetImpl<3, 3, 3, 6, 5, 5, true>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductTet: # of modes / points combo not implemented.");
                } break;
            case 4:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 5: IProductTetImpl<4, 4, 4, 5, 4, 4, true>
                        (in, out); break;
                    case 6: IProductTetImpl<4, 4, 4, 6, 5, 5, true>
                        (in, out); break;
                    case 7: IProductTetImpl<4, 4, 4, 7, 6, 6, true>
                        (in, out); break;
                    case 8: IProductTetImpl<4, 4, 4, 8, 7, 7, true>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductTet: # of modes / points combo not implemented.");
                } break;
            case 5:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 6: IProductTetImpl<5, 5, 5, 6, 5, 5, true>
                        (in, out); break;
                    case 7: IProductTetImpl<5, 5, 5, 7, 6, 6, true>
                        (in, out); break;
                    case 8: IProductTetImpl<5, 5, 5, 8, 7, 7, true>
                        (in, out); break;
                    case 9: IProductTetImpl<5, 5, 5, 9, 8, 8, true>
                        (in, out); break;
                    case 10: IProductTetImpl<5, 5, 5, 10, 9, 9, true>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductTet: # of modes / points combo not implemented.");
                } break;
            case 6:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 7: IProductTetImpl<6, 6, 6, 7, 6, 6, true>
                        (in, out); break;
                    case 8: IProductTetImpl<6, 6, 6, 8, 7, 7, true>
                        (in, out); break;
                    case 9: IProductTetImpl<6, 6, 6, 9, 8, 8, true>
                        (in, out); break;
                    case 10: IProductTetImpl<6, 6, 6, 10, 9, 9, true>
                        (in, out); break;
                    case 11: IProductTetImpl<6, 6, 6, 11, 10, 10, true>
                        (in, out); break;
                    case 12: IProductTetImpl<6, 6, 6, 12, 11, 11, true>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductTet: # of modes / points combo not implemented.");
                } break;
            case 7:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 8: IProductTetImpl<7, 7, 7, 8, 7, 7, true>
                        (in, out); break;
                    case 9: IProductTetImpl<7, 7, 7, 9, 8, 8, true>
                        (in, out); break;
                    case 10: IProductTetImpl<7, 7, 7, 10, 9, 9, true>
                        (in, out); break;
                    case 11: IProductTetImpl<7, 7, 7, 11, 10, 10, true>
                        (in, out); break;
                    case 12: IProductTetImpl<7, 7, 7, 12, 11, 11, true>
                        (in, out); break;
                    case 13: IProductTetImpl<7, 7, 7, 13, 12, 12, true>
                        (in, out); break;
                    case 14: IProductTetImpl<7, 7, 7, 14, 13, 13, true>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductTet: # of modes / points combo not implemented.");
                } break;
            case 8:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 9: IProductTetImpl<8, 8, 8, 9, 8, 8, true>
                        (in, out); break;
                    case 10: IProductTetImpl<8, 8, 8, 10, 9, 9, true>
                        (in, out); break;
                    case 11: IProductTetImpl<8, 8, 8, 11, 10, 10, true>
                        (in, out); break;
                    case 12: IProductTetImpl<8, 8, 8, 12, 11, 11, true>
                        (in, out); break;
                    case 13: IProductTetImpl<8, 8, 8, 13, 12, 12, true>
                        (in, out); break;
                    case 14: IProductTetImpl<8, 8, 8, 14, 13, 13, true>
                        (in, out); break;
                    case 15: IProductTetImpl<8, 8, 8, 15, 14, 14, true>
                        (in, out); break;
                    case 16: IProductTetImpl<8, 8, 8, 16, 15, 15, true>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductTet: # of modes / points combo not implemented.");
                } break;
            default: NEKERROR(ErrorUtil::efatal,
                "IProductTet: # of modes / points combo not implemented.");

            }
        }
        else
        {
            switch(m_basis[0]->GetNumModes())
            {
                case 2: switch(m_basis[0]->GetNumPoints())
                {
                    case 3: IProductTetImpl<2, 2, 2, 3, 2, 2, false>
                        (in, out); break;
                    case 4: IProductTetImpl<2, 2, 2, 4, 3, 3, false>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductTet: # of modes / points combo not implemented.");
                } break;
            case 3:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 4: IProductTetImpl<3, 3, 3, 4, 3, 3, false>
                        (in, out); break;
                    case 5: IProductTetImpl<3, 3, 3, 5, 4, 4, false>
                        (in, out); break;
                    case 6: IProductTetImpl<3, 3, 3, 6, 5, 5, false>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductTet: # of modes / points combo not implemented.");
                } break;
            case 4:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 5: IProductTetImpl<4, 4, 4, 5, 4, 4, false>
                        (in, out); break;
                    case 6: IProductTetImpl<4, 4, 4, 6, 5, 5, false>
                        (in, out); break;
                    case 7: IProductTetImpl<4, 4, 4, 7, 6, 6, false>
                        (in, out); break;
                    case 8: IProductTetImpl<4, 4, 4, 8, 7, 7, false>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductTet: # of modes / points combo not implemented.");
                } break;
            case 5:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 6: IProductTetImpl<5, 5, 5, 6, 5, 5, false>
                        (in, out); break;
                    case 7: IProductTetImpl<5, 5, 5, 7, 6, 6, false>
                        (in, out); break;
                    case 8: IProductTetImpl<5, 5, 5, 8, 7, 7, false>
                        (in, out); break;
                    case 9: IProductTetImpl<5, 5, 5, 9, 8, 8, false>
                        (in, out); break;
                    case 10: IProductTetImpl<5, 5, 5, 10, 9, 9, false>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductTet: # of modes / points combo not implemented.");
                } break;
            case 6:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 7: IProductTetImpl<6, 6, 6, 7, 6, 6, false>
                        (in, out); break;
                    case 8: IProductTetImpl<6, 6, 6, 8, 7, 7, false>
                        (in, out); break;
                    case 9: IProductTetImpl<6, 6, 6, 9, 8, 8, false>
                        (in, out); break;
                    case 10: IProductTetImpl<6, 6, 6, 10, 9, 9, false>
                        (in, out); break;
                    case 11: IProductTetImpl<6, 6, 6, 11, 10, 10, false>
                        (in, out); break;
                    case 12: IProductTetImpl<6, 6, 6, 12, 11, 11, false>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductTet: # of modes / points combo not implemented.");
                } break;
            case 7:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 8: IProductTetImpl<7, 7, 7, 8, 7, 7, false>
                        (in, out); break;
                    case 9: IProductTetImpl<7, 7, 7, 9, 8, 8, false>
                        (in, out); break;
                    case 10: IProductTetImpl<7, 7, 7, 10, 9, 9, false>
                        (in, out); break;
                    case 11: IProductTetImpl<7, 7, 7, 11, 10, 10, false>
                        (in, out); break;
                    case 12: IProductTetImpl<7, 7, 7, 12, 11, 11, false>
                        (in, out); break;
                    case 13: IProductTetImpl<7, 7, 7, 13, 12, 12, false>
                        (in, out); break;
                    case 14: IProductTetImpl<7, 7, 7, 14, 13, 13, false>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductTet: # of modes / points combo not implemented.");
                } break;
            case 8:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 9: IProductTetImpl<8, 8, 8, 9, 8, 8, false>
                        (in, out); break;
                    case 10: IProductTetImpl<8, 8, 8, 10, 9, 9, false>
                        (in, out); break;
                    case 11: IProductTetImpl<8, 8, 8, 11, 10, 10, false>
                        (in, out); break;
                    case 12: IProductTetImpl<8, 8, 8, 12, 11, 11, false>
                        (in, out); break;
                    case 13: IProductTetImpl<8, 8, 8, 13, 12, 12, false>
                        (in, out); break;
                    case 14: IProductTetImpl<8, 8, 8, 14, 13, 13, false>
                        (in, out); break;
                    case 15: IProductTetImpl<8, 8, 8, 15, 14, 14, false>
                        (in, out); break;
                    case 16: IProductTetImpl<8, 8, 8, 16, 15, 15, false>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "IProductTet: # of modes / points combo not implemented.");
                } break;
            default: NEKERROR(ErrorUtil::efatal,
                "IProductTet: # of modes / points combo not implemented.");

            }

        }
    }

    template<int NM0, int NM1, int NM2, int NQ0, int NQ1, int NQ2, bool CORRECT>
    void IProductTetImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &output)
    {
        auto* inptr = input.data();
        auto* outptr = output.data();

        constexpr auto nqTot = NQ0 * NQ1 * NQ2;
        constexpr auto nqBlocks = NQ0 * NQ1 * NQ2 * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        vec_t wsp[NQ1 * NQ2 + NQ2];

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpOut(m_nmTot);
        vec_t* jac_ptr;
        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            if (DEFORMED)
            {
                jac_ptr = &((*this->m_jac)[nqTot*e]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e]);
            }

            // Load and transpose data
            load_interleave(inptr, nqTot, tmpIn);

            IProductTetKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, CORRECT,
                                 false, false, DEFORMED>(
                tmpIn, this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
                this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                wsp, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr += nqBlocks;
            outptr += nmBlocks;
        }
    }
public:

    NekDouble Ndof() final
    {
        return m_nmTot * this->m_nElmt;
    }

private:
    /// Padded basis
    int m_nmTot;
};

} // namespace MatrixFree
} // namespace Nektar

#endif
