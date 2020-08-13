#ifndef MF_PHYSDERIV_H
#define MF_PHYSDERIV_H

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ShapeType.hpp>
#include <LibUtilities/Foundations/Basis.h>

#include "Operator.hpp"
#include "PhysDerivKernels.hpp"

namespace Nektar
{
namespace MatrixFree
{

template<bool DEFORMED = false>
struct PhysDerivQuad : public PhysDeriv, public Helper<2, DEFORMED>
{
    PhysDerivQuad(std::vector<LibUtilities::BasisSharedPtr> basis,
                     int nElmt)
    : PhysDeriv(basis, nElmt),
      Helper<2, DEFORMED>(basis, nElmt),
      m_nmTot(LibUtilities::StdQuadData::getNumberOfCoefficients(
                this->m_nm[0], this->m_nm[1]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<PhysDerivQuad<DEFORMED>>(basis, nElmt);
    }

    static NekDouble FlopsPerElement(
        const int nq0,
        const int nq1)
    {
        int derivTensor = 2 * nq0 *nq1 *nq0 + 2 * nq0 * nq1 *nq1;
        int deriv = nq0 *nq1 * 6;
        return (derivTensor + deriv);
    }

    double GFlops() final
    {
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();

        int flops = this->m_nElmt * PhysDerivQuad::FlopsPerElement(nq0, nq1);
        return flops * 1e-9;
    }

    NekDouble Ndof() final
    {
        return m_nmTot * this->m_nElmt;
    }

    NekDouble NLoads() final
    {
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();

        int t_d0 = nq0 * nq1 * nq0 * 2;
        int t_d1 = nq0 * nq1 * nq1 * 2;
        int physDerivTensor = t_d0 + t_d1;

        int physDeriv = nq0 * nq1 * 2;
        int load_expected = physDerivTensor + physDeriv;

        return this->m_nElmt * load_expected;
    }

    NekDouble NStores() final
    {
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();

        int physDerivTensor = nq0*nq1*2;
        int physDeriv = nq0 * nq1 * 2;

        return this->m_nElmt * (physDerivTensor + physDeriv);
    }

    void operator()(const Array<OneD, const NekDouble> &in,
                          Array<OneD,       NekDouble> &out_d0,
                          Array<OneD,       NekDouble> &out_d1,
                          Array<OneD,       NekDouble> &out_d2) final
    {
        //Only for 3D, but need to implement since its abstract
        boost::ignore_unused(in, out_d0, out_d1, out_d2);
        NEKERROR(ErrorUtil::efatal,
                "Something went horribly wrong... calling 3D op for 2D op");
    }

    void operator()(const Array<OneD, const NekDouble> &in,
                          Array<OneD,       NekDouble> &out_d0,
                          Array<OneD,       NekDouble> &out_d1) final
    {
        // Check preconditions
        ASSERTL0(m_basis[0]->GetNumPoints() == m_basis[1]->GetNumPoints(),
            "MatrixFree op requires homogenous points");

        switch(m_basis[0]->GetNumPoints())
        {
            case 2:
                PhysDerivQuadImpl<2,2>(in, out_d0, out_d1); break;
            case 3:
                PhysDerivQuadImpl<3,3>(in, out_d0, out_d1); break;
            case 4:
                PhysDerivQuadImpl<4,4>(in, out_d0, out_d1); break;
            case 5:
                PhysDerivQuadImpl<5,5>(in, out_d0, out_d1); break;
            case 6:
                PhysDerivQuadImpl<6,6>(in, out_d0, out_d1); break;
            case 7:
                PhysDerivQuadImpl<7,7>(in, out_d0, out_d1); break;
            case 8:
                PhysDerivQuadImpl<8,8>(in, out_d0, out_d1); break;
            case 9:
                PhysDerivQuadImpl<9,9>(in, out_d0, out_d1); break;
            case 10:
                PhysDerivQuadImpl<10,10>(in, out_d0, out_d1); break;
            case 11:
                PhysDerivQuadImpl<11,11>(in, out_d0, out_d1); break;
            case 12:
                PhysDerivQuadImpl<12,12>(in, out_d0, out_d1); break;
            case 13:
                PhysDerivQuadImpl<13,13>(in, out_d0, out_d1); break;
            case 14:
                PhysDerivQuadImpl<14,14>(in, out_d0, out_d1); break;
            case 15:
                PhysDerivQuadImpl<15,15>(in, out_d0, out_d1); break;
            case 16:
                PhysDerivQuadImpl<16,16>(in, out_d0, out_d1); break;
            default: NEKERROR(ErrorUtil::efatal,
                "PhysDerivQuad: # of modes / points combo not implemented.");
        }
    }

    template<int NQ0, int NQ1>
    void PhysDerivQuadImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &out_d0,
              Array<OneD,       NekDouble> &out_d1)
    {
        auto* inptr = &input[0];
        auto* outptr_d0 = &out_d0[0];
        auto* outptr_d1 = &out_d1[0];

        constexpr auto ndf = 4;
        constexpr auto nqTot = NQ0 * NQ1;
        constexpr auto nqBlocks = nqTot * vec_t::width;

        // Get size of derivative factor block
        auto dfSize = ndf;
        if (DEFORMED)
        {
            dfSize *= nqTot;
        }

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpOut_d0(nqTot),
            tmpOut_d1(nqTot);
        const vec_t* df_ptr;
        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            df_ptr = &(this->m_df[e*dfSize]);

            // Load and transpose data
            load_interleave(inptr, nqTot, tmpIn);

            PhysDerivQuadKernel<NQ0, NQ1, DEFORMED>(
                tmpIn,
                this->m_Z[0], this->m_Z[1],
                this->m_D[0], this->m_D[1],
                df_ptr,
                tmpOut_d0, tmpOut_d1);

            // de-interleave and store data
            deinterleave_store(tmpOut_d0, nqTot, outptr_d0);
            deinterleave_store(tmpOut_d1, nqTot, outptr_d1);

            inptr += nqBlocks;
            outptr_d0 += nqBlocks;
            outptr_d1 += nqBlocks;
        }
    }

private:
    int m_nmTot;

};

template<bool DEFORMED=false>
struct PhysDerivTri : public PhysDeriv, public Helper<2,DEFORMED>
{
    PhysDerivTri(std::vector<LibUtilities::BasisSharedPtr> basis,
                    int nElmt)
        : PhysDeriv(basis, nElmt),
          Helper<2, DEFORMED>(basis, nElmt),
          m_nmTot(LibUtilities::StdTriData::getNumberOfCoefficients(
                      this->m_nm[0], this->m_nm[1]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<PhysDerivTri<DEFORMED>>(basis, nElmt);
    }

    static NekDouble FlopsPerElement(
        const int nq0,
        const int nq1)
    {
        int derivTensor = 2 * nq0 *nq1 *nq0 + 2 * nq0 * nq1 *nq1;

        int d0 = nq1 * 2 + nq1 * nq0;
        int d1 = nq1 * nq0 * 4;
        int df = nq1 * nq0 * 6;

        return derivTensor + d0 + d1 + df;
    }

    double GFlops() final
    {
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();

        int flops = m_nElmt * PhysDerivTri::FlopsPerElement(nq0, nq1);
        return flops * 1e-9;
    }

    NekDouble Ndof() final
    {
        return m_nmTot * this->m_nElmt;
    }


    NekDouble NLoads() final
    {
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();

        int t_d0 = nq0 * nq1 * nq0 * 2;
        int t_d1 = nq0 * nq1 * nq1 * 2;
        int physDerivTensor = t_d0 + t_d1;
        int physDeriv = nq1*(1 + nq0*3);

        int load_expected = physDerivTensor + physDeriv;

        return this->m_nElmt * load_expected;
    }

    NekDouble NStores() final
    {
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();

        int physDerivTensor = nq0*nq1*2;
        int physDeriv = nq1*nq0*2;

        return this->m_nElmt * (physDerivTensor + physDeriv);
    }

    void operator()(const Array<OneD, const NekDouble> &in,
                          Array<OneD,       NekDouble> &out_d0,
                          Array<OneD,       NekDouble> &out_d1,
                          Array<OneD,       NekDouble> &out_d2) final
    {
        // Only for 3D, but need to implement since its abstract
        boost::ignore_unused(in, out_d0, out_d1, out_d2);
        NEKERROR(ErrorUtil::efatal,
                "Something went horribly wrong... calling 3D op for 2D element");
    }

    void operator()(const Array<OneD, const NekDouble> &in,
                          Array<OneD,       NekDouble> &out_d0,
                          Array<OneD,       NekDouble> &out_d1)  final
    {
        // Check preconditions
        ASSERTL0(m_basis[0]->GetNumPoints() == (m_basis[1]->GetNumPoints()+1),
            "MatrixFree op requires homogenous points");

        switch(m_basis[0]->GetNumPoints())
        {
            case 2:
                PhysDerivTriImpl<2,1>(in, out_d0, out_d1); break;
            case 3:
                PhysDerivTriImpl<3,2>(in, out_d0, out_d1); break;
            case 4:
                PhysDerivTriImpl<4,3>(in, out_d0, out_d1); break;
            case 5:
                PhysDerivTriImpl<5,4>(in, out_d0, out_d1); break;
            case 6:
                PhysDerivTriImpl<6,5>(in, out_d0, out_d1); break;
            case 7:
                PhysDerivTriImpl<7,6>(in, out_d0, out_d1); break;
            case 8:
                PhysDerivTriImpl<8,7>(in, out_d0, out_d1); break;
            case 9:
                PhysDerivTriImpl<9,8>(in, out_d0, out_d1); break;
            case 10:
                PhysDerivTriImpl<10,9>(in, out_d0, out_d1); break;
            case 11:
                PhysDerivTriImpl<11,10>(in, out_d0, out_d1); break;
            case 12:
                PhysDerivTriImpl<12,11>(in, out_d0, out_d1); break;
            case 13:
                PhysDerivTriImpl<13,12>(in, out_d0, out_d1); break;
            case 14:
                PhysDerivTriImpl<14,13>(in, out_d0, out_d1); break;
            case 15:
                PhysDerivTriImpl<15,14>(in, out_d0, out_d1); break;
            case 16:
                PhysDerivTriImpl<16,15>(in, out_d0, out_d1); break;
            case 17:
                PhysDerivTriImpl<17,16>(in, out_d0, out_d1); break;
            default: NEKERROR(ErrorUtil::efatal,
                "PhysDerivTri: # of modes / points combo not implemented.");
        }
    }

    template<int NQ0, int NQ1>
    void PhysDerivTriImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &out_d0,
              Array<OneD,       NekDouble> &out_d1)
    {
        auto* inptr = &input[0];
        auto* outptr_d0 = &out_d0[0];
        auto* outptr_d1 = &out_d1[0];

        constexpr auto ndf = 4;
        constexpr auto nqTot = NQ0 * NQ1;
        constexpr auto nqBlocks = nqTot * vec_t::width;

        // Get size of derivative factor block
        auto dfSize = ndf;
        if (DEFORMED)
        {
            dfSize *= nqTot;
        }

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpOut_d0(nqTot),
            tmpOut_d1(nqTot);
        const vec_t* df_ptr;
        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            df_ptr = &(this->m_df[e*dfSize]);

            // Load and transpose data
            load_interleave(inptr, nqTot, tmpIn);

            PhysDerivTriKernel<NQ0, NQ1, DEFORMED>(
                tmpIn,
                this->m_Z[0], this->m_Z[1],
                this->m_D[0], this->m_D[1],
                df_ptr,
                tmpOut_d0, tmpOut_d1);

            // de-interleave and store data
            deinterleave_store(tmpOut_d0, nqTot, outptr_d0);
            deinterleave_store(tmpOut_d1, nqTot, outptr_d1);

            inptr += nqBlocks;
            outptr_d0 += nqBlocks;
            outptr_d1 += nqBlocks;
        }
    }
private:
    int m_nmTot;
};

template<bool DEFORMED = false>
struct PhysDerivHex : public PhysDeriv, public Helper<3,DEFORMED>
{
    PhysDerivHex(std::vector<LibUtilities::BasisSharedPtr> basis,
                    int nElmt)
        : PhysDeriv(basis, nElmt),
          Helper<3,DEFORMED>(basis, nElmt),
          m_nmTot(LibUtilities::StdHexData::getNumberOfCoefficients(
                      this->m_nm[0], this->m_nm[1], this->m_nm[2]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<PhysDerivHex<DEFORMED>>(basis, nElmt);
    }

    static NekDouble FlopsPerElement(
        const int nq0,
        const int nq1,
        const int nq2)
    {
        //PhysDerivTensor
        int pdt0 = nq0 * nq1 * nq2 * nq0 * 2;
        int pdt1 = nq2 * nq0 * nq1 * nq1 * 2;
        int pdt2 = nq0 * nq1 * nq2 * nq2 * 2;
        int physDerivTensor = pdt0 + pdt1 + pdt2;

        //PhysDeriv
        int physDeriv = nq2 * nq1 * nq0 * 15;

        return (physDerivTensor + physDeriv);
    }

    double GFlops() final
    {
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        int flops = this->m_nElmt * PhysDerivHex::FlopsPerElement(nq0, nq1, nq2);
        return flops * 1e-9;
    }

    NekDouble Ndof() final
    {
        return m_nmTot * this->m_nElmt;
    }

    NekDouble NLoads() final
    {
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        int load_d0 = nq0 * nq1 * nq2 * nq0 * 2;
        int load_d1 = nq2 * nq0 * nq1 * nq1 * 2;
        int load_d2 = nq0 * nq1 * nq2 * nq2 * 2;;
        int physDerivTensor = load_d0 + load_d1 + load_d2;
        int physDeriv = nq2 * nq1 * nq0 * 3;

        return m_nElmt * (physDeriv + physDerivTensor);

    }

    NekDouble NStores() final
    {
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        int store_d0 = nq0 * nq1 * nq2;
        int store_d1 = nq2 * nq0 * nq1;
        int store_d2 = nq0 * nq1 * nq2;
        int physDerivTensor = store_d0 + store_d1 + store_d2;
        int physDeriv = nq2 * nq1 * nq0 * 3;

        return m_nElmt * (physDeriv + physDerivTensor);
    }

    void operator()(const Array<OneD, const NekDouble> &in,
                                Array<OneD,       NekDouble> &out_d0,
                                Array<OneD,       NekDouble> &out_d1) final
    {
        boost::ignore_unused(in, out_d0, out_d1);
        throw; //Only for 2D, but need to implement since its abstract
    }

    void operator()(const Array<OneD, const NekDouble> &in,
                                  Array<OneD,       NekDouble> &out_d0,
                                  Array<OneD,       NekDouble> &out_d1,
                                  Array<OneD,       NekDouble> &out_d2) final
    {
        // Check preconditions
        ASSERTL0(m_basis[0]->GetNumPoints() == m_basis[1]->GetNumPoints() &&
            m_basis[0]->GetNumPoints() == m_basis[2]->GetNumPoints(),
            "MatrixFree op requires homogenous points");

        switch(m_basis[0]->GetNumPoints())
        {
            case 2:
                PhysDerivHexImpl<2,2,2>(in, out_d0, out_d1, out_d2); break;
            case 3:
                PhysDerivHexImpl<3,3,3>(in, out_d0, out_d1, out_d2); break;
            case 4:
                PhysDerivHexImpl<4,4,4>(in, out_d0, out_d1, out_d2); break;
            case 5:
                PhysDerivHexImpl<5,5,5>(in, out_d0, out_d1, out_d2); break;
            case 6:
                PhysDerivHexImpl<6,6,6>(in, out_d0, out_d1, out_d2); break;
            case 7:
                PhysDerivHexImpl<7,7,7>(in, out_d0, out_d1, out_d2); break;
            case 8:
                PhysDerivHexImpl<8,8,8>(in, out_d0, out_d1, out_d2); break;
            case 9:
                PhysDerivHexImpl<9,9,9>(in, out_d0, out_d1, out_d2); break;
            case 10:
                PhysDerivHexImpl<10,10,10>(in, out_d0, out_d1, out_d2); break;
            case 11:
                PhysDerivHexImpl<11,11,11>(in, out_d0, out_d1, out_d2); break;
            case 12:
                PhysDerivHexImpl<12,12,12>(in, out_d0, out_d1, out_d2); break;
            case 13:
                PhysDerivHexImpl<13,13,13>(in, out_d0, out_d1, out_d2); break;
            case 14:
                PhysDerivHexImpl<14,14,14>(in, out_d0, out_d1, out_d2); break;
            case 15:
                PhysDerivHexImpl<15,15,15>(in, out_d0, out_d1, out_d2); break;
            case 16:
                PhysDerivHexImpl<16,16,16>(in, out_d0, out_d1, out_d2); break;
            default: NEKERROR(ErrorUtil::efatal,
                "PhysDerivHex: # of modes / points combo not implemented.");
        }
    }

    template<int NQ0, int NQ1, int NQ2>
    void PhysDerivHexImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &out_d0,
              Array<OneD,       NekDouble> &out_d1,
              Array<OneD,       NekDouble> &out_d2)
    {
        using namespace tinysimd;
        using vec_t = simd<NekDouble>;

        auto* inptr = &input[0];
        auto* outptr_d0 = &out_d0[0];
        auto* outptr_d1 = &out_d1[0];
        auto* outptr_d2 = &out_d2[0];

        constexpr auto ndf = 9;
        constexpr auto nqTot = NQ0 * NQ1 * NQ2;
        constexpr auto nqBlocks = nqTot * vec_t::width;

        // Get size of derivative factor block
        int dfSize{};
        if (DEFORMED)
        {
            dfSize = ndf*nqTot;
        }
        else
        {
            dfSize = ndf;
        }

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpd0(nqTot),
            tmpd1(nqTot), tmpd2(nqTot);
        const vec_t *df_ptr;
        for (int e = 0; e < this->m_nBlocks; ++e)
        {

            df_ptr = &(this->m_df[e*dfSize]);

            // Load and transpose data
            load_interleave(inptr, nqTot, tmpIn);

            PhysDerivHexKernel<NQ0, NQ1, NQ2, DEFORMED>(
                tmpIn,
                this->m_Z[0], this->m_Z[1], this->m_Z[2],
                this->m_D[0], this->m_D[1], this->m_D[2],
                df_ptr,
                tmpd0, tmpd1, tmpd2);

            // de-interleave and store data
            deinterleave_store(tmpd0, nqTot, outptr_d0);
            deinterleave_store(tmpd1, nqTot, outptr_d1);
            deinterleave_store(tmpd2, nqTot, outptr_d2);

            inptr += nqBlocks;
            outptr_d0 += nqBlocks;
            outptr_d1 += nqBlocks;
            outptr_d2 += nqBlocks;
        }
    }
private:
    int m_nmTot;
};


template<bool DEFORMED = false>
struct PhysDerivPrism : public PhysDeriv, public Helper<3, DEFORMED>
{
    PhysDerivPrism(std::vector<LibUtilities::BasisSharedPtr> basis,
                      int nElmt)
        : PhysDeriv(basis, nElmt),
          Helper<3, DEFORMED>(basis, nElmt),
          m_nmTot(LibUtilities::StdPrismData::getNumberOfCoefficients(
                      this->m_nm[0], this->m_nm[1], this->m_nm[2]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<PhysDerivPrism<DEFORMED>>(basis, nElmt);
    }

    NekDouble Ndof() final
    {
        return m_nmTot * this->m_nElmt;
    }

    static NekDouble FlopsPerElement(
        const int nq0,
        const int nq1,
        const int nq2)
    {
        //PhysDerivTensor
        int pdt0 = nq0 * nq1 * nq2 * nq0 * 2;
        int pdt1 = nq2 * nq0 * nq1 * nq1 * 2;
        int pdt2 = nq0 * nq1 * nq2 * nq2 * 2;
        int physDerivTensor = pdt0 + pdt1 + pdt2;

        //PhysDeriv
        int physDeriv = nq2*(2 + nq1*nq0*(20));

        return (physDeriv + physDerivTensor);
    }

    double GFlops() final
    {
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        int flops = m_nElmt * PhysDerivPrism::FlopsPerElement(nq0, nq1, nq2);
        return flops * 1e-9;
    }

    NekDouble NLoads() final
    {
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        int load_d0 = nq0 * nq1 * nq2 * nq0 * 2;
        int load_d1 = nq2 * nq0 * nq1 * nq1 * 2;
        int load_d2 = nq0 * nq1 * nq2 * nq2 * 2;;
        int physDerivTensor = load_d0 + load_d1 + load_d2;
        int physDeriv = nq2 * (1 + nq1 * nq0 *4);

        return m_nElmt * (physDeriv + physDerivTensor);

    }

    NekDouble NStores() final
    {
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        int store_d0 = nq0 * nq1 * nq2;
        int store_d1 = nq2 * nq0 * nq1;
        int store_d2 = nq0 * nq1 * nq2;
        int physDerivTensor = store_d0 + store_d1 + store_d2;
        int physDeriv = nq2*nq1*nq0*3;

        return m_nElmt * (physDeriv + physDerivTensor);
    }

    void operator()(const Array<OneD, const NekDouble> &in,
                          Array<OneD,       NekDouble> &out_d0,
                          Array<OneD,       NekDouble> &out_d1) final
    {
        // Only for 2D, but need to implement since its abstract
        boost::ignore_unused(in, out_d0, out_d1);
        NEKERROR(ErrorUtil::efatal,
                "Something went horribly wrong... calling 2D op for 3D element");
    }

    void operator()(const Array<OneD, const NekDouble> &in,
                          Array<OneD,       NekDouble> &out_d0,
                          Array<OneD,       NekDouble> &out_d1,
                          Array<OneD,       NekDouble> &out_d2) final
    {
         // Check preconditions
        ASSERTL0(m_basis[0]->GetNumPoints() == m_basis[1]->GetNumPoints() &&
            m_basis[0]->GetNumPoints() == (m_basis[2]->GetNumPoints()+1),
            "MatrixFree op requires homogenous points");

        switch(m_basis[0]->GetNumPoints())
        {
            case 2:
                PhysDerivPrismImpl<2,2,1>(in, out_d0, out_d1, out_d2); break;
            case 3:
                PhysDerivPrismImpl<3,3,2>(in, out_d0, out_d1, out_d2); break;
            case 4:
                PhysDerivPrismImpl<4,4,3>(in, out_d0, out_d1, out_d2); break;
            case 5:
                PhysDerivPrismImpl<5,5,4>(in, out_d0, out_d1, out_d2); break;
            case 6:
                PhysDerivPrismImpl<6,6,5>(in, out_d0, out_d1, out_d2); break;
            case 7:
                PhysDerivPrismImpl<7,7,6>(in, out_d0, out_d1, out_d2); break;
            case 8:
                PhysDerivPrismImpl<8,8,7>(in, out_d0, out_d1, out_d2); break;
            case 9:
                PhysDerivPrismImpl<9,9,8>(in, out_d0, out_d1, out_d2); break;
            case 10:
                PhysDerivPrismImpl<10,10,9>(in, out_d0, out_d1, out_d2); break;
            case 11:
                PhysDerivPrismImpl<11,11,10>(in, out_d0, out_d1, out_d2); break;
            case 12:
                PhysDerivPrismImpl<12,12,11>(in, out_d0, out_d1, out_d2); break;
            case 13:
                PhysDerivPrismImpl<13,13,12>(in, out_d0, out_d1, out_d2); break;
            case 14:
                PhysDerivPrismImpl<14,14,13>(in, out_d0, out_d1, out_d2); break;
            case 15:
                PhysDerivPrismImpl<15,15,14>(in, out_d0, out_d1, out_d2); break;
            case 16:
                PhysDerivPrismImpl<16,16,15>(in, out_d0, out_d1, out_d2); break;
            default: NEKERROR(ErrorUtil::efatal,
                "PhysDerivPrism: # of modes / points combo not implemented.");

        }
    }

    template<int NQ0, int NQ1, int NQ2>
    void PhysDerivPrismImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &out_d0,
              Array<OneD,       NekDouble> &out_d1,
              Array<OneD,       NekDouble> &out_d2)
    {
        const auto* inptr = input.data();
        auto* outptr_d0 = out_d0.data();
        auto* outptr_d1 = out_d1.data();
        auto* outptr_d2 = out_d2.data();

        constexpr auto ndf = 9;
        constexpr auto nqTot = NQ0 * NQ1 * NQ2;
        constexpr auto nqBlocks = nqTot * vec_t::width;

        // Get size of derivative factor block
        auto dfSize = ndf;
        if (DEFORMED)
        {
            dfSize *= nqTot;
        }

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpOut_d0(nqTot),
            tmpOut_d1(nqTot), tmpOut_d2(nqTot);
        const vec_t* df_ptr;
        for (int e = 0; e < this->m_nBlocks; ++e)
        {

            df_ptr = &(this->m_df[e*dfSize]);

            // Load and transpose data
            load_interleave(inptr, nqTot, tmpIn);

            PhysDerivPrismKernel<NQ0, NQ1, NQ2, DEFORMED>(
                tmpIn,
                this->m_Z[0], this->m_Z[1], this->m_Z[2],
                this->m_D[0], this->m_D[1], this->m_D[2],
                df_ptr,
                tmpOut_d0, tmpOut_d1, tmpOut_d2);

            // de-interleave and store data
            deinterleave_store(tmpOut_d0, nqTot, outptr_d0);
            deinterleave_store(tmpOut_d1, nqTot, outptr_d1);
            deinterleave_store(tmpOut_d2, nqTot, outptr_d2);

            inptr += nqBlocks;
            outptr_d0 += nqBlocks;
            outptr_d1 += nqBlocks;
            outptr_d2 += nqBlocks;
        }
    }
private:
    int m_nmTot;

};

template<bool DEFORMED = false>
struct PhysDerivPyr : public PhysDeriv, public Helper<3, DEFORMED>
{
    PhysDerivPyr(std::vector<LibUtilities::BasisSharedPtr> basis,
                      int nElmt)
        : PhysDeriv(basis, nElmt),
          Helper<3, DEFORMED>(basis, nElmt),
          m_nmTot(LibUtilities::StdPyrData::getNumberOfCoefficients(
                      this->m_nm[0], this->m_nm[1], this->m_nm[2]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<PhysDerivPyr<DEFORMED>>(basis, nElmt);
    }

    NekDouble Ndof() final
    {
        return m_nmTot * this->m_nElmt;
    }

    static NekDouble FlopsPerElement(
        const int nq0,
        const int nq1,
        const int nq2)
    {
        //PhysDerivTensor
        int pdt0 = nq0 * nq1 * nq2 * nq0 * 2;
        int pdt1 = nq2 * nq0 * nq1 * nq1 * 2;
        int pdt2 = nq0 * nq1 * nq2 * nq2 * 2;
        int physDerivTensor = pdt0 + pdt1 + pdt2;

        //PhysDeriv
        int physDeriv = nq2*(2 + nq1*nq0*(20));

        return (physDeriv + physDerivTensor);
    }

    double GFlops() final
    {
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        int flops = m_nElmt * PhysDerivPyr::FlopsPerElement(nq0, nq1, nq2);
        return flops * 1e-9;
    }

    NekDouble NLoads() final
    {
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        int load_d0 = nq0 * nq1 * nq2 * nq0 * 2;
        int load_d1 = nq2 * nq0 * nq1 * nq1 * 2;
        int load_d2 = nq0 * nq1 * nq2 * nq2 * 2;;
        int physDerivTensor = load_d0 + load_d1 + load_d2;
        int physDeriv = nq2 * (1 + nq1 * nq0 *4);

        return m_nElmt * (physDeriv + physDerivTensor);

    }

    NekDouble NStores() final
    {
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        int store_d0 = nq0 * nq1 * nq2;
        int store_d1 = nq2 * nq0 * nq1;
        int store_d2 = nq0 * nq1 * nq2;
        int physDerivTensor = store_d0 + store_d1 + store_d2;
        int physDeriv = nq2*nq1*nq0*3;

        return m_nElmt * (physDeriv + physDerivTensor);
    }

    void operator()(const Array<OneD, const NekDouble> &in,
                          Array<OneD,       NekDouble> &out_d0,
                          Array<OneD,       NekDouble> &out_d1) final
    {
        // Only for 2D, but need to implement since its abstract
        boost::ignore_unused(in, out_d0, out_d1);
        NEKERROR(ErrorUtil::efatal,
                "Something went horribly wrong... calling 2D op for 3D element");
    }

    void operator()(const Array<OneD, const NekDouble> &in,
                          Array<OneD,       NekDouble> &out_d0,
                          Array<OneD,       NekDouble> &out_d1,
                          Array<OneD,       NekDouble> &out_d2) final
    {
         // Check preconditions
        ASSERTL0(m_basis[0]->GetNumPoints() == m_basis[1]->GetNumPoints() &&
            m_basis[0]->GetNumPoints() == (m_basis[2]->GetNumPoints()+1),
            "MatrixFree op requires homogenous points");

        switch(m_basis[0]->GetNumPoints())
        {
            case 2:
                PhysDerivPyrImpl<2,2,1>(in, out_d0, out_d1, out_d2); break;
            case 3:
                PhysDerivPyrImpl<3,3,2>(in, out_d0, out_d1, out_d2); break;
            case 4:
                PhysDerivPyrImpl<4,4,3>(in, out_d0, out_d1, out_d2); break;
            case 5:
                PhysDerivPyrImpl<5,5,4>(in, out_d0, out_d1, out_d2); break;
            case 6:
                PhysDerivPyrImpl<6,6,5>(in, out_d0, out_d1, out_d2); break;
            case 7:
                PhysDerivPyrImpl<7,7,6>(in, out_d0, out_d1, out_d2); break;
            case 8:
                PhysDerivPyrImpl<8,8,7>(in, out_d0, out_d1, out_d2); break;
            case 9:
                PhysDerivPyrImpl<9,9,8>(in, out_d0, out_d1, out_d2); break;
            case 10:
                PhysDerivPyrImpl<10,10,9>(in, out_d0, out_d1, out_d2); break;
            case 11:
                PhysDerivPyrImpl<11,11,10>(in, out_d0, out_d1, out_d2); break;
            case 12:
                PhysDerivPyrImpl<12,12,11>(in, out_d0, out_d1, out_d2); break;
            case 13:
                PhysDerivPyrImpl<13,13,12>(in, out_d0, out_d1, out_d2); break;
            case 14:
                PhysDerivPyrImpl<14,14,13>(in, out_d0, out_d1, out_d2); break;
            case 15:
                PhysDerivPyrImpl<15,15,14>(in, out_d0, out_d1, out_d2); break;
            case 16:
                PhysDerivPyrImpl<16,16,15>(in, out_d0, out_d1, out_d2); break;
            default: NEKERROR(ErrorUtil::efatal,
                "PhysDerivPyr: # of modes / points combo not implemented.");

        }
    }

    template<int NQ0, int NQ1, int NQ2>
    void PhysDerivPyrImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &out_d0,
              Array<OneD,       NekDouble> &out_d1,
              Array<OneD,       NekDouble> &out_d2)
    {
        const auto* inptr = input.data();
        auto* outptr_d0 = out_d0.data();
        auto* outptr_d1 = out_d1.data();
        auto* outptr_d2 = out_d2.data();

        constexpr auto ndf = 9;
        constexpr auto nqTot = NQ0 * NQ1 * NQ2;
        constexpr auto nqBlocks = nqTot * vec_t::width;

        // Get size of derivative factor block
        auto dfSize = ndf;
        if (DEFORMED)
        {
            dfSize *= nqTot;
        }

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpOut_d0(nqTot),
            tmpOut_d1(nqTot), tmpOut_d2(nqTot);
        const vec_t* df_ptr;
        for (int e = 0; e < this->m_nBlocks; ++e)
        {

            df_ptr = &(this->m_df[e*dfSize]);

            // Load and transpose data
            load_interleave(inptr, nqTot, tmpIn);

            PhysDerivPyrKernel<NQ0, NQ1, NQ2, DEFORMED>(
                tmpIn,
                this->m_Z[0], this->m_Z[1], this->m_Z[2],
                this->m_D[0], this->m_D[1], this->m_D[2],
                df_ptr,
                tmpOut_d0, tmpOut_d1, tmpOut_d2);

            // de-interleave and store data
            deinterleave_store(tmpOut_d0, nqTot, outptr_d0);
            deinterleave_store(tmpOut_d1, nqTot, outptr_d1);
            deinterleave_store(tmpOut_d2, nqTot, outptr_d2);

            inptr += nqBlocks;
            outptr_d0 += nqBlocks;
            outptr_d1 += nqBlocks;
            outptr_d2 += nqBlocks;
        }
    }
private:
    int m_nmTot;

};

template<bool DEFORMED = false>
struct PhysDerivTet : public PhysDeriv, public Helper<3, DEFORMED>
{
    PhysDerivTet(std::vector<LibUtilities::BasisSharedPtr> basis,
                    int nElmt)
        : PhysDeriv(basis, nElmt),
          Helper<3, DEFORMED>(basis, nElmt),
          m_nmTot(LibUtilities::StdTetData::getNumberOfCoefficients(
                      this->m_nm[0], this->m_nm[1], this->m_nm[2]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<PhysDerivTet<DEFORMED>>(basis, nElmt);
    }

    static NekDouble FlopsPerElement(
        const int nq0,
        const int nq1,
        const int nq2)
    {
        //PhysDerivTensor
        int pdt0 = nq0 * nq1 * nq2 * nq0 * 2;
        int pdt1 = nq2 * nq0 * nq1 * nq1 * 2;
        int pdt2 = nq0 * nq1 * nq2 * nq2 * 2;
        int physDerivTensor = pdt0 + pdt1 + pdt2;

        //PhysDeriv
        int d0 = nq2*(2 + nq1*(3 + nq0));
        int d1 = nq2*(2 + nq1*nq0*5);
        int d2 = nq2*nq1*(2 + nq0*3);
        int df = nq2*nq1*nq0*15;
        int physDeriv = d0 + d1 + d2 + df;

        return (physDerivTensor + physDeriv);
    }

    double GFlops() final
    {
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        int flops = m_nElmt * PhysDerivTet::FlopsPerElement(nq0, nq1, nq2);
        return flops * 1e-9;
    }

    NekDouble NLoads() final
    {
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        //physDerivTensor
        int load_d0 = nq0 * nq1 * nq2 * nq0 * 2;
        int load_d1 = nq2 * nq0 * nq1 * nq1 * 2;
        int load_d2 = nq0 * nq1 * nq2 * nq2 * 2;;
        int physDerivTensor = load_d0 + load_d1 + load_d2;

        int load_dd0 = nq2 * (1 + nq1 * (1+ nq0));
        int load_dd1 = nq2 * (1 + nq1 * nq0 * 3);
        int load_dd2 = nq2 * (nq1 * (1 + nq0*3));
        int load_coef = nq2 * nq1 * nq0 * 3;
        int physDeriv = load_dd0 + load_dd1 + load_dd2 + load_coef;

        return m_nElmt * (physDeriv + physDerivTensor);

    }

    NekDouble NStores() final
    {
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        int store_d0 = nq0 * nq1 * nq2;
        int store_d1 = nq2 * nq0 * nq1;
        int store_d2 = nq0 * nq1 * nq2;
        int physDerivTensor = store_d0 + store_d1 + store_d2;

        int store_dd = nq2 * nq1 * nq0 * 2;
        int store_coef = nq2 * nq1 *nq0 * 3;
        int physDeriv = store_dd * 3 + store_coef;

        return m_nElmt * (physDeriv + physDerivTensor);
    }

    NekDouble Ndof() final
    {
        return m_nmTot * this->m_nElmt;
    }

    void operator()(const Array<OneD, const NekDouble> &in,
                                Array<OneD,       NekDouble> &out_d0,
                                Array<OneD,       NekDouble> &out_d1) final
    {
        // Only for 2D, but need to implement since its abstract
        boost::ignore_unused(in, out_d0, out_d1);
        NEKERROR(ErrorUtil::efatal,
            "Something went horribly wrong... calling 2D op for 3D element");
    }

    void operator()(const Array<OneD, const NekDouble> &in,
                                  Array<OneD,       NekDouble> &out_d0,
                                  Array<OneD,       NekDouble> &out_d1,
                                  Array<OneD,       NekDouble> &out_d2) final
    {
        // Check preconditions
        ASSERTL0(m_basis[0]->GetNumPoints() == (m_basis[1]->GetNumPoints()+1) &&
            m_basis[0]->GetNumPoints() == (m_basis[2]->GetNumPoints()+1),
            "MatrixFree op requires homogenous points");

        switch(m_basis[0]->GetNumPoints())
        {
            case 3:
                PhysDerivTetImpl<3,2,2>(in, out_d0, out_d1, out_d2); break;
            case 4:
                PhysDerivTetImpl<4,3,3>(in, out_d0, out_d1, out_d2); break;
            case 5:
                PhysDerivTetImpl<5,4,4>(in, out_d0, out_d1, out_d2); break;
            case 6:
                PhysDerivTetImpl<6,5,5>(in, out_d0, out_d1, out_d2); break;
            case 7:
                PhysDerivTetImpl<7,6,6>(in, out_d0, out_d1, out_d2); break;
            case 8:
                PhysDerivTetImpl<8,7,7>(in, out_d0, out_d1, out_d2); break;
            case 9:
                PhysDerivTetImpl<9,8,8>(in, out_d0, out_d1, out_d2); break;
            case 10:
                PhysDerivTetImpl<10,9,9>(in, out_d0, out_d1, out_d2); break;
            case 11:
                PhysDerivTetImpl<11,10,10>(in, out_d0, out_d1, out_d2); break;
            case 12:
                PhysDerivTetImpl<12,11,11>(in, out_d0, out_d1, out_d2); break;
            case 13:
                PhysDerivTetImpl<13,12,12>(in, out_d0, out_d1, out_d2); break;
            case 14:
                PhysDerivTetImpl<14,13,13>(in, out_d0, out_d1, out_d2); break;
            case 15:
                PhysDerivTetImpl<15,14,14>(in, out_d0, out_d1, out_d2); break;
            case 16:
                PhysDerivTetImpl<16,15,15>(in, out_d0, out_d1, out_d2); break;
            default: NEKERROR(ErrorUtil::efatal,
                "PhysDerivTet: # of points combo not implemented.");
        }
    }

    template<int NQ0, int NQ1, int NQ2>
    void PhysDerivTetImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &out_d0,
              Array<OneD,       NekDouble> &out_d1,
              Array<OneD,       NekDouble> &out_d2)
    {
        const auto* inptr = input.data();
        auto* outptr_d0 = out_d0.data();
        auto* outptr_d1 = out_d1.data();
        auto* outptr_d2 = out_d2.data();

        constexpr auto ndf = 9;
        constexpr auto nqTot = NQ0 * NQ1 * NQ2;
        constexpr auto nqBlocks = nqTot * vec_t::width;


        // Get size of derivative factor block
        auto dfSize = ndf;
        if (DEFORMED)
        {
            dfSize *= nqTot;
        }

        std::vector<vec_t, allocator<vec_t>> diff0(nqTot), diff1(nqTot),
            diff2(nqTot);

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpOut_d0(nqTot),
            tmpOut_d1(nqTot), tmpOut_d2(nqTot);
        const vec_t* df_ptr;

        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            df_ptr = &(this->m_df[e*dfSize]);

            // Load and transpose data
            load_interleave(inptr, nqTot, tmpIn);

            PhysDerivTetKernel<NQ0, NQ1, NQ2, DEFORMED>(
                tmpIn,
                this->m_Z[0], this->m_Z[1], this->m_Z[2],
                this->m_D[0], this->m_D[1], this->m_D[2],
                df_ptr,
                diff0, diff1, diff2,
                tmpOut_d0, tmpOut_d1, tmpOut_d2);

            // de-interleave and store data
            deinterleave_store(tmpOut_d0, nqTot, outptr_d0);
            deinterleave_store(tmpOut_d1, nqTot, outptr_d1);
            deinterleave_store(tmpOut_d2, nqTot, outptr_d2);

            inptr += nqBlocks;
            outptr_d0 += nqBlocks;
            outptr_d1 += nqBlocks;
            outptr_d2 += nqBlocks;
        }
    }
private:
    int m_nmTot;
};

} // namespace MatrixFree
} // namespace Nektar

#endif
