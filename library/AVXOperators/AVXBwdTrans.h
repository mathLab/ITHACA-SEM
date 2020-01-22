#ifndef AVXBWDTRANS_H
#define AVXBWDTRANS_H

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ShapeType.hpp>
#include <LibUtilities/Foundations/Basis.h>

#include "VecData.hpp"
#include "Operator.hpp"
#include "AVXBwdTransKernels.hpp"

namespace Nektar {
namespace AVX {

template<int VW>
struct AVXBwdTransQuad : public BwdTrans, public AVXHelper<VW, 2>
{
    AVXBwdTransQuad(std::vector<LibUtilities::BasisSharedPtr> basis,
                    int nElmt)
        : BwdTrans(basis, nElmt),
          AVXHelper<VW, 2>(basis, nElmt),
          m_nmTot(LibUtilities::StdQuadData::getNumberOfCoefficients(
                      this->m_nm[0], this->m_nm[1]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<AVXBwdTransQuad<VW>>(basis, nElmt);
    }

    static NekDouble FlopsPerElement(
        const int nm,
        const int nq0,
        const int nq1)
    {
        int loop_i = nq0*nm*nm*2;
        int loop_j = nq1*nq0*nm*2;
        return loop_i + loop_j;
    }

    virtual NekDouble GFlops()
    {

        const int nm = m_basis[0]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();

        int flops = m_nElmt * AVXBwdTransQuad::FlopsPerElement(nm, nq0, nq1);
        return flops * 1e-9;
    }

    virtual NekDouble NLoads()
    {
        const int nm0 = m_basis[0]->GetNumModes();
        const int nm1 = m_basis[1]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();

        int load_iqp = nq0 * nm1 * nm0 * 2;
        int load_jiq = nq1 * nq0 * nm1 * 2;
        int load_expected = load_iqp + load_jiq;

        return m_nElmt * load_expected;

    }

    virtual NekDouble NStores()
    {
        const int nm0 = m_basis[0]->GetNumModes();
        const int nm1 = m_basis[1]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();

        int store_iq = nq0 * nm1;
        int store_ij = nq0 * nq1;
        int store_expected = store_iq + store_ij;

        return m_nElmt * store_expected;
    }

    virtual NekDouble Ndof()
    {
        return m_nmTot * this->m_nElmt;
    }

    virtual void operator()(const Array<OneD, const NekDouble> &in,
                                  Array<OneD,       NekDouble> &out)
    {
        switch(m_basis[0]->GetNumModes())
        {
            case 2:  AVXBwdTransQuadImpl<2 ,2 ,3 ,3 >(in, out); break;
            case 3:  AVXBwdTransQuadImpl<3 ,3 ,4 ,4 >(in, out); break;
            case 4:  AVXBwdTransQuadImpl<4 ,4 ,5 ,5 >(in, out); break;
            case 5:  AVXBwdTransQuadImpl<5 ,5 ,6 ,6 >(in, out); break;
            case 6:  AVXBwdTransQuadImpl<6 ,6 ,7 ,7 >(in, out); break;
            case 7:  AVXBwdTransQuadImpl<7 ,7 ,8 ,8 >(in, out); break;
            case 8:  AVXBwdTransQuadImpl<8 ,8 ,9 ,9 >(in, out); break;
            case 9:  AVXBwdTransQuadImpl<9 ,9 ,10,10>(in, out); break;
            case 10: AVXBwdTransQuadImpl<10,10,11,11>(in, out); break;
            case 11: AVXBwdTransQuadImpl<11,11,12,12>(in, out); break;
        }
    }

    template<int NM0, int NM1, int NQ0, int NQ1>
    void AVXBwdTransQuadImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &output)
    {
        using T = VecData<double, VW>;
        auto *inptr = &input[0];
        auto *outptr = &output[0];

        constexpr int nmTot = NM0 * NM1;
        constexpr int nqTot = NQ0 * NQ1;
        constexpr int nqBlocks = NQ0 * NQ1 * VW;
        const int nmBlocks = m_nmTot * VW;

        T p_sums[NM0 * NQ0]; //Sums over q for each quadpt p_i
        AlignedVector<T> tmpIn(nmTot), tmpOut(nqTot);

        for (int e = 0; e < this->m_nBlocks; e++)
        {
            // Load and transpose data
            T::load_interleave(inptr, nmTot, tmpIn);

            AVXBwdTransQuadKernel<NM0, NM1, NQ0, NQ1, VW>(
                tmpIn, this->m_bdata[0], this->m_bdata[1],
                p_sums, tmpOut);

            // de-interleave and store data
            T::deinterleave_store(tmpOut, nqTot, outptr);

            inptr += nmBlocks;
            outptr += nqBlocks;
        }
    }

private:
    int m_nmTot;
};

template<int VW>
struct AVXBwdTransTri : public BwdTrans, public AVXHelper<VW, 2>
{
    AVXBwdTransTri(std::vector<LibUtilities::BasisSharedPtr> basis,
                   int nElmt)
        : BwdTrans(basis, nElmt),
          AVXHelper<VW, 2>(basis, nElmt),
          m_nmTot(LibUtilities::StdTriData::getNumberOfCoefficients(
                      this->m_nm[0], this->m_nm[1]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<AVXBwdTransTri<VW>>(basis, nElmt);
    }

    static NekDouble FlopsPerElement(
        const int nm,
        const int nq0,
        const int nq1)
    {
        const int pq_loop = nm * (nm+1) / 2;
        const int inner = 2 * nm + 3;

        return nq1 * (2*pq_loop + nq0*inner);
    }

    virtual NekDouble GFlops()
    {
        const int nm = m_basis[0]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();

        int flops = m_nElmt * AVXBwdTransTri::FlopsPerElement(nm, nq0, nq1);

        return flops * 1e-9;
    }

    virtual NekDouble NLoads()
    {
        const int nm0 = m_basis[0]->GetNumModes();
        const int nm1 = m_basis[1]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();

        int load_jpq = nq1 * (nm0 * (nm0 + 1) / 2) * 2;
        int load_jip =  nq1 * nq0 * nm0 * 2;
        int load_corr = nq1 * nq0 * 3;
        int load_expected = load_jpq + load_jip + load_corr;

        return m_nElmt * load_expected;

    }

    virtual NekDouble NStores()
    {
        const int nm0 = m_basis[0]->GetNumModes();
        const int nm1 = m_basis[1]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();

        int store_jp = nq1 * nm0;
        int store_ji = nq1 * nq0;
        int store_expected = store_jp + store_ji;

        return m_nElmt * store_expected;
    }

    virtual NekDouble Ndof()
    {
        return m_nmTot * this->m_nElmt;
    }

    virtual void operator()(const Array<OneD, const NekDouble> &in,
                                  Array<OneD,       NekDouble> &out)
    {
        if (m_basis[0]->GetBasisType() == LibUtilities::eModified_A)
        {
            switch(m_basis[0]->GetNumModes())
            {
                case 2:  AVXBwdTransTriImpl<2 ,2 ,3 ,2 ,true>(in, out); break;
                case 3:  AVXBwdTransTriImpl<3 ,3 ,4 ,3 ,true>(in, out); break;
                case 4:  AVXBwdTransTriImpl<4 ,4 ,5 ,4 ,true>(in, out); break;
                case 5:  AVXBwdTransTriImpl<5 ,5 ,6 ,5 ,true>(in, out); break;
                case 6:  AVXBwdTransTriImpl<6 ,6 ,7 ,6 ,true>(in, out); break;
                case 7:  AVXBwdTransTriImpl<7 ,7 ,8 ,7 ,true>(in, out); break;
                case 8:  AVXBwdTransTriImpl<8 ,8 ,9 ,8 ,true>(in, out); break;
                case 9:  AVXBwdTransTriImpl<9 ,9 ,10,9 ,true>(in, out); break;
                case 10: AVXBwdTransTriImpl<10,10,11,10,true>(in, out); break;
                case 11: AVXBwdTransTriImpl<11,11,12,11,true>(in, out); break;
            }
        }
        else
        {
            switch(m_basis[0]->GetNumModes())
            {
                case 2:  AVXBwdTransTriImpl<2 ,2 ,3 ,2 ,false>(in, out); break;
                case 3:  AVXBwdTransTriImpl<3 ,3 ,4 ,3 ,false>(in, out); break;
                case 4:  AVXBwdTransTriImpl<4 ,4 ,5 ,4 ,false>(in, out); break;
                case 5:  AVXBwdTransTriImpl<5 ,5 ,6 ,5 ,false>(in, out); break;
                case 6:  AVXBwdTransTriImpl<6 ,6 ,7 ,6 ,false>(in, out); break;
                case 7:  AVXBwdTransTriImpl<7 ,7 ,8 ,7 ,false>(in, out); break;
                case 8:  AVXBwdTransTriImpl<8 ,8 ,9 ,8 ,false>(in, out); break;
                case 9:  AVXBwdTransTriImpl<9 ,9 ,10,9 ,false>(in, out); break;
                case 10: AVXBwdTransTriImpl<10,10,11,10,false>(in, out); break;
                case 11: AVXBwdTransTriImpl<11,11,12,11,false>(in, out); break;
            }
        }
    }

    template<int NM0, int NM1, int NQ0, int NQ1, bool CORRECT>
    void AVXBwdTransTriImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &output)
    {
        using T = VecData<double, VW>;
        auto *inptr = &input[0];
        auto *outptr = &output[0];

        constexpr int nqBlocks = NQ0 * NQ1 * VW;
        const int nmBlocks = m_nmTot * VW;

        T q_sums[NM0]; //Sums over q for each p

        for(int e = 0; e < this->m_nBlocks; e++)
        {
            AVXBwdTransTriKernel<NM0, NM1, NQ0, NQ1, VW, CORRECT>(
                inptr,
                this->m_bdata[0], this->m_bdata[1],
                q_sums,
                outptr);

            inptr += nmBlocks;
            outptr += nqBlocks;
        }
    }

private:
    int m_nmTot;
};

template<int VW>
struct AVXBwdTransHex : public BwdTrans, public AVXHelper<VW, 3>
{
    AVXBwdTransHex(std::vector<LibUtilities::BasisSharedPtr> basis,
                   int nElmt )
    : BwdTrans(basis, nElmt),
        AVXHelper<VW, 3>(basis, nElmt),
        m_nmTot(LibUtilities::StdHexData::getNumberOfCoefficients(
                    this->m_nm[0], this->m_nm[1], this->m_nm[2]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<AVXBwdTransHex<VW>>(basis, nElmt);
    }

    static NekDouble FlopsPerElement(
        const int nm,
        const int nq0,
        const int nq1,
        const int nq2)
    {
        int loop_i = nq0 * nm *  nm  * nm  * 2;
        int loop_j = nq1 * nq0 * nm  * nm  * 2;
        int loop_k = nq2 * nq1 * nq0 * nm  * 2;

        return (loop_i + loop_j + loop_k);
    }

    virtual NekDouble GFlops() override
    {
        const int nm = m_basis[0]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        int flops = this->m_nElmt * AVXBwdTransHex::FlopsPerElement(nm,nq0,nq1,nq2);
        return flops * 1e-9;
    }

    virtual NekDouble Ndof()
    {
        return m_nmTot * this->m_nElmt;
    }

    virtual NekDouble NLoads() override
    {
        const int nm0 = m_basis[0]->GetNumModes();
        const int nm1 = m_basis[1]->GetNumModes();
        const int nm2 = m_basis[2]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        int load_irqp = nq0 * nm2 * nm1 * nm0 * 2;
        int load_jirq = nq1 * nq0 * nm2 * nm1 * 2;
        int load_kjir = nq2 * nq1 * nq0 * nm2 * 2;
        int load_expected = load_irqp + load_jirq + load_kjir;

        return load_expected * m_nElmt;
    }

    virtual NekDouble NStores() override
    {
        const int nm0 = m_basis[0]->GetNumModes();
        const int nm1 = m_basis[1]->GetNumModes();
        const int nm2 = m_basis[2]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        int store_irq = nq0 * nm2 * nm1;
        int store_jir = nq1 * nq0 * nm2;
        int store_kji = nq0 * nq1 * nq2;
        int store_expected = store_irq + store_jir + store_kji;

        return store_expected * m_nElmt;
    }

    virtual void operator()(const Array<OneD, const NekDouble> &in,
                                      Array<OneD,       NekDouble> &out)
    {
        switch(m_basis[0]->GetNumModes())
        {
            case 2:  AVXBwdTransHexImpl<2 ,2 ,2 ,3 ,3 ,3 >(in, out); break;
            case 3:  AVXBwdTransHexImpl<3 ,3 ,3 ,4 ,4 ,4 >(in, out); break;
            case 4:  AVXBwdTransHexImpl<4 ,4 ,4 ,5 ,5 ,5 >(in, out); break;
            case 5:  AVXBwdTransHexImpl<5 ,5 ,5 ,6 ,6 ,6 >(in, out); break;
            case 6:  AVXBwdTransHexImpl<6 ,6 ,6 ,7 ,7 ,7 >(in, out); break;
            case 7:  AVXBwdTransHexImpl<7 ,7 ,7 ,8 ,8 ,8 >(in, out); break;
            case 8:  AVXBwdTransHexImpl<8 ,8 ,8 ,9 ,9 ,9 >(in, out); break;
            case 9:  AVXBwdTransHexImpl<9 ,9 ,9 ,10,10,10>(in, out); break;
            case 10: AVXBwdTransHexImpl<10,10,10,11,11,11>(in, out); break;
            case 11: AVXBwdTransHexImpl<11,11,11,12,12,12>(in, out); break;
        }
    }

    template<int NM0, int NM1, int NM2, int NQ0, int NQ1, int NQ2>
    void AVXBwdTransHexImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &output)
    {
        using T = VecData<double, VW>;
        auto *inptr = &input[0];
        auto *outptr = &output[0];

        constexpr int nqBlocks = NQ0 * NQ1 * NQ2 * VW;
        const int nmBlocks = m_nmTot * VW;

        //T sum_irq[NM2 * NQ0 * NM1], sum_jir[NM2 * NQ1 * NQ0];
        T sum_irq[NQ0 * NQ1 * NQ2], sum_jir[NQ0 * NQ1 * NQ2];

        for(int e = 0; e < this->m_nBlocks; e++){
            AVXBwdTransHexKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, VW>(
                inptr,
                this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
                sum_irq, sum_jir,
                outptr);

            inptr += nmBlocks;
            outptr += nqBlocks;
        }
    }

private:
    int m_nmTot;
};

template<int VW>
struct AVXBwdTransTet : public BwdTrans, public AVXHelper<VW, 3>
{
    AVXBwdTransTet(std::vector<LibUtilities::BasisSharedPtr> basis,
                   int nElmt)
    : BwdTrans(basis, nElmt),
        AVXHelper<VW, 3>(basis, nElmt),
        m_nmTot(LibUtilities::StdTetData::getNumberOfCoefficients(
                    this->m_nm[0], this->m_nm[1], this->m_nm[2]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<AVXBwdTransTet<VW>>(basis, nElmt);
    }

    static NekDouble FlopsPerElement(
        const int nm,
        const int nq0,
        const int nq1,
        const int nq2)
    {
        int pqr_loop = 2.0 * nm*(nm+1)*(nm+2) / 6;
        int pq_loop = 2.0 * nm*(nm+1) / 2;
        int i_loop = 2.0 * nm + 8 + 4*(nm - 1);

        return ( nq2 * (pqr_loop + nq1*(pq_loop + nq0*(i_loop) )) );
    }

    virtual NekDouble GFlops() override
    {
        const int nm = m_basis[0]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        int flop_estimate = m_nElmt * AVXBwdTransTet::FlopsPerElement(nm, nq0, nq1, nq2);
        return flop_estimate * 1e-9; //Convert to gigaflops.
    }

    virtual NekDouble Ndof() override
    {
        return m_nmTot * this->m_nElmt;
    }

    virtual NekDouble NLoads() override
    {
        const int nm0 = m_basis[0]->GetNumModes();
        const int nm1 = m_basis[1]->GetNumModes();
        const int nm2 = m_basis[2]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        int load_kpqr = nq2 * (nm0*(nm0+1)/2) * nm1 * 2;
        int load_kjpq = nq2 * nq1 * nm0 * nm1 * 2;
        int load_kjip = nq2 * nq1 * nq0 * nm0 * 2;
        int corr = nq2 * nq1 * nq0 * (2 + nm1*2);
        int load_expected = load_kpqr + load_kjpq + load_kjip + corr;

        return load_expected * m_nElmt;
    }

    virtual NekDouble NStores() override
    {
        const int nm0 = m_basis[0]->GetNumModes();
        const int nm1 = m_basis[1]->GetNumModes();
        const int nm2 = m_basis[2]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        int store_kpq = nq2 * nm0 * nm1;
        int store_kjp = nq2 * nq1 * nm0;
        int store_kji = nq2 * nq1 * nq0;
        int store_expected = store_kpq + store_kjp + store_kji;

        return store_expected * m_nElmt;
    }

    virtual void operator()(const Array<OneD, const NekDouble> &in,
                                    Array<OneD,       NekDouble> &out)
    {
        if (m_basis[0]->GetBasisType() == LibUtilities::eModified_A)
        {
            switch(m_basis[0]->GetNumModes())
            {
                case 2:  AVXBwdTransTetImpl<2 ,2 ,2 ,3 ,2 ,2 ,true>(in, out); break;
                case 3:  AVXBwdTransTetImpl<3 ,3 ,3 ,4 ,3 ,3 ,true>(in, out); break;
                case 4:  AVXBwdTransTetImpl<4 ,4 ,4 ,5 ,4 ,4 ,true>(in, out); break;
                case 5:  AVXBwdTransTetImpl<5 ,5 ,5 ,6 ,5 ,5 ,true>(in, out); break;
                case 6:  AVXBwdTransTetImpl<6 ,6 ,6 ,7 ,6 ,6 ,true>(in, out); break;
                case 7:  AVXBwdTransTetImpl<7 ,7 ,7 ,8 ,7 ,7 ,true>(in, out); break;
                case 8:  AVXBwdTransTetImpl<8 ,8 ,8 ,9 ,8 ,8 ,true>(in, out); break;
                case 9:  AVXBwdTransTetImpl<9 ,9 ,9 ,10,9 ,9 ,true>(in, out); break;
                case 10: AVXBwdTransTetImpl<10,10,10,11,10,10,true>(in, out); break;
                case 11: AVXBwdTransTetImpl<11,11,11,12,11,11,true>(in, out); break;
            }
        }
        else
        {
            switch(m_basis[0]->GetNumModes())
            {
                case 2:  AVXBwdTransTetImpl<2 ,2 ,2 ,3 ,2 ,2 ,false>(in, out); break;
                case 3:  AVXBwdTransTetImpl<3 ,3 ,3 ,4 ,3 ,3 ,false>(in, out); break;
                case 4:  AVXBwdTransTetImpl<4 ,4 ,4 ,5 ,4 ,4 ,false>(in, out); break;
                case 5:  AVXBwdTransTetImpl<5 ,5 ,5 ,6 ,5 ,5 ,false>(in, out); break;
                case 6:  AVXBwdTransTetImpl<6 ,6 ,6 ,7 ,6 ,6 ,false>(in, out); break;
                case 7:  AVXBwdTransTetImpl<7 ,7 ,7 ,8 ,7 ,7 ,false>(in, out); break;
                case 8:  AVXBwdTransTetImpl<8 ,8 ,8 ,9 ,8 ,8 ,false>(in, out); break;
                case 9:  AVXBwdTransTetImpl<9 ,9 ,9 ,10,9 ,9 ,false>(in, out); break;
                case 10: AVXBwdTransTetImpl<10,10,10,11,10,10,false>(in, out); break;
                case 11: AVXBwdTransTetImpl<11,11,11,12,11,11,false>(in, out); break;
            }
        }
    }

    template<int NM0, int NM1, int NM2, int NQ0, int NQ1, int NQ2, bool CORRECT>
    void AVXBwdTransTetImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &output)
    {
        using T = VecData<double, VW>;
        auto *inptr = &input[0];
        auto *outptr = &output[0];

        constexpr int nqBlocks = NQ0 * NQ1 * NQ2 * VW;
        const int nmBlocks = m_nmTot * VW;

        T fpq[NM0 * NM1], fp[NM0];

        for(int e = 0; e < this->m_nBlocks; e++){
            AVXBwdTransTetKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, VW, CORRECT>(
                inptr,
                this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
                fpq, fp,
                outptr);

            inptr += nmBlocks;
            outptr += nqBlocks;
        }
    }
private:
    int m_nmTot;
};

template<int VW>
struct AVXBwdTransPrism : public BwdTrans, public AVXHelper<VW, 3>
{
    AVXBwdTransPrism(std::vector<LibUtilities::BasisSharedPtr> basis,
                     int nElmt)
    : BwdTrans(basis, nElmt),
        AVXHelper<VW, 3>(basis, nElmt),
        m_nmTot(LibUtilities::StdPrismData::getNumberOfCoefficients(
                    this->m_nm[0], this->m_nm[1], this->m_nm[2]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<AVXBwdTransPrism<VW>>(basis, nElmt);
    }

    static NekDouble FlopsPerElement(
        const int nm,
        const int nq0,
        const int nq1,
        const int nq2)
    {
        int loop_pqr = 2 * nm * nm * (nm + 1) / 2.0;
        int loop_pq = nm * nm * 2.0;
        int loop_p = nm * 2.0;
        int q_correction = nm * 6.0;

        return nq2 * (loop_pqr + nq1 *( loop_pq + nq0*(loop_p + q_correction)));
    }

    virtual double GFlops() override
    {
        const int nm = m_basis[0]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        int flops = m_nElmt * AVXBwdTransPrism::FlopsPerElement(nm, nq0, nq1, nq2);

        return flops * 1e-9;

    }

    virtual NekDouble Ndof() override
    {
        return m_nmTot * this->m_nElmt;
    }

    virtual NekDouble NLoads() override
    {
        const int nm0 = m_basis[0]->GetNumModes();
        const int nm1 = m_basis[1]->GetNumModes();
        const int nm2 = m_basis[2]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        int load_kpqr = nq2 * (nm0*(nm0+1)*(nm0+2) / 6) * 2;
        int load_kjpq = nq2 * nq1 * (nm0*(nm0+1) / 2) * 2;
        int load_kjip = nq2 * nq1 * nq0 * nm0 * 2;
        int corr = nq2 * nq1 * nq0 * (8 + 4 *(nm2 - 1));
        int load_expected = load_kpqr + load_kjpq + load_kjip + corr;

        return load_expected * m_nElmt;
    }

    virtual NekDouble NStores() override
    {
        const int nm0 = m_basis[0]->GetNumModes();
        const int nm1 = m_basis[1]->GetNumModes();
        const int nm2 = m_basis[2]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        int store_kpr = nq2 * (nm0*(nm0+1)/2);
        int store_kjp = nq2 * nq1 * nm0;
        int store_kji = nq2 * nq1 * nq0;
        int store_expected = store_kpr + store_kjp + store_kji;

        return store_expected * m_nElmt;
    }

    virtual void operator()(const Array<OneD, const NekDouble> &in,
                                    Array<OneD,       NekDouble> &out)
    {
        if (m_basis[0]->GetBasisType() == LibUtilities::eModified_A)
        {
            switch(m_basis[0]->GetNumModes())
            {
                case 2:  AVXBwdTransPrismImpl<2 ,2 ,2 ,3 ,3 ,2 ,true>(in, out); break;
                case 3:  AVXBwdTransPrismImpl<3 ,3 ,3 ,4 ,4 ,3 ,true>(in, out); break;
                case 4:  AVXBwdTransPrismImpl<4 ,4 ,4 ,5 ,5 ,4 ,true>(in, out); break;
                case 5:  AVXBwdTransPrismImpl<5 ,5 ,5 ,6 ,6 ,5 ,true>(in, out); break;
                case 6:  AVXBwdTransPrismImpl<6 ,6 ,6 ,7 ,7 ,6 ,true>(in, out); break;
                case 7:  AVXBwdTransPrismImpl<7 ,7 ,7 ,8 ,8 ,7 ,true>(in, out); break;
                case 8:  AVXBwdTransPrismImpl<8 ,8 ,8 ,9 ,9 ,8 ,true>(in, out); break;
                case 9:  AVXBwdTransPrismImpl<9 ,9 ,9 ,10,10,9 ,true>(in, out); break;
                case 10: AVXBwdTransPrismImpl<10,10,10,11,11,10,true>(in, out); break;
                case 11: AVXBwdTransPrismImpl<11,11,11,12,12,11,true>(in, out); break;
            }
        }
        else
        {
            switch(m_basis[0]->GetNumModes())
            {
                case 2:  AVXBwdTransPrismImpl<2 ,2 ,2 ,3 ,3 ,2 ,false>(in, out); break;
                case 3:  AVXBwdTransPrismImpl<3 ,3 ,3 ,4 ,4 ,3 ,false>(in, out); break;
                case 4:  AVXBwdTransPrismImpl<4 ,4 ,4 ,5 ,5 ,4 ,false>(in, out); break;
                case 5:  AVXBwdTransPrismImpl<5 ,5 ,5 ,6 ,6 ,5 ,false>(in, out); break;
                case 6:  AVXBwdTransPrismImpl<6 ,6 ,6 ,7 ,7 ,6 ,false>(in, out); break;
                case 7:  AVXBwdTransPrismImpl<7 ,7 ,7 ,8 ,8 ,7 ,false>(in, out); break;
                case 8:  AVXBwdTransPrismImpl<8 ,8 ,8 ,9 ,9 ,8 ,false>(in, out); break;
                case 9:  AVXBwdTransPrismImpl<9 ,9 ,9 ,10,10,9 ,false>(in, out); break;
                case 10: AVXBwdTransPrismImpl<10,10,10,11,11,10,false>(in, out); break;
                case 11: AVXBwdTransPrismImpl<11,11,11,12,12,11,false>(in, out); break;
            }
        }
    }

    template<int NM0, int NM1, int NM2, int NQ0, int NQ1, int NQ2, bool CORRECT>
    void AVXBwdTransPrismImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &output)
    {
        using T = VecData<double, VW>;
        auto *inptr = &input[0];
        auto *outptr = &output[0];

        constexpr int nqBlocks = NQ0 * NQ1 * NQ2 * VW;
        const int nmBlocks = m_nmTot * VW;

        T fpq[NM0 * NM1], fp[NM0];

        for(int e = 0; e < this->m_nBlocks; e++){
            AVXBwdTransPrismKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, VW, CORRECT>(
                inptr,
                this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
                fpq, fp,
                outptr);

            inptr += nmBlocks;
            outptr += nqBlocks;
        }
    }


private:
    int m_nmTot;
};

} // namespace AVX
} // namespace Nektar

#endif
