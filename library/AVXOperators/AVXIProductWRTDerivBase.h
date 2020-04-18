#ifndef NEKTAR_LIBRARY_AVXIPRODUCTWRTDERIVBASE_H
#define NEKTAR_LIBRARY_AVXIPRODUCTWRTDERIVBASE_H

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ShapeType.hpp>
#include <LibUtilities/Foundations/Basis.h>

#include "VecData.hpp"
#include "Operator.hpp"

// #include "AVXIProduct.h"
#include "AVXIProductKernels.hpp"

namespace Nektar
{
namespace AVX
{

template<int VW, bool DEFORMED = false>
struct AVXIProductWRTDerivBaseQuad : public IProductWRTDerivBase, public AVXHelper<VW, 2, DEFORMED>
{
    AVXIProductWRTDerivBaseQuad(std::vector<LibUtilities::BasisSharedPtr> basis,
                   int nElmt)
        : IProductWRTDerivBase(basis, nElmt),
          AVXHelper<VW, 2, DEFORMED>(basis, nElmt),
          m_nmTot(LibUtilities::StdQuadData::getNumberOfCoefficients(
                      this->m_nm[0], this->m_nm[1]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<AVXIProductWRTDerivBaseQuad<VW,DEFORMED>>(basis, nElmt);
    }

    virtual void operator()(const Array<OneD, const NekDouble> &in0,
                            const Array<OneD, const NekDouble> &in1,
                                  Array<OneD,       NekDouble> &out)
    {
        switch(m_basis[0]->GetNumModes())
        {
            case 2:  AVXIProductWRTDerivBaseQuadImpl<2 ,2 ,3 ,3 >(in0, in1, out); break;
            case 3:  AVXIProductWRTDerivBaseQuadImpl<3 ,3 ,4 ,4 >(in0, in1, out); break;
            case 4:  AVXIProductWRTDerivBaseQuadImpl<4 ,4 ,5 ,5 >(in0, in1, out); break;
            case 5:  AVXIProductWRTDerivBaseQuadImpl<5 ,5 ,6 ,6 >(in0, in1, out); break;
            case 6:  AVXIProductWRTDerivBaseQuadImpl<6 ,6 ,7 ,7 >(in0, in1, out); break;
            case 7:  AVXIProductWRTDerivBaseQuadImpl<7 ,7 ,8 ,8 >(in0, in1, out); break;
            case 8:  AVXIProductWRTDerivBaseQuadImpl<8 ,8 ,9 ,9 >(in0, in1, out); break;
            case 9:  AVXIProductWRTDerivBaseQuadImpl<9 ,9 ,10,10>(in0, in1, out); break;
            case 10: AVXIProductWRTDerivBaseQuadImpl<10,10,11,11>(in0, in1, out); break;
            case 11: AVXIProductWRTDerivBaseQuadImpl<11,11,12,12>(in0, in1, out); break;
        }
    }

    template<int NM0, int NM1, int NQ0, int NQ1>
    void AVXIProductWRTDerivBaseQuadImpl(
        const Array<OneD, const NekDouble> &input0,
        const Array<OneD, const NekDouble> &input1,
              Array<OneD,       NekDouble> &output)
    {
        using T = VecData<double, VW>;
        auto *inptr0 = &input0[0];
        auto *inptr1 = &input1[0];
        auto *outptr = &output[0];

        constexpr int ndf = 4;
        constexpr int nqTot = NQ0 * NQ1;
        constexpr int nqBlocks = nqTot * VW;
        const int nmBlocks = m_nmTot * VW;

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

        T sums_j[NQ1]; //Sums over eta0 for each value of eta1;
        AlignedVector<T> tmpIn0(nqTot), tmpIn1(nqTot), tmp0(nqTot), tmp1(nqTot),
            tmpOut(m_nmTot), tmpOut2(m_nmTot);

        const T *df_ptr;
        const T *jac_ptr;
        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            // Jacobian
            jac_ptr = &(this->m_jac[dJSize*e]);

            // Derivative factor
            df_ptr = &(this->m_df[e*dfSize]);

            // Load and transpose data
            T::load_interleave(inptr0, nqTot, tmpIn0);
            T::load_interleave(inptr1, nqTot, tmpIn1);

            // Calculate dx/dxi in[0] + dy/dxi in[1]
            T df0, df1, df2, df3;
            if(!DEFORMED)
            {
                df0 = df_ptr[0];
                df1 = df_ptr[1];
                df2 = df_ptr[2];
                df3 = df_ptr[3];
            }
            for (int i = 0; i < nqTot; ++i)
            {
                tmp0[i] = df0 * tmpIn0[i] + df2 * tmpIn1[i];
                tmp1[i] = df1 * tmpIn0[i] + df3 * tmpIn1[i];
            }

            // IP DB0 B1
            AVXIProductQuadKernel<NM0, NM1, NQ0, NQ1, VW, false, false, DEFORMED>(
                tmp0, this->m_dbdata[0], this->m_bdata[1],
                this->m_w[0], this->m_w[1], jac_ptr,
                sums_j, tmpOut);

            // IP DB1 B0
            AVXIProductQuadKernel<NM0, NM1, NQ0, NQ1, VW, false, true, DEFORMED>(
                tmp1, this->m_bdata[0], this->m_dbdata[1],
                this->m_w[0], this->m_w[1], jac_ptr,
                sums_j, tmpOut);

            // // Add
            // for (int i = 0; i < m_nmTot; ++i)
            // {
            //     tmpOut[i] = tmpOut[i] + tmpOut2[i];
            // }

            // de-interleave and store data
            T::deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr0 += nqBlocks;
            inptr1 += nqBlocks;
            outptr += nmBlocks;
        }
    }
public:
    // static NekDouble FlopsPerElement(
    //     const int nm,
    //     const int nq0,
    //     const int nq1)
    // {
    //     int loop_ji = nm * nq0 * nq1 * 4;
    //     int loop_qj = nm*nm*nq1 * 3;
    //     return ( loop_ji + loop_qj);
    // }

    // virtual NekDouble GFlops() override
    // {
    //     const int nm = m_basis[0]->GetNumModes();
    //     const int nq0 = m_basis[0]->GetNumPoints();
    //     const int nq1 = m_basis[1]->GetNumPoints();

    //     int flops = this->m_nElmt * AVXIProductQuad::FlopsPerElement(nm, nq0, nq1);
    //     return 1e-9 * flops;
    // }

    // virtual NekDouble NStores() override
    // {
    //     const int nm = m_basis[0]->GetNumModes();
    //     const int nq0 = m_basis[0]->GetNumPoints();
    //     const int nq1 = m_basis[1]->GetNumPoints();

    //     int store_pj = nm * nq1;
    //     int store_pq = nm * nm;
    //     int store_expected = store_pj + store_pq;

    //     return this->m_nElmt * store_expected;
    // }

    // virtual NekDouble NLoads() override
    // {
    //     const int nm = m_basis[0]->GetNumModes();
    //     const int nq0 = m_basis[0]->GetNumPoints();
    //     const int nq1 = m_basis[1]->GetNumPoints();

    //     int load_pji = nm * nq1 * nq0 * 3;
    //     int load_pqj =  nm * nm * nq1 * 3;
    //     int load_expected = load_pji + load_pqj;

    //     return this->m_nElmt * load_expected;
    // }

    // virtual NekDouble Ndof() override
    // {
    //     return m_nmTot * this->m_nElmt;
    // }

private:
    int m_nmTot;
};

// template<int VW, bool DEFORMED = false>
// struct AVXIProductTri : public IProduct, public AVXHelper<VW, 2, DEFORMED>
// {
//     AVXIProductTri(std::vector<LibUtilities::BasisSharedPtr> basis,
//                    int nElmt)
//         : IProduct(basis, nElmt),
//           AVXHelper<VW, 2, DEFORMED>(basis, nElmt),
//           m_nmTot(LibUtilities::StdTriData::getNumberOfCoefficients(
//                       this->m_nm[0], this->m_nm[1]))
//     {
//     }

//     static std::shared_ptr<Operator> Create(
//         std::vector<LibUtilities::BasisSharedPtr> basis,
//         int nElmt)
//     {
//         return std::make_shared<AVXIProductTri<VW, DEFORMED>>(basis, nElmt);
//     }

//     virtual void operator()(const Array<OneD, const NekDouble> &in,
//                                   Array<OneD,       NekDouble> &out)
//     {
//         if (m_basis[0]->GetBasisType() == LibUtilities::eModified_A)
//         {
//             switch(m_basis[0]->GetNumModes())
//             {
//                 case 2:  AVXIProductTriImpl<2 ,2 ,3 ,2 ,true>(in, out); break;
//                 case 3:  AVXIProductTriImpl<3 ,3 ,4 ,3 ,true>(in, out); break;
//                 case 4:  AVXIProductTriImpl<4 ,4 ,5 ,4 ,true>(in, out); break;
//                 case 5:  AVXIProductTriImpl<5 ,5 ,6 ,5 ,true>(in, out); break;
//                 case 6:  AVXIProductTriImpl<6 ,6 ,7 ,6 ,true>(in, out); break;
//                 case 7:  AVXIProductTriImpl<7 ,7 ,8 ,7 ,true>(in, out); break;
//                 case 8:  AVXIProductTriImpl<8 ,8 ,9 ,8 ,true>(in, out); break;
//                 case 9:  AVXIProductTriImpl<9 ,9 ,10,9 ,true>(in, out); break;
//                 case 10: AVXIProductTriImpl<10,10,11,10,true>(in, out); break;
//                 case 11: AVXIProductTriImpl<11,11,12,11,true>(in, out); break;
//             }
//         }
//         else
//         {
//             switch(m_basis[0]->GetNumModes())
//             {
//                 case 2:  AVXIProductTriImpl<2 ,2 ,3 ,2 ,false>(in, out); break;
//                 case 3:  AVXIProductTriImpl<3 ,3 ,4 ,3 ,false>(in, out); break;
//                 case 4:  AVXIProductTriImpl<4 ,4 ,5 ,4 ,false>(in, out); break;
//                 case 5:  AVXIProductTriImpl<5 ,5 ,6 ,5 ,false>(in, out); break;
//                 case 6:  AVXIProductTriImpl<6 ,6 ,7 ,6 ,false>(in, out); break;
//                 case 7:  AVXIProductTriImpl<7 ,7 ,8 ,7 ,false>(in, out); break;
//                 case 8:  AVXIProductTriImpl<8 ,8 ,9 ,8 ,false>(in, out); break;
//                 case 9:  AVXIProductTriImpl<9 ,9 ,10,9 ,false>(in, out); break;
//                 case 10: AVXIProductTriImpl<10,10,11,10,false>(in, out); break;
//                 case 11: AVXIProductTriImpl<11,11,12,11,false>(in, out); break;
//             }
//         }
//     }

//     template<int NM0, int NM1, int NQ0, int NQ1, bool CORRECT>
//     void AVXIProductTriImpl(
//         const Array<OneD, const NekDouble> &input,
//               Array<OneD,       NekDouble> &output)
//     {
//         using T = VecData<double, VW>;
//         auto *inptr = &input[0];
//         auto *outptr = &output[0];

//         constexpr int nqBlocks = NQ0 * NQ1 * VW;
//         const int nmBlocks = m_nmTot * VW;

//         T eta0_sums[NQ1]; //Sums over eta0 for each value of eta1;

//         for(int e =0; e < this->m_nBlocks; e++)
//         {

//             VecData<NekDouble, VW> *jac_ptr;
//             if(DEFORMED){
//                 jac_ptr = &(this->m_jac[NQ0*NQ1*e]);
//             }
//             else{
//                 jac_ptr = &(this->m_jac[e]);
//             }

//             AVXIProductTriKernel
//                 <NM0, NM1, NQ0, NQ1, VW, CORRECT,
//                  false, false, DEFORMED>(
//                      inptr, this->m_bdata[0], this->m_bdata[1],
//                      this->m_w[0], this->m_w[1], jac_ptr,
//                      eta0_sums, outptr);

//             inptr += nqBlocks;
//             outptr += nmBlocks;
//         }
//     }
// public:
//     static NekDouble FlopsPerElement(
//         const int nm,
//         const int nq0,
//         const int nq1)
//     {
//         int ploop = nm * 4 * nq0 * nq1;
//         int pqloop = (nm * (nm+1) / 2) * nq1 * 3;

//         int corr_loop;
//         if(DEFORMED){
//             //We can't premultiply the jacobian in the outer loop.
//             corr_loop = nq1 * (1 + nq0*5);
//         }
//         else{
//             corr_loop = nq1 * ( 2 + nq0 * 4);
//         }

//         return (ploop + pqloop + corr_loop);
//     }

//     virtual NekDouble GFlops() override
//     {
//         const int nm = m_basis[0]->GetNumModes();
//         const int nq0 = m_basis[0]->GetNumPoints();
//         const int nq1 = m_basis[1]->GetNumPoints();

//         int flops = m_nElmt * AVXIProductTri::FlopsPerElement(nm, nq0, nq1);
//         return flops * 1e-9;
//     }

//     virtual NekDouble NLoads() override
//     {
//         const int nm0 = m_basis[0]->GetNumModes();
//         const int nq0 = m_basis[0]->GetNumPoints();
//         const int nq1 = m_basis[1]->GetNumPoints();

//         int load_pji = nm0 * nq1 * nq0 * 3;
//         int load_pqj = nm0*(nm0+1) * nq1 *0.5 * 3;
//         int load_corr = 1 + nq1*(2 + nq0*3);
//         int load_expected = load_pji + load_pqj + load_corr;

//         return this->m_nElmt * load_expected;
//     }

//     virtual NekDouble NStores() override
//     {
//         const int nm0 = m_basis[0]->GetNumModes();
//         const int nq0 = m_basis[0]->GetNumPoints();
//         const int nq1 = m_basis[1]->GetNumPoints();

//         int store_pj = nm0*nq1;
//         int store_pq = nm0*(nm0+1) / 2.0;
//         int store_corr = 1;
//         int store_expected = store_pj + store_pq + store_corr;

//         return this->m_nElmt * store_expected;
//     }

//     virtual NekDouble Ndof() override
//     {
//         return m_nmTot * this->m_nElmt;
//     }

// private:
//     int m_nmTot;
// };

// template<int VW, bool DEFORMED = false>
// struct AVXIProductHex : public IProduct, public AVXHelper<VW, 3, DEFORMED>
// {
//     AVXIProductHex(std::vector<LibUtilities::BasisSharedPtr> basis,
//                    int nElmt)
//         : IProduct(basis, nElmt),
//           AVXHelper<VW, 3, DEFORMED>(basis, nElmt),
//           m_nmTot(LibUtilities::StdHexData::getNumberOfCoefficients(
//                       this->m_nm[0], this->m_nm[1], this->m_nm[2]))
//     {
//     }

//     static std::shared_ptr<Operator> Create(
//         std::vector<LibUtilities::BasisSharedPtr> basis,
//         int nElmt)
//     {
//         return std::make_shared<AVXIProductHex<VW, DEFORMED>>(basis, nElmt);
//     }

//     virtual void operator()(const Array<OneD, const NekDouble> &in,
//                                   Array<OneD,       NekDouble> &out)
//     {
//         switch(m_basis[0]->GetNumModes())
//         {
//             case 2:  AVXIProductHexImpl<2 ,2 ,2 ,3 ,3 ,3 >(in, out); break;
//             case 3:  AVXIProductHexImpl<3 ,3 ,3 ,4 ,4 ,4 >(in, out); break;
//             case 4:  AVXIProductHexImpl<4 ,4 ,4 ,5 ,5 ,5 >(in, out); break;
//             case 5:  AVXIProductHexImpl<5 ,5 ,5 ,6 ,6 ,6 >(in, out); break;
//             case 6:  AVXIProductHexImpl<6 ,6 ,6 ,7 ,7 ,7 >(in, out); break;
//             case 7:  AVXIProductHexImpl<7 ,7 ,7 ,8 ,8 ,8 >(in, out); break;
//             case 8:  AVXIProductHexImpl<8 ,8 ,8 ,9 ,9 ,9 >(in, out); break;
//             case 9:  AVXIProductHexImpl<9 ,9 ,9 ,10,10,10>(in, out); break;
//             case 10: AVXIProductHexImpl<10,10,10,11,11,11>(in, out); break;
//             case 11: AVXIProductHexImpl<11,11,11,12,12,12>(in, out); break;
//         }
//     }

//     template<int NM0, int NM1, int NM2, int NQ0, int NQ1, int NQ2>
//     void AVXIProductHexImpl(
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

//         for (int e = 0; e < this->m_nBlocks; ++e)
//         {
//             VecData<NekDouble, VW> *jac_ptr;
//             if(DEFORMED){
//                 jac_ptr = &(this->m_jac[NQ0*NQ1*NQ2*e]);
//             }
//             else{
//                 jac_ptr = &(this->m_jac[e]);
//             }

//             AVXIProductHexKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, VW, false, false, DEFORMED>(
//                 inptr, this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
//                 this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
//                 sums_kj, sums_k,
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
//         int loop_kji = nq2 * nq1 * nq0 * 4;
//         int loop_kj = nq2 * nq1 * 3;
//         int loop_k = nq2 * 3;

//         return nm*(loop_kji + nm*(loop_kj + nm*(loop_k)));
//     }

//     virtual NekDouble GFlops() override
//     {
//         const int nm = m_basis[0]->GetNumModes();
//         const int nq0 = m_basis[0]->GetNumPoints();
//         const int nq1 = m_basis[1]->GetNumPoints();
//         const int nq2 = m_basis[2]->GetNumPoints();

//         int expected = AVXIProductHex::FlopsPerElement(nm, nq0, nq1, nq2);
//         int flops = this->m_nElmt * expected;
//         return flops * 1e-9;
//     }

//     virtual NekDouble Ndof() override
//     {
//         return m_nmTot * this->m_nElmt;
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
//         int load_pqkj = nm0 * nm1 * nq2 * nq1 * 3;
//         int load_pqrk = nm0 * nm1 * nm2 * nq2 * 3;
//         int load_expected = load_pkji + load_pqkj + load_pqrk;


//         return this->m_nElmt * load_expected;
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
//         int store_pqk = nm0 * nm1 * nq2;
//         int store_pqr = nm0 * nm1 * nm2;
//         int store_expected = store_pkj + store_pqk + store_pqr;

//         return this->m_nElmt * store_expected;
//     }

// private:
//     /// Padded basis
//     int m_nmTot;
// };

// template<int VW, bool DEFORMED = false>
// struct AVXIProductPrism : public IProduct, public AVXHelper<VW, 3, DEFORMED>
// {
//     AVXIProductPrism(std::vector<LibUtilities::BasisSharedPtr> basis,
//                      int nElmt)
//         : IProduct(basis, nElmt),
//           AVXHelper<VW, 3, DEFORMED>(basis, nElmt),
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
// struct AVXIProductTet : public IProduct, public AVXHelper<VW, 3, DEFORMED>
// {
//     AVXIProductTet(std::vector<LibUtilities::BasisSharedPtr> basis,
//                    int nElmt)
//         : IProduct(basis, nElmt),
//           AVXHelper<VW, 3, DEFORMED>(basis, nElmt),
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

} // namespace AVX
} // namespace Nektar

#endif
