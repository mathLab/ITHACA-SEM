#ifndef AVXPHYSDERIV_H
#define AVXPHYSDERIV_H

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ShapeType.hpp>
#include <LibUtilities/Foundations/Basis.h>

#include "VecData.hpp"
#include "Operator.hpp"
#include "AVXPhysDerivKernels.hpp"

namespace Nektar
{
namespace AVX
{

template<int VW, bool DEFORMED = false>
struct AVXPhysDerivQuad : public PhysDeriv, public AVXHelper<VW,2, DEFORMED>
{
    AVXPhysDerivQuad(std::vector<LibUtilities::BasisSharedPtr> basis,
                     int nElmt)
    : PhysDeriv(basis, nElmt),
      AVXHelper<VW, 2, DEFORMED>(basis, nElmt),
      m_nmTot(LibUtilities::StdQuadData::getNumberOfCoefficients(
                this->m_nm[0], this->m_nm[1]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<AVXPhysDerivQuad<VW, DEFORMED>>(basis, nElmt);
    }

    static NekDouble FlopsPerElement(
        const int nq0,
        const int nq1)
    {
        int derivTensor = 2 * nq0 *nq1 *nq0 + 2 * nq0 * nq1 *nq1;
        int deriv = nq0 *nq1 * 6;
        return (derivTensor + deriv);
    }

    virtual double GFlops() override
    {
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();

        int flops = this->m_nElmt * AVXPhysDerivQuad::FlopsPerElement(nq0, nq1);
        return flops * 1e-9;
    }

    virtual NekDouble Ndof() override
    {
        return m_nmTot * this->m_nElmt;
    }

    virtual NekDouble NLoads() override
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

    virtual NekDouble NStores() override
    {
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();

        int physDerivTensor = nq0*nq1*2;
        int physDeriv = nq0 * nq1 * 2;

        return this->m_nElmt * (physDerivTensor + physDeriv);
    }

    virtual void operator()(const Array<OneD, const NekDouble> &in,
                                  Array<OneD,       NekDouble> &out_d0,
                                  Array<OneD,       NekDouble> &out_d1,
                                  Array<OneD,       NekDouble> &out_d2) override
    {
        //Only for 3D, but need to implement since its abstract
        boost::ignore_unused(in, out_d0, out_d1, out_d2);
        NEKERROR(ErrorUtil::efatal,
                "Something went horribly wrong... calling 3D op for 2D op");
    }

    virtual void operator()(const Array<OneD, const NekDouble> &in,
                                  Array<OneD,       NekDouble> &out_d0,
                                  Array<OneD,       NekDouble> &out_d1) override
    {

        switch(m_basis[0]->GetNumPoints())
        {
            case 2:
                AVXPhysDerivQuadImpl<2,2>(in, out_d0, out_d1); break;
            case 3:
                AVXPhysDerivQuadImpl<3,3>(in, out_d0, out_d1); break;
            case 4:
                AVXPhysDerivQuadImpl<4,4>(in, out_d0, out_d1); break;
            case 5:
                AVXPhysDerivQuadImpl<5,5>(in, out_d0, out_d1); break;
            case 6:
                AVXPhysDerivQuadImpl<6,6>(in, out_d0, out_d1); break;
            case 7:
                AVXPhysDerivQuadImpl<7,7>(in, out_d0, out_d1); break;
            case 8:
                AVXPhysDerivQuadImpl<8,8>(in, out_d0, out_d1); break;
            case 9:
                AVXPhysDerivQuadImpl<9,9>(in, out_d0, out_d1); break;
            case 10:
                AVXPhysDerivQuadImpl<10,10>(in, out_d0, out_d1); break;
            default: NEKERROR(ErrorUtil::efatal,
                "AVXPhysDerivQuad: # of modes / points combo not implemented.");
        }
    }

    template<int NQ0, int NQ1>
    void AVXPhysDerivQuadImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &out_d0,
              Array<OneD,       NekDouble> &out_d1)
    {
        using T = VecData<NekDouble, VW>;
        auto *inptr = &input[0];
        auto *outptr_d0 = &out_d0[0];
        auto *outptr_d1 = &out_d1[0];

        constexpr int ndf = 4;
        constexpr int nqTot = NQ0 * NQ1;
        constexpr int nqBlocks = nqTot * VW;

        // Get size of derivative factor block
        int dfSize{};
        if(DEFORMED)
        {
            dfSize = ndf*nqTot;
        }
        else
        {
            dfSize = ndf;
        }

        AlignedVector<T> tmpIn(nqTot), tmpOut_d0(nqTot), tmpOut_d1(nqTot);
        const T *df_ptr;
        for (int e = 0; e < this->m_nBlocks; e++)
        {
            df_ptr = &(this->m_df[e*dfSize]);

            // Load and transpose data
            T::load_interleave(inptr, nqTot, tmpIn);

            AVXPhysDerivQuadKernel<NQ0, NQ1, VW, DEFORMED>(
                tmpIn,
                this->m_Z[0], this->m_Z[1],
                this->m_D[0], this->m_D[1],
                df_ptr,
                tmpOut_d0, tmpOut_d1);

            // de-interleave and store data
            T::deinterleave_store(tmpOut_d0, nqTot, outptr_d0);
            T::deinterleave_store(tmpOut_d1, nqTot, outptr_d1);

            inptr += nqBlocks;
            outptr_d0 += nqBlocks;
            outptr_d1 += nqBlocks;

        }
    }

private:
    int m_nmTot;

};

// template<int VW, bool DEFORMED=false>
// struct AVXPhysDerivTri : public PhysDeriv, public AVXHelper<VW,2,DEFORMED>
// {
//     AVXPhysDerivTri(std::vector<LibUtilities::BasisSharedPtr> basis,
//                     int nElmt)
//         : PhysDeriv(basis, nElmt),
//           AVXHelper<VW, 2, DEFORMED>(basis, nElmt),
//           m_nmTot(LibUtilities::StdTriData::getNumberOfCoefficients(
//                       this->m_nm[0], this->m_nm[1]))
//     {
//     }

//     static std::shared_ptr<Operator> Create(
//         std::vector<LibUtilities::BasisSharedPtr> basis,
//         int nElmt)
//     {
//         return std::make_shared<AVXPhysDerivTri<VW,DEFORMED>>(basis, nElmt);
//     }

//     static NekDouble FlopsPerElement(
//         const int nq0,
//         const int nq1)
//     {
//         int derivTensor = 2 * nq0 *nq1 *nq0 + 2 * nq0 * nq1 *nq1;

//         int d0 = nq1 * 2 + nq1 * nq0;
//         int d1 = nq1 * nq0 * 4;
//         int df = nq1 * nq0 * 6;

//         return derivTensor + d0 + d1 + df;
//     }

//     virtual double GFlops() override
//     {
//         const int nq0 = m_basis[0]->GetNumPoints();
//         const int nq1 = m_basis[1]->GetNumPoints();

//         int flops = m_nElmt * AVXPhysDerivTri::FlopsPerElement(nq0, nq1);
//         return flops * 1e-9;
//     }

//     virtual NekDouble Ndof() override
//     {
//         return m_nmTot * this->m_nElmt;
//     }


//     virtual NekDouble NLoads() override
//     {
//         const int nq0 = m_basis[0]->GetNumPoints();
//         const int nq1 = m_basis[1]->GetNumPoints();

//         int t_d0 = nq0 * nq1 * nq0 * 2;
//         int t_d1 = nq0 * nq1 * nq1 * 2;
//         int physDerivTensor = t_d0 + t_d1;
//         int physDeriv = nq1*(1 + nq0*3);

//         int load_expected = physDerivTensor + physDeriv;

//         return this->m_nElmt * load_expected;
//     }

//     virtual NekDouble NStores() override
//     {
//         const int nq0 = m_basis[0]->GetNumPoints();
//         const int nq1 = m_basis[1]->GetNumPoints();

//         int physDerivTensor = nq0*nq1*2;
//         int physDeriv = nq1*nq0*2;

//         return this->m_nElmt * (physDerivTensor + physDeriv);
//     }

//     virtual void operator()(const Array<OneD, const NekDouble> &in,
//                                   Array<OneD,       NekDouble> &out_d0,
//                                   Array<OneD,       NekDouble> &out_d1,
//                                   Array<OneD,       NekDouble> &out_d2) override
//     {
//         throw;//Only for 3D, but need to implement since its abstract
//     }

//     virtual void operator()(const Array<OneD, const NekDouble> &in,
//                                   Array<OneD,       NekDouble> &out_d0,
//                                   Array<OneD,       NekDouble> &out_d1)  override
//     {
//         switch(m_basis[0]->GetNumModes())
//         {
//             case 2:
//                 AVXPhysDerivTriImpl<3,2>(in, out_d0, out_d1);
//                 break;
//             case 3:
//                 AVXPhysDerivTriImpl<4,3>(in, out_d0, out_d1);
//                 break;
//             case 4:
//                 AVXPhysDerivTriImpl<5,4>(in, out_d0, out_d1);
//                 break;
//             case 5:
//                 AVXPhysDerivTriImpl<6,5>(in, out_d0, out_d1);
//                 break;
//             case 6:
//                 AVXPhysDerivTriImpl<7,6>(in, out_d0, out_d1);
//                 break;
//             case 7:
//                 AVXPhysDerivTriImpl<8,7>(in, out_d0, out_d1);
//                 break;
//             case 8:
//                 AVXPhysDerivTriImpl<9,8>(in, out_d0, out_d1);
//                 break;
//             case 9:
//                 AVXPhysDerivTriImpl<10,9>(in, out_d0, out_d1);
//                 break;
//             case 10:
//                 AVXPhysDerivTriImpl<11,10>(in, out_d0, out_d1);
//                 break;
//             case 11:
//                 AVXPhysDerivTriImpl<12,11>(in, out_d0, out_d1);
//                 break;
//             case 12:
//                 AVXPhysDerivTriImpl<13,12>(in, out_d0, out_d1);
//                 break;
//             case 13:
//                 AVXPhysDerivTriImpl<14,13>(in, out_d0, out_d1);
//                 break;
//             case 14:
//                 AVXPhysDerivTriImpl<15,14>(in, out_d0, out_d1);
//                 break;
//             case 15:
//                 AVXPhysDerivTriImpl<16,15>(in, out_d0, out_d1);
//                 break;
//             case 16:
//                 AVXPhysDerivTriImpl<17,16>(in, out_d0, out_d1);
//                 break;
//         }
//     }

//     template<int NQ0, int NQ1>
//     void AVXPhysDerivTriImpl(
//         const Array<OneD, const NekDouble> &input,
//               Array<OneD,       NekDouble> &out_d0,
//               Array<OneD,       NekDouble> &out_d1)
//     {
//         using T = VecData<double, VW>;
//         const NekDouble *inptr = &input[0];
//         NekDouble *outptr_d0 = &out_d0[0];
//         NekDouble *outptr_d1 = &out_d1[0];

//         constexpr int nqBlocks = NQ0 * NQ1 * VW;
//         constexpr int ndf = 4;
//         constexpr int nq = NQ0*NQ1;

//         for(int e = 0; e < this->m_nBlocks; e++)
//         {
//             const VecData<double, VW> *df_ptr;
//             if(DEFORMED){
//                 df_ptr = &(this->m_df[e*ndf*nq]);
//             }
//             else{
//                 df_ptr = &(this->m_df[e*ndf]);
//             }

//             AVXPhysDerivTriKernel<NQ0, NQ1, VW, DEFORMED>(
//                 inptr,
//                 this->m_Z[0], this->m_Z[1],
//                 this->m_D[0], this->m_D[1],
//                 df_ptr,
//                 outptr_d0, outptr_d1);

//             inptr += nqBlocks;
//             outptr_d0 += nqBlocks;
//             outptr_d1 += nqBlocks;
//         }
//     }
// private:
//     int m_nmTot;
// };

// template<int VW, bool DEFORMED = false>
// struct AVXPhysDerivHex : public PhysDeriv, public AVXHelper<VW,3,DEFORMED>
// {
//     AVXPhysDerivHex(std::vector<LibUtilities::BasisSharedPtr> basis,
//                     int nElmt)
//         : PhysDeriv(basis, nElmt),
//           AVXHelper<VW, 3,DEFORMED>(basis, nElmt),
//           m_nmTot(LibUtilities::StdHexData::getNumberOfCoefficients(
//                       this->m_nm[0], this->m_nm[1], this->m_nm[2]))
//     {
//     }

//     static std::shared_ptr<Operator> Create(
//         std::vector<LibUtilities::BasisSharedPtr> basis,
//         int nElmt)
//     {
//         return std::make_shared<AVXPhysDerivHex<VW,DEFORMED>>(basis, nElmt);
//     }

//     static NekDouble FlopsPerElement(
//         const int nq0,
//         const int nq1,
//         const int nq2)
//     {
//         //PhysDerivTensor
//         int pdt0 = nq0 * nq1 * nq2 * nq0 * 2;
//         int pdt1 = nq2 * nq0 * nq1 * nq1 * 2;
//         int pdt2 = nq0 * nq1 * nq2 * nq2 * 2;
//         int physDerivTensor = pdt0 + pdt1 + pdt2;

//         //PhysDeriv
//         int physDeriv = nq2 * nq1 * nq0 * 15;

//         return (physDerivTensor + physDeriv);
//     }

//     virtual double GFlops() override
//     {
//         const int nq0 = m_basis[0]->GetNumPoints();
//         const int nq1 = m_basis[1]->GetNumPoints();
//         const int nq2 = m_basis[2]->GetNumPoints();

//         int flops = this->m_nElmt * AVXPhysDerivHex::FlopsPerElement(nq0, nq1, nq2);
//         return flops * 1e-9;
//     }

//     virtual NekDouble Ndof() override
//     {
//         return m_nmTot * this->m_nElmt;
//     }

//     virtual NekDouble NLoads() override
//     {
//         const int nq0 = m_basis[0]->GetNumPoints();
//         const int nq1 = m_basis[1]->GetNumPoints();
//         const int nq2 = m_basis[2]->GetNumPoints();

//         int load_d0 = nq0 * nq1 * nq2 * nq0 * 2;
//         int load_d1 = nq2 * nq0 * nq1 * nq1 * 2;
//         int load_d2 = nq0 * nq1 * nq2 * nq2 * 2;;
//         int physDerivTensor = load_d0 + load_d1 + load_d2;
//         int physDeriv = nq2 * nq1 * nq0 * 3;

//         return m_nElmt * (physDeriv + physDerivTensor);

//     }

//     virtual NekDouble NStores() override
//     {
//         const int nq0 = m_basis[0]->GetNumPoints();
//         const int nq1 = m_basis[1]->GetNumPoints();
//         const int nq2 = m_basis[2]->GetNumPoints();

//         int store_d0 = nq0 * nq1 * nq2;
//         int store_d1 = nq2 * nq0 * nq1;
//         int store_d2 = nq0 * nq1 * nq2;
//         int physDerivTensor = store_d0 + store_d1 + store_d2;
//         int physDeriv = nq2 * nq1 * nq0 * 3;

//         return m_nElmt * (physDeriv + physDerivTensor);
//     }

//     virtual void operator()(const Array<OneD, const NekDouble> &in,
//                                 Array<OneD,       NekDouble> &out_d0,
//                                 Array<OneD,       NekDouble> &out_d1) override
//     {
//         throw; //Only for 2D, but need to implement since its abstract
//     }

//     virtual void operator()(const Array<OneD, const NekDouble> &in,
//                                   Array<OneD,       NekDouble> &out_d0,
//                                   Array<OneD,       NekDouble> &out_d1,
//                                   Array<OneD,       NekDouble> &out_d2) override
//     {
//         switch(m_basis[0]->GetNumModes())
//         {
//             case 2:
//                 AVXPhysDerivHexImpl<3,3,3>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 3:
//                 AVXPhysDerivHexImpl<4,4,4>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 4:
//                 AVXPhysDerivHexImpl<5,5,5>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 5:
//                 AVXPhysDerivHexImpl<6,6,6>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 6:
//                 AVXPhysDerivHexImpl<7,7,7>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 7:
//                 AVXPhysDerivHexImpl<8,8,8>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 8:
//                 AVXPhysDerivHexImpl<9,9,9>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 9:
//                 AVXPhysDerivHexImpl<10,10,10>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 10:
//                 AVXPhysDerivHexImpl<11,11,11>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 11:
//                 AVXPhysDerivHexImpl<12,12,12>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 12:
//                 AVXPhysDerivHexImpl<13,13,13>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 13:
//                 AVXPhysDerivHexImpl<14,14,14>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 14:
//                 AVXPhysDerivHexImpl<15,15,15>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 15:
//                 AVXPhysDerivHexImpl<16,16,16>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 16:
//                 AVXPhysDerivHexImpl<17,17,17>(in, out_d0, out_d1, out_d2);
//                 break;
//         }
//     }

//     template<int NQ0, int NQ1, int NQ2>
//     void AVXPhysDerivHexImpl(
//         const Array<OneD, const NekDouble> &input,
//               Array<OneD,       NekDouble> &out_d0,
//               Array<OneD,       NekDouble> &out_d1,
//               Array<OneD,       NekDouble> &out_d2)
//     {
//         using T = VecData<double, VW>;
//         const NekDouble *inptr = &input[0];
//         NekDouble *outptr_d0 = &out_d0[0];
//         NekDouble *outptr_d1 = &out_d1[0];
//         NekDouble *outptr_d2 = &out_d2[0];

//         constexpr int nqBlocks = NQ0 * NQ1 * NQ2 * VW;
//         constexpr int ndf = 9;
//         constexpr int nq = NQ0*NQ1*NQ2;

//         const VecData<double, VW> *df_ptr;

//         for(int e = 0; e < this->m_nBlocks; e++){

//             if(DEFORMED){
//                 df_ptr = &(this->m_df[e*ndf*nq]);
//             }
//             else{
//                 df_ptr = &(this->m_df[e*ndf]);
//             }

//             AVXPhysDerivHexKernel<NQ0, NQ1, NQ2, VW, DEFORMED>(
//                 inptr,
//                 this->m_Z[0], this->m_Z[1], this->m_Z[2],
//                 this->m_D[0], this->m_D[1], this->m_D[2],
//                 df_ptr,
//                 outptr_d0, outptr_d1, outptr_d2);

//             inptr += nqBlocks;
//             outptr_d0 += nqBlocks;
//             outptr_d1 += nqBlocks;
//             outptr_d2 += nqBlocks;
//         }
//     }
// private:
//     int m_nmTot;
// };


// template<int VW, bool DEFORMED = false>
// struct AVXPhysDerivPrism : public PhysDeriv, public AVXHelper<VW,3,DEFORMED>
// {
//     AVXPhysDerivPrism(std::vector<LibUtilities::BasisSharedPtr> basis,
//                       int nElmt)
//         : PhysDeriv(basis, nElmt),
//           AVXHelper<VW, 3,DEFORMED>(basis, nElmt),
//           m_nmTot(LibUtilities::StdPrismData::getNumberOfCoefficients(
//                       this->m_nm[0], this->m_nm[1], this->m_nm[2]))
//     {
//     }

//     static std::shared_ptr<Operator> Create(
//         std::vector<LibUtilities::BasisSharedPtr> basis,
//         int nElmt)
//     {
//         return std::make_shared<AVXPhysDerivPrism<VW,DEFORMED>>(basis, nElmt);
//     }

//     virtual NekDouble Ndof() override
//     {
//         return m_nmTot * this->m_nElmt;
//     }

//     static NekDouble FlopsPerElement(
//         const int nq0,
//         const int nq1,
//         const int nq2)
//     {
//         //PhysDerivTensor
//         int pdt0 = nq0 * nq1 * nq2 * nq0 * 2;
//         int pdt1 = nq2 * nq0 * nq1 * nq1 * 2;
//         int pdt2 = nq0 * nq1 * nq2 * nq2 * 2;
//         int physDerivTensor = pdt0 + pdt1 + pdt2;

//         //PhysDeriv
//         int physDeriv = nq2*(2 + nq1*nq0*(20));

//         return (physDeriv + physDerivTensor);
//     }

//     virtual double GFlops() override
//     {
//         const int nq0 = m_basis[0]->GetNumPoints();
//         const int nq1 = m_basis[1]->GetNumPoints();
//         const int nq2 = m_basis[2]->GetNumPoints();

//         int flops = m_nElmt * AVXPhysDerivPrism::FlopsPerElement(nq0, nq1, nq2);
//         return flops * 1e-9;
//     }

//     virtual NekDouble NLoads() override
//     {
//         const int nq0 = m_basis[0]->GetNumPoints();
//         const int nq1 = m_basis[1]->GetNumPoints();
//         const int nq2 = m_basis[2]->GetNumPoints();

//         int load_d0 = nq0 * nq1 * nq2 * nq0 * 2;
//         int load_d1 = nq2 * nq0 * nq1 * nq1 * 2;
//         int load_d2 = nq0 * nq1 * nq2 * nq2 * 2;;
//         int physDerivTensor = load_d0 + load_d1 + load_d2;
//         int physDeriv = nq2 * (1 + nq1 * nq0 *4);

//         return m_nElmt * (physDeriv + physDerivTensor);

//     }

//     virtual NekDouble NStores() override
//     {
//         const int nq0 = m_basis[0]->GetNumPoints();
//         const int nq1 = m_basis[1]->GetNumPoints();
//         const int nq2 = m_basis[2]->GetNumPoints();

//         int store_d0 = nq0 * nq1 * nq2;
//         int store_d1 = nq2 * nq0 * nq1;
//         int store_d2 = nq0 * nq1 * nq2;
//         int physDerivTensor = store_d0 + store_d1 + store_d2;
//         int physDeriv = nq2*nq1*nq0*3;

//         return m_nElmt * (physDeriv + physDerivTensor);
//     }

//     virtual void operator()(const Array<OneD, const NekDouble> &in,
//                                 Array<OneD,       NekDouble> &out_d0,
//                                 Array<OneD,       NekDouble> &out_d1) override
//     {
//         throw; //Only for 2D, but need to implement since its abstract
//     }

//     virtual void operator()(const Array<OneD, const NekDouble> &in,
//                                   Array<OneD,       NekDouble> &out_d0,
//                                   Array<OneD,       NekDouble> &out_d1,
//                                   Array<OneD,       NekDouble> &out_d2) override
//     {
//         switch(m_basis[0]->GetNumModes())
//         {
//             case 2:
//                 AVXPhysDerivPrismImpl<3,3,2>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 3:
//                 AVXPhysDerivPrismImpl<4,4,3>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 4:
//                 AVXPhysDerivPrismImpl<5,5,4>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 5:
//                 AVXPhysDerivPrismImpl<6,6,5>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 6:
//                 AVXPhysDerivPrismImpl<7,7,6>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 7:
//                 AVXPhysDerivPrismImpl<8,8,7>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 8:
//                 AVXPhysDerivPrismImpl<9,9,8>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 9:
//                 AVXPhysDerivPrismImpl<10,10,9>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 10:
//                 AVXPhysDerivPrismImpl<11,11,10>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 11:
//                 AVXPhysDerivPrismImpl<12,12,11>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 12:
//                 AVXPhysDerivPrismImpl<13,13,12>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 13:
//                 AVXPhysDerivPrismImpl<14,14,13>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 14:
//                 AVXPhysDerivPrismImpl<15,15,14>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 15:
//                 AVXPhysDerivPrismImpl<16,16,15>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 16:
//                 AVXPhysDerivPrismImpl<17,17,16>(in, out_d0, out_d1, out_d2);
//                 break;
//         }
//     }

//     template<int NQ0, int NQ1, int NQ2>
//     void AVXPhysDerivPrismImpl(
//         const Array<OneD, const NekDouble> &input,
//               Array<OneD,       NekDouble> &out_d0,
//               Array<OneD,       NekDouble> &out_d1,
//               Array<OneD,       NekDouble> &out_d2)
//     {
//         using T = VecData<double, VW>;
//         const NekDouble *inptr = &input[0];
//         NekDouble *outptr_d0 = &out_d0[0];
//         NekDouble *outptr_d1 = &out_d1[0];
//         NekDouble *outptr_d2 = &out_d2[0];

//         constexpr int nqBlocks = NQ0 * NQ1 * NQ2 * VW;
//         constexpr int ndf = 9;
//         constexpr int nq = NQ0*NQ1*NQ2;

//         const VecData<double, VW> *df_ptr;

//         for(int e = 0; e < this->m_nBlocks; e++){

//             if(DEFORMED){
//                 df_ptr = &(this->m_df[e*ndf*nq]);
//             }
//             else{
//                 df_ptr = &(this->m_df[e*ndf]);
//             }

//             AVXPhysDerivPrismKernel<NQ0, NQ1, NQ2, VW,DEFORMED>(
//                 inptr,
//                 this->m_Z[0], this->m_Z[1], this->m_Z[2],
//                 this->m_D[0], this->m_D[1], this->m_D[2],
//                 df_ptr,
//                 outptr_d0, outptr_d1, outptr_d2);

//             inptr += nqBlocks;
//             outptr_d0 += nqBlocks;
//             outptr_d1 += nqBlocks;
//             outptr_d2 += nqBlocks;
//         }
//     }
// private:
//     int m_nmTot;

// };

// template<int VW, bool DEFORMED = false>
// struct AVXPhysDerivTet : public PhysDeriv, public AVXHelper<VW,3, DEFORMED>
// {
//     AVXPhysDerivTet(std::vector<LibUtilities::BasisSharedPtr> basis,
//                     int nElmt)
//         : PhysDeriv(basis, nElmt),
//           AVXHelper<VW, 3, DEFORMED>(basis, nElmt),
//           m_nmTot(LibUtilities::StdTetData::getNumberOfCoefficients(
//                       this->m_nm[0], this->m_nm[1], this->m_nm[2]))
//     {
//     }

//     static std::shared_ptr<Operator> Create(
//         std::vector<LibUtilities::BasisSharedPtr> basis,
//         int nElmt)
//     {
//         return std::make_shared<AVXPhysDerivTet<VW, DEFORMED>>(basis, nElmt);
//     }

//     static NekDouble FlopsPerElement(
//         const int nq0,
//         const int nq1,
//         const int nq2)
//     {
//         //PhysDerivTensor
//         int pdt0 = nq0 * nq1 * nq2 * nq0 * 2;
//         int pdt1 = nq2 * nq0 * nq1 * nq1 * 2;
//         int pdt2 = nq0 * nq1 * nq2 * nq2 * 2;
//         int physDerivTensor = pdt0 + pdt1 + pdt2;

//         //PhysDeriv
//         int d0 = nq2*(2 + nq1*(3 + nq0));
//         int d1 = nq2*(2 + nq1*nq0*5);
//         int d2 = nq2*nq1*(2 + nq0*3);
//         int df = nq2*nq1*nq0*15;
//         int physDeriv = d0 + d1 + d2 + df;

//         return (physDerivTensor + physDeriv);
//     }

//     virtual double GFlops() override
//     {
//         const int nq0 = m_basis[0]->GetNumPoints();
//         const int nq1 = m_basis[1]->GetNumPoints();
//         const int nq2 = m_basis[2]->GetNumPoints();

//         int flops = m_nElmt * AVXPhysDerivTet::FlopsPerElement(nq0, nq1, nq2);
//         return flops * 1e-9;
//     }

//     virtual NekDouble NLoads() override
//     {
//         const int nq0 = m_basis[0]->GetNumPoints();
//         const int nq1 = m_basis[1]->GetNumPoints();
//         const int nq2 = m_basis[2]->GetNumPoints();

//         //physDerivTensor
//         int load_d0 = nq0 * nq1 * nq2 * nq0 * 2;
//         int load_d1 = nq2 * nq0 * nq1 * nq1 * 2;
//         int load_d2 = nq0 * nq1 * nq2 * nq2 * 2;;
//         int physDerivTensor = load_d0 + load_d1 + load_d2;

//         int load_dd0 = nq2 * (1 + nq1 * (1+ nq0));
//         int load_dd1 = nq2 * (1 + nq1 * nq0 * 3);
//         int load_dd2 = nq2 * (nq1 * (1 + nq0*3));
//         int load_coef = nq2 * nq1 * nq0 * 3;
//         int physDeriv = load_dd0 + load_dd1 + load_dd2 + load_coef;

//         return m_nElmt * (physDeriv + physDerivTensor);

//     }

//     virtual NekDouble NStores() override
//     {
//         const int nq0 = m_basis[0]->GetNumPoints();
//         const int nq1 = m_basis[1]->GetNumPoints();
//         const int nq2 = m_basis[2]->GetNumPoints();

//         int store_d0 = nq0 * nq1 * nq2;
//         int store_d1 = nq2 * nq0 * nq1;
//         int store_d2 = nq0 * nq1 * nq2;
//         int physDerivTensor = store_d0 + store_d1 + store_d2;

//         int store_dd = nq2 * nq1 * nq0 * 2;
//         int store_coef = nq2 * nq1 *nq0 * 3;
//         int physDeriv = store_dd * 3 + store_coef;

//         return m_nElmt * (physDeriv + physDerivTensor);
//     }

//     virtual NekDouble Ndof() override
//     {
//         return m_nmTot * this->m_nElmt;
//     }

//     virtual void operator()(const Array<OneD, const NekDouble> &in,
//                                 Array<OneD,       NekDouble> &out_d0,
//                                 Array<OneD,       NekDouble> &out_d1) override
//     {
//         throw; //Only for 2D, but need to implement since its abstract
//     }

//     virtual void operator()(const Array<OneD, const NekDouble> &in,
//                                   Array<OneD,       NekDouble> &out_d0,
//                                   Array<OneD,       NekDouble> &out_d1,
//                                   Array<OneD,       NekDouble> &out_d2) override
//     {
//         switch(m_basis[0]->GetNumModes())
//         {
//             case 2:
//                 AVXPhysDerivTetImpl<3,2,2>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 3:
//                 AVXPhysDerivTetImpl<4,3,3>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 4:
//                 AVXPhysDerivTetImpl<5,4,4>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 5:
//                 AVXPhysDerivTetImpl<6,5,5>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 6:
//                 AVXPhysDerivTetImpl<7,6,6>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 7:
//                 AVXPhysDerivTetImpl<8,7,7>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 8:
//                 AVXPhysDerivTetImpl<9,8,8>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 9:
//                 AVXPhysDerivTetImpl<10,9,9>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 10:
//                 AVXPhysDerivTetImpl<11,10,10>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 11:
//                 AVXPhysDerivTetImpl<12,11,11>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 12:
//                 AVXPhysDerivTetImpl<13,12,12>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 13:
//                 AVXPhysDerivTetImpl<14,13,13>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 14:
//                 AVXPhysDerivTetImpl<15,14,14>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 15:
//                 AVXPhysDerivTetImpl<16,15,15>(in, out_d0, out_d1, out_d2);
//                 break;
//             case 16:
//                 AVXPhysDerivTetImpl<17,16,16>(in, out_d0, out_d1, out_d2);
//                 break;
//         }
//     }

//     template<int NQ0, int NQ1, int NQ2>
//     void AVXPhysDerivTetImpl(
//         const Array<OneD, const NekDouble> &input,
//               Array<OneD,       NekDouble> &out_d0,
//               Array<OneD,       NekDouble> &out_d1,
//               Array<OneD,       NekDouble> &out_d2)
//     {
//         using T = VecData<double, VW>;
//         const NekDouble *inptr = &input[0];
//         NekDouble *outptr_d0 = &out_d0[0];
//         NekDouble *outptr_d1 = &out_d1[0];
//         NekDouble *outptr_d2 = &out_d2[0];


//         constexpr int nqBlocks = NQ0 * NQ1 * NQ2 * VW;
//         constexpr int ndf = 9;
//         constexpr int nq = NQ0 * NQ1 * NQ2;

//         AlignedVector<NekDouble> diff0_arr(nqBlocks), diff1_arr(nqBlocks), diff2_arr(nqBlocks);
//         NekDouble *diff0 = diff0_arr.data();
//         NekDouble *diff1 = diff1_arr.data();
//         NekDouble *diff2 = diff2_arr.data();

//         const VecData<double, VW> *df_ptr;

//         for(int e = 0; e < this->m_nBlocks; e++){

//             if(DEFORMED){
//                 df_ptr = &(this->m_df[e*ndf*nq]);
//             }
//             else{
//                 df_ptr = &(this->m_df[e*ndf]);
//             }
//             AVXPhysDerivTetKernel<NQ0, NQ1, NQ2, VW,DEFORMED>(
//                 inptr,
//                 this->m_Z[0], this->m_Z[1], this->m_Z[2],
//                 this->m_D[0], this->m_D[1], this->m_D[2],
//                 df_ptr,
//                 diff0, diff1, diff2,
//                 outptr_d0, outptr_d1, outptr_d2);

//             inptr += nqBlocks;
//             outptr_d0 += nqBlocks;
//             outptr_d1 += nqBlocks;
//             outptr_d2 += nqBlocks;
//         }
//     }
// private:
//     int m_nmTot;
// };

} // namespace AVX
} // namespace Nektar

#endif
