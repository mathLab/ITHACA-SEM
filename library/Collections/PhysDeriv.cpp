///////////////////////////////////////////////////////////////////////////////
//
// File: PhysDeriv.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: PhysDeriv operator implementations
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <Collections/Operator.h>
#include <Collections/Collection.h>

using namespace std;

namespace Nektar {
namespace Collections {

using LibUtilities::eSegment;
using LibUtilities::eQuadrilateral;
using LibUtilities::eTriangle;
using LibUtilities::eHexahedron;
using LibUtilities::eTetrahedron;
using LibUtilities::ePrism;
using LibUtilities::ePyramid;

/**
 * @brief Phys deriv operator using standard matrix approach
 */
class PhysDeriv_StdMat : public Operator
{
    public:
        OPERATOR_CREATE(PhysDeriv_StdMat)

        virtual ~PhysDeriv_StdMat()
        {
        }

        virtual void operator()(
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output0,
                      Array<OneD,       NekDouble> &output1,
                      Array<OneD,       NekDouble> &output2,
                      Array<OneD,       NekDouble> &wsp)
        {
            int nPhys = m_stdExp->GetTotPoints();
            int ntot = m_numElmt*nPhys;
            Array<OneD, NekDouble> tmp0,tmp1,tmp2;
            Array<OneD, Array<OneD, NekDouble> > Diff(3);
            Array<OneD, Array<OneD, NekDouble> > out(3);
            out[0] = output0;  out[1] = output1;    out[2] = output2;

            for(int i = 0; i < m_dim; ++i)
            {
                Diff[i] = wsp + i*ntot;
            }

            // calculate local derivatives
            for(int i = 0; i < m_dim; ++i)
            {
                Blas::Dgemm('N', 'N', m_derivMat[i]->GetRows(), m_numElmt,
                            m_derivMat[i]->GetColumns(), 1.0,
                            m_derivMat[i]->GetRawPtr(),
                            m_derivMat[i]->GetRows(), input.get(), nPhys,
                            0.0, &Diff[i][0],nPhys);
            }

            // calculate full derivative
            for(int i = 0; i < m_coordim; ++i)
            {
                Vmath::Zero(ntot,out[i],1);
                for(int j = 0; j < m_dim; ++j)
                {
                    Vmath::Vvtvp (ntot, m_derivFac[i*m_dim+j], 1,
                                        Diff[j],               1,
                                        out[i],                1,
                                        out[i],                1);
                }
            }
        }

        virtual void operator()(
                      int                           dir,
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output,
                      Array<OneD,       NekDouble> &wsp)
        {
            int nPhys = m_stdExp->GetTotPoints();
            int ntot = m_numElmt*nPhys;
            Array<OneD, NekDouble> tmp0,tmp1,tmp2;
            Array<OneD, Array<OneD, NekDouble> > Diff(3);

            for(int i = 0; i < m_dim; ++i)
            {
                Diff[i] = wsp + i*ntot;
            }

            // calculate local derivatives
            for(int i = 0; i < m_dim; ++i)
            {
                Blas::Dgemm('N', 'N', m_derivMat[i]->GetRows(), m_numElmt,
                            m_derivMat[i]->GetColumns(), 1.0,
                            m_derivMat[i]->GetRawPtr(),
                            m_derivMat[i]->GetRows(), input.get(), nPhys,
                            0.0, &Diff[i][0],nPhys);
            }

            // calculate full derivative
            Vmath::Zero(ntot,output,1);
            for(int j = 0; j < m_dim; ++j)
            {
                Vmath::Vvtvp (ntot, m_derivFac[dir*m_dim+j], 1,
                                    Diff[j],               1,
                                    output,                1,
                                    output,                1);
            }
        }

    protected:
        Array<OneD, DNekMatSharedPtr>   m_derivMat;
        Array<TwoD, const NekDouble>    m_derivFac;
        int                             m_dim;
        int                             m_coordim;

    private:
        PhysDeriv_StdMat(
                vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                CoalescedGeomDataSharedPtr                pGeomData)
            : Operator(pCollExp, pGeomData)
        {
            int nqtot = 1;
            LibUtilities::PointsKeyVector PtsKey = m_stdExp->GetPointsKeys();
            m_dim = PtsKey.size();
            m_coordim = pCollExp[0]->GetCoordim();

            for(int i = 0; i < m_dim; ++i)
            {
                nqtot *= PtsKey[i].GetNumPoints();
            }
            // set up a PhysDeriv StdMat.
            m_derivMat = Array<OneD, DNekMatSharedPtr>(m_dim);
            for(int i = 0; i < m_dim; ++i)
            {
                Array<OneD, NekDouble> tmp(nqtot),tmp1(nqtot);
                m_derivMat[i] = MemoryManager<DNekMat>
                                            ::AllocateSharedPtr(nqtot,nqtot);
                for(int j = 0; j < nqtot; ++j)
                {
                    Vmath::Zero(nqtot,tmp,1);
                    tmp[j] = 1.0;
                    m_stdExp->PhysDeriv(i,tmp,tmp1);
                    Vmath::Vcopy(nqtot, &tmp1[0], 1,
                                 &(m_derivMat[i]->GetPtr())[0] + j*nqtot, 1);
                }
            }
            m_derivFac = pGeomData->GetDerivFactors(pCollExp);
            m_wspSize = 3*nqtot*m_numElmt;
        }
};

/// Factory initialisation for the PhysDeriv_StdMat operators
OperatorKey PhysDeriv_StdMat::m_typeArr[] =
{
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eSegment,       ePhysDeriv, eStdMat, false),
        PhysDeriv_StdMat::create, "PhysDeriv_StdMat_Seg"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTriangle,      ePhysDeriv, eStdMat, false),
        PhysDeriv_StdMat::create, "PhysDeriv_StdMat_Tri"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTriangle,      ePhysDeriv, eStdMat, true),
        PhysDeriv_StdMat::create, "PhysDeriv_StdMat_NodalTri"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eQuadrilateral, ePhysDeriv, eStdMat, false),
        PhysDeriv_StdMat::create, "PhysDeriv_StdMat_Quad"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTetrahedron,   ePhysDeriv, eStdMat, false),
        PhysDeriv_StdMat::create, "PhysDeriv_StdMat_Tet"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTetrahedron,   ePhysDeriv, eStdMat, true),
        PhysDeriv_StdMat::create, "PhysDeriv_StdMat_NodalTet"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePyramid,       ePhysDeriv, eStdMat, false),
        PhysDeriv_StdMat::create, "PhysDeriv_StdMat_Pyr"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePrism,         ePhysDeriv, eStdMat, false),
        PhysDeriv_StdMat::create, "PhysDeriv_StdMat_Prism"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePrism,         ePhysDeriv, eStdMat, true),
        PhysDeriv_StdMat::create, "PhysDeriv_StdMat_NodalPrism"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eHexahedron,    ePhysDeriv, eStdMat, false),
        PhysDeriv_StdMat::create, "PhysDeriv_StdMat_Hex"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePyramid, ePhysDeriv, eSumFac, false),
        PhysDeriv_StdMat::create, "PhysDeriv_SumFac_Pyr")
};


/**
 * @brief Phys deriv operator using element-wise operation
 */
class PhysDeriv_IterPerExp : public Operator
{
    public:
        OPERATOR_CREATE(PhysDeriv_IterPerExp)

        virtual ~PhysDeriv_IterPerExp()
        {
        }

        virtual void operator()(
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output0,
                      Array<OneD,       NekDouble> &output1,
                      Array<OneD,       NekDouble> &output2,
                      Array<OneD,       NekDouble> &wsp)
        {
            int nPhys = m_stdExp->GetTotPoints();
            int ntot = m_numElmt*nPhys;
            Array<OneD, NekDouble> tmp0,tmp1,tmp2;
            Array<OneD, Array<OneD, NekDouble> > Diff(3);
            Array<OneD, Array<OneD, NekDouble> > out(3);
            out[0] = output0;  out[1] = output1;  out[2] = output2;

            for(int i = 0; i < m_dim; ++i)
            {
                Diff[i] = wsp + i*ntot;
            }

            // calculate local derivatives
            for (int i = 0; i < m_numElmt; ++i)
            {
                m_stdExp->PhysDeriv(input + i*nPhys,
                                    tmp0 = Diff[0] + i*nPhys,
                                    tmp1 = Diff[1] + i*nPhys,
                                    tmp2 = Diff[2] + i*nPhys);
            }

            // calculate full derivative
            for(int i = 0; i < m_coordim; ++i)
            {
                Vmath::Vmul(ntot,m_derivFac[i*m_dim],1,Diff[0],1,out[i],1);
                for(int j = 1; j < m_dim; ++j)
                {
                    Vmath::Vvtvp (ntot, m_derivFac[i*m_dim+j], 1,
                                        Diff[j],               1,
                                        out[i],                1,
                                        out[i],                1);
                }
            }
        }

        virtual void operator()(
                      int                           dir,
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output,
                      Array<OneD,       NekDouble> &wsp)
        {
            int nPhys = m_stdExp->GetTotPoints();
            int ntot = m_numElmt*nPhys;
            Array<OneD, NekDouble> tmp0,tmp1,tmp2;
            Array<OneD, Array<OneD, NekDouble> > Diff(3);

            for(int i = 0; i < m_dim; ++i)
            {
                Diff[i] = wsp + i*ntot;
            }

            // calculate local derivatives
            for (int i = 0; i < m_numElmt; ++i)
            {
                m_stdExp->PhysDeriv(input + i*nPhys,
                                    tmp0 = Diff[0] + i*nPhys,
                                    tmp1 = Diff[1] + i*nPhys,
                                    tmp2 = Diff[2] + i*nPhys);
            }

            // calculate full derivative
            Vmath::Vmul(ntot,m_derivFac[dir*m_dim],1,Diff[0],1,output,1);
            for(int j = 1; j < m_dim; ++j)
            {
                Vmath::Vvtvp (ntot, m_derivFac[dir*m_dim+j], 1,
                                    Diff[j],               1,
                                    output,                1,
                                    output,                1);
            }
        }

    protected:
        Array<TwoD, const NekDouble>    m_derivFac;
        int                             m_dim;
        int                             m_coordim;

    private:
        PhysDeriv_IterPerExp(
                vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                CoalescedGeomDataSharedPtr                pGeomData)
            : Operator(pCollExp, pGeomData)
        {
            int nqtot = 1;
            LibUtilities::PointsKeyVector PtsKey = m_stdExp->GetPointsKeys();
            m_dim = PtsKey.size();
            m_coordim = pCollExp[0]->GetCoordim();

            for(int i = 0; i < m_dim; ++i)
            {
                nqtot *= PtsKey[i].GetNumPoints();
            }
            m_derivFac = pGeomData->GetDerivFactors(pCollExp);
            m_wspSize = 3*nqtot*m_numElmt;
        }
};

/// Factory initialisation for the PhysDeriv_IterPerExp operators
OperatorKey PhysDeriv_IterPerExp::m_typeArr[] =
{
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eSegment,       ePhysDeriv, eIterPerExp,false),
        PhysDeriv_IterPerExp::create, "PhysDeriv_IterPerExp_Seg"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTriangle,      ePhysDeriv, eIterPerExp,false),
        PhysDeriv_IterPerExp::create, "PhysDeriv_IterPerExp_Tri"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTriangle,      ePhysDeriv, eIterPerExp,true),
        PhysDeriv_IterPerExp::create, "PhysDeriv_IterPerExp_NodalTri"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eQuadrilateral, ePhysDeriv, eIterPerExp,false),
        PhysDeriv_IterPerExp::create, "PhysDeriv_IterPerExp_Quad"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTetrahedron,   ePhysDeriv, eIterPerExp,false),
        PhysDeriv_IterPerExp::create, "PhysDeriv_IterPerExp_Tet"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTetrahedron,   ePhysDeriv, eIterPerExp,true),
        PhysDeriv_IterPerExp::create, "PhysDeriv_IterPerExp_NodalTet"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePyramid,       ePhysDeriv, eIterPerExp,false),
        PhysDeriv_IterPerExp::create, "PhysDeriv_IterPerExp_Pyr"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePrism,         ePhysDeriv, eIterPerExp,false),
        PhysDeriv_IterPerExp::create, "PhysDeriv_IterPerExp_Prism"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePrism,         ePhysDeriv, eIterPerExp,true),
        PhysDeriv_IterPerExp::create, "PhysDeriv_IterPerExp_NodalPrism"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eHexahedron,    ePhysDeriv, eIterPerExp,false),
        PhysDeriv_IterPerExp::create, "PhysDeriv_IterPerExp_Hex")
};


/**
 * @brief Phys deriv operator using original LocalRegions implementation.
 */
class PhysDeriv_NoCollection : public Operator
{
    public:
        OPERATOR_CREATE(PhysDeriv_NoCollection)

        virtual ~PhysDeriv_NoCollection()
        {
        }

        virtual void operator()(
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output0,
                      Array<OneD,       NekDouble> &output1,
                      Array<OneD,       NekDouble> &output2,
                      Array<OneD,       NekDouble> &wsp)
        {
            boost::ignore_unused(wsp);

            const int nPhys   = m_expList[0]->GetTotPoints();
            Array<OneD, NekDouble> tmp0,tmp1,tmp2;

            // calculate local derivatives
            switch (m_expList[0]->GetShapeDimension())
            {
                case 1:
                {
                    for (int i = 0; i < m_numElmt; ++i)
                    {
                        m_expList[i]->PhysDeriv(input + i*nPhys,
                                        tmp0 = output0 + i*nPhys);
                    }
                    break;
                }
                case 2:
                {
                    for (int i = 0; i < m_numElmt; ++i)
                    {
                        m_expList[i]->PhysDeriv(input + i*nPhys,
                                        tmp0 = output0 + i*nPhys,
                                        tmp1 = output1 + i*nPhys);
                    }
                    break;
                }
                case 3:
                {
                    for (int i = 0; i < m_numElmt; ++i)
                    {
                        m_expList[i]->PhysDeriv(input + i*nPhys,
                                        tmp0 = output0 + i*nPhys,
                                        tmp1 = output1 + i*nPhys,
                                        tmp2 = output2 + i*nPhys);
                    }
                    break;
                }
                default:
                    ASSERTL0(false, "Unknown dimension.");
            }
        }

        virtual void operator()(
                      int                           dir,
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output,
                      Array<OneD,       NekDouble> &wsp)
        {
            boost::ignore_unused(wsp);

            const int nPhys   = m_expList[0]->GetTotPoints();
            Array<OneD, NekDouble> tmp;

            // calculate local derivatives
            for (int i = 0; i < m_numElmt; ++i)
            {
                m_expList[i]->PhysDeriv(dir, input + i*nPhys,
                                             tmp = output + i*nPhys);
            }
        }

    protected:
        vector<StdRegions::StdExpansionSharedPtr> m_expList;

    private:
        PhysDeriv_NoCollection(
                vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                CoalescedGeomDataSharedPtr                pGeomData)
            : Operator(pCollExp,pGeomData)
        {
            m_expList = pCollExp;
        }
};

/// Factory initialisation for the PhysDeriv_NoCollection operators
OperatorKey PhysDeriv_NoCollection::m_typeArr[] =
{
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eSegment,       ePhysDeriv, eNoCollection,false),
        PhysDeriv_NoCollection::create, "PhysDeriv_NoCollection_Seg"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTriangle,      ePhysDeriv, eNoCollection,false),
        PhysDeriv_NoCollection::create, "PhysDeriv_NoCollection_Tri"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTriangle,      ePhysDeriv, eNoCollection,true),
        PhysDeriv_NoCollection::create, "PhysDeriv_NoCollection_NodalTri"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eQuadrilateral, ePhysDeriv, eNoCollection,false),
        PhysDeriv_NoCollection::create, "PhysDeriv_NoCollection_Quad"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTetrahedron,   ePhysDeriv, eNoCollection,false),
        PhysDeriv_NoCollection::create, "PhysDeriv_NoCollection_Tet"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTetrahedron,   ePhysDeriv, eNoCollection,true),
        PhysDeriv_NoCollection::create, "PhysDeriv_NoCollection_NodalTet"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePyramid,       ePhysDeriv, eNoCollection,false),
        PhysDeriv_NoCollection::create, "PhysDeriv_NoCollection_Pyr"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePrism,         ePhysDeriv, eNoCollection,false),
        PhysDeriv_NoCollection::create, "PhysDeriv_NoCollection_Prism"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePrism,         ePhysDeriv, eNoCollection,true),
        PhysDeriv_NoCollection::create, "PhysDeriv_NoCollection_NodalPrism"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eHexahedron,    ePhysDeriv, eNoCollection,false),
        PhysDeriv_NoCollection::create, "PhysDeriv_NoCollection_Hex")
};


/**
 * @brief Phys deriv operator using sum-factorisation (Segment)
 */
class PhysDeriv_SumFac_Seg : public Operator
{
    public:
        OPERATOR_CREATE(PhysDeriv_SumFac_Seg)

        virtual ~PhysDeriv_SumFac_Seg()
        {
        }

        virtual void operator()(
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output0,
                      Array<OneD,       NekDouble> &output1,
                      Array<OneD,       NekDouble> &output2,
                      Array<OneD,       NekDouble> &wsp)
        {
            const int nqcol   = m_nquad0*m_numElmt;

            ASSERTL1(wsp.size() == m_wspSize,
                     "Incorrect workspace size");
            ASSERTL1(input.size() >= nqcol,
                     "Incorrect input size");

            Array<OneD, NekDouble> diff0(nqcol, wsp);

            Blas::Dgemm('N', 'N', m_nquad0, m_numElmt,
                        m_nquad0, 1.0, m_Deriv0, m_nquad0,
                        input.get(), m_nquad0, 0.0,
                        diff0.get(), m_nquad0);

            Vmath::Vmul  (nqcol, m_derivFac[0], 1, diff0, 1, output0, 1);

            if (m_coordim == 2)
            {
                Vmath::Vmul  (nqcol, m_derivFac[1], 1, diff0, 1, output1, 1);
            }
            else if (m_coordim == 3)
            {
                Vmath::Vmul  (nqcol, m_derivFac[1], 1, diff0, 1, output1, 1);
                Vmath::Vmul  (nqcol, m_derivFac[2], 1, diff0, 1, output2, 1);
            }
        }

        virtual void operator()(
                      int                           dir,
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output,
                      Array<OneD,       NekDouble> &wsp)
        {
            const int nqcol   = m_nquad0*m_numElmt;

            ASSERTL1(wsp.size() == m_wspSize,
                     "Incorrect workspace size");
            ASSERTL1(input.size() >= nqcol,
                     "Incorrect input size");

            Array<OneD, NekDouble> diff0(nqcol, wsp);

            Blas::Dgemm('N', 'N', m_nquad0, m_numElmt,
                        m_nquad0, 1.0, m_Deriv0, m_nquad0,
                        input.get(), m_nquad0, 0.0,
                        diff0.get(), m_nquad0);

            Vmath::Vmul(nqcol, m_derivFac[dir], 1, diff0, 1, output, 1);
        }

    protected:
        int                             m_coordim;
        const int                       m_nquad0;
        Array<TwoD, const NekDouble>    m_derivFac;
        NekDouble                      *m_Deriv0;

    private:
        PhysDeriv_SumFac_Seg(
                vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                CoalescedGeomDataSharedPtr                pGeomData)
            : Operator (pCollExp, pGeomData),
              m_nquad0 (m_stdExp->GetNumPoints(0))
        {
            LibUtilities::PointsKeyVector PtsKey = m_stdExp->GetPointsKeys();
            m_coordim = pCollExp[0]->GetCoordim();

            m_derivFac = pGeomData->GetDerivFactors(pCollExp);

            m_Deriv0 = &((m_stdExp->GetBasis(0)->GetD())->GetPtr())[0];
            m_wspSize = m_nquad0*m_numElmt;
        }

};

/// Factory initialisation for the PhysDeriv_SumFac_Seg operators
OperatorKey PhysDeriv_SumFac_Seg::m_type = GetOperatorFactory().
    RegisterCreatorFunction(
        OperatorKey(eSegment, ePhysDeriv, eSumFac,false),
        PhysDeriv_SumFac_Seg::create, "PhysDeriv_SumFac_Seg");



/**
 * @brief Phys deriv operator using sum-factorisation (Quad)
 */
class PhysDeriv_SumFac_Quad : public Operator
{
    public:
        OPERATOR_CREATE(PhysDeriv_SumFac_Quad)

        virtual ~PhysDeriv_SumFac_Quad()
        {
        }

        virtual void operator()(
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output0,
                      Array<OneD,       NekDouble> &output1,
                      Array<OneD,       NekDouble> &output2,
                      Array<OneD,       NekDouble> &wsp)
        {
            const int nqtot   = m_nquad0 * m_nquad1;
            const int nqcol   = nqtot*m_numElmt;

            ASSERTL1(wsp.size() == m_wspSize,
                     "Incorrect workspace size");
            ASSERTL1(input.size() >= nqcol,
                     "Incorrect input size");

            Array<OneD, NekDouble> diff0(nqcol, wsp             );
            Array<OneD, NekDouble> diff1(nqcol, wsp    +   nqcol);

            Blas::Dgemm('N', 'N', m_nquad0, m_nquad1*m_numElmt,
                        m_nquad0, 1.0, m_Deriv0, m_nquad0,
                        input.get(), m_nquad0, 0.0,
                        diff0.get(), m_nquad0);

            int cnt = 0;
            for (int i = 0; i < m_numElmt; ++i, cnt += nqtot)
            {
                Blas::Dgemm('N', 'T', m_nquad0, m_nquad1, m_nquad1, 1.0,
                            input.get() + cnt, m_nquad0,
                            m_Deriv1, m_nquad1, 0.0,
                            diff1.get() + cnt, m_nquad0);
            }

            Vmath::Vmul  (nqcol, m_derivFac[0], 1, diff0, 1, output0, 1);
            Vmath::Vvtvp (nqcol, m_derivFac[1], 1, diff1, 1, output0, 1,
                                                             output0, 1);
            Vmath::Vmul  (nqcol, m_derivFac[2], 1, diff0, 1, output1, 1);
            Vmath::Vvtvp (nqcol, m_derivFac[3], 1, diff1, 1, output1, 1,
                                                             output1, 1);

            if (m_coordim == 3)
            {
                Vmath::Vmul  (nqcol, m_derivFac[4], 1, diff0, 1, output2, 1);
                Vmath::Vvtvp (nqcol, m_derivFac[5], 1, diff1, 1, output2, 1,
                                                                 output2, 1);
            }
        }

        virtual void operator()(
                      int                           dir,
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output,
                      Array<OneD,       NekDouble> &wsp)
        {
            const int nqtot   = m_nquad0 * m_nquad1;
            const int nqcol   = nqtot*m_numElmt;

            ASSERTL1(wsp.size() == m_wspSize,
                     "Incorrect workspace size");
            ASSERTL1(input.size() >= nqcol,
                     "Incorrect input size");

            Array<OneD, NekDouble> diff0(nqcol, wsp             );
            Array<OneD, NekDouble> diff1(nqcol, wsp    +   nqcol);

            Blas::Dgemm('N', 'N', m_nquad0, m_nquad1*m_numElmt,
                        m_nquad0, 1.0, m_Deriv0, m_nquad0,
                        input.get(), m_nquad0, 0.0,
                        diff0.get(), m_nquad0);

            int cnt = 0;
            for (int i = 0; i < m_numElmt; ++i, cnt += nqtot)
            {
                Blas::Dgemm('N', 'T', m_nquad0, m_nquad1, m_nquad1, 1.0,
                            input.get() + cnt, m_nquad0,
                            m_Deriv1, m_nquad1, 0.0,
                            diff1.get() + cnt, m_nquad0);
            }

            Vmath::Vmul  (nqcol, m_derivFac[2*dir]  , 1, diff0, 1, output, 1);
            Vmath::Vvtvp (nqcol, m_derivFac[2*dir+1], 1, diff1, 1, output, 1,
                                                                   output, 1);
        }

    protected:
        int                             m_coordim;
        const int                       m_nquad0;
        const int                       m_nquad1;
        Array<TwoD, const NekDouble>    m_derivFac;
        NekDouble                      *m_Deriv0;
        NekDouble                      *m_Deriv1;

    private:
        PhysDeriv_SumFac_Quad(
                vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                CoalescedGeomDataSharedPtr                pGeomData)
            : Operator (pCollExp, pGeomData),
              m_nquad0 (m_stdExp->GetNumPoints(0)),
              m_nquad1 (m_stdExp->GetNumPoints(1))
        {
            LibUtilities::PointsKeyVector PtsKey = m_stdExp->GetPointsKeys();
            m_coordim = pCollExp[0]->GetCoordim();

            m_derivFac = pGeomData->GetDerivFactors(pCollExp);

            m_Deriv0 = &((m_stdExp->GetBasis(0)->GetD())->GetPtr())[0];
            m_Deriv1 = &((m_stdExp->GetBasis(1)->GetD())->GetPtr())[0];
            m_wspSize = 2 * m_nquad0*m_nquad1*m_numElmt;
        }
};

/// Factory initialisation for the PhysDeriv_SumFac_Quad operators
OperatorKey PhysDeriv_SumFac_Quad::m_type = GetOperatorFactory().
    RegisterCreatorFunction(
        OperatorKey(eQuadrilateral, ePhysDeriv, eSumFac, false),
        PhysDeriv_SumFac_Quad::create, "PhysDeriv_SumFac_Quad");


/**
 * @brief Phys deriv operator using sum-factorisation (Tri)
 */
class PhysDeriv_SumFac_Tri : public Operator
{
    public:
        OPERATOR_CREATE(PhysDeriv_SumFac_Tri)

        virtual ~PhysDeriv_SumFac_Tri()
        {
        }

        virtual void operator()(
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output0,
                      Array<OneD,       NekDouble> &output1,
                      Array<OneD,       NekDouble> &output2,
                      Array<OneD,       NekDouble> &wsp)
        {
            const int nqtot   = m_nquad0 * m_nquad1;
            const int nqcol   = nqtot*m_numElmt;

            ASSERTL1(wsp.size() == m_wspSize,
                     "Incorrect workspace size");
            ASSERTL1(input.size() >= nqcol,
                     "Incorrect input size");

            Array<OneD, NekDouble> diff0(nqcol, wsp             );
            Array<OneD, NekDouble> diff1(nqcol, wsp    +   nqcol);

            // Tensor Product Derivative
            Blas::Dgemm('N', 'N', m_nquad0, m_nquad1*m_numElmt,
                        m_nquad0, 1.0, m_Deriv0, m_nquad0,
                        input.get(), m_nquad0, 0.0,
                        diff0.get(), m_nquad0);

            int cnt = 0;
            for (int i = 0; i < m_numElmt; ++i, cnt += nqtot)
            {
                // scale diff0 by geometric factor: 2/(1-z1)
                Vmath::Vmul(nqtot,&m_fac1[0],1,diff0.get()+cnt,1,
                            diff0.get()+cnt,1);

                Blas::Dgemm('N', 'T', m_nquad0, m_nquad1, m_nquad1, 1.0,
                            input.get() + cnt, m_nquad0,
                            m_Deriv1, m_nquad1, 0.0,
                            diff1.get() + cnt, m_nquad0);

                // add to diff1 by diff0 scaled by: (1_z0)/(1-z1)
                Vmath::Vvtvp(nqtot,m_fac0.get(),1,diff0.get()+cnt,1,
                             diff1.get()+cnt,1,diff1.get()+cnt,1);
            }


            Vmath::Vmul  (nqcol, m_derivFac[0], 1, diff0, 1, output0, 1);
            Vmath::Vvtvp (nqcol, m_derivFac[1], 1, diff1, 1,
                          output0, 1, output0, 1);
            Vmath::Vmul  (nqcol, m_derivFac[2], 1, diff0, 1, output1, 1);
            Vmath::Vvtvp (nqcol, m_derivFac[3], 1, diff1, 1,
                          output1, 1, output1, 1);

            if (m_coordim == 3)
            {
                Vmath::Vmul  (nqcol, m_derivFac[4], 1, diff0, 1,
                              output2, 1);
                Vmath::Vvtvp (nqcol, m_derivFac[5], 1, diff1, 1,
                              output2, 1, output2, 1);
            }
        }

        virtual void operator()(
                      int                           dir,
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output,
                      Array<OneD,       NekDouble> &wsp)
        {
            const int nqtot   = m_nquad0 * m_nquad1;
            const int nqcol   = nqtot*m_numElmt;

            ASSERTL1(wsp.size() == m_wspSize,
                     "Incorrect workspace size");
            ASSERTL1(input.size() >= nqcol,
                     "Incorrect input size");

            Array<OneD, NekDouble> diff0(nqcol, wsp             );
            Array<OneD, NekDouble> diff1(nqcol, wsp    +   nqcol);

            // Tensor Product Derivative
            Blas::Dgemm('N', 'N', m_nquad0, m_nquad1*m_numElmt,
                        m_nquad0, 1.0, m_Deriv0, m_nquad0,
                        input.get(), m_nquad0, 0.0,
                        diff0.get(), m_nquad0);

            int cnt = 0;
            for (int i = 0; i < m_numElmt; ++i, cnt += nqtot)
            {
                // scale diff0 by geometric factor: 2/(1-z1)
                Vmath::Vmul(nqtot,&m_fac1[0],1,diff0.get()+cnt,1,
                            diff0.get()+cnt,1);

                Blas::Dgemm('N', 'T', m_nquad0, m_nquad1, m_nquad1, 1.0,
                            input.get() + cnt, m_nquad0,
                            m_Deriv1, m_nquad1, 0.0,
                            diff1.get() + cnt, m_nquad0);

                // add to diff1 by diff0 scaled by: (1_z0)/(1-z1)
                Vmath::Vvtvp(nqtot,m_fac0.get(),1,diff0.get()+cnt,1,
                             diff1.get()+cnt,1,diff1.get()+cnt,1);
            }


            Vmath::Vmul  (nqcol, m_derivFac[2*dir]  , 1, diff0, 1, output, 1);
            Vmath::Vvtvp (nqcol, m_derivFac[2*dir+1], 1, diff1, 1, output, 1,
                                                                   output, 1);
        }

    protected:
        int                             m_coordim;
        const int                       m_nquad0;
        const int                       m_nquad1;
        Array<TwoD, const NekDouble>    m_derivFac;
        NekDouble                      *m_Deriv0;
        NekDouble                      *m_Deriv1;
        Array<OneD, NekDouble>          m_fac0;
        Array<OneD, NekDouble>          m_fac1;

    private:
        PhysDeriv_SumFac_Tri(
                vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                CoalescedGeomDataSharedPtr                pGeomData)
            : Operator (pCollExp, pGeomData),
              m_nquad0 (m_stdExp->GetNumPoints(0)),
              m_nquad1 (m_stdExp->GetNumPoints(1))
        {
            LibUtilities::PointsKeyVector PtsKey = m_stdExp->GetPointsKeys();
            m_coordim = pCollExp[0]->GetCoordim();

            m_derivFac = pGeomData->GetDerivFactors(pCollExp);

            const Array<OneD, const NekDouble>& z0
                                            = m_stdExp->GetBasis(0)->GetZ();
            const Array<OneD, const NekDouble>& z1
                                            = m_stdExp->GetBasis(1)->GetZ();
            m_fac0 = Array<OneD, NekDouble>(m_nquad0*m_nquad1);
            // set up geometric factor: 0.5*(1+z0)
            for (int i = 0; i < m_nquad0; ++i)
            {
                for(int j = 0; j < m_nquad1; ++j)
                {
                    m_fac0[i+j*m_nquad0] = 0.5*(1+z0[i]);
                }
            }

            m_fac1 = Array<OneD, NekDouble>(m_nquad0*m_nquad1);
            // set up geometric factor: 2/(1-z1)
            for (int i = 0; i < m_nquad0; ++i)
            {
                for(int j = 0; j < m_nquad1; ++j)
                {
                    m_fac1[i+j*m_nquad0] = 2.0/(1-z1[j]);
                }
            }


            m_Deriv0 = &((m_stdExp->GetBasis(0)->GetD())->GetPtr())[0];
            m_Deriv1 = &((m_stdExp->GetBasis(1)->GetD())->GetPtr())[0];
            m_wspSize = 2 * m_nquad0*m_nquad1*m_numElmt;
        }
};

/// Factory initialisation for the PhysDeriv_SumFac_Tri operators
OperatorKey PhysDeriv_SumFac_Tri::m_typeArr[] =
{
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTriangle, ePhysDeriv, eSumFac,false),
        PhysDeriv_SumFac_Tri::create, "PhysDeriv_SumFac_Tri"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTriangle, ePhysDeriv, eSumFac,true),
        PhysDeriv_SumFac_Tri::create, "PhysDeriv_SumFac_NodalTri")
};


/**
 * @brief Phys deriv operator using sum-factorisation (Hex)
 */
class PhysDeriv_SumFac_Hex : public Operator
{
    public:
        OPERATOR_CREATE(PhysDeriv_SumFac_Hex)

        virtual ~PhysDeriv_SumFac_Hex()
        {
        }

        virtual void operator()(
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output0,
                      Array<OneD,       NekDouble> &output1,
                      Array<OneD,       NekDouble> &output2,
                      Array<OneD,       NekDouble> &wsp)
        {
            int nPhys = m_stdExp->GetTotPoints();
            int ntot = m_numElmt*nPhys;
            Array<OneD, NekDouble> tmp0,tmp1,tmp2;
            Array<OneD, Array<OneD, NekDouble> > Diff(3);
            Array<OneD, Array<OneD, NekDouble> > out(3);
            out[0] = output0;  out[1] = output1;    out[2] = output2;

            for(int i = 0; i < 3; ++i)
            {
                Diff[i] = wsp + i*ntot;
            }

            Blas::Dgemm('N','N', m_nquad0,m_nquad1*m_nquad2*m_numElmt,
                        m_nquad0,1.0, m_Deriv0,m_nquad0,&input[0],
                        m_nquad0,0.0,&Diff[0][0],m_nquad0);

            for(int  i = 0; i < m_numElmt; ++i)
            {
                for (int j = 0; j < m_nquad2; ++j)
                {
                    Blas::Dgemm('N', 'T', m_nquad0, m_nquad1, m_nquad1,
                                1.0, &input[i*nPhys+j*m_nquad0*m_nquad1],
                                m_nquad0, m_Deriv1, m_nquad1, 0.0,
                                &Diff[1][i*nPhys+j*m_nquad0*m_nquad1],
                                m_nquad0);
                }

                Blas::Dgemm('N','T',m_nquad0*m_nquad1,m_nquad2,m_nquad2,
                            1.0, &input[i*nPhys],m_nquad0*m_nquad1,
                            m_Deriv2,m_nquad2, 0.0,&Diff[2][i*nPhys],
                            m_nquad0*m_nquad1);
            }

            // calculate full derivative
            for(int i = 0; i < m_coordim; ++i)
            {
                Vmath::Vmul(ntot,m_derivFac[i*3],1,Diff[0],1,out[i],1);
                for(int j = 1; j < 3; ++j)
                {
                    Vmath::Vvtvp (ntot, m_derivFac[i*3+j], 1,
                                        Diff[j],               1,
                                        out[i],                1,
                                        out[i],                1);
                }
            }
        }

        virtual void operator()(
                      int                           dir,
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output,
                      Array<OneD,       NekDouble> &wsp)
        {
            int nPhys = m_stdExp->GetTotPoints();
            int ntot = m_numElmt*nPhys;
            Array<OneD, NekDouble> tmp0,tmp1,tmp2;
            Array<OneD, Array<OneD, NekDouble> > Diff(3);

            for(int i = 0; i < 3; ++i)
            {
                Diff[i] = wsp + i*ntot;
            }

            Blas::Dgemm('N','N', m_nquad0,m_nquad1*m_nquad2*m_numElmt,
                        m_nquad0,1.0, m_Deriv0,m_nquad0,&input[0],
                        m_nquad0,0.0,&Diff[0][0],m_nquad0);

            for(int  i = 0; i < m_numElmt; ++i)
            {
                for (int j = 0; j < m_nquad2; ++j)
                {
                    Blas::Dgemm('N', 'T', m_nquad0, m_nquad1, m_nquad1,
                                1.0, &input[i*nPhys+j*m_nquad0*m_nquad1],
                                m_nquad0, m_Deriv1, m_nquad1, 0.0,
                                &Diff[1][i*nPhys+j*m_nquad0*m_nquad1],
                                m_nquad0);
                }

                Blas::Dgemm('N','T',m_nquad0*m_nquad1,m_nquad2,m_nquad2,
                            1.0, &input[i*nPhys],m_nquad0*m_nquad1,
                            m_Deriv2,m_nquad2, 0.0,&Diff[2][i*nPhys],
                            m_nquad0*m_nquad1);
            }

            // calculate full derivative
            Vmath::Vmul(ntot,m_derivFac[dir*3],1,Diff[0],1,output,1);
            for(int j = 1; j < 3; ++j)
            {
                Vmath::Vvtvp (ntot, m_derivFac[dir*3+j], 1,
                                    Diff[j],               1,
                                    output,                1,
                                    output,                1);
            }
        }

    protected:
        Array<TwoD, const NekDouble>    m_derivFac;
        int                             m_coordim;
        const int                       m_nquad0;
        const int                       m_nquad1;
        const int                       m_nquad2;
        NekDouble                      *m_Deriv0;
        NekDouble                      *m_Deriv1;
        NekDouble                      *m_Deriv2;

    private:
        PhysDeriv_SumFac_Hex(
                vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                CoalescedGeomDataSharedPtr                pGeomData)
            : Operator(pCollExp, pGeomData),
              m_nquad0  (m_stdExp->GetNumPoints(0)),
              m_nquad1  (m_stdExp->GetNumPoints(1)),
              m_nquad2  (m_stdExp->GetNumPoints(2))
        {
            LibUtilities::PointsKeyVector PtsKey = m_stdExp->GetPointsKeys();

            m_coordim = pCollExp[0]->GetCoordim();

            m_derivFac = pGeomData->GetDerivFactors(pCollExp);

            m_Deriv0 = &((m_stdExp->GetBasis(0)->GetD())->GetPtr())[0];
            m_Deriv1 = &((m_stdExp->GetBasis(1)->GetD())->GetPtr())[0];
            m_Deriv2 = &((m_stdExp->GetBasis(2)->GetD())->GetPtr())[0];

            m_wspSize = 3*m_nquad0*m_nquad1*m_nquad2*m_numElmt;
        }
};

/// Factory initialisation for the PhysDeriv_SumFac_Hex operators
OperatorKey PhysDeriv_SumFac_Hex::m_typeArr[] =
{
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eHexahedron, ePhysDeriv, eSumFac, false),
        PhysDeriv_SumFac_Hex::create, "PhysDeriv_SumFac_Hex")
};


/**
 * @brief Phys deriv operator using sum-factorisation (Tet)
 */
class PhysDeriv_SumFac_Tet : public Operator
{
    public:
        OPERATOR_CREATE(PhysDeriv_SumFac_Tet)

        virtual ~PhysDeriv_SumFac_Tet()
        {
        }

        virtual void operator()(
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output0,
                      Array<OneD,       NekDouble> &output1,
                      Array<OneD,       NekDouble> &output2,
                      Array<OneD,       NekDouble> &wsp)
        {
            int nPhys = m_stdExp->GetTotPoints();
            int ntot = m_numElmt*nPhys;
            Array<OneD, NekDouble> tmp0,tmp1,tmp2;
            Array<OneD, Array<OneD, NekDouble> > Diff(3);
            Array<OneD, Array<OneD, NekDouble> > out(3);
            out[0] = output0;  out[1] = output1;    out[2] = output2;

            for(int i = 0; i < 3; ++i)
            {
                Diff[i] = wsp + i*ntot;
            }

            // dEta0
            Blas::Dgemm('N','N', m_nquad0,m_nquad1*m_nquad2*m_numElmt,
                        m_nquad0,1.0, m_Deriv0,m_nquad0,&input[0],
                        m_nquad0,0.0,&Diff[0][0],m_nquad0);

            // dEta2
            for(int  i = 0; i < m_numElmt; ++i)
            {
                Blas::Dgemm('N','T',m_nquad0*m_nquad1,m_nquad2,m_nquad2,
                            1.0, &input[i*nPhys],m_nquad0*m_nquad1,
                            m_Deriv2,m_nquad2, 0.0,&Diff[2][i*nPhys],
                            m_nquad0*m_nquad1);
            }

            for(int  i = 0; i < m_numElmt; ++i)
            {

                // dEta1
                for (int j = 0; j < m_nquad2; ++j)
                {
                    Blas::Dgemm('N', 'T', m_nquad0, m_nquad1, m_nquad1,
                                1.0, &input[i*nPhys+j*m_nquad0*m_nquad1],
                                m_nquad0, m_Deriv1, m_nquad1, 0.0,
                                &Diff[1][i*nPhys+j*m_nquad0*m_nquad1],
                                m_nquad0);
                }

                // dxi2 = (1 + eta_1)/(1 -eta_2)*dEta1 + dEta2
                Vmath::Vvtvp(nPhys, m_fac3.get(),            1,
                                    Diff[1].get() + i*nPhys, 1,
                                    Diff[2].get() + i*nPhys, 1,
                                    Diff[2].get() + i*nPhys, 1);

                // dxi1 =  2/(1 - eta_2) dEta1
                Vmath::Vmul(nPhys,  m_fac2.get(),            1,
                                    Diff[1].get() + i*nPhys, 1,
                                    Diff[1].get() + i*nPhys, 1);

                // dxi1 = 2.0(1+eta_0)/((1-eta_1)(1-eta_2)) dEta0 + dxi1
                Vmath::Vvtvp(nPhys, m_fac1.get(),            1,
                                    Diff[0].get() + i*nPhys, 1,
                                    Diff[1].get() + i*nPhys, 1,
                                    Diff[1].get() + i*nPhys, 1);

                // dxi2 = 2.0(1+eta_0)/((1-eta_1)(1-eta_2)) dEta0 + dxi2
                Vmath::Vvtvp(nPhys, m_fac1.get(),            1,
                                    Diff[0].get() + i*nPhys, 1,
                                    Diff[2].get() + i*nPhys, 1,
                                    Diff[2].get() + i*nPhys, 1);

                // dxi0 = 4.0/((1-eta_1)(1-eta_2)) dEta0
                Vmath::Vmul(nPhys,  m_fac0.get(),            1,
                                    Diff[0].get() + i*nPhys, 1,
                                    Diff[0].get() + i*nPhys, 1);

            }

            // calculate full derivative
            for(int i = 0; i < m_coordim; ++i)
            {
                Vmath::Vmul(ntot,m_derivFac[i*3],1,Diff[0],1,out[i],1);
                for(int j = 1; j < 3; ++j)
                {
                    Vmath::Vvtvp (ntot, m_derivFac[i*3+j], 1,
                                        Diff[j], 1, out[i], 1, out[i], 1);
                }
            }
        }

        virtual void operator()(
                      int                           dir,
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output,
                      Array<OneD,       NekDouble> &wsp)
        {
            int nPhys = m_stdExp->GetTotPoints();
            int ntot = m_numElmt*nPhys;
            Array<OneD, NekDouble> tmp0,tmp1,tmp2;
            Array<OneD, Array<OneD, NekDouble> > Diff(3);

            for(int i = 0; i < 3; ++i)
            {
                Diff[i] = wsp + i*ntot;
            }

            // dEta0
            Blas::Dgemm('N','N', m_nquad0,m_nquad1*m_nquad2*m_numElmt,
                        m_nquad0,1.0, m_Deriv0,m_nquad0,&input[0],
                        m_nquad0,0.0,&Diff[0][0],m_nquad0);

            // dEta2
            for(int  i = 0; i < m_numElmt; ++i)
            {
                Blas::Dgemm('N','T',m_nquad0*m_nquad1,m_nquad2,m_nquad2,
                            1.0, &input[i*nPhys],m_nquad0*m_nquad1,
                            m_Deriv2,m_nquad2, 0.0,&Diff[2][i*nPhys],
                            m_nquad0*m_nquad1);
            }

            for(int  i = 0; i < m_numElmt; ++i)
            {

                // dEta1
                for (int j = 0; j < m_nquad2; ++j)
                {
                    Blas::Dgemm('N', 'T', m_nquad0, m_nquad1, m_nquad1,
                                1.0, &input[i*nPhys+j*m_nquad0*m_nquad1],
                                m_nquad0, m_Deriv1, m_nquad1, 0.0,
                                &Diff[1][i*nPhys+j*m_nquad0*m_nquad1],
                                m_nquad0);
                }

                // dxi2 = (1 + eta_1)/(1 -eta_2)*dEta1 + dEta2
                Vmath::Vvtvp(nPhys, m_fac3.get(),            1,
                                    Diff[1].get() + i*nPhys, 1,
                                    Diff[2].get() + i*nPhys, 1,
                                    Diff[2].get() + i*nPhys, 1);

                // dxi1 =  2/(1 - eta_2) dEta1
                Vmath::Vmul(nPhys,  m_fac2.get(),            1,
                                    Diff[1].get() + i*nPhys, 1,
                                    Diff[1].get() + i*nPhys, 1);

                // dxi1 = 2.0(1+eta_0)/((1-eta_1)(1-eta_2)) dEta0 + dxi1
                Vmath::Vvtvp(nPhys, m_fac1.get(),            1,
                                    Diff[0].get() + i*nPhys, 1,
                                    Diff[1].get() + i*nPhys, 1,
                                    Diff[1].get() + i*nPhys, 1);

                // dxi2 = 2.0(1+eta_0)/((1-eta_1)(1-eta_2)) dEta0 + dxi2
                Vmath::Vvtvp(nPhys, m_fac1.get(),            1,
                                    Diff[0].get() + i*nPhys, 1,
                                    Diff[2].get() + i*nPhys, 1,
                                    Diff[2].get() + i*nPhys, 1);

                // dxi0 = 4.0/((1-eta_1)(1-eta_2)) dEta0
                Vmath::Vmul(nPhys,  m_fac0.get(),            1,
                                    Diff[0].get() + i*nPhys, 1,
                                    Diff[0].get() + i*nPhys, 1);

            }

            // calculate full derivative
            Vmath::Vmul(ntot,m_derivFac[dir*3],1,Diff[0],1,output,1);
            for(int j = 1; j < 3; ++j)
            {
                Vmath::Vvtvp (ntot, m_derivFac[dir*3+j], 1,
                                    Diff[j], 1, output, 1, output, 1);
            }
        }

    protected:
        Array<TwoD, const NekDouble>    m_derivFac;
        int                             m_coordim;
        const int                       m_nquad0;
        const int                       m_nquad1;
        const int                       m_nquad2;
        NekDouble                      *m_Deriv0;
        NekDouble                      *m_Deriv1;
        NekDouble                      *m_Deriv2;
        Array<OneD, NekDouble>          m_fac0;
        Array<OneD, NekDouble>          m_fac1;
        Array<OneD, NekDouble>          m_fac2;
        Array<OneD, NekDouble>          m_fac3;

    private:
        PhysDeriv_SumFac_Tet(
                vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                CoalescedGeomDataSharedPtr                pGeomData)
            : Operator(pCollExp, pGeomData),
              m_nquad0  (m_stdExp->GetNumPoints(0)),
              m_nquad1  (m_stdExp->GetNumPoints(1)),
              m_nquad2  (m_stdExp->GetNumPoints(2))
        {
            LibUtilities::PointsKeyVector PtsKey = m_stdExp->GetPointsKeys();

            m_coordim = pCollExp[0]->GetCoordim();

            m_derivFac = pGeomData->GetDerivFactors(pCollExp);

            m_Deriv0 = &((m_stdExp->GetBasis(0)->GetD())->GetPtr())[0];
            m_Deriv1 = &((m_stdExp->GetBasis(1)->GetD())->GetPtr())[0];
            m_Deriv2 = &((m_stdExp->GetBasis(2)->GetD())->GetPtr())[0];

            m_wspSize = 3*m_nquad0*m_nquad1*m_nquad2*m_numElmt;

            const Array<OneD, const NekDouble>& z0
                                            = m_stdExp->GetBasis(0)->GetZ();
            const Array<OneD, const NekDouble>& z1
                                            = m_stdExp->GetBasis(1)->GetZ();
            const Array<OneD, const NekDouble>& z2
                                            = m_stdExp->GetBasis(2)->GetZ();

            m_fac0 = Array<OneD, NekDouble>(m_nquad0*m_nquad1*m_nquad2);
            m_fac1 = Array<OneD, NekDouble>(m_nquad0*m_nquad1*m_nquad2);
            m_fac2 = Array<OneD, NekDouble>(m_nquad0*m_nquad1*m_nquad2);
            m_fac3 = Array<OneD, NekDouble>(m_nquad0*m_nquad1*m_nquad2);
            // calculate 2.0/((1-eta_1)(1-eta_2))
            for (int i = 0; i < m_nquad0; ++i)
            {
                for(int j = 0; j < m_nquad1; ++j)
                {
                    for(int k = 0; k < m_nquad2; ++k)
                    {

                        m_fac0[i + j*m_nquad0 + k*m_nquad0*m_nquad1]
                               = 4.0/((1-z1[j])*(1-z2[k]));
                        m_fac1[i + j*m_nquad0 + k*m_nquad0*m_nquad1]
                               = 2.0*(1+z0[i])/((1-z1[j])*(1-z2[k]));
                        m_fac2[i + j*m_nquad0 + k*m_nquad0*m_nquad1]
                               = 2.0/(1-z2[k]);
                        m_fac3[i + j*m_nquad0 + k*m_nquad0*m_nquad1]
                               = (1+z1[j])/(1-z2[k]);
                    }
                }
            }

        }
};

/// Factory initialisation for the PhysDeriv_SumFac_Tet operators
OperatorKey PhysDeriv_SumFac_Tet::m_typeArr[] =
{
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTetrahedron, ePhysDeriv, eSumFac, false),
        PhysDeriv_SumFac_Tet::create, "PhysDeriv_SumFac_Tet")
};


/**
 * @brief Phys deriv operator using sum-factorisation (Prism)
 */
class PhysDeriv_SumFac_Prism : public Operator
{
    public:
        OPERATOR_CREATE(PhysDeriv_SumFac_Prism)

        virtual ~PhysDeriv_SumFac_Prism()
        {
        }

        virtual void operator()(
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output0,
                      Array<OneD,       NekDouble> &output1,
                      Array<OneD,       NekDouble> &output2,
                      Array<OneD,       NekDouble> &wsp)
        {
            int nPhys = m_stdExp->GetTotPoints();
            int ntot = m_numElmt*nPhys;
            Array<OneD, NekDouble> tmp0,tmp1,tmp2;
            Array<OneD, Array<OneD, NekDouble> > Diff(3);
            Array<OneD, Array<OneD, NekDouble> > out(3);
            out[0] = output0; out[1] = output1; out[2] = output2;

            for(int i = 0; i < 3; ++i)
            {
                Diff[i] = wsp + i*ntot;
            }

            // dEta0
            Blas::Dgemm('N','N', m_nquad0,m_nquad1*m_nquad2*m_numElmt,
                        m_nquad0,1.0, m_Deriv0,m_nquad0,&input[0],
                        m_nquad0,0.0,&Diff[0][0],m_nquad0);

            int cnt = 0;
            for(int  i = 0; i < m_numElmt; ++i)
            {

                // dEta 1
                for (int j = 0; j < m_nquad2; ++j)
                {
                    Blas::Dgemm('N', 'T', m_nquad0, m_nquad1, m_nquad1,
                                1.0, &input[i*nPhys+j*m_nquad0*m_nquad1],
                                m_nquad0, m_Deriv1, m_nquad1, 0.0,
                                &Diff[1][i*nPhys+j*m_nquad0*m_nquad1],
                                m_nquad0);
                }

                // dEta 2
                Blas::Dgemm('N','T',m_nquad0*m_nquad1,m_nquad2,m_nquad2,
                            1.0, &input[i*nPhys],m_nquad0*m_nquad1,
                            m_Deriv2,m_nquad2, 0.0,&Diff[2][i*nPhys],
                            m_nquad0*m_nquad1);

                // dxi0 = 2/(1-eta_2) d Eta_0
                Vmath::Vmul(nPhys,&m_fac0[0],1,Diff[0].get()+cnt,1,
                            Diff[0].get()+cnt,1);

                // dxi2 = (1+eta0)/(1-eta_2) d Eta_0 + d/dEta2;
                Vmath::Vvtvp(nPhys,&m_fac1[0],1,Diff[0].get()+cnt,1,
                             Diff[2].get()+cnt,1,Diff[2].get()+cnt,1);
                cnt += nPhys;
            }

            // calculate full derivative
            for(int i = 0; i < m_coordim; ++i)
            {
                Vmath::Vmul(ntot,m_derivFac[i*3],1,Diff[0],1,out[i],1);
                for(int j = 1; j < 3; ++j)
                {
                    Vmath::Vvtvp (ntot, m_derivFac[i*3+j], 1,
                                        Diff[j], 1, out[i], 1, out[i], 1);
                }
            }
        }

        virtual void operator()(
                      int                           dir,
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output,
                      Array<OneD,       NekDouble> &wsp)
        {
            int nPhys = m_stdExp->GetTotPoints();
            int ntot = m_numElmt*nPhys;
            Array<OneD, NekDouble> tmp0,tmp1,tmp2;
            Array<OneD, Array<OneD, NekDouble> > Diff(3);

            for(int i = 0; i < 3; ++i)
            {
                Diff[i] = wsp + i*ntot;
            }

            // dEta0
            Blas::Dgemm('N','N', m_nquad0,m_nquad1*m_nquad2*m_numElmt,
                        m_nquad0,1.0, m_Deriv0,m_nquad0,&input[0],
                        m_nquad0,0.0,&Diff[0][0],m_nquad0);

            int cnt = 0;
            for(int  i = 0; i < m_numElmt; ++i)
            {

                // dEta 1
                for (int j = 0; j < m_nquad2; ++j)
                {
                    Blas::Dgemm('N', 'T', m_nquad0, m_nquad1, m_nquad1,
                                1.0, &input[i*nPhys+j*m_nquad0*m_nquad1],
                                m_nquad0, m_Deriv1, m_nquad1, 0.0,
                                &Diff[1][i*nPhys+j*m_nquad0*m_nquad1],
                                m_nquad0);
                }

                // dEta 2
                Blas::Dgemm('N','T',m_nquad0*m_nquad1,m_nquad2,m_nquad2,
                            1.0, &input[i*nPhys],m_nquad0*m_nquad1,
                            m_Deriv2,m_nquad2, 0.0,&Diff[2][i*nPhys],
                            m_nquad0*m_nquad1);

                // dxi0 = 2/(1-eta_2) d Eta_0
                Vmath::Vmul(nPhys,&m_fac0[0],1,Diff[0].get()+cnt,1,
                            Diff[0].get()+cnt,1);

                // dxi2 = (1+eta0)/(1-eta_2) d Eta_0 + d/dEta2;
                Vmath::Vvtvp(nPhys,&m_fac1[0],1,Diff[0].get()+cnt,1,
                             Diff[2].get()+cnt,1,Diff[2].get()+cnt,1);
                cnt += nPhys;
            }

            // calculate full derivative
            Vmath::Vmul(ntot,m_derivFac[dir*3],1,Diff[0],1,output,1);
            for(int j = 1; j < 3; ++j)
            {
                Vmath::Vvtvp (ntot, m_derivFac[dir*3+j], 1,
                                    Diff[j], 1, output, 1, output, 1);
            }
        }

    protected:
        Array<TwoD, const NekDouble>    m_derivFac;
        int                             m_coordim;
        const int                       m_nquad0;
        const int                       m_nquad1;
        const int                       m_nquad2;
        NekDouble                      *m_Deriv0;
        NekDouble                      *m_Deriv1;
        NekDouble                      *m_Deriv2;
        Array<OneD, NekDouble>          m_fac0;
        Array<OneD, NekDouble>          m_fac1;

    private:
        PhysDeriv_SumFac_Prism(
                vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                CoalescedGeomDataSharedPtr                pGeomData)
            : Operator(pCollExp, pGeomData),
              m_nquad0  (m_stdExp->GetNumPoints(0)),
              m_nquad1  (m_stdExp->GetNumPoints(1)),
              m_nquad2  (m_stdExp->GetNumPoints(2))
        {
            LibUtilities::PointsKeyVector PtsKey = m_stdExp->GetPointsKeys();

            m_coordim = pCollExp[0]->GetCoordim();

            m_derivFac = pGeomData->GetDerivFactors(pCollExp);

            const Array<OneD, const NekDouble>& z0
                                            = m_stdExp->GetBasis(0)->GetZ();
            const Array<OneD, const NekDouble>& z2
                                            = m_stdExp->GetBasis(2)->GetZ();
            m_fac0 = Array<OneD, NekDouble>(m_nquad0*m_nquad1*m_nquad2);
            m_fac1 = Array<OneD, NekDouble>(m_nquad0*m_nquad1*m_nquad2);
            for (int i = 0; i < m_nquad0; ++i)
            {
                for(int j = 0; j < m_nquad1; ++j)
                {
                    for(int k = 0; k < m_nquad2; ++k)
                    {
                        m_fac0[i+j*m_nquad0 + k*m_nquad0*m_nquad1] =
                            2.0/(1-z2[k]);
                        m_fac1[i+j*m_nquad0 + k*m_nquad0*m_nquad1] =
                            0.5*(1+z0[i]);
                    }
                }
            }



            m_Deriv0 = &((m_stdExp->GetBasis(0)->GetD())->GetPtr())[0];
            m_Deriv1 = &((m_stdExp->GetBasis(1)->GetD())->GetPtr())[0];
            m_Deriv2 = &((m_stdExp->GetBasis(2)->GetD())->GetPtr())[0];

            m_wspSize = 3*m_nquad0*m_nquad1*m_nquad2*m_numElmt;
        }
};

/// Factory initialisation for the PhysDeriv_SumFac_Prism operators
OperatorKey PhysDeriv_SumFac_Prism::m_typeArr[] = {
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePrism, ePhysDeriv, eSumFac, false),
        PhysDeriv_SumFac_Prism::create, "PhysDeriv_SumFac_Prism")
};


/**
 * @brief Phys deriv operator using sum-factorisation (Pyramid)
 */
class PhysDeriv_SumFac_Pyr : public Operator
{
    public:
        OPERATOR_CREATE(PhysDeriv_SumFac_Pyr)

        virtual ~PhysDeriv_SumFac_Pyr()
        {
        }

        virtual void operator()(
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output0,
                      Array<OneD,       NekDouble> &output1,
                      Array<OneD,       NekDouble> &output2,
                      Array<OneD,       NekDouble> &wsp)
        {
            int nPhys = m_stdExp->GetTotPoints();
            int ntot = m_numElmt*nPhys;
            Array<OneD, NekDouble> tmp0,tmp1,tmp2;
            Array<OneD, Array<OneD, NekDouble> > Diff(3);
            Array<OneD, Array<OneD, NekDouble> > out(3);
            out[0] = output0; out[1] = output1; out[2] = output2;

            for(int i = 0; i < 3; ++i)
            {
                Diff[i] = wsp + i*ntot;
            }

            // dEta0
            Blas::Dgemm('N','N', m_nquad0,m_nquad1*m_nquad2*m_numElmt,
                        m_nquad0,1.0, m_Deriv0,m_nquad0,&input[0],
                        m_nquad0,0.0,&Diff[0][0],m_nquad0);

            int cnt = 0;
            for(int  i = 0; i < m_numElmt; ++i)
            {

                // dEta 1
                for (int j = 0; j < m_nquad2; ++j)
                {
                    Blas::Dgemm('N', 'T', m_nquad0, m_nquad1, m_nquad1,
                                1.0, &input[i*nPhys+j*m_nquad0*m_nquad1],
                                m_nquad0, m_Deriv1, m_nquad1, 0.0,
                                &Diff[1][i*nPhys+j*m_nquad0*m_nquad1],
                                m_nquad0);
                }

                // dEta 2
                Blas::Dgemm('N','T',m_nquad0*m_nquad1,m_nquad2,m_nquad2,
                            1.0, &input[i*nPhys],m_nquad0*m_nquad1,
                            m_Deriv2,m_nquad2, 0.0,&Diff[2][i*nPhys],
                            m_nquad0*m_nquad1);

                // dxi0 = 2/(1-eta_2) d Eta_0
                Vmath::Vmul(nPhys,&m_fac0[0],1,Diff[0].get()+cnt,1,
                            Diff[0].get()+cnt,1);

                // dxi1 = 2/(1-eta_2) d Eta_1
                Vmath::Vmul(nPhys,&m_fac0[0],1,Diff[1].get()+cnt,1,
                            Diff[1].get()+cnt,1);

                // dxi2 = (1+eta0)/(1-eta_2) d Eta_0 + d/dEta2;
                Vmath::Vvtvp(nPhys,&m_fac1[0],1,Diff[0].get()+cnt,1,
                             Diff[2].get()+cnt,1,Diff[2].get()+cnt,1);
                // dxi2 += (1+eta1)/(1-eta_2) d Eta_1
                Vmath::Vvtvp(nPhys,&m_fac2[0],1,Diff[1].get()+cnt,1,
                             Diff[2].get()+cnt,1,Diff[2].get()+cnt,1);
                cnt += nPhys;
            }

            // calculate full derivative
            for(int i = 0; i < m_coordim; ++i)
            {
                Vmath::Vmul(ntot,m_derivFac[i*3],1,Diff[0],1,out[i],1);
                for(int j = 1; j < 3; ++j)
                {
                    Vmath::Vvtvp (ntot, m_derivFac[i*3+j], 1,
                                        Diff[j], 1, out[i], 1, out[i], 1);
                }
            }
        }

        virtual void operator()(
                      int                           dir,
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output,
                      Array<OneD,       NekDouble> &wsp)
        {
            int nPhys = m_stdExp->GetTotPoints();
            int ntot = m_numElmt*nPhys;
            Array<OneD, NekDouble> tmp0,tmp1,tmp2;
            Array<OneD, Array<OneD, NekDouble> > Diff(3);

            for(int i = 0; i < 3; ++i)
            {
                Diff[i] = wsp + i*ntot;
            }

            // dEta0
            Blas::Dgemm('N','N', m_nquad0,m_nquad1*m_nquad2*m_numElmt,
                        m_nquad0,1.0, m_Deriv0,m_nquad0,&input[0],
                        m_nquad0,0.0,&Diff[0][0],m_nquad0);

            int cnt = 0;
            for(int  i = 0; i < m_numElmt; ++i)
            {
                // dEta 1
                for (int j = 0; j < m_nquad2; ++j)
                {
                    Blas::Dgemm('N', 'T', m_nquad0, m_nquad1, m_nquad1,
                                1.0, &input[i*nPhys+j*m_nquad0*m_nquad1],
                                m_nquad0, m_Deriv1, m_nquad1, 0.0,
                                &Diff[1][i*nPhys+j*m_nquad0*m_nquad1],
                                m_nquad0);
                }

                // dEta 2
                Blas::Dgemm('N','T',m_nquad0*m_nquad1,m_nquad2,m_nquad2,
                            1.0, &input[i*nPhys],m_nquad0*m_nquad1,
                            m_Deriv2,m_nquad2, 0.0,&Diff[2][i*nPhys],
                            m_nquad0*m_nquad1);

                // dxi0 = 2/(1-eta_2) d Eta_0
                Vmath::Vmul(nPhys,&m_fac0[0],1,Diff[0].get()+cnt,1,
                            Diff[0].get()+cnt,1);

                // dxi1 = 2/(1-eta_2) d Eta_1
                Vmath::Vmul(nPhys,&m_fac0[0],1,Diff[1].get()+cnt,1,
                            Diff[1].get()+cnt,1);

                // dxi2 = (1+eta0)/(1-eta_2) d Eta_0 + d/dEta2;
                Vmath::Vvtvp(nPhys,&m_fac1[0],1,Diff[0].get()+cnt,1,
                             Diff[2].get()+cnt,1,Diff[2].get()+cnt,1);
                // dxi2 = (1+eta1)/(1-eta_2) d Eta_1 + d/dEta2;
                Vmath::Vvtvp(nPhys,&m_fac2[0],1,Diff[1].get()+cnt,1,
                             Diff[2].get()+cnt,1,Diff[2].get()+cnt,1);
                cnt += nPhys;
            }

            // calculate full derivative
            Vmath::Vmul(ntot,m_derivFac[dir*3],1,Diff[0],1,output,1);
            for(int j = 1; j < 3; ++j)
            {
                Vmath::Vvtvp (ntot, m_derivFac[dir*3+j], 1,
                                    Diff[j], 1, output, 1, output, 1);
            }
        }

    protected:
        Array<TwoD, const NekDouble>    m_derivFac;
        int                             m_coordim;
        const int                       m_nquad0;
        const int                       m_nquad1;
        const int                       m_nquad2;
        NekDouble                      *m_Deriv0;
        NekDouble                      *m_Deriv1;
        NekDouble                      *m_Deriv2;
        Array<OneD, NekDouble>          m_fac0;
        Array<OneD, NekDouble>          m_fac1;
        Array<OneD, NekDouble>          m_fac2;

    private:
        PhysDeriv_SumFac_Pyr(
                vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                CoalescedGeomDataSharedPtr                pGeomData)
            : Operator(pCollExp, pGeomData),
              m_nquad0  (m_stdExp->GetNumPoints(0)),
              m_nquad1  (m_stdExp->GetNumPoints(1)),
              m_nquad2  (m_stdExp->GetNumPoints(2))
        {
            LibUtilities::PointsKeyVector PtsKey = m_stdExp->GetPointsKeys();

            m_coordim = pCollExp[0]->GetCoordim();

            m_derivFac = pGeomData->GetDerivFactors(pCollExp);

            const Array<OneD, const NekDouble>& z0
                                            = m_stdExp->GetBasis(0)->GetZ();
            const Array<OneD, const NekDouble>& z1
                                            = m_stdExp->GetBasis(1)->GetZ();
            const Array<OneD, const NekDouble>& z2
                                            = m_stdExp->GetBasis(2)->GetZ();
            m_fac0 = Array<OneD, NekDouble>(m_nquad0*m_nquad1*m_nquad2);
            m_fac1 = Array<OneD, NekDouble>(m_nquad0*m_nquad1*m_nquad2);
            m_fac2 = Array<OneD, NekDouble>(m_nquad0*m_nquad1*m_nquad2);

            int nq0_nq1 = m_nquad0*m_nquad1;
            for (int i = 0; i < m_nquad0; ++i)
            {
                for(int j = 0; j < m_nquad1; ++j)
                {
                    int ifac = i+j*m_nquad0;
                    for(int k = 0; k < m_nquad2; ++k)
                    {
                        m_fac0[ifac + k*nq0_nq1] =
                            2.0/(1-z2[k]);
                        m_fac1[ifac + k*nq0_nq1] =
                            0.5*(1+z0[i]);
                        m_fac2[ifac + k*nq0_nq1] =
                            0.5*(1+z1[j]);
                    }
                }
            }

            m_Deriv0 = &((m_stdExp->GetBasis(0)->GetD())->GetPtr())[0];
            m_Deriv1 = &((m_stdExp->GetBasis(1)->GetD())->GetPtr())[0];
            m_Deriv2 = &((m_stdExp->GetBasis(2)->GetD())->GetPtr())[0];

            m_wspSize = 3*m_nquad0*m_nquad1*m_nquad2*m_numElmt;
        }
};

/// Factory initialisation for the PhysDeriv_SumFac_Pyr operators
OperatorKey PhysDeriv_SumFac_Pyr::m_typeArr[] = {
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePyramid, ePhysDeriv, eSumFac, false),
        PhysDeriv_SumFac_Pyr::create, "PhysDeriv_SumFac_Pyr")
};

}
}
