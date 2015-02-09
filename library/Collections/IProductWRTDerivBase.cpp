///////////////////////////////////////////////////////////////////////////////
//
// File: IProductWRTDerivBase.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
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
// Description: IProductWRTDerivBase operator implementations
//
///////////////////////////////////////////////////////////////////////////////

#include <loki/Singleton.h>
#include <Collections/Operator.h>
#include <Collections/Collection.h>
#include <Collections/IProduct.h>

namespace Nektar {
namespace Collections {

class IProductWRTDerivBase_IterPerExp : public Operator
{
public:
    IProductWRTDerivBase_IterPerExp(vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                         CoalescedGeomDataSharedPtr GeomData)
        : Operator(pCollExp, GeomData)
    {
        LibUtilities::PointsKeyVector PtsKey = m_stdExp->GetPointsKeys();
        m_dim = PtsKey.size();
        m_coordim = m_stdExp->GetCoordim();

        int nqtot  = m_stdExp->GetTotPoints();

        m_derivFac = GeomData->GetDerivFactors(pCollExp);
        m_jac = GeomData->GetJac(pCollExp);
        m_wspSize = m_dim*nqtot*m_numElmt;
    }

    virtual void operator()(const Array<OneD, const NekDouble> &entry0,
                            Array<OneD, NekDouble> &entry1,
                            Array<OneD, NekDouble> &entry2,
                            Array<OneD, NekDouble> &entry3,
                            Array<OneD, NekDouble> &wsp)
    {
        unsigned int nPhys  = m_stdExp->GetTotPoints();
        unsigned int ntot   = m_numElmt*nPhys;
        unsigned int nmodes = m_stdExp->GetNcoeffs();
        unsigned int nmax   = max(ntot,m_numElmt*nmodes);
        Array<OneD, Array<OneD, const NekDouble> > in(3);
        Array<OneD, NekDouble> output, tmp1;
        Array<OneD, Array<OneD, NekDouble> > tmp(3);

        in[0] = entry0; in[1] = entry1; in[2] = entry2;

        output = (m_coordim == 3)? entry3: (m_coordim == 2)?
            entry2: entry1;

        for(int i = 0; i < m_dim; ++i)
        {
            tmp[i] = wsp + i*nmax;
        }

        // calculate dx/dxi in[0] + dy/dxi in[2] + dz/dxi in[3]
        for(int i = 0; i < m_dim; ++i)
        {
            Vmath::Vmul (ntot,m_derivFac[i],1, in[0],1,
                         tmp[i],1);
            for(int j = 1; j < m_coordim; ++j)
            {
                Vmath::Vvtvp (ntot,m_derivFac[i +j*m_dim],1,
                              in[j],1, tmp[i], 1, tmp[i],1);
            }
        }

        // calculate Iproduct WRT Std Deriv
        // first component
        Vmath::Vmul(ntot,m_jac,1,tmp[0],1,tmp[0],1);
        for(int n = 0; n < m_numElmt; ++n)
        {
            m_stdExp->IProductWRTDerivBase(0,tmp[0]+n*nPhys,
                                           tmp1 = output + n*nmodes);
        }

        // other components
        for(int i = 1; i < m_dim; ++i)
        {
            // multiply by Jacobian
            Vmath::Vmul(ntot,m_jac,1,tmp[i],1,tmp[i],1);
            for(int n = 0; n < m_numElmt; ++n)
            {
                m_stdExp->IProductWRTDerivBase(i,tmp[i]+n*nPhys,tmp[0]);
                Vmath::Vadd(nmodes,tmp[0],1,output+n*nmodes,1,
                            tmp1 = output+n*nmodes,1);
            }
        }
    }

    OPERATOR_CREATE(IProductWRTDerivBase_IterPerExp)

    Array<TwoD, const NekDouble>  m_derivFac;
    Array<OneD, const NekDouble> m_jac;
    int m_dim;
    int m_coordim;
};

OperatorKey IProductWRTDerivBase_IterPerExp::m_typeArr[] = {
    GetOperatorFactory().RegisterCreatorFunction(
          OperatorKey(LibUtilities::eSegment,
                      eIProductWRTDerivBase, eIterPerExp,false),
          IProductWRTDerivBase_IterPerExp::create,
          "IProductWRTDerivBase_IterPerExp_Seg"),
    GetOperatorFactory().RegisterCreatorFunction(
          OperatorKey(LibUtilities::eTriangle,
                      eIProductWRTDerivBase, eIterPerExp,false),
          IProductWRTDerivBase_IterPerExp::create,
          "IProductWRTDerivBase_IterPerExp_Tri"),
    GetOperatorFactory().RegisterCreatorFunction(
          OperatorKey(LibUtilities::eTriangle,
                      eIProductWRTDerivBase, eIterPerExp,true),
          IProductWRTDerivBase_IterPerExp::create,
          "IProductWRTDerivBase_IterPerExp_NodalTri"),
    GetOperatorFactory().RegisterCreatorFunction(
          OperatorKey(LibUtilities::eQuadrilateral,
                      eIProductWRTDerivBase, eIterPerExp,false),
          IProductWRTDerivBase_IterPerExp::create,
          "IProductWRTDerivBase_IterPerExp_Quad"),
    GetOperatorFactory().RegisterCreatorFunction(
          OperatorKey(LibUtilities::eTetrahedron,
                      eIProductWRTDerivBase, eIterPerExp,false),
          IProductWRTDerivBase_IterPerExp::create,
          "IProductWRTDerivBase_IterPerExp_Tet"),
    GetOperatorFactory().RegisterCreatorFunction(
          OperatorKey(LibUtilities::eTetrahedron,
                      eIProductWRTDerivBase, eIterPerExp,true),
          IProductWRTDerivBase_IterPerExp::create,
          "IProductWRTDerivBase_IterPerExp_NodalTet"),
    GetOperatorFactory().RegisterCreatorFunction(
          OperatorKey(LibUtilities::ePyramid,
                      eIProductWRTDerivBase, eIterPerExp,false),
          IProductWRTDerivBase_IterPerExp::create,
          "IProductWRTDerivBase_IterPerExp_Pyr"),
    GetOperatorFactory().RegisterCreatorFunction(
          OperatorKey(LibUtilities::ePrism,
                      eIProductWRTDerivBase, eIterPerExp,false),
          IProductWRTDerivBase_IterPerExp::create,
          "IProductWRTDerivBase_IterPerExp_Prism"),
    GetOperatorFactory().RegisterCreatorFunction(
          OperatorKey(LibUtilities::ePrism,
                      eIProductWRTDerivBase, eIterPerExp,true),
          IProductWRTDerivBase_IterPerExp::create,
          "IProductWRTDerivBase_IterPerExp_NodalPrism"),
    GetOperatorFactory().RegisterCreatorFunction(
          OperatorKey(LibUtilities::eHexahedron,
                      eIProductWRTDerivBase, eIterPerExp,false),
          IProductWRTDerivBase_IterPerExp::create,
          "IProductWRTDerivBase_IterPerExp_Hex")
};

class IProductWRTDerivBase_StdMat : public Operator
{
public:
    IProductWRTDerivBase_StdMat(vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                                 CoalescedGeomDataSharedPtr GeomData)
        : Operator(pCollExp,GeomData)
    {
        LibUtilities::PointsKeyVector PtsKey = m_stdExp->GetPointsKeys();
        m_dim = PtsKey.size();
        m_coordim = m_stdExp->GetCoordim();

        int nqtot  = m_stdExp->GetTotPoints();
        int nmodes = m_stdExp->GetNcoeffs();

        // set up a IProductWRTDerivBase StdMat.
        m_iProdWRTStdDBase = Array<OneD, DNekMatSharedPtr>(m_dim);
        for(int i = 0; i < m_dim; ++i)
        {
            Array<OneD, NekDouble> tmp(nqtot),tmp1(nmodes);
            m_iProdWRTStdDBase[i] = MemoryManager<DNekMat>::AllocateSharedPtr(nmodes,nqtot);
            for(int j = 0; j < nqtot; ++j)
            {
                Vmath::Zero(nqtot,tmp,1);
                tmp[j] = 1.0;
                m_stdExp->IProductWRTDerivBase(i,tmp,tmp1);
                Vmath::Vcopy(nmodes,&tmp1[0],1,
                             &(m_iProdWRTStdDBase[i]->GetPtr())[0]+j*nmodes,1);
            }
        }
        m_derivFac = GeomData->GetDerivFactors(pCollExp);
        m_jac      = GeomData->GetJac(pCollExp);
        m_wspSize = m_dim*nqtot*m_numElmt;
    }

    virtual void operator()(const Array<OneD, const NekDouble> &entry0,
                            Array<OneD,       NekDouble> &entry1,
                            Array<OneD,       NekDouble> &entry2,
                            Array<OneD,       NekDouble> &entry3,
                            Array<OneD,       NekDouble> &wsp)
    {
        int nPhys = m_stdExp->GetTotPoints();
        int ntot = m_numElmt*nPhys;
        int nmodes = m_stdExp->GetNcoeffs();
        Array<OneD, Array<OneD, const NekDouble> > in(3);
        Array<OneD, NekDouble> output;
        Array<OneD, Array<OneD, NekDouble> > tmp(3);

        in[0] = entry0; in[1] = entry1;
        in[2] = entry2;

        output = (m_coordim == 3)? entry3: (m_coordim == 2)?
            entry2: entry1;

        for(int i = 0; i < m_dim; ++i)
        {
            tmp[i] = wsp + i*ntot;
        }

        // calculate dx/dxi in[0] + dy/dxi in[1] + dz/dxi in[2]
        for(int i = 0; i < m_dim; ++i)
        {
            Vmath::Vmul (ntot,m_derivFac[i],1, in[0],1,
                         tmp[i],1);
            for(int j = 1; j < m_coordim; ++j)
            {
                Vmath::Vvtvp (ntot,m_derivFac[i +j*m_dim],1,
                              in[j],1, tmp[i], 1, tmp[i],1);
            }
        }

        // calculate Iproduct WRT Std Deriv

        // First component
        Vmath::Vmul(ntot,m_jac,1,tmp[0],1,tmp[0],1);
        Blas::Dgemm('N', 'N', m_iProdWRTStdDBase[0]->GetRows(),
                    m_numElmt,m_iProdWRTStdDBase[0]->GetColumns(),
                    1.0, m_iProdWRTStdDBase[0]->GetRawPtr(),
                    m_iProdWRTStdDBase[0]->GetRows(),
                    tmp[0].get(), nPhys, 0.0,
                    output.get(), nmodes);

        // Other components
        for(int i = 1; i < m_dim; ++i)
        {
            Vmath::Vmul(ntot,m_jac,1,tmp[i],1,tmp[i],1);
            Blas::Dgemm('N', 'N', m_iProdWRTStdDBase[i]->GetRows(),
                        m_numElmt,m_iProdWRTStdDBase[i]->GetColumns(),
                        1.0, m_iProdWRTStdDBase[i]->GetRawPtr(),
                        m_iProdWRTStdDBase[i]->GetRows(),
                        tmp[i].get(), nPhys, 1.0,
                        output.get(), nmodes);
        }
    }

    OPERATOR_CREATE(IProductWRTDerivBase_StdMat)

    Array<OneD, DNekMatSharedPtr> m_iProdWRTStdDBase;
    Array<TwoD, const NekDouble>  m_derivFac;
    Array<OneD, const NekDouble>  m_jac;
    int m_dim;
    int m_coordim;
};

OperatorKey IProductWRTDerivBase_StdMat::m_typeArr[] = {
    GetOperatorFactory().RegisterCreatorFunction(
          OperatorKey(LibUtilities::eSegment,
                      eIProductWRTDerivBase, eStdMat,false),
          IProductWRTDerivBase_StdMat::create,
          "IProductWRTDerivBase_StdMat_Seg"),
    GetOperatorFactory().RegisterCreatorFunction(
          OperatorKey(LibUtilities::eTriangle,
                      eIProductWRTDerivBase, eStdMat,false),
          IProductWRTDerivBase_StdMat::create,
          "IProductWRTDerivBase_StdMat_Tri"),
    GetOperatorFactory().RegisterCreatorFunction(
          OperatorKey(LibUtilities::eTriangle,
                      eIProductWRTDerivBase, eStdMat,true),
          IProductWRTDerivBase_StdMat::create,
          "IProductWRTDerivBase_StdMat_NodalTri"),
    GetOperatorFactory().RegisterCreatorFunction(
          OperatorKey(LibUtilities::eTriangle,
                      eIProductWRTDerivBase, eSumFac,true),
          IProductWRTDerivBase_StdMat::create,
          "IProductWRTDerivBase_SumFac_NodalTri"),
    GetOperatorFactory().RegisterCreatorFunction(
          OperatorKey(LibUtilities::eQuadrilateral,
                      eIProductWRTDerivBase, eStdMat,false),
          IProductWRTDerivBase_StdMat::create,
          "IProductWRTDerivBase_StdMat_Quad"),
    GetOperatorFactory().RegisterCreatorFunction(
          OperatorKey(LibUtilities::eTetrahedron,
                      eIProductWRTDerivBase, eStdMat,false),
          IProductWRTDerivBase_StdMat::create,
          "IProductWRTDerivBase_StdMat_Tet"),
    GetOperatorFactory().RegisterCreatorFunction(
          OperatorKey(LibUtilities::eTetrahedron,
                      eIProductWRTDerivBase, eStdMat,true),
          IProductWRTDerivBase_StdMat::create,
          "IProductWRTDerivBase_StdMat_NodalTet"),
    GetOperatorFactory().RegisterCreatorFunction(
          OperatorKey(LibUtilities::eTetrahedron,
                      eIProductWRTDerivBase, eSumFac,true),
          IProductWRTDerivBase_StdMat::create,
          "IProductWRTDerivBase_SumFac_NodalTet"),
    GetOperatorFactory().RegisterCreatorFunction(
          OperatorKey(LibUtilities::ePyramid,
                      eIProductWRTDerivBase, eStdMat,false),
          IProductWRTDerivBase_StdMat::create,
          "IProductWRTDerivBase_StdMat_Pyr"),
    GetOperatorFactory().RegisterCreatorFunction(
          OperatorKey(LibUtilities::ePyramid,
                      eIProductWRTDerivBase, eSumFac,false),
          IProductWRTDerivBase_StdMat::create,
          "IProductWRTDerivBase_SumFac_Pyr"),
    GetOperatorFactory().RegisterCreatorFunction(
          OperatorKey(LibUtilities::ePrism,
                      eIProductWRTDerivBase, eStdMat,false),
          IProductWRTDerivBase_StdMat::create,
          "IProductWRTDerivBase_StdMat_Prism"),
    GetOperatorFactory().RegisterCreatorFunction(
          OperatorKey(LibUtilities::ePrism,
                      eIProductWRTDerivBase, eStdMat,true),
          IProductWRTDerivBase_StdMat::create,
          "IProductWRTDerivBase_StdMat_NodalPrism"),
    GetOperatorFactory().RegisterCreatorFunction(
          OperatorKey(LibUtilities::ePrism,
                      eIProductWRTDerivBase, eSumFac,true),
          IProductWRTDerivBase_StdMat::create,
          "IProductWRTDerivBase_SumFac_NodalPrism"),
    GetOperatorFactory().RegisterCreatorFunction(
          OperatorKey(LibUtilities::eHexahedron,
                      eIProductWRTDerivBase, eStdMat,false),
          IProductWRTDerivBase_StdMat::create,
          "IProductWRTDerivBase_StdMat_Hex")
};

class IProductWRTDerivBase_NoCollection : public Operator
{
public:
    IProductWRTDerivBase_NoCollection(
            vector<StdRegions::StdExpansionSharedPtr> pCollExp,
            CoalescedGeomDataSharedPtr GeomData)
        : Operator(pCollExp,GeomData)
    {
        m_expList = pCollExp;
        m_dim = pCollExp[0]->GetNumBases();
        m_coordim = pCollExp[0]->GetCoordim();
    }

    virtual void operator()(const Array<OneD, const NekDouble> &entry0,
                            Array<OneD, NekDouble> &entry1,
                            Array<OneD, NekDouble> &entry2,
                            Array<OneD, NekDouble> &entry3,
                            Array<OneD, NekDouble> &wsp)
    {
        unsigned int nmodes = m_expList[0]->GetNcoeffs();
        unsigned int nPhys  = m_expList[0]->GetTotPoints();
        Array<OneD, NekDouble> tmp(nmodes),tmp1;

        Array<OneD, Array<OneD, const NekDouble> > in(3);
        Array<OneD, NekDouble> output;
        in[0] = entry0; in[1] = entry1; in[2] = entry2;

        output = (m_coordim == 3)? entry3: (m_coordim == 2)?
            entry2: entry1;

        for(int n = 0; n < m_numElmt; ++n)
        {
            m_expList[n]->IProductWRTDerivBase(0, in[0] + n * nPhys,
                                                  tmp1 = output + n * nmodes);
        }

        for(int i = 1; i < m_dim; ++i)
        {
            for(int n = 0; n < m_numElmt; ++n)
            {
                m_expList[n]->IProductWRTDerivBase(i,in[i]+n*nPhys,tmp);

                Vmath::Vadd(nmodes,tmp,1,output+n*nmodes,1,
                            tmp1 = output+n*nmodes,1);
            }
        }
    }

    OPERATOR_CREATE(IProductWRTDerivBase_NoCollection)

    int m_dim;
    int m_coordim;
    vector<StdRegions::StdExpansionSharedPtr> m_expList;
};

OperatorKey IProductWRTDerivBase_NoCollection::m_typeArr[] = {
    GetOperatorFactory().RegisterCreatorFunction(
          OperatorKey(LibUtilities::eSegment,
                      eIProductWRTDerivBase, eNoCollection,false),
          IProductWRTDerivBase_NoCollection::create,
          "IProductWRTDerivBase_NoCollection_Seg"),
    GetOperatorFactory().RegisterCreatorFunction(
          OperatorKey(LibUtilities::eTriangle,
                      eIProductWRTDerivBase, eNoCollection,false),
          IProductWRTDerivBase_NoCollection::create,
          "IProductWRTDerivBase_NoCollection_Tri"),
    GetOperatorFactory().RegisterCreatorFunction(
          OperatorKey(LibUtilities::eTriangle,
                      eIProductWRTDerivBase, eNoCollection,true),
          IProductWRTDerivBase_NoCollection::create,
          "IProductWRTDerivBase_NoCollection_NodalTri"),
    GetOperatorFactory().RegisterCreatorFunction(
          OperatorKey(LibUtilities::eQuadrilateral,
                      eIProductWRTDerivBase, eNoCollection,false),
          IProductWRTDerivBase_NoCollection::create,
          "IProductWRTDerivBase_NoCollection_Quad"),
    GetOperatorFactory().RegisterCreatorFunction(
          OperatorKey(LibUtilities::eTetrahedron,
                      eIProductWRTDerivBase, eNoCollection,false),
          IProductWRTDerivBase_NoCollection::create,
          "IProductWRTDerivBase_NoCollection_Tet"),
    GetOperatorFactory().RegisterCreatorFunction(
          OperatorKey(LibUtilities::eTetrahedron,
                      eIProductWRTDerivBase, eNoCollection,true),
          IProductWRTDerivBase_NoCollection::create,
          "IProductWRTDerivBase_NoCollection_NodalTet"),
    GetOperatorFactory().RegisterCreatorFunction(
          OperatorKey(LibUtilities::ePyramid,
                      eIProductWRTDerivBase, eNoCollection,false),
          IProductWRTDerivBase_NoCollection::create,
          "IProductWRTDerivBase_NoCollection_Pyr"),
    GetOperatorFactory().RegisterCreatorFunction(
          OperatorKey(LibUtilities::ePrism,
                      eIProductWRTDerivBase, eNoCollection,false),
          IProductWRTDerivBase_NoCollection::create,
          "IProductWRTDerivBase_NoCollection_Prism"),
    GetOperatorFactory().RegisterCreatorFunction(
          OperatorKey(LibUtilities::ePrism,
                      eIProductWRTDerivBase, eNoCollection,true),
          IProductWRTDerivBase_NoCollection::create,
          "IProductWRTDerivBase_NoCollection_NodalPrism"),
    GetOperatorFactory().RegisterCreatorFunction(
          OperatorKey(LibUtilities::eHexahedron,
                      eIProductWRTDerivBase, eNoCollection,false),
          IProductWRTDerivBase_NoCollection::create,
          "IProductWRTDerivBase_NoCollection_Hex")
};


/*
 * ----------------------------------------------------------
 * IProductWRTDerivBase operators
 * ----------------------------------------------------------
 */
class IProductWRTDerivBase_SumFac_Seg : public Operator
{
public:
    IProductWRTDerivBase_SumFac_Seg(vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                                    CoalescedGeomDataSharedPtr GeomData)
        : Operator  (pCollExp, GeomData),
          m_nquad0  (m_stdExp->GetNumPoints(0)),
          m_nmodes0 (m_stdExp->GetBasisNumModes(0)),
          m_derbase0(m_stdExp->GetBasis(0)->GetDbdata())
    {
        m_wspSize = m_numElmt*m_nquad0;
        m_derivFac = GeomData->GetDerivFactors(pCollExp);
        m_jac = GeomData->GetJacWithStdWeights(pCollExp);
    }

    virtual void operator()(const Array<OneD, const NekDouble> &input,
                            Array<OneD,       NekDouble> &output,
                            Array<OneD,       NekDouble> &output1,
                            Array<OneD,       NekDouble> &output2,
                            Array<OneD,       NekDouble> &wsp)
    {

        Vmath::Vmul(m_numElmt*m_nquad0,m_jac,1,input,1,wsp,1);
        Vmath::Vmul(m_numElmt*m_nquad0,&m_derivFac[0][0],1,&wsp[0],1,&wsp[0],1);

        // out = B0*in;
        Blas::Dgemm('T','N', m_nmodes0,m_numElmt,m_nquad0,1.0, m_derbase0.get(), m_nquad0,
                    &wsp[0], m_nquad0, 0.0,&output[0], m_nmodes0);
    }

    OPERATOR_CREATE(IProductWRTDerivBase_SumFac_Seg)

    protected:
    const int  m_nquad0;
    const int  m_nmodes0;
    Array<OneD, const NekDouble> m_jac;
    Array<OneD, const NekDouble> m_derbase0;
    Array<TwoD, const NekDouble> m_derivFac;
};

OperatorKey IProductWRTDerivBase_SumFac_Seg::m_type = GetOperatorFactory().
    RegisterCreatorFunction(OperatorKey(LibUtilities::eSegment, eIProductWRTDerivBase,
                                        eSumFac,false),
                            IProductWRTDerivBase_SumFac_Seg::create,
                            "IProductWRTDerivBase_SumFac_Seg");


/*
 * ----------------------------------------------------------
 * IProductWRTDerivBase operator
 * ----------------------------------------------------------
 */
class IProductWRTDerivBase_SumFac_Quad : public Operator
{
public:
    IProductWRTDerivBase_SumFac_Quad(vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                                     CoalescedGeomDataSharedPtr GeomData)
        : Operator(pCollExp, GeomData),
          m_nquad0  (m_stdExp->GetNumPoints(0)),
          m_nquad1  (m_stdExp->GetNumPoints(1)),
          m_nmodes0 (m_stdExp->GetBasisNumModes(0)),
          m_nmodes1 (m_stdExp->GetBasisNumModes(1)),
          m_colldir0(m_stdExp->GetBasis(0)->Collocation()),
          m_colldir1(m_stdExp->GetBasis(1)->Collocation()),
          m_base0   (m_stdExp->GetBasis(0)->GetBdata()),
          m_base1   (m_stdExp->GetBasis(1)->GetBdata()),
          m_derbase0   (m_stdExp->GetBasis(0)->GetDbdata()),
          m_derbase1   (m_stdExp->GetBasis(1)->GetDbdata())
    {
        LibUtilities::PointsKeyVector PtsKey = m_stdExp->GetPointsKeys();
        m_coordim = m_stdExp->GetCoordim();

        m_derivFac = GeomData->GetDerivFactors(pCollExp);
        m_jac      = GeomData->GetJacWithStdWeights(pCollExp);
        m_wspSize  = 4*m_numElmt*(max(m_nquad0*m_nquad1,m_nmodes0*m_nmodes1));
    }

    virtual void operator()(const Array<OneD, const NekDouble> &entry0,
                            Array<OneD, NekDouble> &entry1,
                            Array<OneD, NekDouble> &entry2,
                            Array<OneD, NekDouble> &entry3,
                            Array<OneD, NekDouble> &wsp)
    {
        unsigned int nPhys  = m_stdExp->GetTotPoints();
        unsigned int ntot   = m_numElmt*nPhys;
        unsigned int nmodes = m_stdExp->GetNcoeffs();
        unsigned int nmax   = max(ntot,m_numElmt*nmodes);
        Array<OneD, Array<OneD, const NekDouble> > in(3);
        Array<OneD, NekDouble> output, wsp1;
        Array<OneD, Array<OneD, NekDouble> > tmp(2);

        in[0] = entry0; in[1] = entry1; in[2] = entry2;

        output = (m_coordim == 2)? entry2: entry3;

        tmp[0] = wsp; tmp[1] = wsp + nmax;
        wsp1   = wsp + 2*nmax;

        // calculate dx/dxi in[0] + dy/dxi in[1]
        for(int i = 0; i < 2; ++i)
        {
            Vmath::Vmul (ntot,m_derivFac[i],1, in[0],1,
                         tmp[i],1);
            for(int j = 1; j < m_coordim; ++j)
            {
                Vmath::Vvtvp (ntot,m_derivFac[i +j*2],1,
                              in[j],1, tmp[i], 1, tmp[i],1);
            }
        }

        // Iproduct wrt derivative of base 0
        QuadIProduct(false, m_colldir1,m_numElmt,
                     m_nquad0,   m_nquad1,
                     m_nmodes0,  m_nmodes1,
                     m_derbase0, m_base1,
                     m_jac, tmp[0], output, wsp1);

        // Iproduct wrt derivative of base 1
        QuadIProduct(m_colldir0, false, m_numElmt,
                     m_nquad0,   m_nquad1,
                     m_nmodes0,  m_nmodes1,
                     m_base0, m_derbase1,
                     m_jac, tmp[1],  tmp[0], wsp1);

        Vmath::Vadd(m_numElmt*nmodes,tmp[0],1,output,1,output,1);
    }

    OPERATOR_CREATE(IProductWRTDerivBase_SumFac_Quad)

    const int  m_nquad0;
    const int  m_nquad1;
    const int  m_nmodes0;
    const int  m_nmodes1;
    const bool m_colldir0;
    const bool m_colldir1;
    int m_coordim;
    Array<TwoD, const NekDouble> m_derivFac;
    Array<OneD, const NekDouble> m_jac;
    Array<OneD, const NekDouble> m_base0;
    Array<OneD, const NekDouble> m_base1;
    Array<OneD, const NekDouble> m_derbase0;
    Array<OneD, const NekDouble> m_derbase1;
};

OperatorKey IProductWRTDerivBase_SumFac_Quad::m_type =
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(LibUtilities::eQuadrilateral,
                    eIProductWRTDerivBase, eSumFac,false),
        IProductWRTDerivBase_SumFac_Quad::create,
        "IProductWRTDerivBase_IterPerExp_Quad");


/*
 * ----------------------------------------------------------
 * IProductWRTDerivBase operator
 * ----------------------------------------------------------
 */
class IProductWRTDerivBase_SumFac_Tri : public Operator
{
public:
    IProductWRTDerivBase_SumFac_Tri(vector<StdRegions::StdExpansionSharedPtr> pCollExp,
           CoalescedGeomDataSharedPtr GeomData)
        : Operator(pCollExp, GeomData),
          m_nquad0  (m_stdExp->GetNumPoints(0)),
          m_nquad1  (m_stdExp->GetNumPoints(1)),
          m_nmodes0 (m_stdExp->GetBasisNumModes(0)),
          m_nmodes1 (m_stdExp->GetBasisNumModes(1)),
          m_colldir0(m_stdExp->GetBasis(0)->Collocation()),
          m_colldir1(m_stdExp->GetBasis(1)->Collocation()),
          m_base0   (m_stdExp->GetBasis(0)->GetBdata()),
          m_base1   (m_stdExp->GetBasis(1)->GetBdata()),
          m_derbase0   (m_stdExp->GetBasis(0)->GetDbdata()),
          m_derbase1   (m_stdExp->GetBasis(1)->GetDbdata())
    {
        LibUtilities::PointsKeyVector PtsKey = m_stdExp->GetPointsKeys();
        m_coordim = m_stdExp->GetCoordim();

        m_derivFac = GeomData->GetDerivFactors(pCollExp);
        m_jac      = GeomData->GetJacWithStdWeights(pCollExp);
        m_wspSize  = 4*m_numElmt*(max(m_nquad0*m_nquad1,m_nmodes0*m_nmodes1));

        if(m_stdExp->GetBasis(0)->GetBasisType() == LibUtilities::eModified_A)
        {
            m_sortTopVertex = true;
        }
        else
        {
            m_sortTopVertex = false;
        }

        const Array<OneD, const NekDouble>& z0 = m_stdExp->GetBasis(0)->GetZ();
        const Array<OneD, const NekDouble>& z1 = m_stdExp->GetBasis(1)->GetZ();

        m_fac0 = Array<OneD, NekDouble>(m_nquad0*m_nquad1);
        // set up geometric factor: 2/(1-z1)
        for (int i = 0; i < m_nquad0; ++i)
        {
            for(int j = 0; j < m_nquad1; ++j)
            {
                m_fac0[i+j*m_nquad0] = 2.0/(1-z1[j]);
            }
        }

        m_fac1 = Array<OneD, NekDouble>(m_nquad0*m_nquad1);
        // set up geometric factor: (1+z0)/(1-z1)
        for (int i = 0; i < m_nquad0; ++i)
        {
            for(int j = 0; j < m_nquad1; ++j)
            {
                m_fac1[i+j*m_nquad0] = (1+z0[i])/(1-z1[j]);
            }
        }
    }

    virtual void operator()(const Array<OneD, const NekDouble> &entry0,
                            Array<OneD, NekDouble> &entry1,
                            Array<OneD, NekDouble> &entry2,
                            Array<OneD, NekDouble> &entry3,
                            Array<OneD, NekDouble> &wsp)
    {
        unsigned int nPhys  = m_stdExp->GetTotPoints();
        unsigned int ntot   = m_numElmt*nPhys;
        unsigned int nmodes = m_stdExp->GetNcoeffs();
        unsigned int nmax   = max(ntot,m_numElmt*nmodes);
        Array<OneD, Array<OneD, const NekDouble> > in(3);
        Array<OneD, NekDouble> output, wsp1;
        Array<OneD, Array<OneD, NekDouble> > tmp(2);

        in[0] = entry0; in[1] = entry1; in[2] = entry2;

        output = (m_coordim == 2)? entry2: entry3;

        tmp[0] = wsp; tmp[1] = wsp + nmax;
        wsp1   = wsp + 2*nmax;


        // calculate (dphi/dx,in[0]) = ((dphi/dxi_0 dxi_0/dx + dphi/dxi_1 dxi_1/dx),in[0])
        //     +     (dphi/dy,in[1]) = ((dphi/dxi_0 dxi_0/dy + dphi/dxi_1 dxi_1/dy),in[1])
        //
        // Note dphi/dxi_0  =
        //             dphi/deta_0 deta_0/dxi_0 = dphi/deta_0 2/(1-eta_1)
        //
        //      dphi/dxi_1  =
        //             dphi/deta_1 deta_1/dxi_1 + dphi/deta_1 deta_1/dxi_1 =
        //             dphi/deta_0 (1+eta_0)/(1-eta_1) + dphi/deta_1
        //
        // and so the full inner products are
        //
        // (dphi/dx,in[0]) + (dphi/dy,in[1])
        //    = (dphi/deta_0, ((2/(1-eta_1) (dxi_0/dx in[0] + dxi_0/dy in[1])
        //            + (1_eta_0)/(1-eta_1) (dxi_1/dx in[0] + dxi_1/dy in[1]))
        //    + (dphi/deta_1, (dxi_1/dx in[0] + dxi_1/dy in[1]))

        for(int i = 0; i < 2; ++i)
        {
            Vmath::Vmul (ntot,m_derivFac[i],1, in[0],1, tmp[i],1);

            for(int j = 1; j < m_coordim; ++j)
            {
                Vmath::Vvtvp (ntot,m_derivFac[i +j*2],1,
                              in[j],1, tmp[i], 1, tmp[i],1);
            }
        }

        // Multiply by factor: 2/(1-z1)
        for (int i = 0; i < m_numElmt; ++i)
        {
            // scale tmp[0] by geometric factor: 2/(1-z1)
            Vmath::Vmul(nPhys,&m_fac0[0],1,tmp[0].get()+i*nPhys,1,
                        tmp[0].get()+i*nPhys,1);

            // scale tmp[1] by geometric factor (1+z0)/(1-z1)
            Vmath::Vvtvp(nPhys,&m_fac1[0],1,tmp[1].get()+i*nPhys,1,
                         tmp[0].get()+i*nPhys,1,tmp[0].get()+i*nPhys,1);
        }

        // Iproduct wrt derivative of base 0
        TriIProduct(m_sortTopVertex, m_numElmt, m_nquad0, m_nquad1,
                    m_nmodes0,  m_nmodes1, m_derbase0, m_base1,
                    m_jac, tmp[0], output, wsp1);

        // Iproduct wrt derivative of base 1
        TriIProduct(m_sortTopVertex, m_numElmt, m_nquad0, m_nquad1,
                    m_nmodes0,  m_nmodes1, m_base0, m_derbase1,
                    m_jac, tmp[1], tmp[0], wsp1);

        Vmath::Vadd(m_numElmt*nmodes,tmp[0],1,output,1,output,1);
    }

    OPERATOR_CREATE(IProductWRTDerivBase_SumFac_Tri)

    const int  m_nquad0;
    const int  m_nquad1;
    const int  m_nmodes0;
    const int  m_nmodes1;
    const bool m_colldir0;
    const bool m_colldir1;
    int m_coordim;
    Array<TwoD, const NekDouble> m_derivFac;
    Array<OneD, const NekDouble> m_jac;
    Array<OneD, const NekDouble> m_base0;
    Array<OneD, const NekDouble> m_base1;
    Array<OneD, const NekDouble> m_derbase0;
    Array<OneD, const NekDouble> m_derbase1;
    Array<OneD, NekDouble> m_fac0;
    Array<OneD, NekDouble> m_fac1;
    bool m_sortTopVertex;
};

OperatorKey IProductWRTDerivBase_SumFac_Tri::m_type =
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(LibUtilities::eTriangle,
                    eIProductWRTDerivBase, eSumFac,false),
        IProductWRTDerivBase_SumFac_Tri::create,
        "IProductWRTDerivBase_IterPerExp_Tri");

/*
 * ----------------------------------------------------------
 * IProductWRTDerivBase operators
 * ----------------------------------------------------------
 */
class IProductWRTDerivBase_SumFac_Hex : public Operator
{
public:
    IProductWRTDerivBase_SumFac_Hex(vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                               CoalescedGeomDataSharedPtr GeomData)
        : Operator  (pCollExp, GeomData),
          m_nquad0  (m_stdExp->GetNumPoints(0)),
          m_nquad1  (m_stdExp->GetNumPoints(1)),
          m_nquad2  (m_stdExp->GetNumPoints(2)),
          m_nmodes0 (m_stdExp->GetBasisNumModes(0)),
          m_nmodes1 (m_stdExp->GetBasisNumModes(1)),
          m_nmodes2 (m_stdExp->GetBasisNumModes(2)),
          m_colldir0(m_stdExp->GetBasis(0)->Collocation()),
          m_colldir1(m_stdExp->GetBasis(1)->Collocation()),
          m_colldir2(m_stdExp->GetBasis(2)->Collocation()),
          m_base0    (m_stdExp->GetBasis(0)->GetBdata()),
          m_base1    (m_stdExp->GetBasis(1)->GetBdata()),
          m_base2    (m_stdExp->GetBasis(2)->GetBdata()),
          m_derbase0    (m_stdExp->GetBasis(0)->GetDbdata()),
          m_derbase1    (m_stdExp->GetBasis(1)->GetDbdata()),
          m_derbase2    (m_stdExp->GetBasis(2)->GetDbdata())

    {
        m_jac     = GeomData->GetJacWithStdWeights(pCollExp);
        m_wspSize = 6*m_numElmt*(max(m_nquad0*m_nquad1*m_nquad2,m_nmodes0*m_nmodes1*m_nmodes2));
        m_derivFac = GeomData->GetDerivFactors(pCollExp);
    }

    virtual void operator()(const Array<OneD, const NekDouble> &entry0,
                            Array<OneD, NekDouble> &entry1,
                            Array<OneD, NekDouble> &entry2,
                            Array<OneD, NekDouble> &entry3,
                            Array<OneD, NekDouble> &wsp)
    {
        unsigned int nPhys  = m_stdExp->GetTotPoints();
        unsigned int ntot   = m_numElmt*nPhys;
        unsigned int nmodes = m_stdExp->GetNcoeffs();
        unsigned int nmax  = max(ntot,m_numElmt*nmodes);
        Array<OneD, Array<OneD, const NekDouble> > in(3);
        Array<OneD, NekDouble> output, wsp1;
        Array<OneD, Array<OneD, NekDouble> > tmp(3);

        in[0] = entry0; in[1] = entry1;
        in[2] = entry2;

        output =  entry3;

        for(int i = 0; i < 3; ++i)
        {
            tmp[i] = wsp + i*nmax;
        }

        // calculate dx/dxi in[0] + dy/dxi in[1] + dz/dxi in[2]
        for(int i = 0; i < 3; ++i)
        {
            Vmath::Vmul (ntot,m_derivFac[i],1, in[0],1,
                         tmp[i],1);
            for(int j = 1; j < 3; ++j)
            {
                Vmath::Vvtvp (ntot,m_derivFac[i+3*j],1,
                              in[j],1, tmp[i], 1, tmp[i],1);
            }
        }

        wsp1   = wsp + 3*nmax;

        // calculate Iproduct WRT Std Deriv
        HexIProduct(false,m_colldir1,m_colldir2, m_numElmt,
                    m_nquad0,   m_nquad1,  m_nquad2,
                    m_nmodes0,  m_nmodes1, m_nmodes2,
                    m_derbase0, m_base1,   m_base2,
                    m_jac,tmp[0],output,wsp1);

        HexIProduct(m_colldir0,false,m_colldir2, m_numElmt,
                    m_nquad0,  m_nquad1,   m_nquad2,
                    m_nmodes0, m_nmodes1,  m_nmodes2,
                    m_base0,   m_derbase1, m_base2,
                    m_jac,tmp[1],tmp[0],wsp1);
        Vmath::Vadd(m_numElmt*nmodes,tmp[0],1,output,1,output,1);

        HexIProduct(m_colldir0,m_colldir1,false, m_numElmt,
                    m_nquad0,  m_nquad1,  m_nquad2,
                    m_nmodes0, m_nmodes1, m_nmodes2,
                    m_base0,   m_base1,   m_derbase2,
                    m_jac,tmp[2],tmp[0],wsp1);
        Vmath::Vadd(m_numElmt*nmodes,tmp[0],1,output,1,output,1);
    }

    OPERATOR_CREATE(IProductWRTDerivBase_SumFac_Hex)

    protected:
    const int  m_nquad0;
    const int  m_nquad1;
    const int  m_nquad2;
    const int  m_nmodes0;
    const int  m_nmodes1;
    const int  m_nmodes2;
    const bool m_colldir0;
    const bool m_colldir1;
    const bool m_colldir2;
    Array<OneD, const NekDouble> m_jac;
    Array<OneD, const NekDouble> m_base0;
    Array<OneD, const NekDouble> m_base1;
    Array<OneD, const NekDouble> m_base2;
    Array<OneD, const NekDouble> m_derbase0;
    Array<OneD, const NekDouble> m_derbase1;
    Array<OneD, const NekDouble> m_derbase2;
    Array<TwoD, const NekDouble> m_derivFac;
};

OperatorKey IProductWRTDerivBase_SumFac_Hex::m_type = GetOperatorFactory().
    RegisterCreatorFunction(OperatorKey(LibUtilities::eHexahedron,
                                        eIProductWRTDerivBase, eSumFac,false),
                            IProductWRTDerivBase_SumFac_Hex::create,
                            "IProductWRTDerivBase_SumFac_Hex");


/*
 * ----------------------------------------------------------
 * IProductWRTDerivBase operators
 * ----------------------------------------------------------
 */
class IProductWRTDerivBase_SumFac_Tet : public Operator
{
public:
    IProductWRTDerivBase_SumFac_Tet(vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                                    CoalescedGeomDataSharedPtr GeomData)
        : Operator  (pCollExp, GeomData),
          m_nquad0  (m_stdExp->GetNumPoints(0)),
          m_nquad1  (m_stdExp->GetNumPoints(1)),
          m_nquad2  (m_stdExp->GetNumPoints(2)),
          m_nmodes0 (m_stdExp->GetBasisNumModes(0)),
          m_nmodes1 (m_stdExp->GetBasisNumModes(1)),
          m_nmodes2 (m_stdExp->GetBasisNumModes(2)),
          m_base0   (m_stdExp->GetBasis(0)->GetBdata()),
          m_base1   (m_stdExp->GetBasis(1)->GetBdata()),
          m_base2   (m_stdExp->GetBasis(2)->GetBdata()),
          m_derbase0(m_stdExp->GetBasis(0)->GetDbdata()),
          m_derbase1(m_stdExp->GetBasis(1)->GetDbdata()),
          m_derbase2(m_stdExp->GetBasis(2)->GetDbdata())

    {
        m_jac     = GeomData->GetJacWithStdWeights(pCollExp);
        m_wspSize = 6*m_numElmt*(max(m_nquad0*m_nquad1*m_nquad2,m_nmodes0*m_nmodes1*m_nmodes2));
        m_derivFac = GeomData->GetDerivFactors(pCollExp);


        const Array<OneD, const NekDouble>& z0 = m_stdExp->GetBasis(0)->GetZ();
        const Array<OneD, const NekDouble>& z1 = m_stdExp->GetBasis(1)->GetZ();
        const Array<OneD, const NekDouble>& z2 = m_stdExp->GetBasis(2)->GetZ();

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

                    m_fac0[i+j*m_nquad0+k*m_nquad0*m_nquad1] = 4.0/((1-z1[j])*(1-z2[k]));
                    m_fac1[i+j*m_nquad0+k*m_nquad0*m_nquad1] = (1+z0[i])*0.5;
                    m_fac2[i+j*m_nquad0+k*m_nquad0*m_nquad1] = 2.0/(1-z2[k]);
                    m_fac3[i+j*m_nquad0+k*m_nquad0*m_nquad1] = (1+z1[j])*0.5;
                }
            }
        }

        if(m_stdExp->GetBasis(0)->GetBasisType() == LibUtilities::eModified_A)
        {
            m_sortTopEdge = true;
        }
        else
        {
            m_sortTopEdge = false;
        }

    }

    virtual void operator()(const Array<OneD, const NekDouble> &entry0,
                            Array<OneD, NekDouble> &entry1,
                            Array<OneD, NekDouble> &entry2,
                            Array<OneD, NekDouble> &entry3,
                            Array<OneD, NekDouble> &wsp)
    {
        unsigned int nPhys  = m_stdExp->GetTotPoints();
        unsigned int ntot   = m_numElmt*nPhys;
        unsigned int nmodes = m_stdExp->GetNcoeffs();
        unsigned int nmax  = max(ntot,m_numElmt*nmodes);
        Array<OneD, Array<OneD, const NekDouble> > in(3);
        Array<OneD, NekDouble> output, wsp1;
        Array<OneD, Array<OneD, NekDouble> > tmp(3);

        in[0] = entry0; in[1] = entry1;
        in[2] = entry2;

        output =  entry3;

        for(int i = 0; i < 3; ++i)
        {
            tmp[i] = wsp + i*nmax;
        }


        // calculate (dphi/dx,in[0]) = ((dphi/dxi_0 dxi_0/dx + dphi/dxi_1 dxi_1/dx + dphi/dxi_2 dxi_2/dx),in[0])
        //     +     (dphi/dy,in[1]) = ((dphi/dxi_0 dxi_0/dy + dphi/dxi_1 dxi_1/dy + dphi/dxi_2 dxi_2/dy),in[1])
        //     +     (dphi/dz,in[2]) = ((dphi/dxi_0 dxi_0/dz + dphi/dxi_1 dxi_1/dz + dphi/dxi_2 dxi_2/dz),in[1])
        //
        // Note dphi/dxi_0  =
        //             dphi/deta_0 4/((1-eta_1)(1-eta2))
        //
        //      dphi/dxi_1  =
        //             dphi/deta_0 2(1+eta_0)/((1-eta_1)(1-eta_2)) + dphi/deta_1 2/(1-eta_2)
        //
        //      dphi/dxi_2  =
        //             dphi/deta_0 2(1+eta_0)/((1-eta_1)(1-eta_2)) + dphi/deta_1 (1+eta_1)/(1-eta_2)  + dphi/deta_2
        //
        // and so the full inner products are
        //
        // (dphi/dx,in[0]) + (dphi/dy,in[1]) + (dphi/dz,in[2])
        //    = (dphi/deta_0, fac0 (tmp0 + fac1(tmp1 + tmp2)))
        //    + (dphi/deta_1, fac2 (tmp1 + fac3 tmp2))
        //    + (dphi/deta_2, tmp2)
        //
        // tmp0 = (dxi_0/dx in[0] + dxi_0/dy in[1] + dxi_0/dz in[2])
        // tmp1 = (dxi_1/dx in[0] + dxi_1/dy in[1] + dxi_1/dz in[2])
        // tmp2 = (dxi_2/dx in[0] + dxi_2/dy in[1] + dxi_2/dz in[2])

        // fac0 = 4/((1-eta_1)(1-eta2))
        // fac1 = (1+eta_0)/2
        // fac2 = 2/(1-eta_2)
        // fac3 = (1+eta_1)/2

        for(int i = 0; i < 3; ++i)
        {
            Vmath::Vmul (ntot,m_derivFac[i],1, in[0],1,tmp[i],1);
            for(int j = 1; j < 3; ++j)
            {
                Vmath::Vvtvp (ntot,m_derivFac[i+3*j],1,
                              in[j],1, tmp[i], 1, tmp[i],1);
            }
        }

        wsp1   = wsp + 3*nmax;

        // Sort into eta factors
        for (int i = 0; i < m_numElmt; ++i)
        {
            // add tmp[1] + tmp[2]
            Vmath::Vadd(nPhys,tmp[1].get()+i*nPhys,1,tmp[2].get()+i*nPhys,1,wsp1.get(),1);

            // scale wsp1 by fac1 and add to tmp0
            Vmath::Vvtvp(nPhys,&m_fac1[0],1,wsp1.get(),1,
                         tmp[0].get()+i*nPhys,1,tmp[0].get()+i*nPhys,1);

            // scale tmp[0] by fac0
            Vmath::Vmul(nPhys,&m_fac0[0],1,tmp[0].get()+i*nPhys,1,
                        tmp[0].get()+i*nPhys,1);

            // scale tmp[2] by fac3 and add to tmp1
            Vmath::Vvtvp(nPhys,&m_fac3[0],1,tmp[2].get()+i*nPhys,1,
                         tmp[1].get()+i*nPhys,1,tmp[1].get()+i*nPhys,1);

            // scale tmp[1] by fac2
            Vmath::Vmul(nPhys,&m_fac2[0],1,tmp[1].get()+i*nPhys,1,
                        tmp[1].get()+i*nPhys,1);
        }


        // calculate Iproduct WRT Std Deriv
        TetIProduct(m_sortTopEdge, m_numElmt,
                    m_nquad0,   m_nquad1,  m_nquad2,
                    m_nmodes0,  m_nmodes1, m_nmodes2,
                    m_derbase0, m_base1,   m_base2,
                    m_jac,tmp[0],output,wsp1);

        TetIProduct(m_sortTopEdge, m_numElmt,
                    m_nquad0,  m_nquad1,   m_nquad2,
                    m_nmodes0, m_nmodes1,  m_nmodes2,
                    m_base0,   m_derbase1, m_base2,
                    m_jac,tmp[1],tmp[0],wsp1);
        Vmath::Vadd(m_numElmt*nmodes,tmp[0],1,output,1,output,1);

        TetIProduct(m_sortTopEdge, m_numElmt,
                    m_nquad0,  m_nquad1,  m_nquad2,
                    m_nmodes0, m_nmodes1, m_nmodes2,
                    m_base0,   m_base1,   m_derbase2,
                    m_jac,tmp[2],tmp[0],wsp1);
        Vmath::Vadd(m_numElmt*nmodes,tmp[0],1,output,1,output,1);
    }

    OPERATOR_CREATE(IProductWRTDerivBase_SumFac_Tet)

    protected:
    const int  m_nquad0;
    const int  m_nquad1;
    const int  m_nquad2;
    const int  m_nmodes0;
    const int  m_nmodes1;
    const int  m_nmodes2;
    Array<OneD, const NekDouble> m_jac;
    Array<OneD, const NekDouble> m_base0;
    Array<OneD, const NekDouble> m_base1;
    Array<OneD, const NekDouble> m_base2;
    Array<OneD, const NekDouble> m_derbase0;
    Array<OneD, const NekDouble> m_derbase1;
    Array<OneD, const NekDouble> m_derbase2;
    Array<TwoD, const NekDouble> m_derivFac;
    Array<OneD, NekDouble> m_fac0;
    Array<OneD, NekDouble> m_fac1;
    Array<OneD, NekDouble> m_fac2;
    Array<OneD, NekDouble> m_fac3;
    bool m_sortTopEdge;
};

OperatorKey IProductWRTDerivBase_SumFac_Tet::m_type = GetOperatorFactory().
    RegisterCreatorFunction(OperatorKey(LibUtilities::eTetrahedron,
                                        eIProductWRTDerivBase, eSumFac,false),
                            IProductWRTDerivBase_SumFac_Tet::create,
                            "IProductWRTDerivBase_SumFac_Tet");


/*
 * ----------------------------------------------------------
 * IProductWRTDerivBase operators
 * ----------------------------------------------------------
 */
class IProductWRTDerivBase_SumFac_Prism : public Operator
{
public:
    IProductWRTDerivBase_SumFac_Prism(vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                                      CoalescedGeomDataSharedPtr GeomData)
        : Operator  (pCollExp, GeomData),
          m_nquad0  (m_stdExp->GetNumPoints(0)),
          m_nquad1  (m_stdExp->GetNumPoints(1)),
          m_nquad2  (m_stdExp->GetNumPoints(2)),
          m_nmodes0 (m_stdExp->GetBasisNumModes(0)),
          m_nmodes1 (m_stdExp->GetBasisNumModes(1)),
          m_nmodes2 (m_stdExp->GetBasisNumModes(2)),
          m_base0    (m_stdExp->GetBasis(0)->GetBdata()),
          m_base1    (m_stdExp->GetBasis(1)->GetBdata()),
          m_base2    (m_stdExp->GetBasis(2)->GetBdata()),
          m_derbase0    (m_stdExp->GetBasis(0)->GetDbdata()),
          m_derbase1    (m_stdExp->GetBasis(1)->GetDbdata()),
          m_derbase2    (m_stdExp->GetBasis(2)->GetDbdata())

    {
        m_jac      = GeomData->GetJacWithStdWeights(pCollExp);
        m_wspSize  = 6*m_numElmt*(max(m_nquad0*m_nquad1*m_nquad2,m_nmodes0*m_nmodes1*m_nmodes2));
        m_derivFac = GeomData->GetDerivFactors(pCollExp);

        if(m_stdExp->GetBasis(0)->GetBasisType() == LibUtilities::eModified_A)
        {
            m_sortTopVertex = true;
        }
        else
        {
            m_sortTopVertex = false;
        }

        const Array<OneD, const NekDouble>& z0 = m_stdExp->GetBasis(0)->GetZ();
        const Array<OneD, const NekDouble>& z2 = m_stdExp->GetBasis(2)->GetZ();

        m_fac0 = Array<OneD, NekDouble>(m_nquad0*m_nquad1*m_nquad2);
        m_fac1 = Array<OneD, NekDouble>(m_nquad0*m_nquad1*m_nquad2);

        for (int i = 0; i < m_nquad0; ++i)
        {
            for(int j = 0; j < m_nquad1; ++j)
            {
                for(int k = 0; k < m_nquad2; ++k)
                {
                    // set up geometric factor: 2/(1-z1)
                    m_fac0[i+j*m_nquad0+k*m_nquad0*m_nquad1] = 2.0/(1-z2[k]);
                    // set up geometric factor: (1+z0)/(1-z1)
                    m_fac1[i+j*m_nquad0+k*m_nquad0*m_nquad1] = (1+z0[i])/(1-z2[k]);

                }
            }
        }
    }

    virtual void operator()(const Array<OneD, const NekDouble> &entry0,
                            Array<OneD, NekDouble> &entry1,
                            Array<OneD, NekDouble> &entry2,
                            Array<OneD, NekDouble> &entry3,
                            Array<OneD, NekDouble> &wsp)
    {
        unsigned int nPhys  = m_stdExp->GetTotPoints();
        unsigned int ntot   = m_numElmt*nPhys;
        unsigned int nmodes = m_stdExp->GetNcoeffs();
        unsigned int nmax  = max(ntot,m_numElmt*nmodes);
        Array<OneD, Array<OneD, const NekDouble> > in(3);
        Array<OneD, NekDouble> output, wsp1;
        Array<OneD, Array<OneD, NekDouble> > tmp(3);

        in[0] = entry0; in[1] = entry1;
        in[2] = entry2;

        output =  entry3;

        for(int i = 0; i < 3; ++i)
        {
            tmp[i] = wsp + i*nmax;
        }

        // calculate (dphi/dx,in[0]) = ((dphi/dxi_0 dxi_0/dx + dphi/dxi_1 dxi_1/dx),in[0])
        //     +     (dphi/dy,in[1]) = ((dphi/dxi_0 dxi_0/dy + dphi/dxi_1 dxi_1/dy),in[1])
        //     +     (dphi/dz,in[2]) = ((dphi/dxi_0 dxi_0/dz + dphi/dxi_1 dxi_1/dz),in[2])
        //
        // Note dphi/dxi_0  =
        //             dphi/deta_0 deta_0/dxi_0 = dphi/deta_0 2/(1-eta_2)
        //
        //      dphi/dxi_2  =
        //             dphi/deta_0 deta_0/dxi_2 + dphi/deta_2 deta_2/dxi_2 =
        //             dphi/deta_0 (1+eta_0)/(1-eta_2) + dphi/deta_2
        //
        // and so the full inner products are
        //
        // (dphi/dx,in[0]) + (dphi/dy,in[1]) + (dphi/dz,in[2])
        //    = (dphi/deta_0, ((2/(1-eta_2) (dxi_0/dx in[0] + dxi_0/dy in[1] + dxi_0/dz in[2]   )
        //            + (1_eta_0)/(1-eta_2) (dxi_2/dx in[0] + dxi_2/dy in[1] + dxi_2/dz in[2] ))
        //    + (dphi/deta_1, (dxi_1/dx in[0] + dxi_1/dy in[1] + dxi_1/dz in[2]))
        //    + (dphi/deta_2, (dxi_2/dx in[0] + dxi_2/dy in[1] + dxi_2/dz in[2]))



        for(int i = 0; i < 3; ++i)
        {
            Vmath::Vmul (ntot,m_derivFac[i],1, in[0],1,
                         tmp[i],1);
            for(int j = 1; j < 3; ++j)
            {
                Vmath::Vvtvp (ntot,m_derivFac[i+3*j],1,
                              in[j],1, tmp[i], 1, tmp[i],1);
            }
        }
        wsp1   = wsp + 3*nmax;

        // Sort into eta factors
        for (int i = 0; i < m_numElmt; ++i)
        {
            // scale tmp[0] by fac0
            Vmath::Vmul(nPhys,&m_fac0[0],1,tmp[0].get()+i*nPhys,1,
                        tmp[0].get()+i*nPhys,1);

            // scale tmp[2] by fac1 and add to tmp0
            Vmath::Vvtvp(nPhys,&m_fac1[0],1,tmp[2].get()+i*nPhys,1,
                         tmp[0].get()+i*nPhys,1,tmp[0].get()+i*nPhys,1);
        }

        // calculate Iproduct WRT Std Deriv
        PrismIProduct(m_sortTopVertex, m_numElmt,
                    m_nquad0,   m_nquad1,  m_nquad2,
                    m_nmodes0,  m_nmodes1, m_nmodes2,
                    m_derbase0, m_base1,   m_base2,
                    m_jac,tmp[0],output,wsp1);

        PrismIProduct(m_sortTopVertex, m_numElmt,
                    m_nquad0,  m_nquad1,   m_nquad2,
                    m_nmodes0, m_nmodes1,  m_nmodes2,
                    m_base0,   m_derbase1, m_base2,
                    m_jac,tmp[1],tmp[0],wsp1);
        Vmath::Vadd(m_numElmt*nmodes,tmp[0],1,output,1,output,1);

        PrismIProduct(m_sortTopVertex, m_numElmt,
                    m_nquad0,  m_nquad1,  m_nquad2,
                    m_nmodes0, m_nmodes1, m_nmodes2,
                    m_base0,   m_base1,   m_derbase2,
                    m_jac,tmp[2],tmp[0],wsp1);
        Vmath::Vadd(m_numElmt*nmodes,tmp[0],1,output,1,output,1);
    }

    OPERATOR_CREATE(IProductWRTDerivBase_SumFac_Prism)

    protected:
    const int  m_nquad0;
    const int  m_nquad1;
    const int  m_nquad2;
    const int  m_nmodes0;
    const int  m_nmodes1;
    const int  m_nmodes2;
    Array<OneD, const NekDouble> m_jac;
    Array<OneD, const NekDouble> m_base0;
    Array<OneD, const NekDouble> m_base1;
    Array<OneD, const NekDouble> m_base2;
    Array<OneD, const NekDouble> m_derbase0;
    Array<OneD, const NekDouble> m_derbase1;
    Array<OneD, const NekDouble> m_derbase2;
    Array<TwoD, const NekDouble> m_derivFac;
    Array<OneD, NekDouble> m_fac0;
    Array<OneD, NekDouble> m_fac1;
    bool m_sortTopVertex;
};

OperatorKey IProductWRTDerivBase_SumFac_Prism::m_type = GetOperatorFactory().
    RegisterCreatorFunction(OperatorKey(LibUtilities::ePrism,
                                        eIProductWRTDerivBase, eSumFac,false),
                            IProductWRTDerivBase_SumFac_Prism::create,
                            "IProductWRTDerivBase_SumFac_Prism");

}
}
