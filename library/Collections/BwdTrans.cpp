///////////////////////////////////////////////////////////////////////////////
//
// File: BwdTrans.cpp
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
// Description: BwdTrans operator implementations
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
 * @brief Backward transform operator using standard matrix approach.
 */
class BwdTrans_StdMat : public Operator
{
    public:
        OPERATOR_CREATE(BwdTrans_StdMat)

        virtual ~BwdTrans_StdMat()
        {
        }

        virtual void operator()(
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output,
                      Array<OneD,       NekDouble> &output1,
                      Array<OneD,       NekDouble> &output2,
                      Array<OneD,       NekDouble> &wsp)
        {
            boost::ignore_unused(output1, output2, wsp);
            Blas::Dgemm('N', 'N', m_mat->GetRows(), m_numElmt,
                        m_mat->GetColumns(), 1.0, m_mat->GetRawPtr(),
                        m_mat->GetRows(), input.get(), m_stdExp->GetNcoeffs(),
                        0.0, output.get(), m_stdExp->GetTotPoints());
        }

        virtual void operator()(
                      int                           dir,
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output,
                      Array<OneD,       NekDouble> &wsp)
        {
            boost::ignore_unused(dir, input, output, wsp);
            ASSERTL0(false, "Not valid for this operator.");
        }

    protected:
        DNekMatSharedPtr m_mat;

    private:
        BwdTrans_StdMat(
                vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                CoalescedGeomDataSharedPtr                pGeomData)
            : Operator(pCollExp, pGeomData)
        {
            StdRegions::StdMatrixKey  key(StdRegions::eBwdTrans,
                                          m_stdExp->DetShapeType(), *m_stdExp);
            m_mat = m_stdExp->GetStdMatrix(key);
        }
};

/// Factory initialisation for the BwdTrans_StdMat operators
OperatorKey BwdTrans_StdMat::m_typeArr[] = {
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eSegment,       eBwdTrans, eStdMat,false),
        BwdTrans_StdMat::create, "BwdTrans_StdMat_Seg"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTriangle,      eBwdTrans, eStdMat,false),
        BwdTrans_StdMat::create, "BwdTrans_StdMat_Tri"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTriangle,      eBwdTrans, eStdMat,true),
        BwdTrans_StdMat::create, "BwdTrans_StdMat_NodalTri"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eQuadrilateral, eBwdTrans, eStdMat,false),
        BwdTrans_StdMat::create, "BwdTrans_StdMat_Quad"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTetrahedron,   eBwdTrans, eStdMat,false),
        BwdTrans_StdMat::create, "BwdTrans_StdMat_Tet"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTetrahedron,   eBwdTrans, eStdMat,true),
        BwdTrans_StdMat::create, "BwdTrans_StdMat_NodalTet"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePyramid,       eBwdTrans, eStdMat,false),
        BwdTrans_StdMat::create, "BwdTrans_StdMat_Pyr"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePrism,         eBwdTrans, eStdMat,false),
        BwdTrans_StdMat::create, "BwdTrans_StdMat_Prism"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePrism,         eBwdTrans, eStdMat,true),
        BwdTrans_StdMat::create, "BwdTrans_StdMat_NodalPrism"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eHexahedron,    eBwdTrans, eStdMat,false),
        BwdTrans_StdMat::create, "BwdTrans_StdMat_Hex"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePyramid, eBwdTrans, eSumFac, false),
        BwdTrans_StdMat::create, "BwdTrans_SumFac_Pyr")
};



/**
 * @brief Backward transform operator using default StdRegions operator
 */
class BwdTrans_IterPerExp : public Operator
{
    public:
        OPERATOR_CREATE(BwdTrans_IterPerExp)

        virtual ~BwdTrans_IterPerExp()
        {
        }

        virtual void operator()(
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output,
                      Array<OneD,       NekDouble> &output1,
                      Array<OneD,       NekDouble> &output2,
                      Array<OneD,       NekDouble> &wsp)
        {
            boost::ignore_unused(output1, output2, wsp);

            const int nCoeffs = m_stdExp->GetNcoeffs();
            const int nPhys   = m_stdExp->GetTotPoints();
            Array<OneD, NekDouble> tmp;

            for (int i = 0; i < m_numElmt; ++i)
            {
                m_stdExp->BwdTrans(input + i*nCoeffs, tmp = output + i*nPhys);
            }
        }

        virtual void operator()(
                      int                           dir,
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output,
                      Array<OneD,       NekDouble> &wsp)
        {
            boost::ignore_unused(dir, input, output, wsp);
            ASSERTL0(false, "Not valid for this operator.");
        }

    private:
        BwdTrans_IterPerExp(
                vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                CoalescedGeomDataSharedPtr                pGeomData)
            : Operator(pCollExp, pGeomData)
        {
        }
};

/// Factory initialisation for the BwdTrans_IterPerExp operators
OperatorKey BwdTrans_IterPerExp::m_typeArr[] = {
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eSegment,       eBwdTrans, eIterPerExp,false),
        BwdTrans_IterPerExp::create, "BwdTrans_IterPerExp_Seg"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTriangle,      eBwdTrans, eIterPerExp,false),
        BwdTrans_IterPerExp::create, "BwdTrans_IterPerExp_Tri"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTriangle,      eBwdTrans, eIterPerExp,true),
        BwdTrans_IterPerExp::create, "BwdTrans_IterPerExp_NodalTri"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eQuadrilateral, eBwdTrans, eIterPerExp,false),
        BwdTrans_IterPerExp::create, "BwdTrans_IterPerExp_Quad"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTetrahedron,   eBwdTrans, eIterPerExp,false),
        BwdTrans_IterPerExp::create, "BwdTrans_IterPerExp_Tet"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTetrahedron,   eBwdTrans, eIterPerExp,true),
        BwdTrans_IterPerExp::create, "BwdTrans_IterPerExp_NodalTet"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePyramid,       eBwdTrans, eIterPerExp,false),
        BwdTrans_IterPerExp::create, "BwdTrans_IterPerExp_Pyr"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePrism,         eBwdTrans, eIterPerExp,false),
        BwdTrans_IterPerExp::create, "BwdTrans_IterPerExp_Prism"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePrism,         eBwdTrans, eIterPerExp,true),
        BwdTrans_IterPerExp::create, "BwdTrans_IterPerExp_NodalPrism"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eHexahedron,    eBwdTrans, eIterPerExp,false),
        BwdTrans_IterPerExp::create, "BwdTrans_IterPerExp_Hex"),
};


/**
 * @brief Backward transform operator using LocalRegions implementation.
 */
class BwdTrans_NoCollection : public Operator
{
    public:
        OPERATOR_CREATE(BwdTrans_NoCollection)

        virtual ~BwdTrans_NoCollection()
        {
        }

        virtual void operator()(
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output,
                      Array<OneD,       NekDouble> &output1,
                      Array<OneD,       NekDouble> &output2,
                      Array<OneD,       NekDouble> &wsp)
        {
            boost::ignore_unused(output1, output2, wsp);

            const int nCoeffs = m_expList[0]->GetNcoeffs();
            const int nPhys   = m_expList[0]->GetTotPoints();
            Array<OneD, NekDouble> tmp;

            for (int i = 0; i < m_numElmt; ++i)
            {
                m_expList[i]->BwdTrans(input + i*nCoeffs,
                                       tmp = output + i*nPhys);
            }
        }

        virtual void operator()(
                      int                           dir,
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output,
                      Array<OneD,       NekDouble> &wsp)
        {
            boost::ignore_unused(dir, input, output, wsp);
            ASSERTL0(false, "Not valid for this operator.");
        }

    protected:
        vector<StdRegions::StdExpansionSharedPtr> m_expList;

    private:
        BwdTrans_NoCollection(
                vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                CoalescedGeomDataSharedPtr                pGeomData)
            : Operator(pCollExp, pGeomData)
        {
            m_expList = pCollExp;
        }
};

/// Factory initialisation for the BwdTrans_NoCollection operators
OperatorKey BwdTrans_NoCollection::m_typeArr[] = {
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eSegment,       eBwdTrans, eNoCollection,false),
        BwdTrans_NoCollection::create, "BwdTrans_NoCollection_Seg"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTriangle,      eBwdTrans, eNoCollection,false),
        BwdTrans_NoCollection::create, "BwdTrans_NoCollection_Tri"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTriangle,      eBwdTrans, eNoCollection,true),
        BwdTrans_NoCollection::create, "BwdTrans_NoCollection_NodalTri"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eQuadrilateral, eBwdTrans, eNoCollection,false),
        BwdTrans_NoCollection::create, "BwdTrans_NoCollection_Quad"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTetrahedron,   eBwdTrans, eNoCollection,false),
        BwdTrans_NoCollection::create, "BwdTrans_NoCollection_Tet"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTetrahedron,   eBwdTrans, eNoCollection,true),
        BwdTrans_NoCollection::create, "BwdTrans_NoCollection_NodalTet"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePyramid,       eBwdTrans, eNoCollection,false),
        BwdTrans_NoCollection::create, "BwdTrans_NoCollection_Pyr"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePrism,         eBwdTrans, eNoCollection,false),
        BwdTrans_NoCollection::create, "BwdTrans_NoCollection_Prism"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePrism,         eBwdTrans, eNoCollection,true),
        BwdTrans_NoCollection::create, "BwdTrans_NoCollection_NodalPrism"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eHexahedron,    eBwdTrans, eNoCollection,false),
        BwdTrans_NoCollection::create, "BwdTrans_NoCollection_Hex"),
};


/**
 * @brief Backward transform operator using sum-factorisation (Segment)
 */
class BwdTrans_SumFac_Seg : public Operator
{
    public:
        OPERATOR_CREATE(BwdTrans_SumFac_Seg)

        virtual ~BwdTrans_SumFac_Seg()
        {
        }

        virtual void operator()(
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output,
                      Array<OneD,       NekDouble> &output1,
                      Array<OneD,       NekDouble> &output2,
                      Array<OneD,       NekDouble> &wsp)
        {
            boost::ignore_unused(output1, output2, wsp);
            if(m_colldir0)
            {
                Vmath::Vcopy(m_numElmt*m_nmodes0,input.get(),1,output.get(),1);
            }
            else
            {
                // out = B0*in;
                Blas::Dgemm('N','N', m_nquad0, m_numElmt, m_nmodes0,
                            1.0, m_base0.get(), m_nquad0,
                                 &input[0],     m_nmodes0, 0.0,
                                 &output[0],    m_nquad0);
            }
        }

        virtual void operator()(
                      int                           dir,
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output,
                      Array<OneD,       NekDouble> &wsp)
        {
            boost::ignore_unused(dir, input, output, wsp);
            ASSERTL0(false, "Not valid for this operator.");
        }


    protected:
        const int                       m_nquad0;
        const int                       m_nmodes0;
        const bool                      m_colldir0;
        Array<OneD, const NekDouble>    m_base0;

    private:
        BwdTrans_SumFac_Seg(
                vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                CoalescedGeomDataSharedPtr                pGeomData)
            : Operator  (pCollExp, pGeomData),
              m_nquad0  (m_stdExp->GetNumPoints(0)),
              m_nmodes0 (m_stdExp->GetBasisNumModes(0)),
              m_colldir0(m_stdExp->GetBasis(0)->Collocation()),
              m_base0   (m_stdExp->GetBasis(0)->GetBdata())
        {
            m_wspSize = 0;
        }
};

/// Factory initialisation for the BwdTrans_SumFac_Seg operator
OperatorKey BwdTrans_SumFac_Seg::m_type = GetOperatorFactory().
    RegisterCreatorFunction(
        OperatorKey(eSegment, eBwdTrans, eSumFac, false),
        BwdTrans_SumFac_Seg::create, "BwdTrans_SumFac_Seg");



/**
 * @brief Backward transform operator using sum-factorisation (Quad)
 */
class BwdTrans_SumFac_Quad : public Operator
{
    public:
        OPERATOR_CREATE(BwdTrans_SumFac_Quad)

        virtual ~BwdTrans_SumFac_Quad()
        {
        }

        virtual void operator()(
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output,
                      Array<OneD,       NekDouble> &output1,
                      Array<OneD,       NekDouble> &output2,
                      Array<OneD,       NekDouble> &wsp)
        {
            boost::ignore_unused(output1, output2);

            int i = 0;
            if(m_colldir0 && m_colldir1)
            {
                Vmath::Vcopy(m_numElmt*m_nmodes0*m_nmodes1,
                             input.get(), 1, output.get(), 1);
            }
            else if(m_colldir0)
            {
                Array<OneD, const NekDouble> base1
                        = m_stdExp->GetBasis(1)->GetBdata();

                for(i = 0; i < m_numElmt; ++i)
                {
                    Blas::Dgemm('N', 'T', m_nquad0, m_nquad1, m_nmodes1,
                                1.0, &input[i*m_nquad0*m_nmodes1], m_nquad0,
                                     base1.get(), m_nquad1, 0.0,
                                     &output[i*m_nquad0*m_nquad1], m_nquad0);
                }
            }
            else if(m_colldir1)
            {
                Blas::Dgemm('N', 'N', m_nquad0, m_nmodes1*m_numElmt, m_nmodes0,
                            1.0, m_base0.get(), m_nquad0,
                            &input[0],  m_nmodes0, 0.0,
                            &output[0], m_nquad0);
            }
            else
            {
                ASSERTL1(wsp.size() == m_wspSize,
                         "Incorrect workspace size");

                // Those two calls correpsond to the operation
                // out = B0*in*Transpose(B1);
                Blas::Dgemm('N', 'N', m_nquad0, m_nmodes1*m_numElmt, m_nmodes0,
                            1.0, m_base0.get(), m_nquad0,
                            &input[0], m_nmodes0, 0.0,
                            &wsp[0],   m_nquad0);

                for(i = 0; i < m_numElmt; ++i)
                {
                    Blas::Dgemm('N','T', m_nquad0, m_nquad1, m_nmodes1,
                                1.0, &wsp[i*m_nquad0*m_nmodes1], m_nquad0,
                                m_base1.get(), m_nquad1, 0.0,
                                &output[i*m_nquad0*m_nquad1],    m_nquad0);
                }
            }
        }

        virtual void operator()(
                      int                           dir,
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output,
                      Array<OneD,       NekDouble> &wsp)
        {
            boost::ignore_unused(dir, input, output, wsp);
            ASSERTL0(false, "Not valid for this operator.");
        }

    protected:
        const int                       m_nquad0;
        const int                       m_nquad1;
        const int                       m_nmodes0;
        const int                       m_nmodes1;
        const bool                      m_colldir0;
        const bool                      m_colldir1;
        Array<OneD, const NekDouble>    m_base0;
        Array<OneD, const NekDouble>    m_base1;

    private:
        BwdTrans_SumFac_Quad(
                vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                CoalescedGeomDataSharedPtr                pGeomData)
            : Operator  (pCollExp, pGeomData),
              m_nquad0  (m_stdExp->GetNumPoints(0)),
              m_nquad1  (m_stdExp->GetNumPoints(1)),
              m_nmodes0 (m_stdExp->GetBasisNumModes(0)),
              m_nmodes1 (m_stdExp->GetBasisNumModes(1)),
              m_colldir0(m_stdExp->GetBasis(0)->Collocation()),
              m_colldir1(m_stdExp->GetBasis(1)->Collocation()),
              m_base0   (m_stdExp->GetBasis(0)->GetBdata()),
              m_base1   (m_stdExp->GetBasis(1)->GetBdata())
        {
            m_wspSize = m_nquad0*m_nmodes1*m_numElmt;
        }
};

/// Factory initialisation for the BwdTrans_SumFac_Quad operator
OperatorKey BwdTrans_SumFac_Quad::m_type = GetOperatorFactory().
    RegisterCreatorFunction(
        OperatorKey(eQuadrilateral, eBwdTrans, eSumFac,false),
        BwdTrans_SumFac_Quad::create, "BwdTrans_SumFac_Quad");


/**
 * @brief Backward transform operator using sum-factorisation (Tri)
 */
class BwdTrans_SumFac_Tri : public Operator
{
    public:
        OPERATOR_CREATE(BwdTrans_SumFac_Tri)

        virtual ~BwdTrans_SumFac_Tri()
        {
        }

        virtual void operator()(
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output,
                      Array<OneD,       NekDouble> &output1,
                      Array<OneD,       NekDouble> &output2,
                      Array<OneD,       NekDouble> &wsp)
        {
            boost::ignore_unused(output1, output2);

            ASSERTL1(wsp.size() == m_wspSize,
                     "Incorrect workspace size");

            int ncoeffs = m_stdExp->GetNcoeffs();
            int i       = 0;
            int mode    = 0;

            for (i = mode = 0; i < m_nmodes0; ++i)
            {
                Blas::Dgemm('N', 'N', m_nquad1, m_numElmt, m_nmodes1-i,
                            1.0, m_base1.get()+mode*m_nquad1, m_nquad1,
                            &input[0]+mode, ncoeffs, 0.0,
                            &wsp[i*m_nquad1*m_numElmt],       m_nquad1);
                mode += m_nmodes1-i;
            }

            // fix for modified basis by splitting top vertex mode
            if(m_sortTopVertex)
            {
                for(i = 0; i < m_numElmt; ++i)
                {
                    Blas::Daxpy(m_nquad1, input[1+i*ncoeffs],
                                m_base1.get() + m_nquad1, 1,
                                &wsp[m_nquad1*m_numElmt] + i*m_nquad1, 1);
                }

            }

            Blas::Dgemm('N', 'T', m_nquad0, m_nquad1*m_numElmt, m_nmodes0,
                        1.0, m_base0.get(), m_nquad0,
                        &wsp[0], m_nquad1*m_numElmt, 0.0,
                        &output[0], m_nquad0);
        }

        virtual void operator()(
                      int                           dir,
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output,
                      Array<OneD,       NekDouble> &wsp)
        {
            boost::ignore_unused(dir, input, output, wsp);
            ASSERTL0(false, "Not valid for this operator.");
        }


    protected:
        const int                       m_nquad0;
        const int                       m_nquad1;
        const int                       m_nmodes0;
        const int                       m_nmodes1;
        Array<OneD, const NekDouble>    m_base0;
        Array<OneD, const NekDouble>    m_base1;
        bool                            m_sortTopVertex;

    private:
        BwdTrans_SumFac_Tri(
                vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                CoalescedGeomDataSharedPtr                pGeomData)
            : Operator  (pCollExp, pGeomData),
              m_nquad0  (m_stdExp->GetNumPoints(0)),
              m_nquad1  (m_stdExp->GetNumPoints(1)),
              m_nmodes0 (m_stdExp->GetBasisNumModes(0)),
              m_nmodes1 (m_stdExp->GetBasisNumModes(1)),
              m_base0   (m_stdExp->GetBasis(0)->GetBdata()),
              m_base1   (m_stdExp->GetBasis(1)->GetBdata())
        {
            m_wspSize = m_nquad0 * m_nmodes1 * m_numElmt;
            if(m_stdExp->GetBasis(0)->GetBasisType()
                    == LibUtilities::eModified_A)
            {
                m_sortTopVertex = true;
            }
            else
            {
                m_sortTopVertex = false;
            }
        }

};

/// Factory initialisation for the BwdTrans_SumFac_Tri operator
OperatorKey BwdTrans_SumFac_Tri::m_type = GetOperatorFactory().
    RegisterCreatorFunction(
        OperatorKey(eTriangle, eBwdTrans, eSumFac,false),
        BwdTrans_SumFac_Tri::create, "BwdTrans_SumFac_Tri");


/// Backward transform operator using sum-factorisation (Hex)
class BwdTrans_SumFac_Hex : public Operator
{
    public:
        OPERATOR_CREATE(BwdTrans_SumFac_Hex)

        virtual ~BwdTrans_SumFac_Hex()
        {
        }

        virtual void operator()(
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output,
                      Array<OneD,       NekDouble> &output1,
                      Array<OneD,       NekDouble> &output2,
                      Array<OneD,       NekDouble> &wsp)
        {
            boost::ignore_unused(output1, output2);

            if(m_colldir0 && m_colldir1 && m_colldir2)
            {
                Vmath::Vcopy(m_numElmt * m_nmodes0 * m_nmodes1 * m_nmodes2,
                             input.get(),  1,
                             output.get(), 1);
            }
            else
            {
                ASSERTL1(wsp.size() == m_wspSize,
                         "Incorrect workspace size");

                // Assign second half of workspace for 2nd DGEMM operation.
                int totmodes  = m_nmodes0*m_nmodes1*m_nmodes2;

                Array<OneD, NekDouble> wsp2
                                = wsp + m_nmodes0*m_nmodes1*m_nquad2*m_numElmt;

                //loop over elements  and do bwd trans wrt c
                for(int n = 0; n < m_numElmt; ++n)
                {
                    Blas::Dgemm('N', 'T', m_nquad2, m_nmodes0*m_nmodes1,
                                m_nmodes2, 1.0, m_base2.get(), m_nquad2,
                                &input[n*totmodes], m_nmodes0*m_nmodes1, 0.0,
                                &wsp[n*m_nquad2], m_nquad2*m_numElmt);
                }

                // trans wrt b
                Blas::Dgemm('N', 'T', m_nquad1, m_nquad2*m_numElmt*m_nmodes0,
                            m_nmodes1, 1.0, m_base1.get(), m_nquad1,
                            wsp.get(), m_nquad2*m_numElmt*m_nmodes0, 0.0,
                            wsp2.get(), m_nquad1);

                // trans wrt a
                Blas::Dgemm('N', 'T', m_nquad0, m_nquad1*m_nquad2*m_numElmt,
                            m_nmodes0, 1.0, m_base0.get(), m_nquad0,
                            wsp2.get(), m_nquad1*m_nquad2*m_numElmt, 0.0,
                            output.get(), m_nquad0);
            }
        }

        virtual void operator()(
                      int                           dir,
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output,
                      Array<OneD,       NekDouble> &wsp)
        {
            boost::ignore_unused(dir, input, output, wsp);
            ASSERTL0(false, "Not valid for this operator.");
        }

    protected:
        const int                       m_nquad0;
        const int                       m_nquad1;
        const int                       m_nquad2;
        const int                       m_nmodes0;
        const int                       m_nmodes1;
        const int                       m_nmodes2;
        Array<OneD, const NekDouble>    m_base0;
        Array<OneD, const NekDouble>    m_base1;
        Array<OneD, const NekDouble>    m_base2;
        const bool                      m_colldir0;
        const bool                      m_colldir1;
        const bool                      m_colldir2;

    private:
        BwdTrans_SumFac_Hex(
                vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                CoalescedGeomDataSharedPtr                pGeomData)
            : Operator  (pCollExp, pGeomData),
              m_nquad0  (pCollExp[0]->GetNumPoints(0)),
              m_nquad1  (pCollExp[0]->GetNumPoints(1)),
              m_nquad2  (pCollExp[0]->GetNumPoints(2)),
              m_nmodes0 (pCollExp[0]->GetBasisNumModes(0)),
              m_nmodes1 (pCollExp[0]->GetBasisNumModes(1)),
              m_nmodes2 (pCollExp[0]->GetBasisNumModes(2)),
              m_base0   (pCollExp[0]->GetBasis(0)->GetBdata()),
              m_base1   (pCollExp[0]->GetBasis(1)->GetBdata()),
              m_base2   (pCollExp[0]->GetBasis(2)->GetBdata()),
              m_colldir0(pCollExp[0]->GetBasis(0)->Collocation()),
              m_colldir1(pCollExp[0]->GetBasis(1)->Collocation()),
              m_colldir2(pCollExp[0]->GetBasis(2)->Collocation())
        {
            m_wspSize =  m_numElmt*m_nmodes0*(m_nmodes1*m_nquad2 +
                                              m_nquad1*m_nquad2);
        }
};

/// Factory initialisation for the BwdTrans_SumFac_Hex operator
OperatorKey BwdTrans_SumFac_Hex::m_type = GetOperatorFactory().
    RegisterCreatorFunction(
        OperatorKey(eHexahedron, eBwdTrans, eSumFac,false),
        BwdTrans_SumFac_Hex::create, "BwdTrans_SumFac_Hex");


/**
 * @brief Backward transform operator using sum-factorisation (Tet)
 */
class BwdTrans_SumFac_Tet : public Operator
{
    public:
        OPERATOR_CREATE(BwdTrans_SumFac_Tet)

        virtual ~BwdTrans_SumFac_Tet()
        {
        }

        virtual void operator()(
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output,
                      Array<OneD,       NekDouble> &output1,
                      Array<OneD,       NekDouble> &output2,
                      Array<OneD,       NekDouble> &wsp)
        {
            boost::ignore_unused(output1, output2);

            ASSERTL1(wsp.size() == m_wspSize,
                     "Incorrect workspace size");

            Array<OneD, NekDouble > tmp  = wsp;
            Array<OneD, NekDouble > tmp1 = tmp
                    + m_numElmt*m_nquad2*m_nmodes0*(2*m_nmodes1-m_nmodes0+1)/2;

            int i       = 0;
            int j       = 0;
            int mode    = 0;
            int mode1   = 0;
            int cnt     = 0;
            int ncoeffs = m_stdExp->GetNcoeffs();

            // Perform summation over '2' direction
            for(i = 0; i < m_nmodes0; ++i)
            {
                for(j = 0; j < m_nmodes1-i; ++j, ++cnt)
                {
                    Blas::Dgemm('N', 'N', m_nquad2, m_numElmt, m_nmodes2-i-j,
                                1.0, m_base2.get()+mode*m_nquad2, m_nquad2,
                                input.get()+ mode1,  ncoeffs, 0.0,
                                tmp.get() + cnt*m_nquad2*m_numElmt, m_nquad2);
                    mode  += m_nmodes2-i-j;
                    mode1 += m_nmodes2-i-j;
                }

                //increment mode in case m_nmodes1!=m_nmodes2
                mode +=  (m_nmodes2-m_nmodes1)*(m_nmodes2-m_nmodes1+1)/2;
            }

            // vertex mode - currently (1+c)/2 x (1-b)/2 x (1-a)/2
            // component is evaluated
            if(m_sortTopEdge)
            {
                for(i = 0; i < m_numElmt; ++i)
                {
                    // top singular vertex
                    // (1+c)/2 x (1+b)/2 x (1-a)/2 component
                    Blas::Daxpy(m_nquad2, input[1+i*ncoeffs],
                                m_base2.get() + m_nquad2, 1,
                                &tmp[m_nquad2*m_numElmt] + i*m_nquad2, 1);

                    // top singular vertex
                    // (1+c)/2 x (1-b)/2 x (1+a)/2 component
                    Blas::Daxpy(m_nquad2, input[1+i*ncoeffs],
                                m_base2.get() + m_nquad2, 1,
                                &tmp[m_nmodes1*m_nquad2*m_numElmt]
                                     + i*m_nquad2, 1);
                }
            }

            // Perform summation over '1' direction
            mode = 0;
            for(i = 0; i < m_nmodes0; ++i)
            {
                Blas::Dgemm('N', 'T', m_nquad1, m_nquad2*m_numElmt, m_nmodes1-i,
                            1.0, m_base1.get()+mode*m_nquad1, m_nquad1,
                            tmp.get() + mode*m_nquad2*m_numElmt,
                            m_nquad2*m_numElmt,
                            0.0, tmp1.get() + i*m_nquad1*m_nquad2*m_numElmt,
                            m_nquad1);
                mode  += m_nmodes1-i;
            }

            // fix for modified basis by adding additional split of
            // top and base singular vertex modes as well as singular
            // edge
            if(m_sortTopEdge)
            {
                // this could probably be a dgemv or higher if we
                // made a specialised m_base1[m_nuqad1] array
                // containing multiply copies
                for(i = 0; i < m_numElmt; ++i)
                {
                    // sort out singular vertices and singular
                    // edge components with (1+b)/2 (1+a)/2 form
                    for(int j = 0; j < m_nquad2; ++j)
                    {
                        Blas::Daxpy(m_nquad1,
                                    tmp[m_nquad2*m_numElmt + i*m_nquad2 + j],
                                    m_base1.get() + m_nquad1,1,
                                    &tmp1[m_nquad1*m_nquad2*m_numElmt]
                                        + i*m_nquad1*m_nquad2 + j*m_nquad1, 1);
                    }
                }
            }

            // Perform summation over '0' direction
            Blas::Dgemm('N', 'T', m_nquad0, m_nquad1*m_nquad2*m_numElmt,
                        m_nmodes0, 1.0, m_base0.get(), m_nquad0,
                        tmp1.get(), m_nquad1*m_nquad2*m_numElmt, 0.0,
                        output.get(), m_nquad0);

        }

        virtual void operator()(
                      int                           dir,
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output,
                      Array<OneD,       NekDouble> &wsp)
        {
            boost::ignore_unused(dir, input, output, wsp);
            ASSERTL0(false, "Not valid for this operator.");
        }

    protected:
        const int                       m_nquad0;
        const int                       m_nquad1;
        const int                       m_nquad2;
        const int                       m_nmodes0;
        const int                       m_nmodes1;
        const int                       m_nmodes2;
        Array<OneD, const NekDouble>    m_base0;
        Array<OneD, const NekDouble>    m_base1;
        Array<OneD, const NekDouble>    m_base2;
        bool                            m_sortTopEdge;

    private:
        BwdTrans_SumFac_Tet(
                vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                CoalescedGeomDataSharedPtr                pGeomData)
            : Operator  (pCollExp, pGeomData),
              m_nquad0  (m_stdExp->GetNumPoints(0)),
              m_nquad1  (m_stdExp->GetNumPoints(1)),
              m_nquad2  (m_stdExp->GetNumPoints(2)),
              m_nmodes0 (m_stdExp->GetBasisNumModes(0)),
              m_nmodes1 (m_stdExp->GetBasisNumModes(1)),
              m_nmodes2 (m_stdExp->GetBasisNumModes(2)),
              m_base0   (m_stdExp->GetBasis(0)->GetBdata()),
              m_base1   (m_stdExp->GetBasis(1)->GetBdata()),
              m_base2   (m_stdExp->GetBasis(2)->GetBdata())
        {
            m_wspSize = m_numElmt
                          * (m_nquad2*m_nmodes0*(2*m_nmodes1-m_nmodes0+1)/2
                      + m_nquad2*m_nquad1*m_nmodes0);

            if(m_stdExp->GetBasis(0)->GetBasisType()
                    == LibUtilities::eModified_A)
            {
                m_sortTopEdge = true;
            }
            else
            {
                m_sortTopEdge = false;
            }
        }
};

/// Factory initialisation for the BwdTrans_SumFac_Tet operator
OperatorKey BwdTrans_SumFac_Tet::m_type = GetOperatorFactory().
    RegisterCreatorFunction(
        OperatorKey(eTetrahedron, eBwdTrans, eSumFac,false),
        BwdTrans_SumFac_Tet::create, "BwdTrans_SumFac_Tet");


/**
 * @brief Backward transform operator using sum-factorisation (Prism)
 */
class BwdTrans_SumFac_Prism : public Operator
{
    public:
        OPERATOR_CREATE(BwdTrans_SumFac_Prism)

        virtual ~BwdTrans_SumFac_Prism()
        {
        }

        virtual void operator()(
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output,
                      Array<OneD,       NekDouble> &output1,
                      Array<OneD,       NekDouble> &output2,
                      Array<OneD,       NekDouble> &wsp)
        {
            boost::ignore_unused(output1, output2);

            ASSERTL1(wsp.size() == m_wspSize,
                    "Incorrect workspace size");

            // Assign second half of workspace for 2nd DGEMM operation.
            int totmodes  = m_stdExp->GetNcoeffs();

            Array<OneD, NekDouble> wsp2
                    = wsp + m_nmodes0*m_nmodes1*m_nquad2*m_numElmt;

            Vmath::Zero(m_nmodes0*m_nmodes1*m_nquad2*m_numElmt, wsp, 1);
            int i     = 0;
            int j     = 0;
            int mode  = 0;
            int mode1 = 0;
            int cnt   = 0;
            for (i = mode = mode1 = 0; i < m_nmodes0; ++i)
            {
                cnt = i * m_nquad2 * m_numElmt;
                for (j = 0; j < m_nmodes1; ++j)
                {
                    Blas::Dgemm('N', 'N', m_nquad2, m_numElmt, m_nmodes2-i,
                                1.0,  m_base2.get()+mode*m_nquad2, m_nquad2,
                                input.get()+mode1, totmodes, 0.0,
                                &wsp[j*m_nquad2*m_numElmt*m_nmodes0+ cnt],
                                m_nquad2);
                    mode1 += m_nmodes2-i;
                }
                mode += m_nmodes2-i;
            }

            // fix for modified basis by splitting top vertex mode
            if(m_sortTopVertex)
            {
                for(j = 0; j < m_nmodes1; ++j)
                {
                    for(i = 0; i < m_numElmt; ++i)
                    {
                        Blas::Daxpy(m_nquad2,input[1+i*totmodes+j*m_nmodes2],
                                    m_base2.get()+m_nquad2,1,
                                    &wsp[j*m_nquad2*m_numElmt*m_nmodes0 +
                                         m_nquad2*m_numElmt]+i*m_nquad2,1);
                    }
                }
                // Believe this could be made into a m_nmodes1
                // dgemv if we made an array of m_numElmt copies
                // of m_base2[m_quad2] (which are of size
                // m_nquad2.
            }

            // Perform summation over '1' direction
            Blas::Dgemm('N', 'T', m_nquad1, m_nquad2*m_numElmt*m_nmodes0,
                        m_nmodes1,1.0, m_base1.get(),  m_nquad1,
                        wsp.get(), m_nquad2*m_numElmt*m_nmodes0,
                        0.0, wsp2.get(), m_nquad1);


            // Perform summation over '0' direction
            Blas::Dgemm('N', 'T', m_nquad0, m_nquad1*m_nquad2*m_numElmt,
                        m_nmodes0, 1.0, m_base0.get(),  m_nquad0,
                        wsp2.get(), m_nquad1*m_nquad2*m_numElmt,
                        0.0, output.get(), m_nquad0);
        }

        virtual void operator()(
                      int                           dir,
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output,
                      Array<OneD,       NekDouble> &wsp)
        {
            boost::ignore_unused(dir, input, output, wsp);
            ASSERTL0(false, "Not valid for this operator.");
        }


    protected:
        const int                       m_nquad0;
        const int                       m_nquad1;
        const int                       m_nquad2;
        const int                       m_nmodes0;
        const int                       m_nmodes1;
        const int                       m_nmodes2;
        Array<OneD, const NekDouble>    m_base0;
        Array<OneD, const NekDouble>    m_base1;
        Array<OneD, const NekDouble>    m_base2;
        bool                            m_sortTopVertex;

    private:
        BwdTrans_SumFac_Prism(
                vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                CoalescedGeomDataSharedPtr                pGeomData)
            : Operator  (pCollExp, pGeomData),
              m_nquad0  (m_stdExp->GetNumPoints(0)),
              m_nquad1  (m_stdExp->GetNumPoints(1)),
              m_nquad2  (m_stdExp->GetNumPoints(2)),
              m_nmodes0 (m_stdExp->GetBasisNumModes(0)),
              m_nmodes1 (m_stdExp->GetBasisNumModes(1)),
              m_nmodes2 (m_stdExp->GetBasisNumModes(2)),
              m_base0   (m_stdExp->GetBasis(0)->GetBdata()),
              m_base1   (m_stdExp->GetBasis(1)->GetBdata()),
              m_base2   (m_stdExp->GetBasis(2)->GetBdata())
        {
            m_wspSize = m_numElmt*m_nmodes0*(m_nmodes1*m_nquad2
                                                + m_nquad1*m_nquad2);

            if(m_stdExp->GetBasis(0)->GetBasisType()
                    == LibUtilities::eModified_A)
            {
                m_sortTopVertex = true;
            }
            else
            {
                m_sortTopVertex = false;
            }

        }
};

/// Factory initialisation for the BwdTrans_SumFac_Prism operator
OperatorKey BwdTrans_SumFac_Prism::m_type = GetOperatorFactory().
    RegisterCreatorFunction(
        OperatorKey(ePrism, eBwdTrans, eSumFac,false),
        BwdTrans_SumFac_Prism::create, "BwdTrans_SumFac_Prism");

/**
 * @brief Backward transform operator using sum-factorisation (Pyr)
 */
class BwdTrans_SumFac_Pyr : public Operator
{
    public:
        OPERATOR_CREATE(BwdTrans_SumFac_Pyr)

        virtual ~BwdTrans_SumFac_Pyr()
        {
        }

        virtual void operator()(
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output,
                      Array<OneD,       NekDouble> &output1,
                      Array<OneD,       NekDouble> &output2,
                      Array<OneD,       NekDouble> &wsp)
        {
            boost::ignore_unused(output1, output2);

            ASSERTL1(wsp.size() == m_wspSize,
                    "Incorrect workspace size");

            // Assign second half of workspace for 2nd DGEMM operation.
            int totmodes  = m_stdExp->GetNcoeffs();

            Array<OneD, NekDouble> wsp2
                    = wsp + m_nmodes0*m_nmodes1*m_nquad2*m_numElmt;

            Vmath::Zero(m_nmodes0*m_nmodes1*m_nquad2*m_numElmt, wsp, 1);
            int i     = 0;
            int j     = 0;
            int mode  = 0;
            int mode1 = 0;
            int cnt   = 0;
            for (i = 0; i < m_nmodes0; ++i)
            {
                for (j = 0; j < m_nmodes1; ++j, ++cnt)
                {
                    int ijmax = max(i,j);
                    Blas::Dgemm('N', 'N', m_nquad2, m_numElmt, m_nmodes2-ijmax,
                                1.0,  m_base2.get()+mode*m_nquad2, m_nquad2,
                                input.get()+mode1, totmodes, 0.0,
                                wsp.get() + cnt*m_nquad2*m_numElmt, m_nquad2);
                    mode  += m_nmodes2-ijmax;
                    mode1 += m_nmodes2-ijmax;
                }

                //increment mode in case order1!=order2
                for(j = m_nmodes1; j < m_nmodes2-i; ++j)
                {
                    int ijmax = max(i,j);
                    mode += m_nmodes2-ijmax;
                }
            }

            // vertex mode - currently (1+c)/2 x (1-b)/2 x (1-a)/2
            // component is evaluated
            if(m_sortTopVertex)
            {
                for(i = 0; i < m_numElmt; ++i)
                {
                    // top singular vertex
                    // (1+c)/2 x (1+b)/2 x (1-a)/2 component
                    Blas::Daxpy(m_nquad2, input[1+i*totmodes],
                                m_base2.get() + m_nquad2, 1,
                                &wsp[m_nquad2*m_numElmt] + i*m_nquad2, 1);

                    // top singular vertex
                    // (1+c)/2 x (1-b)/2 x (1+a)/2 component
                    Blas::Daxpy(m_nquad2, input[1+i*totmodes],
                                m_base2.get() + m_nquad2, 1,
                                &wsp[m_nmodes1*m_nquad2*m_numElmt]
                                + i*m_nquad2, 1);

                    // top singular vertex
                    // (1+c)/2 x (1+b)/2 x (1+a)/2 component
                    Blas::Daxpy(m_nquad2, input[1+i*totmodes],
                                m_base2.get() + m_nquad2, 1,
                                &wsp[(m_nmodes1+1)*m_nquad2*m_numElmt]
                                + i*m_nquad2, 1);

                }
            }

            // Perform summation over '1' direction
            mode = 0;
            for(i = 0; i < m_nmodes0; ++i)
            {
                Blas::Dgemm('N', 'T', m_nquad1, m_nquad2*m_numElmt, m_nmodes1,
                            1.0, m_base1.get(), m_nquad1,
                            wsp.get() + mode*m_nquad2*m_numElmt,
                            m_nquad2*m_numElmt,
                            0.0, wsp2.get() + i*m_nquad1*m_nquad2*m_numElmt,
                            m_nquad1);
                mode += m_nmodes1;
            }

            // Perform summation over '0' direction
            Blas::Dgemm('N', 'T', m_nquad0, m_nquad1*m_nquad2*m_numElmt,
                        m_nmodes0, 1.0, m_base0.get(),  m_nquad0,
                        wsp2.get(), m_nquad1*m_nquad2*m_numElmt,
                        0.0, output.get(), m_nquad0);
        }

        virtual void operator()(
                      int                           dir,
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output,
                      Array<OneD,       NekDouble> &wsp)
        {
            boost::ignore_unused(dir, input, output, wsp);
            ASSERTL0(false, "Not valid for this operator.");
        }


    protected:
        const int                       m_nquad0;
        const int                       m_nquad1;
        const int                       m_nquad2;
        const int                       m_nmodes0;
        const int                       m_nmodes1;
        const int                       m_nmodes2;
        Array<OneD, const NekDouble>    m_base0;
        Array<OneD, const NekDouble>    m_base1;
        Array<OneD, const NekDouble>    m_base2;
        bool                            m_sortTopVertex;

    private:
        BwdTrans_SumFac_Pyr(
                vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                CoalescedGeomDataSharedPtr                pGeomData)
            : Operator  (pCollExp, pGeomData),
              m_nquad0  (m_stdExp->GetNumPoints(0)),
              m_nquad1  (m_stdExp->GetNumPoints(1)),
              m_nquad2  (m_stdExp->GetNumPoints(2)),
              m_nmodes0 (m_stdExp->GetBasisNumModes(0)),
              m_nmodes1 (m_stdExp->GetBasisNumModes(1)),
              m_nmodes2 (m_stdExp->GetBasisNumModes(2)),
              m_base0   (m_stdExp->GetBasis(0)->GetBdata()),
              m_base1   (m_stdExp->GetBasis(1)->GetBdata()),
              m_base2   (m_stdExp->GetBasis(2)->GetBdata())
        {
            m_wspSize = m_numElmt*m_nmodes0*m_nquad2*(m_nmodes1 + m_nquad1);

            if(m_stdExp->GetBasis(0)->GetBasisType()
                    == LibUtilities::eModified_A)
            {
                m_sortTopVertex = true;
            }
            else
            {
                m_sortTopVertex = false;
            }

        }
};

/// Factory initialisation for the BwdTrans_SumFac_Pyr operator
OperatorKey BwdTrans_SumFac_Pyr::m_type = GetOperatorFactory().
    RegisterCreatorFunction(
        OperatorKey(ePyramid, eBwdTrans, eSumFac,false),
        BwdTrans_SumFac_Pyr::create, "BwdTrans_SumFac_Pyr");

}

}
