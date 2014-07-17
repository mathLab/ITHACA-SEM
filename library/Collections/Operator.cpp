///////////////////////////////////////////////////////////////////////////////
//
// File: Operator.h
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
// Description: Operator top class definition
//
///////////////////////////////////////////////////////////////////////////////

#include <loki/Singleton.h>
#include <Collections/Operator.h>

namespace Nektar {
namespace Collections {
    bool operator< (OperatorKey const &p1, OperatorKey const &p2)
    {
        if (boost::get<0>(p1) < boost::get<0>(p2))
        {
            return true;
        }
        if (boost::get<0>(p1) > boost::get<0>(p2))
        {
            return false;
        }
        if (boost::get<1>(p1) < boost::get<1>(p2))
        {
            return true;
        }
        if (boost::get<1>(p1) > boost::get<1>(p2))
        {
            return false;
        }
        if (boost::get<2>(p1) < boost::get<2>(p2))
        {
            return true;
        }
        if (boost::get<2>(p1) > boost::get<2>(p2))
        {
            return false;
        }

        return false;
    }

    std::ostream &operator<<(std::ostream &os, OperatorKey const &p)
    {
        os << boost::get<0>(p) << " "
           << boost::get<1>(p) << " "
           << boost::get<2>(p) << endl;
        return os;
    }

    OperatorFactory& GetOperatorFactory()
    {
        typedef Loki::SingletonHolder<OperatorFactory,
                                      Loki::CreateUsingNew,
                                      Loki::NoDestroy > Type;
        return Type::Instance();
    }

    /*
     * ----------------------------------------------------------
     * BwdTrans operators
     * ----------------------------------------------------------
     */

    class BwdTrans_LocMat : public Operator
    {
    public:
        BwdTrans_LocMat(StdRegions::StdExpansionSharedPtr pExp,
                        vector<SpatialDomains::GeometrySharedPtr> pGeom)
            : Operator(pExp, pGeom),
              m_key(StdRegions::eBwdTrans, pExp->DetShapeType(), *pExp)
        {
            m_mat = m_stdExp->GetStdMatrix(m_key);
        }

        virtual void operator()(
            const Array<OneD, const NekDouble> &input,
                  Array<OneD,       NekDouble> &output,
                  Array<OneD,       NekDouble> &wsp)
        {
            Blas::Dgemm('N', 'N', m_mat->GetRows(), m_numElmt,
                        m_mat->GetColumns(), 1.0, m_mat->GetRawPtr(),
                        m_mat->GetRows(), input.get(), m_stdExp->GetNcoeffs(),
                        0.0, output.get(), m_stdExp->GetTotPoints());
        }

        OPERATOR_CREATE(BwdTrans_LocMat)

        StdRegions::StdMatrixKey m_key;
        DNekMatSharedPtr m_mat;
    };

    OperatorKey BwdTrans_LocMat::m_typeArr[] =
    {
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::eSegment, eBwdTrans, eLocMat),
            BwdTrans_LocMat::create, "BwdTrans_LocMat_Seg"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::eTriangle, eBwdTrans, eLocMat),
            BwdTrans_LocMat::create, "BwdTrans_LocMat_Tri"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::eQuadrilateral, eBwdTrans, eLocMat),
            BwdTrans_LocMat::create, "BwdTrans_LocMat_Quad"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::eTetrahedron, eBwdTrans, eLocMat),
            BwdTrans_LocMat::create, "BwdTrans_LocMat_Tet"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::ePyramid, eBwdTrans, eLocMat),
            BwdTrans_LocMat::create, "BwdTrans_LocMat_Pyr"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::ePrism, eBwdTrans, eLocMat),
            BwdTrans_LocMat::create, "BwdTrans_LocMat_Prism"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::eHexahedron, eBwdTrans, eLocMat),
            BwdTrans_LocMat::create, "BwdTrans_LocMat_Hex"),
    };

    class BwdTrans_IterPerExp : public Operator
    {
    public:
        BwdTrans_IterPerExp(StdRegions::StdExpansionSharedPtr pExp,
                            vector<SpatialDomains::GeometrySharedPtr> pGeom)
            : Operator(pExp, pGeom)
        {
        }

        virtual void operator()(
            const Array<OneD, const NekDouble> &input,
                  Array<OneD,       NekDouble> &output,
                  Array<OneD,       NekDouble> &wsp)
        {
            const int nCoeffs = m_stdExp->GetNcoeffs();
            const int nPhys   = m_stdExp->GetTotPoints();
            Array<OneD, NekDouble> tmp;

            for (int i = 0; i < m_numElmt; ++i)
            {
                m_stdExp->BwdTrans(input + i*nCoeffs, tmp = output + i*nPhys);
            }
        }

        OPERATOR_CREATE(BwdTrans_IterPerExp)
    };

    OperatorKey BwdTrans_IterPerExp::m_typeArr[] =
    {
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::eSegment, eBwdTrans, eIterPerExp),
            BwdTrans_IterPerExp::create, "BwdTrans_IterPerExp_Seg"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::eTriangle, eBwdTrans, eIterPerExp),
            BwdTrans_IterPerExp::create, "BwdTrans_IterPerExp_Tri"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::eQuadrilateral, eBwdTrans, eIterPerExp),
            BwdTrans_IterPerExp::create, "BwdTrans_IterPerExp_Quad"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::eTetrahedron, eBwdTrans, eIterPerExp),
            BwdTrans_IterPerExp::create, "BwdTrans_IterPerExp_Tet"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::ePyramid, eBwdTrans, eIterPerExp),
            BwdTrans_IterPerExp::create, "BwdTrans_IterPerExp_Pyr"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::ePrism, eBwdTrans, eIterPerExp),
            BwdTrans_IterPerExp::create, "BwdTrans_IterPerExp_Prism"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::eHexahedron, eBwdTrans, eIterPerExp),
            BwdTrans_IterPerExp::create, "BwdTrans_IterPerExp_Hex"),
    };

    /*
     * ----------------------------------------------------------
     * PhysDeriv operators
     * ----------------------------------------------------------
     */

    class PhysDeriv_IterPerExp : public Operator
    {
    public:
        PhysDeriv_IterPerExp(StdRegions::StdExpansionSharedPtr pExp,
                             vector<SpatialDomains::GeometrySharedPtr> pGeom)
            : Operator(pExp, pGeom)
        {
        }

        virtual void operator()(
            const Array<OneD, const NekDouble> &input,
                  Array<OneD,       NekDouble> &output,
                  Array<OneD,       NekDouble> &wsp)
        {
            const int nPhys = m_stdExp->GetTotPoints();
            Array<OneD, NekDouble> tmp;

            for (int i = 0; i < m_numElmt; ++i)
            {
                m_stdExp->PhysDeriv(input + i*nPhys, tmp = output + i*nPhys);
            }
        }

        OPERATOR_CREATE(PhysDeriv_IterPerExp)
    };

    OperatorKey PhysDeriv_IterPerExp::m_typeArr[] =
    {
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::eQuadrilateral, ePhysDeriv, eIterPerExp),
            PhysDeriv_IterPerExp::create, "PhysDeriv_IterPerExp_Quad"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::eTriangle, ePhysDeriv, eIterPerExp),
            PhysDeriv_IterPerExp::create, "PhysDeriv_IterPerExp_Tri")
    };

    /*
     * ----------------------------------------------------------
     * IProductWRTBase operators
     * ----------------------------------------------------------
     */

    class IProductWRTBase_IterPerExp : public Operator
    {
    public:
        IProductWRTBase_IterPerExp(
            StdRegions::StdExpansionSharedPtr pExp,
            vector<SpatialDomains::GeometrySharedPtr> pGeom)
                : Operator(pExp, pGeom)
        {
        }

        virtual void operator()(
            const Array<OneD, const NekDouble> &input,
                  Array<OneD,       NekDouble> &output,
                  Array<OneD,       NekDouble> &wsp)
        {
            const int nCoeffs = m_stdExp->GetNcoeffs();
            const int nPhys   = m_stdExp->GetTotPoints();
            Array<OneD, NekDouble> tmp;

            for (int i = 0; i < m_numElmt; ++i)
            {
                m_stdExp->IProductWRTBase(
                    input + i*nPhys, tmp = output + i*nCoeffs);
            }
        }

        OPERATOR_CREATE(IProductWRTBase_IterPerExp)
    };

    OperatorKey IProductWRTBase_IterPerExp::m_typeArr[] =
    {
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::eQuadrilateral, eIProductWRTBase, eIterPerExp),
            IProductWRTBase_IterPerExp::create, "IProductWRTBase_IterPerExp_Quad"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::eTriangle, eIProductWRTBase, eIterPerExp),
            IProductWRTBase_IterPerExp::create, "IProductWRTBase_IterPerExp_Tri")
    };

    /*
     * ----------------------------------------------------------
     * FwdTrans operators
     * ----------------------------------------------------------
     */

    class FwdTrans_IterPerExp : public Operator
    {
    public:
        FwdTrans_IterPerExp(StdRegions::StdExpansionSharedPtr pExp,
                            vector<SpatialDomains::GeometrySharedPtr> pGeom)
            : Operator(pExp, pGeom)
        {
        }

        virtual void operator()(
            const Array<OneD, const NekDouble> &input,
                  Array<OneD,       NekDouble> &output,
                  Array<OneD,       NekDouble> &wsp)
        {
            const int nCoeffs = m_stdExp->GetNcoeffs();
            const int nPhys   = m_stdExp->GetTotPoints();
            Array<OneD, NekDouble> tmp;

            for (int i = 0; i < m_numElmt; ++i)
            {
                m_stdExp->FwdTrans(input + i*nCoeffs, tmp = output + i*nPhys);
            }
        }

        OPERATOR_CREATE(FwdTrans_IterPerExp)
    };

    OperatorKey FwdTrans_IterPerExp::m_typeArr[] =
    {
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::eQuadrilateral, eFwdTrans, eIterPerExp),
            FwdTrans_IterPerExp::create, "FwdTrans_IterPerExp_Quad"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::eTriangle, eFwdTrans, eIterPerExp),
            FwdTrans_IterPerExp::create, "FwdTrans_IterPerExp_Tri")
    };
}
}
