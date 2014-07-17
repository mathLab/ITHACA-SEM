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
              m_key(StdRegions::eBwdTrans, pExp->GetShapeType(), *pExp)
        {

        }

        virtual void operator()(
            const Array<OneD, const NekDouble> &input,
                  Array<OneD,       NekDouble> &output,
                  Array<OneD,       NekDouble> &wsp)
        {
            DNekMatSharedPtr mat = m_stdExp->GetStdMatrix(m_key);

            Blas::Dgemm('N', 'N', mat->GetRows(), m_numElmt, mat->GetColumns(),
                        1.0, mat->GetRawPtr(), mat->GetRows(), input.get(),
                        m_stdExp->GetNcoeffs(), 0.0, output.get(),
                        m_stdExp->GetTotPoints());
        }

        OPERATOR_CREATE(BwdTrans_LocMat)

        StdRegions::StdMatrixKey m_key;
    };

    OperatorKey BwdTrans_LocMat::m_typeArr[] =
    {
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::eQuadrilateral, eBwdTrans, eLocMat),
            BwdTrans_LocMat::create, "BwdTrans_LocMat_Quad"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::eTriangle, eBwdTrans, eLocMat),
            BwdTrans_LocMat::create, "BwdTrans_LocMat_Tri")
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
            OperatorKey(LibUtilities::eQuadrilateral, eBwdTrans, eIterPerExp),
            BwdTrans_IterPerExp::create, "BwdTrans_IterPerExp_Quad"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::eTriangle, eBwdTrans, eIterPerExp),
            BwdTrans_IterPerExp::create, "BwdTrans_IterPerExp_Tri")
    };
}
}
