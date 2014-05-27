///////////////////////////////////////////////////////////////////////////////
//
// File: QuadOperators.cpp
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
// Description: Operators specific to the quadrilateral
//
///////////////////////////////////////////////////////////////////////////////

#include <Collections/Operator.h>

namespace Nektar {
namespace Collections {

    class QuadBwdTransLocMat : public Operator
    {
    public:
        QuadBwdTransLocMat(StdRegions::StdExpansionSharedPtr pExp,
                           vector<SpatialDomains::GeometrySharedPtr> pGeom)
            : Operator(pExp, pGeom)
        {
        }

        virtual void operator()(const Array<OneD, const NekDouble> &input,
                                      Array<OneD,       NekDouble> &output)
        {
            StdRegions::StdMatrixKey key(StdRegions::eBwdTrans,
                                         LibUtilities::eQuadrilateral,
                                         *m_stdExp);
            DNekMatSharedPtr mat = m_stdExp->GetStdMatrix(key);

            Blas::Dgemm('N', 'N', mat->GetRows(), m_numElmt, mat->GetColumns(),
                        1.0, mat->GetRawPtr(), mat->GetRows(), input.get(),
                        m_stdExp->GetNcoeffs(), 0.0, output.get(),
                        m_stdExp->GetTotPoints());
        }

        static OperatorKey m_type;
        static OperatorSharedPtr create(
            StdRegions::StdExpansionSharedPtr pExp,
            vector<SpatialDomains::GeometrySharedPtr> pGeom)
        {
            return MemoryManager<QuadBwdTransLocMat>
                ::AllocateSharedPtr(pExp, pGeom);
        }
    };

    OperatorKey QuadBwdTransLocMat::m_type = GetOperatorFactory().
        RegisterCreatorFunction(
            OperatorKey(LibUtilities::eQuadrilateral, eBwdTrans, eLocMat),
            QuadBwdTransLocMat::create, "QuadBwdTransLocMat");
}
}
