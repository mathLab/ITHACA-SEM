///////////////////////////////////////////////////////////////////////////////
//
// File: Collection.h
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
// Description: Collection top class definition
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBRARY_COLLECTIONS_COLLECTION_H
#define NEKTAR_LIBRARY_COLLECTIONS_COLLECTION_H

#include <Collections/Operator.h>

#include <boost/unordered_map.hpp>

namespace Nektar {
namespace Collections {
    /**
     * @brief Collection
     */
    class Collection
    {
    public:
        Collection(StdRegions::StdExpansionSharedPtr pExp,
                   vector<SpatialDomains::GeometrySharedPtr> pGeom)
            : m_stdExp(pExp), m_geom(pGeom)
        {
            OperatorKey bwdLocMat(
                LibUtilities::eQuadrilateral, eBwdTrans, eLocMat);
            OperatorKey derivSumFac(
                LibUtilities::eQuadrilateral, ePhysDeriv, eSumFac);
            m_ops[eBwdTrans] = GetOperatorFactory().CreateInstance(
                bwdLocMat, pExp, pGeom);
            m_ops[ePhysDeriv] = GetOperatorFactory().CreateInstance(
                derivSumFac, pExp, pGeom);
        }

        void ApplyOperator(
            const OperatorType                 &op,
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray)
        {
            Array<OneD, NekDouble> wsp(m_ops[op]->GetWspSize());
            (*m_ops[op])(inarray, outarray, wsp);
        }

    protected:
        StdRegions::StdExpansionSharedPtr                     m_stdExp;
        vector<SpatialDomains::GeometrySharedPtr>             m_geom;
        boost::unordered_map<OperatorType, OperatorSharedPtr> m_ops;
    };
}
}

#endif
