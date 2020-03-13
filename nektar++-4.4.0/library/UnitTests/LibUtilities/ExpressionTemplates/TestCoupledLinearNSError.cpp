///////////////////////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_USE_EXPRESSION_TEMPLATES
#define NEKTAR_USE_EXPRESSION_TEMPLATES
#endif //NEKTAR_USE_EXPRESSION_TEMPLATES

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>
#include <UnitTests/LibUtilities/ExpressionTemplates/CountedObjectExpression.h>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <ExpressionTemplates/Node.hpp>
#include <ExpressionTemplates/RemoveAllUnecessaryTemporaries.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/type_traits.hpp>
#include <LibUtilities/LinearAlgebra/MatrixSize.hpp>

// Tests the error in CoupledLinearNS.cpp:963
// DNekScalBlkMatSharedPtr  m_Cinv; 
// NekVector< NekDouble > F_int
// DNekScalBlkMatSharedPtr  m_D_int; 
// NekVector<  NekDouble > F_p
// DNekScalBlkMatSharedPtr m_Btilde;
// NekVector< NekDouble > F_bnd
// F_int = (*m_Cinv)*(F_int + Transpose(*m_D_int)*F_p - Transpose(*m_Btilde)*F_bnd);

namespace Nektar
{
    namespace UnitTests
    {
        // Start my testing individual portions to help isolate the problem.
        BOOST_AUTO_TEST_CASE(TestScalBlkMatTimeVector)
        {
            // 2x2 block matrix, with submatrices also 2x2
            double data_buf[] = {1, 2, 3, 4};
            DNekMatSharedPtr m0(new DNekMat(2,2, data_buf));
            DNekMatSharedPtr m1(new DNekMat(2,2, data_buf));
            DNekMatSharedPtr m2(new DNekMat(2,2, data_buf));
            DNekMatSharedPtr m3(new DNekMat(2,2, data_buf));

            DNekScalMatSharedPtr s0(new DNekScalMat(2.0, m0));
            DNekScalMatSharedPtr s1(new DNekScalMat(3.0, m1));
            DNekScalMatSharedPtr s2(new DNekScalMat(4.0, m2));
            DNekScalMatSharedPtr s3(new DNekScalMat(5.0, m3));

            DNekScalBlkMatSharedPtr m(new DNekScalBlkMat(2,2, 2, 2));
            m->SetBlock(0, 0, s0);
            m->SetBlock(0, 1, s1);
            m->SetBlock(1, 0, s2);
            m->SetBlock(1, 1, s3);

            double v_buf[] = {1, 2, 3, 4, 5, 6, 7, 8};
            NekVector<NekDouble> v(8, v_buf);

            NekVector<NekDouble> result = (*m)*v;
        }

        // The multiplication compiles, so I assume the problem is in a tree transformation.
        BOOST_AUTO_TEST_CASE(TestScaleBlkMatTermWithTransformation)
        {
            double data_buf[] = {1, 2, 3, 4};
            DNekMatSharedPtr m0(new DNekMat(2,2, data_buf));
            DNekMatSharedPtr m1(new DNekMat(2,2, data_buf));
            DNekMatSharedPtr m2(new DNekMat(2,2, data_buf));
            DNekMatSharedPtr m3(new DNekMat(2,2, data_buf));

            DNekScalMatSharedPtr s0(new DNekScalMat(2.0, m0));
            DNekScalMatSharedPtr s1(new DNekScalMat(3.0, m1));
            DNekScalMatSharedPtr s2(new DNekScalMat(4.0, m2));
            DNekScalMatSharedPtr s3(new DNekScalMat(5.0, m3));

            DNekScalBlkMatSharedPtr m(new DNekScalBlkMat(2,2, 2, 2));
            m->SetBlock(0, 0, s0);
            m->SetBlock(0, 1, s1);
            m->SetBlock(1, 0, s2);
            m->SetBlock(1, 1, s3);

            double rhs_buf[] = {1, 2, 3, 4};
            NekVector<NekDouble> rhs(4, rhs_buf);

            double lhs_buf[] = {1, 2, 3, 4};
            NekVector<NekDouble> lhs(4, lhs_buf);

            NekVector<NekDouble> result = lhs + (*m)*rhs;

            DNekMat result2 = (*m)*((*m)*(*m));

            // This fails.  It appears that it is a failure of the commutative property
            // for matrix/vector multiplication.
            NekVector<NekDouble> result3 = (*m)*(lhs+rhs);

            typedef expt::Node<expt::Node<Nektar::NekMatrix<Nektar::NekMatrix<Nektar::NekMatrix<Nektar::NekDouble,Nektar::StandardMatrixTag>,Nektar::ScaledMatrixTag>,Nektar::BlockMatrixTag>,void,void>,expt::MultiplyOp,expt::Node<expt::Node<Nektar::NekVector<Nektar::NekDouble>,void,void>,expt::AddOp,expt::Node<Nektar::NekVector<Nektar::NekDouble>,void,void> > > TreeType;
            typedef expt::CreateVectorC<int, boost::mpl::int_, 3>::type Indices;

            typedef expt::RemoveUnecessaryTemporaries<TreeType>::TransformedNodeType OptimizedParseTree;

            BOOST_MPL_ASSERT(( boost::is_same<expt::impl::RemoveUnecessaryTemporariesInternal<TreeType, Indices, 0>::RightNode0Type, expt::Node<expt::Node<Nektar::NekVector<Nektar::NekDouble>,void,void>,expt::AddOp,expt::Node<Nektar::NekVector<Nektar::NekDouble>,void,void> > > ));
            BOOST_MPL_ASSERT(( boost::is_same<expt::impl::RemoveUnecessaryTemporariesInternal<TreeType, Indices, 0>::Tree0Type, TreeType > ));
            BOOST_MPL_ASSERT(( boost::is_same<expt::impl::RemoveUnecessaryTemporariesInternal<TreeType, Indices, 0>::Tree1Type, TreeType > ));
            BOOST_MPL_ASSERT(( boost::is_same<expt::impl::RemoveUnecessaryTemporariesInternal<TreeType, Indices, 0>::Tree2Type, TreeType > ));
            BOOST_MPL_ASSERT(( boost::is_same<expt::impl::RemoveUnecessaryTemporariesInternal<TreeType, Indices, 0>::Tree3Type, TreeType > ));
            BOOST_MPL_ASSERT(( boost::is_same<TreeType, OptimizedParseTree> ));

            NekVector<NekDouble> result1 = (*m)*(lhs + (*m)*rhs - (*m)*rhs);
        }
        
    }
}
