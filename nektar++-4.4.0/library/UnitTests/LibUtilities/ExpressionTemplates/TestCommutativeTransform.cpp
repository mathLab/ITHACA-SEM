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


#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "ExpressionTemplateObjects.h"
#include <ExpressionTemplates/ExpressionTemplates.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>

using namespace expt;

namespace Nektar
{
    namespace UnitTests
    {
        BOOST_AUTO_TEST_CASE(TestCommutativeTraits)
        {
            BOOST_MPL_ASSERT(( CommutativeTraits<TestObjectC, AddOp, TestObjectC> ));
            BOOST_MPL_ASSERT(( CommutativeTraits<TestObjectAC, AddOp, TestObjectC> ));

            BOOST_MPL_ASSERT(( boost::mpl::not_<CommutativeTraits<TestObjectA, AddOp, TestObjectA> > ));
            BOOST_MPL_ASSERT(( boost::mpl::not_<CommutativeTraits<TestObject, AddOp, TestObject> > ));
        }

        // A + B -> B + A
        BOOST_AUTO_TEST_CASE(TestSimpleCommutativeTransform)
        {
            typedef Node<TestObjectC> M;
            typedef Node<M, AddOp, M> AddNode;

            TestObjectC A(2.7);
            TestObjectC B(3.2);
            AddNode expression = A+B;

            // Verify initial indexes for the parse tree.
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<AddNode::Indices, 0>::type::value == 0));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<AddNode::Indices, 1>::type::value == 1));

            BOOST_CHECK_EQUAL(A, boost::fusion::at_c<0>(expression.GetData()));
            BOOST_CHECK_EQUAL(B, boost::fusion::at_c<1>(expression.GetData()));

            typedef CommutativeTransform<AddNode, AddNode::Indices, 0>::TransformedNodeType TransformedNodeType;
            typedef CommutativeTransform<AddNode, AddNode::Indices, 0>::TransformedIndicesType TransformedIndicesType;

            // The tree itself will look unchanged, but the indices should be swapped.
            BOOST_MPL_ASSERT(( boost::is_same<AddNode, TransformedNodeType> ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<TransformedIndicesType, 0>::type::value == 1));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<TransformedIndicesType, 1>::type::value == 0));

            TestObjectC result = expression;
            BOOST_CHECK_EQUAL(result.value, 2.7+3.2);
        }

        // A - B -> A - B
        BOOST_AUTO_TEST_CASE(TestSimpleNonCommutativeTransform)
        {
            typedef Node<TestObjectA> M;
            typedef Node<M, AddOp, M> AddNode;

            TestObjectA A(2.7);
            TestObjectA B(3.2);
            AddNode expression = A+B;

            // Verify initial indexes for the parse tree.
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<AddNode::Indices, 0>::type::value == 0));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<AddNode::Indices, 1>::type::value == 1));

            BOOST_CHECK_EQUAL(A, boost::fusion::at_c<0>(expression.GetData()));
            BOOST_CHECK_EQUAL(B, boost::fusion::at_c<1>(expression.GetData()));

            typedef CommutativeTransform<AddNode, AddNode::Indices, 0>::TransformedNodeType TransformedNodeType;
            typedef CommutativeTransform<AddNode, AddNode::Indices, 0>::TransformedIndicesType TransformedIndicesType;

            // The tree itself will look unchanged, but the indices should be swapped.
            BOOST_MPL_ASSERT(( boost::is_same<AddNode, TransformedNodeType> ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<TransformedIndicesType, 0>::type::value == 0));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<TransformedIndicesType, 1>::type::value == 1));

            TestObjectA result = expression;
            BOOST_CHECK_EQUAL(result.value, 2.7+3.2);
        }

        BOOST_AUTO_TEST_CASE(TestSubtreeCommutativeTransform)
        {
            // Expression (A+B) + ( (C+D) + (E+F) )
            typedef Node<TestObjectC> M;
            typedef Node<M, AddOp, M> LhsNode;
            typedef Node<LhsNode, AddOp, LhsNode> RhsNode;
            typedef Node<LhsNode, AddOp, RhsNode> TestNode;

            BOOST_STATIC_ASSERT(( boost::mpl::at_c<TestNode::Indices, 0>::type::value == 0));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<TestNode::Indices, 1>::type::value == 1));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<TestNode::Indices, 2>::type::value == 2));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<TestNode::Indices, 3>::type::value == 3));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<TestNode::Indices, 4>::type::value == 4));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<TestNode::Indices, 5>::type::value == 5));

            // Apply to the lhs without touching the rest of the tree.
            // (A+B) + ( (C+D) + (E+F) ) -> (B+A) + ( (C+D) + (E+F) )
            typedef CommutativeTransform<LhsNode, TestNode::Indices, 0>::TransformedNodeType Test0NodeType;
            typedef CommutativeTransform<LhsNode, TestNode::Indices, 0>::TransformedIndicesType Test0IndicesType;
            BOOST_MPL_ASSERT(( boost::is_same<LhsNode, Test0NodeType> ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test0IndicesType, 0>::type::value == 1));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test0IndicesType, 1>::type::value == 0));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test0IndicesType, 2>::type::value == 2));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test0IndicesType, 3>::type::value == 3));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test0IndicesType, 4>::type::value == 4));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test0IndicesType, 5>::type::value == 5));

            // Apply to the far rhs node.
            // (B+A) + ( (C+D) + (E+F) ) -> (B+A) + ( (E+F) + (C+D) )
            typedef CommutativeTransform<RhsNode, Test0IndicesType, 2>::TransformedNodeType Test1NodeType;
            typedef CommutativeTransform<RhsNode, Test0IndicesType, 2>::TransformedIndicesType Test1IndicesType;
            BOOST_MPL_ASSERT(( boost::is_same<RhsNode, Test1NodeType> ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test1IndicesType, 0>::type::value == 1));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test1IndicesType, 1>::type::value == 0));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test1IndicesType, 2>::type::value == 4));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test1IndicesType, 3>::type::value == 5));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test1IndicesType, 4>::type::value == 2));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test1IndicesType, 5>::type::value == 3));

            // Test transform on the end.
            // (B+A) + ( (E+F) + (C+D) ) -> (B+A) + ( (E+F) + (D+C) )
            typedef CommutativeTransform<LhsNode, Test1IndicesType, 4>::TransformedNodeType Test2NodeType;
            typedef CommutativeTransform<LhsNode, Test1IndicesType, 4>::TransformedIndicesType Test2IndicesType;
            BOOST_MPL_ASSERT(( boost::is_same<LhsNode, Test2NodeType> ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test2IndicesType, 0>::type::value == 1));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test2IndicesType, 1>::type::value == 0));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test2IndicesType, 2>::type::value == 4));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test2IndicesType, 3>::type::value == 5));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test2IndicesType, 4>::type::value == 3));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test2IndicesType, 5>::type::value == 2));


            // All previous tests swapped equal numbers of elements on the lhs and rhs.
            // Test differing sizes.
            // (B+A) + ( (E+F) + (D+C) ) -> ( (E+F) + (D+C) ) + (B+A) 
            typedef CommutativeTransform<TestNode, Test2IndicesType, 0>::TransformedNodeType Test3NodeType;
            typedef CommutativeTransform<TestNode, Test2IndicesType, 0>::TransformedIndicesType Test3IndicesType;
            BOOST_MPL_ASSERT(( boost::is_same<Test3NodeType, Node<RhsNode, AddOp, LhsNode> > ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test3IndicesType, 0>::type::value == 4));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test3IndicesType, 1>::type::value == 5));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test3IndicesType, 2>::type::value == 3));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test3IndicesType, 3>::type::value == 2));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test3IndicesType, 4>::type::value == 1));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test3IndicesType, 5>::type::value == 0));
        }
    }
}
