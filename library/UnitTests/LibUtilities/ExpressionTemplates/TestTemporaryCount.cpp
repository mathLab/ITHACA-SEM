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

#include <UnitTests/LibUtilities/ExpressionTemplates/CountedObjectExpression.h>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <ExpressionTemplates/ExpressionTemplates.hpp>
#include <ExpressionTemplates/Node.hpp>
#include <ExpressionTemplates/RemoveAllUnecessaryTemporaries.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/type_traits.hpp>
#include <LibUtilities/LinearAlgebra/MatrixSize.hpp>

namespace Nektar
{
    namespace UnitTests
    {
        BOOST_AUTO_TEST_CASE(TestConstantCount)
        {
            typedef NekMatrix<double> Matrix;

            typedef expt::Node<Matrix> ConstantNode;

            BOOST_STATIC_ASSERT( expt::TemporaryCount<ConstantNode>::Value == 0 );
        }

        BOOST_AUTO_TEST_CASE(TestUnaryCount)
        {
            typedef NekMatrix<double> Matrix;

            typedef expt::Node<Matrix> ConstantNode;
            typedef expt::Node<ConstantNode, expt::NegateOp> UnaryNode;

            BOOST_STATIC_ASSERT( expt::TemporaryCount<UnaryNode>::Value == 0 );
        }

        BOOST_AUTO_TEST_CASE(TestBinaryCount)
        {
            typedef NekMatrix<double> Matrix;

            typedef expt::Node<Matrix> ConstantNode;
            typedef expt::Node<ConstantNode, expt::MultiplyOp, ConstantNode> BinaryNode;

            BOOST_STATIC_ASSERT( expt::TemporaryCount<BinaryNode>::Value == 0 );

            typedef expt::Node<BinaryNode, expt::MultiplyOp, ConstantNode> LeftTree;
            BOOST_STATIC_ASSERT( expt::TemporaryCount<LeftTree>::Value == 0 );

            typedef expt::Node<ConstantNode, expt::MultiplyOp, BinaryNode> RightTree;
            BOOST_STATIC_ASSERT( expt::TemporaryCount<RightTree>::Value == 1 );
        }

        // Tests for RemoveUnecessaryTemporaries
        BOOST_AUTO_TEST_CASE(TestRemoveUnecessaryTemporariesConstantTree)
        {
            typedef NekMatrix<double> Matrix;
            typedef expt::Node<Matrix> ConstantNode;
            typedef ConstantNode::Indices Indices;

            typedef expt::RemoveUnecessaryTemporaries<ConstantNode>::TransformedNodeType ResultNodeType;
            typedef expt::RemoveUnecessaryTemporaries<ConstantNode>::TransformedIndicesType TransformedIndicesType;

            BOOST_STATIC_ASSERT(( boost::mpl::size<Indices>::type::value == 2 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<TransformedIndicesType, 0>::type::value == 0 ));
            
            BOOST_MPL_ASSERT(( boost::is_same<ResultNodeType, ConstantNode> ));
        }

        BOOST_AUTO_TEST_CASE(TestRemoveUnecessaryTemporariesUnaryOneLevelTree)
        {
            typedef NekMatrix<double> Matrix;

            typedef expt::Node<Matrix> ConstantNode;
            typedef expt::Node<ConstantNode, expt::NegateOp> UnaryNode;

            BOOST_STATIC_ASSERT( expt::TemporaryCount<UnaryNode>::Value == 0 );
        }

        BOOST_AUTO_TEST_CASE(TestRemoveUnecessaryTemporariesBinaryOneLevelTree)
        {
            typedef NekMatrix<double> Matrix;

            typedef expt::Node<Matrix> ConstantNode;
            typedef expt::Node<ConstantNode, expt::MultiplyOp, ConstantNode> BinaryNode;

            typedef BinaryNode::Indices Indices;

            typedef expt::RemoveUnecessaryTemporaries<BinaryNode>::TransformedNodeType ResultNodeType;
            typedef expt::RemoveUnecessaryTemporaries<BinaryNode>::TransformedIndicesType TransformedIndicesType;

            BOOST_STATIC_ASSERT(( boost::mpl::size<TransformedIndicesType>::type::value == 2 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<TransformedIndicesType, 0>::type::value == 0 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<TransformedIndicesType, 1>::type::value == 1 ));
            
            BOOST_MPL_ASSERT(( boost::is_same<ResultNodeType, BinaryNode> ));

        }

        BOOST_AUTO_TEST_CASE(TestRemoveUnecessaryTemporariesBinaryTreeTest4)
        {
            typedef NekMatrix<double> Matrix;

            typedef expt::Node<Matrix> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> BinaryNode;
            typedef expt::Node<ConstantNode, expt::AddOp, BinaryNode> NodeType;
            typedef NodeType::Indices Indices;

            typedef expt::Node<BinaryNode, expt::AddOp, ConstantNode> ExpectedNodeType;

            typedef expt::RemoveUnecessaryTemporaries<NodeType>::TransformedNodeType ResultNodeType;
            typedef expt::RemoveUnecessaryTemporaries<NodeType>::TransformedIndicesType TransformedIndicesType;

            BOOST_STATIC_ASSERT(( boost::mpl::size<TransformedIndicesType>::type::value == 3 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<TransformedIndicesType, 0>::type::value == 0 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<TransformedIndicesType, 1>::type::value == 1 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<TransformedIndicesType, 2>::type::value == 2 ));
            
            BOOST_MPL_ASSERT(( boost::is_same<ExpectedNodeType, ResultNodeType> ));

        }

        BOOST_AUTO_TEST_CASE(TestRemoveUnecessaryTemporariesBinaryTreeTest5)
        {
            typedef NekMatrix<double> Matrix;

            typedef expt::Node<Matrix> ConstantNode;
            typedef expt::Node<ConstantNode, expt::MultiplyOp, ConstantNode> BinaryNode;
            typedef expt::Node<ConstantNode, expt::AddOp, BinaryNode> NodeType;
            typedef NodeType::Indices Indices;

            typedef expt::Node<BinaryNode, expt::AddOp, ConstantNode> ExpectedNodeType;

            typedef expt::RemoveUnecessaryTemporaries<NodeType>::TransformedNodeType ResultNodeType;
            typedef expt::RemoveUnecessaryTemporaries<NodeType>::TransformedIndicesType TransformedIndicesType;

            typedef expt::impl::RemoveUnecessaryTemporariesInternal<NodeType, Indices, 0>::Tree0Type Tree0Type;
            typedef expt::impl::RemoveUnecessaryTemporariesInternal<NodeType, Indices, 0>::Tree1Type Tree1Type;
            typedef expt::impl::RemoveUnecessaryTemporariesInternal<NodeType, Indices, 0>::Tree2Type Tree2Type;

            BOOST_MPL_ASSERT(( boost::is_same<Tree0Type, NodeType> ));
            BOOST_MPL_ASSERT(( boost::is_same<Tree1Type, NodeType> ));
            BOOST_MPL_ASSERT(( boost::is_same<Tree2Type, NodeType> ));

            BOOST_STATIC_ASSERT(( boost::mpl::size<TransformedIndicesType>::type::value == 3 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<TransformedIndicesType, 0>::type::value == 1 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<TransformedIndicesType, 1>::type::value == 2 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<TransformedIndicesType, 2>::type::value == 0 ));
            
            BOOST_MPL_ASSERT(( boost::is_same<ExpectedNodeType, ResultNodeType> ));

        }

        BOOST_AUTO_TEST_CASE(TestRemoveUnecessaryTemporariesBinaryTreeTest5_1)
        {
            //// (A+B) + (BC) -> ((BC) + A) + B
            //typedef NekMatrix<double> Matrix;

            //typedef expt::Node<Matrix> ConstantNode;
            //typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> AddNode;
            //typedef expt::Node<ConstantNode, expt::MultiplyOp, ConstantNode> MultiplyNode;

            //typedef expt::Node<AddNode, expt::AddOp, MultiplyNode> NodeType;
            //typedef NodeType::Indices Indices;

            //typedef expt::Node<BinaryNode, expt::AddOp, ConstantNode> ExpectedNodeType;

            //typedef expt::RemoveUnecessaryTemporaries<NodeType, Indices, 0>::TransformedNodeType ResultNodeType;
            //typedef expt::RemoveUnecessaryTemporaries<NodeType, Indices, 0>::TransformedIndicesType TransformedIndicesType;

            //typedef expt::RemoveUnecessaryTemporaries<NodeType, Indices, 0>::Tree0Type Tree0Type;
            //typedef expt::RemoveUnecessaryTemporaries<NodeType, Indices, 0>::Tree1Type Tree1Type;
            //typedef expt::RemoveUnecessaryTemporaries<NodeType, Indices, 0>::Tree2Type Tree2Type;

            //BOOST_MPL_ASSERT(( boost::is_same<Tree0Type, NodeType> ));
            //BOOST_MPL_ASSERT(( boost::is_same<Tree1Type, NodeType> ));
            //BOOST_MPL_ASSERT(( boost::is_same<Tree2Type, NodeType> ));

            //BOOST_STATIC_ASSERT(( boost::mpl::size<TransformedIndicesType>::type::value == 3 ));
            //BOOST_STATIC_ASSERT(( boost::mpl::at_c<TransformedIndicesType, 0>::type::value == 1 ));
            //BOOST_STATIC_ASSERT(( boost::mpl::at_c<TransformedIndicesType, 1>::type::value == 2 ));
            //BOOST_STATIC_ASSERT(( boost::mpl::at_c<TransformedIndicesType, 2>::type::value == 0 ));
            //
            //BOOST_MPL_ASSERT(( boost::is_same<ExpectedNodeType, ResultNodeType> ));

        }

        BOOST_AUTO_TEST_CASE(TestRemoveUnecessaryTemporariesBinaryTreeTest6)
        {
            typedef NekMatrix<double> Matrix;

            // (A+B) + (C+DE) -> ((A+B) + DE) + C
            typedef expt::Node<Matrix> ConstantNode;
            typedef expt::Node<ConstantNode, expt::MultiplyOp, ConstantNode> DENode;
            typedef expt::Node<ConstantNode, expt::AddOp, DENode> CDENode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> ABNode;
            typedef expt::Node<ABNode, expt::AddOp, CDENode> NodeType;
            typedef NodeType::Indices Indices;

            typedef expt::Node<DENode, expt::AddOp, ConstantNode> B0;
            typedef expt::Node<B0, expt::AddOp, ConstantNode> B1;
            typedef expt::Node<B1, expt::AddOp, ConstantNode> ExpectedNodeType;
            

            // Right node goes from (C+DE) -> (DE+C) after optimizing rhs.
            typedef expt::Node<DENode, expt::AddOp, ConstantNode> ExpectedRightNode0Type;
            typedef expt::impl::RemoveUnecessaryTemporariesInternal<NodeType, Indices, 0>::RightNode0Type RightNode0Type;
            BOOST_MPL_ASSERT(( boost::is_same<ExpectedRightNode0Type, RightNode0Type> ));
            
            typedef expt::Node<ABNode, expt::AddOp, ExpectedRightNode0Type> ExpectedTree0Type;
            typedef expt::impl::RemoveUnecessaryTemporariesInternal<NodeType, Indices, 0>::Tree0Type Tree0Type;
            BOOST_MPL_ASSERT(( boost::is_same<ExpectedTree0Type, Tree0Type> ));

            // Tree goes from (A+B) + (DE+C) -> ((A+B)+DE) + C after the associative transform.
            typedef expt::impl::RemoveUnecessaryTemporariesInternal<NodeType, Indices, 0>::Tree1Type Tree1Type;
            typedef Tree1Type::Right Tree1TypeRight;
            BOOST_MPL_ASSERT(( boost::is_same<ConstantNode, Tree1TypeRight> ));

            typedef Tree1Type::Left Tree1TypeLeft;
            typedef expt::Node<ABNode, expt::AddOp, DENode> ExpectedTree1TypeLeft;
            BOOST_MPL_ASSERT(( boost::is_same<ExpectedTree1TypeLeft, Tree1TypeLeft> ));

            // Optimzing the lhs causes no change
            typedef expt::impl::RemoveUnecessaryTemporariesInternal<NodeType, Indices, 0>::LeftNode2Type LeftNode2Type;
            typedef expt::Node<ABNode, expt::AddOp, DENode> ExpectedTree2TypeLeft;
            BOOST_MPL_ASSERT(( boost::is_same<ExpectedTree2TypeLeft, LeftNode2Type> ));


            typedef expt::RemoveUnecessaryTemporaries<NodeType>::TransformedNodeType ResultNodeType;
            typedef expt::RemoveUnecessaryTemporaries<NodeType>::TransformedIndicesType TransformedIndicesType;

            BOOST_STATIC_ASSERT(( boost::mpl::size<TransformedIndicesType>::type::value == 5 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<TransformedIndicesType, 0>::type::value == 3 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<TransformedIndicesType, 1>::type::value == 4 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<TransformedIndicesType, 2>::type::value == 1 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<TransformedIndicesType, 3>::type::value == 0 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<TransformedIndicesType, 4>::type::value == 2 ));
            
            BOOST_MPL_ASSERT(( boost::is_same<ExpectedNodeType, ResultNodeType> ));

        }
    }
}
