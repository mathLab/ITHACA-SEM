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
#include <UnitTests/LibUtilities/ExpressionTemplates/ExpressionTemplateObjects.h>

#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <ExpressionTemplates/ExpressionTemplates.hpp>
#include <ExpressionTemplates/Node.hpp>
#include <ExpressionTemplates/RemoveAllUnecessaryTemporaries.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/type_traits.hpp>
#include <LibUtilities/LinearAlgebra/MatrixSize.hpp>

using namespace expt;

namespace Nektar
{
    namespace UnitTests
    {
        typedef NekMatrix<double> Matrix;
        BOOST_AUTO_TEST_CASE(TestConstantAndUnaryNodes)
        {
            typedef Node<TestObjectA> ConstantNode;
            typedef Node<ConstantNode, NegateOp> UnaryNode;

            typedef AssociativeTransform<ConstantNode>::TransformedNodeType ConstantResultNode;
            typedef AssociativeTransform<UnaryNode>::TransformedNodeType UnaryResultNode;

            BOOST_MPL_ASSERT(( boost::is_same<ConstantResultNode, ConstantNode> ));
            BOOST_MPL_ASSERT(( boost::is_same<UnaryResultNode, UnaryNode> ));
        }

        // A+B -> A+B
        BOOST_AUTO_TEST_CASE(TestAssociativeTransformAPlusB)
        {
            typedef Node<TestObjectA> LeafNode;

            typedef Node<LeafNode, AddOp, LeafNode > NodeType;

            typedef AssociativeTransform<NodeType>::TransformedNodeType ResultType;
            typedef InverseAssociativeTransform<NodeType>::TransformedNodeType InverseResultType;

            BOOST_MPL_ASSERT(( boost::is_same<NodeType, ResultType>));
            BOOST_MPL_ASSERT(( boost::is_same<NodeType, InverseResultType>));
        }

        // (A+B)+C -> (A+B)+C
        // A+(B+C) -> (A+B)+C
        BOOST_AUTO_TEST_CASE(TestAssociativeTransformThreeNodes)
        {
            typedef Node<TestObjectA> LeafNode;
            typedef Node<LeafNode, AddOp, LeafNode> AddNode;
            typedef Node<AddNode, AddOp, LeafNode> Test1Node;
            typedef Node<LeafNode, AddOp, AddNode> Test2Node;

            typedef AssociativeTransform<Test1Node>::TransformedNodeType Test1ResultNode;
            typedef AssociativeTransform<Test2Node>::TransformedNodeType Test2ResultNode;

            typedef InverseAssociativeTransform<Test1Node>::TransformedNodeType Test1InverseResultNode;
            typedef InverseAssociativeTransform<Test2Node>::TransformedNodeType Test2InverseResultNode;

            typedef Test1Node ExpectedResult1Node;
            typedef Test1Node ExpectedResult2Node;

            typedef Test2Node ExpectedInverseResult1Node;
            typedef Test2Node ExpectedInverseResult2Node;

            BOOST_MPL_ASSERT(( boost::is_same<Test1ResultNode, ExpectedResult1Node> ));
            BOOST_MPL_ASSERT(( boost::is_same<Test2ResultNode, ExpectedResult2Node> ));
            BOOST_MPL_ASSERT(( boost::is_same<Test1InverseResultNode, ExpectedInverseResult1Node> ));
            BOOST_MPL_ASSERT(( boost::is_same<Test2InverseResultNode, ExpectedInverseResult2Node> ));
        }

        BOOST_AUTO_TEST_CASE(TestAssociativeTransformConstantAndUnaryNodes)
        {
            typedef Node<TestObjectA> T;
            BOOST_MPL_ASSERT(( boost::is_same<AssociativeTransform<T>::TransformedNodeType, T >));
            BOOST_MPL_ASSERT(( boost::is_same<InverseAssociativeTransform<T>::TransformedNodeType, T >));

            typedef Node<T, NegateOp> U;
            BOOST_MPL_ASSERT(( boost::is_same<AssociativeTransform<U>::TransformedNodeType, U >));
            BOOST_MPL_ASSERT(( boost::is_same<InverseAssociativeTransform<T>::TransformedNodeType, T >));
        }

        BOOST_AUTO_TEST_CASE(TestAssociativeTransformAPlusBPlusCPlusD)
        {
            typedef NekMatrix<double> Matrix;

            // (A+B) + (C+D) -> ((A+B)+C)+D
            typedef Node<Matrix> ConstantNode;
            typedef Node<ConstantNode, AddOp, ConstantNode> BinaryNode;
            typedef Node<BinaryNode, AddOp, BinaryNode> T0;
            
            typedef Node<BinaryNode, AddOp, ConstantNode> ExpectedLhs;
            typedef Node<ExpectedLhs, AddOp, ConstantNode> T1;

            
            BOOST_MPL_ASSERT(( boost::is_same<AssociativeTransform<T0>::TransformedNodeType, T1 >));
            BOOST_MPL_ASSERT(( boost::is_same<AssociativeTransform<T1>::TransformedNodeType, T1 >));

            // ((A+B)+C)+D -> (A+B) + (C+D)
            BOOST_MPL_ASSERT(( boost::is_same<InverseAssociativeTransform<T1>::TransformedNodeType, T0 >));
            
            // (A+B) + (C+D) -> (A + (B + (C+D))
            typedef Node<ConstantNode, AddOp, BinaryNode> ExpectedRhs;
            typedef Node<ConstantNode, AddOp, ExpectedRhs> T2;
            BOOST_MPL_ASSERT(( boost::is_same<InverseAssociativeTransform<T0>::TransformedNodeType, T2 >));
        }

        // (A+B)+C -> (A+B)+C
        // A+(B+C) -> A+(B+C)
        BOOST_AUTO_TEST_CASE(TestNonAssociativeOperator)
        {
            typedef Node<TestObject> LeafNode;
            typedef Node<LeafNode, AddOp, LeafNode> AddNode;
            typedef Node<AddNode, AddOp, LeafNode> Test1Node;
            typedef Node<LeafNode, AddOp, AddNode> Test2Node;

            typedef AssociativeTransform<Test1Node>::TransformedNodeType Test1ResultNode;
            typedef AssociativeTransform<Test2Node>::TransformedNodeType Test2ResultNode;

            typedef InverseAssociativeTransform<Test1Node>::TransformedNodeType Test1InverseResultNode;
            typedef InverseAssociativeTransform<Test2Node>::TransformedNodeType Test2InverseResultNode;

            typedef Test1Node ExpectedResult1Node;
            typedef Test2Node ExpectedResult2Node;

            typedef Test1Node ExpectedInverseResult1Node;
            typedef Test2Node ExpectedInverseResult2Node;

            BOOST_MPL_ASSERT(( boost::is_same<Test1ResultNode, ExpectedResult1Node> ));
            BOOST_MPL_ASSERT(( boost::is_same<Test2ResultNode, ExpectedResult2Node> ));
            BOOST_MPL_ASSERT(( boost::is_same<Test1InverseResultNode, ExpectedInverseResult1Node> ));
            BOOST_MPL_ASSERT(( boost::is_same<Test2InverseResultNode, ExpectedInverseResult2Node> ));
        }
    }
}
