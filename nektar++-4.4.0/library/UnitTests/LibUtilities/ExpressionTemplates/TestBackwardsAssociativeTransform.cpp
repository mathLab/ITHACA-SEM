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

#include <ExpressionTemplates/ExpressionTemplates.hpp>
#include <LibUtilities/LinearAlgebra/NekVector.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>

using namespace expt;

namespace Nektar
{
    namespace UnitTests
    {
        BOOST_AUTO_TEST_CASE(TestConstantNode)
        {
            typedef Node<NekVector<double> > LeafNode;
            typedef BackwardInverseTransform<LeafNode, LeafNode::Indices, 0>::TransformedNodeType Type;
            BOOST_MPL_ASSERT(( boost::is_same<LeafNode, Type> ));
        }

        struct OperatorWithNoInverse
        {
            template<typename T, typename R>
            struct ResultType
            {
                typedef T type;
            };
        };

        BOOST_AUTO_TEST_CASE(Test120)
        {
            typedef Node<NekVector<double> > LeafNode;
            typedef Node<LeafNode, NegateOp> NegateNode;
            typedef Node<LeafNode, AddOp, NegateNode> Test0;
            typedef Node<LeafNode, SubtractOp, NegateNode> Test1;
            typedef Node<LeafNode, OperatorWithNoInverse, LeafNode> Test2;

            typedef BackwardInverseTransform<Test0, Test0::Indices, 0>::TransformedNodeType Test0Result;
            typedef BackwardInverseTransform<Test1, Test1::Indices, 0>::TransformedNodeType Test1Result;
            typedef BackwardInverseTransform<Test2, Test2::Indices, 0>::TransformedNodeType Test2Result;
            typedef BackwardInverseTransform<LeafNode, LeafNode::Indices, 0>::TransformedNodeType Test3Result;
            typedef BackwardInverseTransform<NegateNode, NegateNode::Indices, 0>::TransformedNodeType Test4Result;
            
            typedef Node<LeafNode, SubtractOp, LeafNode> ExpectedTest0Result;
            typedef Node<LeafNode, AddOp, LeafNode> ExpectedTest1Result;
            typedef Test2 ExpectedTest2Result;
            typedef LeafNode ExpectedTest3Result;
            typedef NegateNode ExpectedTest4Result;

            BOOST_MPL_ASSERT(( boost::is_same<ExpectedTest0Result, Test0Result> ));
            BOOST_MPL_ASSERT(( boost::is_same<ExpectedTest1Result, Test1Result> ));
            BOOST_MPL_ASSERT(( boost::is_same<ExpectedTest2Result, Test2Result> ));
            BOOST_MPL_ASSERT(( boost::is_same<ExpectedTest3Result, Test3Result> ));
            BOOST_MPL_ASSERT(( boost::is_same<ExpectedTest4Result, Test4Result> ));
            
        }
    }
}
