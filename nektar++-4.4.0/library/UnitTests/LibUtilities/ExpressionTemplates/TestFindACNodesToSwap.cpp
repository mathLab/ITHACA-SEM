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

#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <LibUtilities/LinearAlgebra/NekVector.hpp>
#include <ExpressionTemplates/ExpressionTemplates.hpp>

using namespace expt;

namespace Nektar
{
    namespace UnitTests
    {
        BOOST_AUTO_TEST_CASE(TestFindingLeftNode)
        {
            typedef Node<NekVector<double> > V;
            typedef Node<V, AddOp, V> AddNode;
            typedef Node<V, MultiplyOp, V> MultiplyNode;
            typedef Node<V, NegateOp> NegateNode;

            typedef AddNode Test1;
            typedef MultiplyNode Test2;
            typedef Node<NegateNode, AddOp, AddNode> Test3;
            typedef Node<AddNode, AddOp, AddNode> Test4;
            typedef Node<MultiplyNode, AddOp, AddNode> Test5;

            typedef FindFurthestLeftNode<Test1, AddOp>::Type Test1Result;
            typedef FindFurthestLeftNode<Test2, MultiplyOp>::Type Test2Result;
            typedef FindFurthestLeftNode<Test3, AddOp>::Type Test3Result;
            typedef FindFurthestLeftNode<Test4, AddOp>::Type Test4Result;
            typedef FindFurthestLeftNode<Test5, AddOp>::Type Test5Result;

            typedef V ExpectedTest1Result;
            typedef V ExpectedTest2Result;
            typedef NegateNode ExpectedTest3Result;
            typedef V ExpectedTest4Result;
            typedef MultiplyNode ExpectedTest5Result;

            BOOST_MPL_ASSERT(( boost::is_same<ExpectedTest1Result, Test1Result> ));
            BOOST_MPL_ASSERT(( boost::is_same<ExpectedTest2Result, Test2Result> ));
            BOOST_MPL_ASSERT(( boost::is_same<ExpectedTest3Result, Test3Result> ));
            BOOST_MPL_ASSERT(( boost::is_same<ExpectedTest4Result, Test4Result> ));
            BOOST_MPL_ASSERT(( boost::is_same<ExpectedTest5Result, Test5Result> ));
        }

        BOOST_AUTO_TEST_CASE(TestFindingLargestNode)
        {
            typedef Node<NekVector<double> > V;
            typedef Node<V, AddOp, V> AddNode;
            typedef Node<V, MultiplyOp, V> MultiplyNode;
            typedef Node<V, NegateOp> NegateNode;

            typedef AddNode Test1;
            typedef MultiplyNode Test2;
            typedef Node<NegateNode, AddOp, AddNode> Test3;
            typedef Node<AddNode, AddOp, AddNode> Test4;
            typedef Node<MultiplyNode, AddOp, AddNode> Test5;
            typedef Node<AddNode, AddOp, MultiplyNode> Test6;

            typedef Node<V, AddOp, MultiplyNode> T1;
            typedef Node<T1, AddOp, V> Lhs;
            typedef Node<V, AddOp, MultiplyNode> T2;
            typedef Node<V, MultiplyOp, T2> Rhs;
            typedef Node<Lhs, AddOp, Rhs> Test7;

            typedef FindNodeWithMostTemporaries<Test1, AddOp>::Type Test1Result;
            static const unsigned int test1IndexResult = FindNodeWithMostTemporaries<Test1, AddOp>::Index;
            BOOST_STATIC_ASSERT(( FindNodeWithMostTemporaries<Test1, AddOp>::lhsTempCount == 0 ));
            BOOST_STATIC_ASSERT(( FindNodeWithMostTemporaries<Test1, AddOp>::LhsIndex == 0 ));
            BOOST_STATIC_ASSERT(( FindNodeWithMostTemporaries<Test1, AddOp>::RhsIndex == 1 ));
            BOOST_STATIC_ASSERT(( FindNodeWithMostTemporaries<Test1, AddOp>::lhsTempCount == 0 ));
            BOOST_STATIC_ASSERT(( FindNodeWithMostTemporaries<Test1, AddOp>::rhsTempCount == 0 ));

            typedef FindNodeWithMostTemporaries<Test2, MultiplyOp>::Type Test2Result;
            static const unsigned int test2IndexResult = FindNodeWithMostTemporaries<Test2, MultiplyOp>::Index;
            typedef FindNodeWithMostTemporaries<Test3, AddOp>::Type Test3Result;
            static const unsigned int test3IndexResult = FindNodeWithMostTemporaries<Test3, AddOp>::Index;
            typedef FindNodeWithMostTemporaries<Test4, AddOp>::Type Test4Result;
            static const unsigned int test4IndexResult = FindNodeWithMostTemporaries<Test4, AddOp>::Index;
            typedef FindNodeWithMostTemporaries<Test5, AddOp>::Type Test5Result;
            static const unsigned int test5IndexResult = FindNodeWithMostTemporaries<Test5, AddOp>::Index;
            typedef FindNodeWithMostTemporaries<Test6, AddOp>::Type Test6Result;
            static const unsigned int test6IndexResult = FindNodeWithMostTemporaries<Test6, AddOp>::Index;
            typedef FindNodeWithMostTemporaries<Test7, AddOp>::Type Test7Result;
            static const unsigned int test7IndexResult = FindNodeWithMostTemporaries<Test7, AddOp>::Index;

            BOOST_STATIC_ASSERT(( FindNodeWithMostTemporaries<Test6, AddOp>::rhsStartingIndex == 2 ));
            BOOST_STATIC_ASSERT(( FindNodeWithMostTemporaries<Test6, AddOp>::RhsIndex == 2 ));
            BOOST_STATIC_ASSERT(( FindNodeWithMostTemporaries<Test6, AddOp>::LhsIndex == 0 ));
            BOOST_STATIC_ASSERT(( FindNodeWithMostTemporaries<Test6, AddOp>::lhsTempCount == 0 ));
            BOOST_STATIC_ASSERT(( FindNodeWithMostTemporaries<Test6, AddOp>::rhsTempCount == 1 ));

            typedef V ExpectedTest1Result;
            static const unsigned int expectedTest1IndexResult = 0;
            typedef V ExpectedTest2Result;
            static const unsigned int expectedTest2IndexResult = 0;
            typedef NegateNode ExpectedTest3Result;
            static const unsigned int expectedTest3IndexResult = 0;
            typedef V ExpectedTest4Result;
            static const unsigned int expectedTest4IndexResult = 0;
            typedef MultiplyNode ExpectedTest5Result;
            static const unsigned int expectedTest5IndexResult = 0;
            typedef MultiplyNode ExpectedTest6Result;
            static const unsigned int expectedTest6IndexResult = 2;
            typedef Rhs ExpectedTest7Result;
            static const unsigned int expectedTest7IndexResult = 4;

            BOOST_MPL_ASSERT(( boost::is_same<ExpectedTest1Result, Test1Result> ));
            BOOST_STATIC_ASSERT(( test1IndexResult == expectedTest1IndexResult ));
            BOOST_MPL_ASSERT(( boost::is_same<ExpectedTest2Result, Test2Result> ));
            BOOST_STATIC_ASSERT(( test2IndexResult == expectedTest2IndexResult ));
            BOOST_MPL_ASSERT(( boost::is_same<ExpectedTest3Result, Test3Result> ));
            BOOST_STATIC_ASSERT(( test3IndexResult == expectedTest3IndexResult ));
            BOOST_MPL_ASSERT(( boost::is_same<ExpectedTest4Result, Test4Result> ));
            BOOST_STATIC_ASSERT(( test4IndexResult == expectedTest4IndexResult ));
            BOOST_MPL_ASSERT(( boost::is_same<ExpectedTest5Result, Test5Result> ));
            BOOST_STATIC_ASSERT(( test5IndexResult == expectedTest5IndexResult ));
            BOOST_MPL_ASSERT(( boost::is_same<ExpectedTest6Result, Test6Result> ));
            BOOST_STATIC_ASSERT(( test6IndexResult == expectedTest6IndexResult ));
            BOOST_MPL_ASSERT(( boost::is_same<ExpectedTest7Result, Test7Result> ));
            BOOST_STATIC_ASSERT(( test7IndexResult == expectedTest7IndexResult ));
            
        }

        BOOST_AUTO_TEST_CASE(TestIndexSwap)
        {
            typedef CreateVectorC<int, boost::mpl::int_, 10>::type Indices;
            
            typedef impl::SwapNodeIndices<Indices, 0, 2, 4, 7>::Type Test0;
            typedef impl::SwapNodeIndices<Indices, 0, 2, 4, 7>::Range1 Range1;
            

            typedef boost::mpl::vector_c<int, 4, 5, 6, 2, 3, 0, 1, 7, 8, 9> ExpectedTest0Result;

            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ExpectedTest0Result, 0>::type::value == boost::mpl::at_c<Test0, 0>::type::value));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ExpectedTest0Result, 1>::type::value == boost::mpl::at_c<Test0, 1>::type::value));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ExpectedTest0Result, 2>::type::value == boost::mpl::at_c<Test0, 2>::type::value));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ExpectedTest0Result, 3>::type::value == boost::mpl::at_c<Test0, 3>::type::value));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ExpectedTest0Result, 4>::type::value == boost::mpl::at_c<Test0, 4>::type::value));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ExpectedTest0Result, 5>::type::value == boost::mpl::at_c<Test0, 5>::type::value));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ExpectedTest0Result, 6>::type::value == boost::mpl::at_c<Test0, 6>::type::value));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ExpectedTest0Result, 7>::type::value == boost::mpl::at_c<Test0, 7>::type::value));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ExpectedTest0Result, 8>::type::value == boost::mpl::at_c<Test0, 8>::type::value));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ExpectedTest0Result, 9>::type::value == boost::mpl::at_c<Test0, 9>::type::value));

        }

        BOOST_AUTO_TEST_CASE(TestNodeSwap)
        {
            typedef Node<NekVector<double> > V;
            typedef Node<V, AddOp, V> AddNode;
            typedef Node<V, MultiplyOp, V> MultiplyNode;
            typedef Node<V, NegateOp> NegateNode;

            typedef impl::SwapNodes<AddNode, V, NegateNode, 0>::Type Test1Result;
            typedef Node<NegateNode, AddOp, V> ExpectedTest1Result;
        
            typedef impl::SwapNodes<AddNode, V, NegateNode, 1>::Type Test2Result;
            typedef Node<V, AddOp, NegateNode> ExpectedTest2Result;

            typedef Node<V, AddOp, AddNode> Lhs;
            typedef Node<AddNode, AddOp, AddNode> Rhs;
            typedef Node<Lhs, AddNode, Rhs> Test3;
            typedef impl::SwapNodes<Test3, AddNode, NegateNode, 1>::Type Test3Result;
            typedef Node<Node<V, AddOp, NegateNode>, AddNode, Rhs> ExpectedTest3Result;

            typedef impl::SwapNodes<Test3, Lhs, NegateNode, 0>::Type Test4Result;
            typedef Node<NegateNode, AddNode, Rhs> ExpectedTest4Result;

            BOOST_MPL_ASSERT(( boost::is_same<Test1Result, ExpectedTest1Result> ));
            BOOST_MPL_ASSERT(( boost::is_same<Test2Result, ExpectedTest2Result> ));
            BOOST_MPL_ASSERT(( boost::is_same<Test3Result, ExpectedTest3Result> ));
            BOOST_MPL_ASSERT(( boost::is_same<Test4Result, ExpectedTest4Result> ));
        }

        BOOST_AUTO_TEST_CASE(TestSwapACNodesIfNeeded)
        {
            typedef Node<NekVector<double> > V;
            typedef Node<V, MultiplyOp, V> M;
            typedef Node<V, AddOp, V> A;
            typedef Node<Node<A, AddOp, M>, AddOp, V> Test1;

            typedef SwapACNodesIfNeeded<Test1, Test1::Indices, 0>::TransformedNodeType Test1NodeType;
            typedef SwapACNodesIfNeeded<Test1, Test1::Indices, 0>::LhsNode Test1LhsNodeType;
            static const unsigned int TestLeftIndex = SwapACNodesIfNeeded<Test1, Test1::Indices, 0>::LhsIndex;
            typedef SwapACNodesIfNeeded<Test1, Test1::Indices, 0>::MostTemps Test1RhsNodeType;
            static const unsigned int TestRhsIndex = SwapACNodesIfNeeded<Test1, Test1::Indices, 0>::RhsIndex;

            typedef FindFurthestLeftNode<Test1, AddOp>::Type FurthestLeft;
            typedef SwapACNodesIfNeeded<Test1, Test1::Indices, 0>::TransformedIndicesType Test1IndicesType;

            typedef Node<Node<Node<M, AddOp, V>, AddOp, V>, AddOp, V> ExpectedTest1NodeType;

            BOOST_MPL_ASSERT(( boost::is_same<Test1NodeType, ExpectedTest1NodeType> ));
            BOOST_MPL_ASSERT(( boost::is_same<V, Test1LhsNodeType> ));
            BOOST_MPL_ASSERT(( boost::is_same<V, FurthestLeft> ));
            BOOST_STATIC_ASSERT(( TestLeftIndex == 0 ));
            BOOST_MPL_ASSERT(( boost::is_same<M, Test1RhsNodeType> ));
            BOOST_STATIC_ASSERT(( TestRhsIndex == 2 ));       
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test1IndicesType, 0>::type::value == 2));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test1IndicesType, 1>::type::value == 3));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test1IndicesType, 2>::type::value == 1));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test1IndicesType, 3>::type::value == 0));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test1IndicesType, 4>::type::value == 4));
            
        }

        BOOST_AUTO_TEST_CASE(TestSwapACNodesIfNeeded1)
        {
            // ( A + (BC+D) ) + E
            typedef Node<NekVector<double> > V;
            typedef Node<V, MultiplyOp, V> M;
            typedef Node<M, AddOp, V> BCDNode;
            typedef Node<V, AddOp, BCDNode> ABCDNode;
            typedef Node<ABCDNode, AddOp, V> Test1;

            typedef SwapACNodesIfNeeded<Test1, Test1::Indices, 0>::TransformedNodeType Test1NodeType;
            typedef SwapACNodesIfNeeded<Test1, Test1::Indices, 0>::TransformedIndicesType Test1IndicesType;

            typedef SwapACNodesIfNeeded<Test1, Test1::Indices, 0>::LhsNode Test1LhsNodeType;
            static const unsigned int TestLeftIndex = SwapACNodesIfNeeded<Test1, Test1::Indices, 0>::LhsIndex;
            typedef SwapACNodesIfNeeded<Test1, Test1::Indices, 0>::MostTemps Test1RhsNodeType;
            static const unsigned int TestRhsIndex = SwapACNodesIfNeeded<Test1, Test1::Indices, 0>::RhsIndex;

            typedef FindFurthestLeftNode<Test1, AddOp>::Type FurthestLeft;
            BOOST_MPL_ASSERT(( boost::is_same<FurthestLeft, V> ));
            BOOST_MPL_ASSERT(( boost::is_same<Test1LhsNodeType, V> ));
            BOOST_MPL_ASSERT(( boost::is_same<Test1RhsNodeType, M> ));
            BOOST_STATIC_ASSERT(( TestLeftIndex == 0 ));
            BOOST_STATIC_ASSERT(( TestRhsIndex == 1 ));


            // Expected is (BC + (A + D)) + E
            typedef Node<Node<M, AddOp, Node<V, AddOp, V> >, AddOp, V> ExpectedTest1NodeType;

            BOOST_MPL_ASSERT(( boost::is_same<Test1NodeType, ExpectedTest1NodeType> ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test1IndicesType, 0>::type::value == 1));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test1IndicesType, 1>::type::value == 2));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test1IndicesType, 2>::type::value == 0));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test1IndicesType, 3>::type::value == 3));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test1IndicesType, 4>::type::value == 4));
            
        }

        BOOST_AUTO_TEST_CASE(TestSwapACNodesIfNeeded2)
        {
            // ((((BC * (D+E) + A) + F) + G
            typedef Node<NekVector<double> > V;
            typedef Node<V, MultiplyOp, V> M;
            typedef Node<V, AddOp, V> A;
            typedef Node<M, MultiplyOp, A> BCDENode;
            typedef Node<BCDENode, AddOp, V> T0;
            typedef Node<T0, AddOp, V> T1;
            typedef Node<T1, AddOp, V> Test1;

            typedef SwapACNodesIfNeeded<Test1, Test1::Indices, 0>::TransformedNodeType Test1NodeType;
            typedef SwapACNodesIfNeeded<Test1, Test1::Indices, 0>::TransformedIndicesType Test1IndicesType;

            // Expected is ((((D+E)*C)*B)+A)+F
            typedef Node<A, MultiplyOp, V> A0;
            typedef Node<A0, MultiplyOp, V> A1;
            typedef Node<A1, AddOp, V> A2;
            typedef Node<A2, AddOp, V> A3;
            typedef Node<A3, AddOp, V> ExpectedTest1NodeType;

            //typedef Test1 ExpectedTest1NodeType;

            BOOST_MPL_ASSERT(( boost::is_same<Test1NodeType, ExpectedTest1NodeType> ));
            ////BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test1IndicesType, 0>::type::value == 1));
            ////BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test1IndicesType, 1>::type::value == 2));
            ////BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test1IndicesType, 2>::type::value == 0));
            ////BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test1IndicesType, 3>::type::value == 3));
            ////BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test1IndicesType, 4>::type::value == 4));
            
        }

        BOOST_AUTO_TEST_CASE(TestSwapWithUnary)
        {
            typedef Node<NekVector<double> > V;
            typedef Node<NekMatrix<double> > M;
            typedef Node<Node<M, MultiplyOp, V>, NegateOp, void> LhsNodeType;
            typedef Node<LhsNodeType, AddOp, V> NodeType;

            typedef SwapACNodesIfNeeded<NodeType, NodeType::Indices, 0>::TransformedNodeType Test1ResultNodeType;
            typedef SwapACNodesIfNeeded<NodeType, NodeType::Indices, 0>::TransformedIndicesType Test1ResultIndicesType;

            BOOST_MPL_ASSERT(( boost::is_same<NodeType, Test1ResultNodeType> ));
        }
    }
}

