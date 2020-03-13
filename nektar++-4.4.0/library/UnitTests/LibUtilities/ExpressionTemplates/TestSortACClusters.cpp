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

        BOOST_AUTO_TEST_CASE(TestSortAssociativeCommutativeClusterAPlusB)
        {
        //    // A+B -> A+B
        //    typedef NekMatrix<double> Matrix;
        //    typedef Node<Matrix> ConstantNode;
        //    typedef Node<ConstantNode, AddOp, ConstantNode> T0;
        //    typedef T0::Indices IndicesType;

        //    typedef SortAssociativeCommutativeCluster<T0, AddOp, IndicesType, 0>::TransformedNodeType ResultNodeType;
        //    typedef SortAssociativeCommutativeCluster<T0, AddOp, IndicesType, 0>::TransformedIndicesType IndicesType;

        //    BOOST_MPL_ASSERT(( boost::is_same<ResultNodeType, T0> ));

        //    BOOST_STATIC_ASSERT(( boost::mpl::size<IndicesType>::type::value == 2 ));
        //    BOOST_STATIC_ASSERT(( boost::mpl::at_c<IndicesType, 0>::type::value == 0 ));
        //    BOOST_STATIC_ASSERT(( boost::mpl::at_c<IndicesType, 1>::type::value == 1 ));
        }

        BOOST_AUTO_TEST_CASE(TestSortAssociativeCommutativeClusterAPlusBPlusC)
        {
            // A+B+C -> A+B+C
//            typedef NekMatrix<double> Matrix;
//            typedef Node<Matrix> ConstantNode;
//            typedef Node<ConstantNode, AddOp, ConstantNode> BinaryNode;
//            typedef Node<BinaryNode, AddOp, ConstantNode> T0;
//            typedef T0::Indices IndicesType;
//
//            typedef SortAssociativeCommutativeCluster<T0, AddOp, IndicesType, 0>::TransformedNodeType ResultNodeType;
//            typedef SortAssociativeCommutativeCluster<T0, AddOp, IndicesType, 0>::TransformedIndicesType IndicesType;
//
//            BOOST_MPL_ASSERT(( boost::is_same<ResultNodeType, T0> ));
//
//            BOOST_STATIC_ASSERT(( boost::mpl::size<IndicesType>::type::value == 3 ));
//            BOOST_STATIC_ASSERT(( boost::mpl::at_c<IndicesType, 0>::type::value == 0 ));
//            BOOST_STATIC_ASSERT(( boost::mpl::at_c<IndicesType, 1>::type::value == 1 ));
//            BOOST_STATIC_ASSERT(( boost::mpl::at_c<IndicesType, 2>::type::value == 2 ));
        }

        BOOST_AUTO_TEST_CASE(TestSortAssociativeCommutativeClusterAPlusBTimesC)
        {
            // A+BC -> BC+A
//            typedef NekMatrix<double> Matrix;
//            typedef Node<Matrix> ConstantNode;
//            typedef Node<ConstantNode, AddOp, ConstantNode> BinaryNode;
//            typedef Node<ConstantNode, AddOp, BinaryNode> T0;
//            typedef Node<BinaryNode, AddOp, ConstantNode> T1;
//            typedef T0::Indices IndicesType;
//
//            BOOST_STATIC_ASSERT(( SortAssociativeCommutativeCluster<T0, AddOp, IndicesType, 0>::Id == 1 )); 
//            typedef SortAssociativeCommutativeCluster<T0, AddOp, IndicesType, 0>::TransformedNodeType ResultNodeType;
//            typedef SortAssociativeCommutativeCluster<T0, AddOp, IndicesType, 0>::TransformedIndicesType ResultIndicesType;
//
//            BOOST_MPL_ASSERT(( boost::is_same<ResultNodeType, T1> ));
//
//            BOOST_STATIC_ASSERT(( boost::mpl::size<ResultIndicesType>::type::value == 3 ));
//            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultIndicesType, 0>::type::value == 1 ));
//            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultIndicesType, 1>::type::value == 2 ));
//            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultIndicesType, 2>::type::value == 0 ));
        }

        BOOST_AUTO_TEST_CASE(TestSortAssociativeCommutativeClusterAPlusBPlusCTimesD)
        {
            // A+B+CD -> CD+A+B
//            typedef NekMatrix<double> Matrix;
//            typedef Node<Matrix> ConstantNode;
//            typedef Node<ConstantNode, AddOp, ConstantNode> AddNode;
//            typedef Node<ConstantNode, MultiplyOp, ConstantNode> MultiplyNode;
//            typedef Node<AddNode, AddOp, MultiplyNode> T0;
//            typedef T0::Indices T0Indices;
//
//            typedef Node<MultiplyNode, AddOp, ConstantNode> Temp1;
//            typedef Node<Temp1, AddOp, ConstantNode> T1;
//
//            typedef SortAssociativeCommutativeCluster<T0, AddOp, T0Indices, 0>::TransformedNodeType ResultNodeType;
//            typedef SortAssociativeCommutativeCluster<T0, AddOp, T0Indices, 0>::TransformedIndicesType IndicesType;
//
//            BOOST_MPL_ASSERT(( boost::is_same<SortAssociativeCommutativeCluster<T0, AddOp, IndicesType, 0>::Tree0, T0> ));
//
//            // Check first inverse associative transform.
//            // A+B+CD -> A+(B+CD)
//            typedef Node<ConstantNode, AddOp, MultiplyNode> B0;
//            typedef Node<ConstantNode, AddOp, B0> B1;
//            BOOST_MPL_ASSERT(( boost::is_same<SortAssociativeCommutativeCluster<T0, AddOp, IndicesType, 0>::Tree1, B1> ));
//
//            // Check the commutative transform.
//            // A+(B+CD) -> A + (CD+B)
//            typedef Node<MultiplyNode, AddOp, ConstantNode> C0;
//            typedef Node<ConstantNode, AddOp, C0> C1;
//            
//            BOOST_MPL_ASSERT(( boost::is_same<SortAssociativeCommutativeCluster<T0, AddOp, IndicesType, 0>::Tree2, C1> ));
//
//            // Second associative to complete the transfomration.
//            // A + (CB+B) -> (A + CB) + B
//            typedef Node<ConstantNode, AddOp, MultiplyNode> D0;
//            typedef Node<D0, AddOp, ConstantNode> D1;
//            BOOST_MPL_ASSERT(( boost::is_same<SortAssociativeCommutativeCluster<T0, AddOp, IndicesType, 0>::Tree3, D1> ));
//
//            BOOST_MPL_ASSERT(( boost::is_same<ResultNodeType, T1> ));
//
//            BOOST_STATIC_ASSERT(( boost::mpl::size<IndicesType>::type::value == 4 ));
//            BOOST_STATIC_ASSERT(( boost::mpl::at_c<IndicesType, 0>::type::value == 2 ));
//            BOOST_STATIC_ASSERT(( boost::mpl::at_c<IndicesType, 1>::type::value == 3 ));
//            BOOST_STATIC_ASSERT(( boost::mpl::at_c<IndicesType, 2>::type::value == 0 ));
//            BOOST_STATIC_ASSERT(( boost::mpl::at_c<IndicesType, 3>::type::value == 1 ));
        }
    }
}