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

        //BOOST_AUTO_TEST_CASE(TestAssociativeTransform)
        //{
        //    typedef NekMatrix<double> Matrix;

        //    {
        //        // Base case tests.
        //        typedef Node<Matrix> T;
        //        BOOST_MPL_ASSERT(( boost::is_same<typename AssociativeTreeTransform<T>::TransformedNodeType, T >));

        //        typedef Node<Matrix, expt::NegateOp> U;
        //        BOOST_MPL_ASSERT(( boost::is_same<typename AssociativeTreeTransform<U>::TransformedNodeType, U >));
        //    }
        //    {
        //        // A+B -> A+B
        //        typedef Node<Node<Matrix>, expt::AddOp, Node<Matrix> > T;

        //        BOOST_MPL_ASSERT(( boost::is_same<typename AssociativeTreeTransform<T>::LeftChildType, Node<Matrix> >));
        //        BOOST_MPL_ASSERT(( boost::is_same<typename AssociativeTreeTransform<T>::RightChildType, Node<Matrix> >));
        //        BOOST_MPL_ASSERT(( boost::is_same<typename AssociativeTreeTransform<T>::TransformedRightChildType, Node<Matrix> >));

        //        BOOST_MPL_ASSERT(( boost::is_same<typename AssociativeTreeTransform<T>::UpdatedLeftChildType, Node<Matrix> >));
        //        BOOST_MPL_ASSERT(( boost::is_same<typename AssociativeTreeTransform<T>::FinalRightChildType, Node<Matrix> >));

        //        BOOST_MPL_ASSERT(( boost::is_same<typename AssociativeTreeTransform<T>::TransformedNodeType, T> ));
        //    }

        //    {
        //        // A + (B+C) -> (A+B) + C
        //        typedef Node<Node<Matrix>, expt::AddOp, Node<Matrix> > R;
        //        typedef Node<Matrix> L;
        //        typedef Node<L, expt::AddOp, R> Expression;

        //        typedef Node<R, expt::AddOp, L> ExpectedType;

        //        BOOST_MPL_ASSERT(( boost::is_same<typename AssociativeTreeTransform<Expression>::TransformedNodeType, ExpectedType> ));
        //    }

        //    {
        //        // (A+B) + (C+D) -> ((A+B)+C)+D
        //        typedef Node<Node<Matrix>, expt::AddOp, Node<Matrix> > Binary;
        //        typedef Node<Binary, expt::AddOp, Binary> Expression;

        //        typedef Node<Binary, expt::AddOp, Node<Matrix> > FirstTerm;
        //        typedef Node<FirstTerm, expt::AddOp, Node<Matrix> > ExpectedType;

        //        BOOST_MPL_ASSERT(( boost::is_same<typename AssociativeTreeTransform<Expression>::TransformedNodeType, ExpectedType> ));
        //    }

        //    {
        //        // (A*B) + (C+D) -> ((A*B)+C)+D
        //        typedef Node<Node<Matrix>, expt::MultiplyOp, Node<Matrix> > MultiplicationBinary;
        //        typedef Node<Node<Matrix>, expt::AddOp, Node<Matrix> > AddBinary;
        //        typedef Node<MultiplicationBinary, expt::AddOp, AddBinary> Expression;

        //        typedef Node<MultiplicationBinary, expt::AddOp, Node<Matrix> > FirstTerm;
        //        typedef Node<FirstTerm, expt::AddOp, Node<Matrix> > ExpectedType;

        //        BOOST_MPL_ASSERT(( boost::is_same<typename AssociativeTreeTransform<Expression>::TransformedNodeType, ExpectedType> ));
        //    }

        //    {
        //        // (A+B) + (C*D) -> (A+B) + (C*D)
        //        typedef Node<Node<Matrix>, expt::MultiplyOp, Node<Matrix> > MultiplicationBinary;
        //        typedef Node<Node<Matrix>, expt::AddOp, Node<Matrix> > AddBinary;
        //        typedef Node<AddBinary, expt::AddOp, MultiplicationBinary> Expression;

        //        typedef Expression ExpectedType;

        //        BOOST_MPL_ASSERT(( boost::is_same<typename AssociativeTreeTransform<Expression>::TransformedNodeType, ExpectedType> ));
        //    }
        //}

        BOOST_AUTO_TEST_CASE(TestAssociativeTransformConstantAndUnaryNodes)
        {
            typedef NekMatrix<double> Matrix;

            typedef expt::Node<Matrix> T;
            BOOST_MPL_ASSERT(( boost::is_same<expt::AssociativeTransform<T>::TransformedNodeType, T >));
            BOOST_MPL_ASSERT(( boost::is_same<expt::InverseAssociativeTransform<T>::TransformedNodeType, T >));

            typedef expt::Node<Matrix, expt::NegateOp> U;
            BOOST_MPL_ASSERT(( boost::is_same<expt::AssociativeTransform<U>::TransformedNodeType, U >));
            BOOST_MPL_ASSERT(( boost::is_same<expt::InverseAssociativeTransform<T>::TransformedNodeType, T >));
        }


        BOOST_AUTO_TEST_CASE(TestAssociativeTransformAPlusB)
        {
            typedef NekMatrix<double> Matrix;

            // A+B -> A+B
            typedef expt::Node<expt::Node<Matrix>, expt::AddOp, expt::Node<Matrix> > T;

            BOOST_MPL_ASSERT(( boost::is_same<expt::AssociativeTransform<T>::TransformedNodeType, T >));
            BOOST_MPL_ASSERT(( boost::is_same<expt::InverseAssociativeTransform<T>::TransformedNodeType, T >));
        }

        BOOST_AUTO_TEST_CASE(TestAssociativeTransformAPlusBPlusC)
        {
            typedef NekMatrix<double> Matrix;

            // A + (B+C) -> (A+B) + C
            typedef expt::Node<Matrix> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> BinaryNode;
            typedef expt::Node<ConstantNode, expt::AddOp, BinaryNode> T0;
            typedef expt::Node<BinaryNode, expt::AddOp, ConstantNode> T1;

            BOOST_MPL_ASSERT(( boost::is_same<expt::AssociativeTransform<T0>::TransformedNodeType, T1 >));
            BOOST_MPL_ASSERT(( boost::is_same<expt::AssociativeTransform<T1>::TransformedNodeType, T1 >));

            BOOST_MPL_ASSERT(( boost::is_same<expt::InverseAssociativeTransform<T1>::TransformedNodeType, T0 >));
            BOOST_MPL_ASSERT(( boost::is_same<expt::InverseAssociativeTransform<T0>::TransformedNodeType, T0 >));
        }

        BOOST_AUTO_TEST_CASE(TestAssociativeTransformAPlusBPlusCPlusD)
        {
            typedef NekMatrix<double> Matrix;

            // (A+B) + (C+D) -> ((A+B)+C)+D
            typedef expt::Node<Matrix> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> BinaryNode;
            typedef expt::Node<BinaryNode, expt::AddOp, BinaryNode> T0;
            
            typedef expt::Node<BinaryNode, expt::AddOp, ConstantNode> ExpectedLhs;
            typedef expt::Node<ExpectedLhs, expt::AddOp, ConstantNode> T1;

            
            BOOST_MPL_ASSERT(( boost::is_same<expt::AssociativeTransform<T0>::TransformedNodeType, T1 >));
            BOOST_MPL_ASSERT(( boost::is_same<expt::AssociativeTransform<T1>::TransformedNodeType, T1 >));

            // ((A+B)+C)+D -> (A+B) + (C+D)
            BOOST_MPL_ASSERT(( boost::is_same<expt::InverseAssociativeTransform<T1>::TransformedNodeType, T0 >));
            
            // (A+B) + (C+D) -> (A + (B + (C+D))
            typedef expt::Node<ConstantNode, expt::AddOp, BinaryNode> ExpectedRhs;
            typedef expt::Node<ConstantNode, expt::AddOp, ExpectedRhs> T2;
            BOOST_MPL_ASSERT(( boost::is_same<expt::InverseAssociativeTransform<T0>::TransformedNodeType, T2 >));
        }


            //{
            //    // (A*B) + (C+D) -> ((A*B)+C)+D
            //    typedef expt::Node<expt::Node<Matrix>, expt::MultiplyOp, expt::Node<Matrix> > MultiplicationBinary;
            //    typedef expt::Node<expt::Node<Matrix>, expt::AddOp, expt::Node<Matrix> > AddBinary;
            //    typedef expt::Node<MultiplicationBinary, expt::AddOp, AddBinary> Expression;

            //    typedef expt::Node<MultiplicationBinary, expt::AddOp, expt::Node<Matrix> > FirstTerm;
            //    typedef expt::Node<FirstTerm, expt::AddOp, expt::Node<Matrix> > ExpectedType;

            //    BOOST_MPL_ASSERT(( boost::is_same<typename AssociativeTreeTransform<Expression>::TransformedNodeType, ExpectedType> ));
            //}

            //{
            //    // (A+B) + (C*D) -> (A+B) + (C*D)
            //    typedef expt::Node<expt::Node<Matrix>, expt::MultiplyOp, expt::Node<Matrix> > MultiplicationBinary;
            //    typedef expt::Node<expt::Node<Matrix>, expt::AddOp, expt::Node<Matrix> > AddBinary;
            //    typedef expt::Node<AddBinary, expt::AddOp, MultiplicationBinary> Expression;

            //    typedef Expression ExpectedType;

            //    BOOST_MPL_ASSERT(( boost::is_same<typename AssociativeTreeTransform<Expression>::TransformedNodeType, ExpectedType> ));
            //}
        //}

        BOOST_AUTO_TEST_CASE(TestSortAssociativeCommutativeClusterAPlusB)
        {
        //    // A+B -> A+B
        //    typedef NekMatrix<double> Matrix;
        //    typedef expt::Node<Matrix> ConstantNode;
        //    typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> T0;
        //    typedef T0::Indices IndicesType;

        //    typedef SortAssociativeCommutativeCluster<T0, expt::AddOp, IndicesType, 0>::TransformedNodeType ResultNodeType;
        //    typedef SortAssociativeCommutativeCluster<T0, expt::AddOp, IndicesType, 0>::TransformedIndicesType IndicesType;

        //    BOOST_MPL_ASSERT(( boost::is_same<ResultNodeType, T0> ));

        //    BOOST_STATIC_ASSERT(( boost::mpl::size<IndicesType>::type::value == 2 ));
        //    BOOST_STATIC_ASSERT(( boost::mpl::at_c<IndicesType, 0>::type::value == 0 ));
        //    BOOST_STATIC_ASSERT(( boost::mpl::at_c<IndicesType, 1>::type::value == 1 ));
        }

        BOOST_AUTO_TEST_CASE(TestSortAssociativeCommutativeClusterAPlusBPlusC)
        {
            // A+B+C -> A+B+C
//            typedef NekMatrix<double> Matrix;
//            typedef expt::Node<Matrix> ConstantNode;
//            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> BinaryNode;
//            typedef expt::Node<BinaryNode, expt::AddOp, ConstantNode> T0;
//            typedef T0::Indices IndicesType;
//
//            typedef SortAssociativeCommutativeCluster<T0, expt::AddOp, IndicesType, 0>::TransformedNodeType ResultNodeType;
//            typedef SortAssociativeCommutativeCluster<T0, expt::AddOp, IndicesType, 0>::TransformedIndicesType IndicesType;
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
//            typedef expt::Node<Matrix> ConstantNode;
//            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> BinaryNode;
//            typedef expt::Node<ConstantNode, expt::AddOp, BinaryNode> T0;
//            typedef expt::Node<BinaryNode, expt::AddOp, ConstantNode> T1;
//            typedef T0::Indices IndicesType;
//
//            BOOST_STATIC_ASSERT(( SortAssociativeCommutativeCluster<T0, expt::AddOp, IndicesType, 0>::Id == 1 )); 
//            typedef SortAssociativeCommutativeCluster<T0, expt::AddOp, IndicesType, 0>::TransformedNodeType ResultNodeType;
//            typedef SortAssociativeCommutativeCluster<T0, expt::AddOp, IndicesType, 0>::TransformedIndicesType ResultIndicesType;
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
//            typedef expt::Node<Matrix> ConstantNode;
//            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> AddNode;
//            typedef expt::Node<ConstantNode, expt::MultiplyOp, ConstantNode> MultiplyNode;
//            typedef expt::Node<AddNode, expt::AddOp, MultiplyNode> T0;
//            typedef T0::Indices T0Indices;
//
//            typedef expt::Node<MultiplyNode, expt::AddOp, ConstantNode> Temp1;
//            typedef expt::Node<Temp1, expt::AddOp, ConstantNode> T1;
//
//            typedef SortAssociativeCommutativeCluster<T0, expt::AddOp, T0Indices, 0>::TransformedNodeType ResultNodeType;
//            typedef SortAssociativeCommutativeCluster<T0, expt::AddOp, T0Indices, 0>::TransformedIndicesType IndicesType;
//
//            BOOST_MPL_ASSERT(( boost::is_same<SortAssociativeCommutativeCluster<T0, expt::AddOp, IndicesType, 0>::Tree0, T0> ));
//
//            // Check first inverse associative transform.
//            // A+B+CD -> A+(B+CD)
//            typedef expt::Node<ConstantNode, expt::AddOp, MultiplyNode> B0;
//            typedef expt::Node<ConstantNode, expt::AddOp, B0> B1;
//            BOOST_MPL_ASSERT(( boost::is_same<SortAssociativeCommutativeCluster<T0, expt::AddOp, IndicesType, 0>::Tree1, B1> ));
//
//            // Check the commutative transform.
//            // A+(B+CD) -> A + (CD+B)
//            typedef expt::Node<MultiplyNode, expt::AddOp, ConstantNode> C0;
//            typedef expt::Node<ConstantNode, expt::AddOp, C0> C1;
//            
//            BOOST_MPL_ASSERT(( boost::is_same<SortAssociativeCommutativeCluster<T0, expt::AddOp, IndicesType, 0>::Tree2, C1> ));
//
//            // Second associative to complete the transfomration.
//            // A + (CB+B) -> (A + CB) + B
//            typedef expt::Node<ConstantNode, expt::AddOp, MultiplyNode> D0;
//            typedef expt::Node<D0, expt::AddOp, ConstantNode> D1;
//            BOOST_MPL_ASSERT(( boost::is_same<SortAssociativeCommutativeCluster<T0, expt::AddOp, IndicesType, 0>::Tree3, D1> ));
//
//            BOOST_MPL_ASSERT(( boost::is_same<ResultNodeType, T1> ));
//
//            BOOST_STATIC_ASSERT(( boost::mpl::size<IndicesType>::type::value == 4 ));
//            BOOST_STATIC_ASSERT(( boost::mpl::at_c<IndicesType, 0>::type::value == 2 ));
//            BOOST_STATIC_ASSERT(( boost::mpl::at_c<IndicesType, 1>::type::value == 3 ));
//            BOOST_STATIC_ASSERT(( boost::mpl::at_c<IndicesType, 2>::type::value == 0 ));
//            BOOST_STATIC_ASSERT(( boost::mpl::at_c<IndicesType, 3>::type::value == 1 ));
        }

        BOOST_AUTO_TEST_CASE(TestNegateOp)
        {
            NekVector<double> v(3);
            
            //NekVector<double> result = -v;
        }
    }
}
