/////////////////////////////////////////////////////////////////////////////////
////
//// For more information, please see: http://www.nektar.info
////
//// The MIT License
////
//// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//// Department of Aeronautics, Imperial College London (UK), and Scientific
//// Computing and Imaging Institute, University of Utah (USA).
////
//// License for the specific language governing rights and limitations under
//// Permission is hereby granted, free of charge, to any person obtaining a
//// copy of this software and associated documentation files (the "Software"),
//// to deal in the Software without restriction, including without limitation
//// the rights to use, copy, modify, merge, publish, distribute, sublicense,
//// and/or sell copies of the Software, and to permit persons to whom the
//// Software is furnished to do so, subject to the following conditions:
////
//// The above copyright notice and this permission notice shall be included
//// in all copies or substantial portions of the Software.
////
//// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//// DEALINGS IN THE SOFTWARE.
////
/////////////////////////////////////////////////////////////////////////////////
//
//#ifndef NEKTAR_USE_EXPRESSION_TEMPLATES
//#define NEKTAR_USE_EXPRESSION_TEMPLATES
//#endif //NEKTAR_USE_EXPRESSION_TEMPLATES
//
//#include <boost/test/auto_unit_test.hpp>
//#include <boost/test/test_case_template.hpp>
//#include <boost/test/floating_point_comparison.hpp>
//#include <boost/test/unit_test.hpp>
//
//#include <UnitTests/LibUtilities/ExpressionTemplates/CountedObjectExpression.h>
//#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
//#include <ExpressionTemplates/ExpressionTemplates.hpp>
//#include <ExpressionTemplates/Node.hpp>
//#include <ExpressionTemplates/RemoveAllUnecessaryTemporaries.hpp>
//#include <boost/mpl/assert.hpp>
//#include <boost/type_traits.hpp>
//#include <LibUtilities/LinearAlgebra/MatrixSize.hpp>
//
//namespace Nektar
//{
//    namespace UnitTests
//    {
//        BOOST_AUTO_TEST_CASE(TestUpdateLeftChildIfNeeded)
//        {
//            typedef NekMatrix<double> Matrix;
//
//            typedef expt::Node<double> T1;
//
//            // Constant node - no change
//            BOOST_MPL_ASSERT(( boost::is_same<typename UpdateLeftChildIfNeeded<T1, void, void>::UpdatedLeftChild, T1> ));
//            BOOST_MPL_ASSERT(( boost::is_same<typename UpdateLeftChildIfNeeded<T1, void, void>::UpdatedRightChild, void> ));
//
//            // Unary Node - no change.
//            BOOST_MPL_ASSERT(( boost::is_same<typename UpdateLeftChildIfNeeded<T1, expt::NegateOp, void>::UpdatedLeftChild, T1> ));
//            BOOST_MPL_ASSERT(( boost::is_same<typename UpdateLeftChildIfNeeded<T1, expt::NegateOp, void>::UpdatedRightChild, void> ));
//
//            // Binary node with constant children - no change
//            typedef expt::Node<double> L;
//            typedef expt::Node<int> R;
//            BOOST_MPL_ASSERT(( boost::is_same<typename UpdateLeftChildIfNeeded<L, expt::AddOp, R>::UpdatedLeftChild, L> ));
//            BOOST_MPL_ASSERT(( boost::is_same<typename UpdateLeftChildIfNeeded<L, expt::AddOp, R>::UpdatedRightChild, R> ));
//
//            {
//                // Deeper tree with no change.
//                typedef expt::Node<expt::Node<Matrix>, expt::AddOp, expt::Node<Matrix> > L;
//                typedef expt::Node<Matrix> R;
//
//                BOOST_MPL_ASSERT(( boost::is_same<typename UpdateLeftChildIfNeeded<L, expt::AddOp, R>::UpdatedLeftChild, L> ));
//                BOOST_MPL_ASSERT(( boost::is_same<typename UpdateLeftChildIfNeeded<L, expt::AddOp, R>::UpdatedRightChild, R> ));
//            }
//
//            {
//                // Child moves from the right.
//                typedef expt::Node<Matrix> L;
//                typedef expt::Node<expt::Node<Matrix>, expt::AddOp, expt::Node<Matrix> > R;
//                
//                typedef expt::Node<expt::Node<Matrix>, expt::AddOp, expt::Node<Matrix> > ExpectedL;
//                typedef expt::Node<Matrix> ExpectedR;
//
//                BOOST_MPL_ASSERT(( boost::is_same<typename UpdateLeftChildIfNeeded<L, expt::AddOp, R>::UpdatedLeftChild, ExpectedL> ));
//                BOOST_MPL_ASSERT(( boost::is_same<typename UpdateLeftChildIfNeeded<L, expt::AddOp, R>::UpdatedRightChild, ExpectedR> ));
//            }
//
//            {
//                // Child doesn't move because operators are not the same.
//                typedef expt::Node<Matrix> L;
//                typedef expt::Node<expt::Node<Matrix>, expt::AddOp, expt::Node<Matrix> > R;
//                
//                typedef expt::Node<expt::Node<Matrix>, expt::AddOp, expt::Node<Matrix> > ExpectedR;
//                typedef expt::Node<Matrix> ExpectedL;
//
//                BOOST_MPL_ASSERT(( boost::is_same<typename UpdateLeftChildIfNeeded<L, expt::SubtractOp, R>::UpdatedLeftChild, ExpectedL> ));
//                BOOST_MPL_ASSERT(( boost::is_same<typename UpdateLeftChildIfNeeded<L, expt::SubtractOp, R>::UpdatedRightChild, ExpectedR> ));
//            }
//
//            {
//                // Child doesn't move because operators are not associative
//                typedef expt::Node<Matrix> L;
//                typedef expt::Node<expt::Node<Matrix>, expt::SubtractOp, expt::Node<Matrix> > R;
//                
//                typedef expt::Node<expt::Node<Matrix>, expt::SubtractOp, expt::Node<Matrix> > ExpectedR;
//                typedef expt::Node<Matrix> ExpectedL;
//
//                BOOST_MPL_ASSERT(( boost::is_same<typename UpdateLeftChildIfNeeded<L, expt::SubtractOp, R>::UpdatedLeftChild, ExpectedL> ));
//                BOOST_MPL_ASSERT(( boost::is_same<typename UpdateLeftChildIfNeeded<L, expt::SubtractOp, R>::UpdatedRightChild, ExpectedR> ));
//            }
//        }
//
//        BOOST_AUTO_TEST_CASE(TestAssociativeTransform)
//        {
//            typedef NekMatrix<double> Matrix;
//
//            {
//                // Base case tests.
//                typedef Node<Matrix> T;
//                BOOST_MPL_ASSERT(( boost::is_same<typename AssociativeTreeTransform<T>::TransformedNodeType, T >));
//
//                typedef Node<Matrix, expt::NegateOp> U;
//                BOOST_MPL_ASSERT(( boost::is_same<typename AssociativeTreeTransform<U>::TransformedNodeType, U >));
//            }
//            {
//                // A+B -> A+B
//                typedef Node<Node<Matrix>, expt::AddOp, Node<Matrix> > T;
//
//                BOOST_MPL_ASSERT(( boost::is_same<typename AssociativeTreeTransform<T>::LeftChildType, Node<Matrix> >));
//                BOOST_MPL_ASSERT(( boost::is_same<typename AssociativeTreeTransform<T>::RightChildType, Node<Matrix> >));
//                BOOST_MPL_ASSERT(( boost::is_same<typename AssociativeTreeTransform<T>::TransformedRightChildType, Node<Matrix> >));
//
//                BOOST_MPL_ASSERT(( boost::is_same<typename AssociativeTreeTransform<T>::UpdatedLeftChildType, Node<Matrix> >));
//                BOOST_MPL_ASSERT(( boost::is_same<typename AssociativeTreeTransform<T>::FinalRightChildType, Node<Matrix> >));
//
//                BOOST_MPL_ASSERT(( boost::is_same<typename AssociativeTreeTransform<T>::TransformedNodeType, T> ));
//            }
//
//            {
//                // A + (B+C) -> (A+B) + C
//                typedef Node<Node<Matrix>, expt::AddOp, Node<Matrix> > R;
//                typedef Node<Matrix> L;
//                typedef Node<L, expt::AddOp, R> Expression;
//
//                typedef Node<R, expt::AddOp, L> ExpectedType;
//
//                BOOST_MPL_ASSERT(( boost::is_same<typename AssociativeTreeTransform<Expression>::TransformedNodeType, ExpectedType> ));
//            }
//
//            {
//                // (A+B) + (C+D) -> ((A+B)+C)+D
//                typedef Node<Node<Matrix>, expt::AddOp, Node<Matrix> > Binary;
//                typedef Node<Binary, expt::AddOp, Binary> Expression;
//
//                typedef Node<Binary, expt::AddOp, Node<Matrix> > FirstTerm;
//                typedef Node<FirstTerm, expt::AddOp, Node<Matrix> > ExpectedType;
//
//                BOOST_MPL_ASSERT(( boost::is_same<typename AssociativeTreeTransform<Expression>::TransformedNodeType, ExpectedType> ));
//            }
//
//            {
//                // (A*B) + (C+D) -> ((A*B)+C)+D
//                typedef Node<Node<Matrix>, expt::MultiplyOp, Node<Matrix> > MultiplicationBinary;
//                typedef Node<Node<Matrix>, expt::AddOp, Node<Matrix> > AddBinary;
//                typedef Node<MultiplicationBinary, expt::AddOp, AddBinary> Expression;
//
//                typedef Node<MultiplicationBinary, expt::AddOp, Node<Matrix> > FirstTerm;
//                typedef Node<FirstTerm, expt::AddOp, Node<Matrix> > ExpectedType;
//
//                BOOST_MPL_ASSERT(( boost::is_same<typename AssociativeTreeTransform<Expression>::TransformedNodeType, ExpectedType> ));
//            }
//
//            {
//                // (A+B) + (C*D) -> (A+B) + (C*D)
//                typedef Node<Node<Matrix>, expt::MultiplyOp, Node<Matrix> > MultiplicationBinary;
//                typedef Node<Node<Matrix>, expt::AddOp, Node<Matrix> > AddBinary;
//                typedef Node<AddBinary, expt::AddOp, MultiplicationBinary> Expression;
//
//                typedef Expression ExpectedType;
//
//                BOOST_MPL_ASSERT(( boost::is_same<typename AssociativeTreeTransform<Expression>::TransformedNodeType, ExpectedType> ));
//            }
//        }
//    }
//}
