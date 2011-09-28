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
#include <math.h>

namespace Nektar
{
    namespace UnitTests
    {
        BOOST_AUTO_TEST_CASE(BaseCases)
        {
            typedef NekMatrix<double> Matrix;

            // A -> A
            typedef expt::Node<Matrix> ConstantNode;
            typedef expt::RemoveUnecessaryTemporaries<ConstantNode, ConstantNode::Indices, 0>::TransformedNodeType ExpectedConstantNode;

            BOOST_MPL_ASSERT(( boost::is_same<ConstantNode, ExpectedConstantNode> ));

            // -A -> -A
            typedef expt::Node<ConstantNode, expt::NegateOp> UnaryNode;
            typedef expt::RemoveUnecessaryTemporaries<UnaryNode, UnaryNode::Indices, 0>::TransformedNodeType ExpectedUnaryNode;
            BOOST_MPL_ASSERT(( boost::is_same<UnaryNode, ExpectedUnaryNode> ));

            // A+B -> A+B
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> BinaryAdd;
            typedef expt::RemoveUnecessaryTemporaries<BinaryAdd, BinaryAdd::Indices, 0>::TransformedNodeType ExpectedBinaryAddNode;
            typedef expt::RemoveUnecessaryTemporaries<BinaryAdd, BinaryAdd::Indices, 0>::TransformedIndicesType ExpectedBinaryAddIndices;

            BOOST_MPL_ASSERT(( boost::is_same<BinaryAdd, ExpectedBinaryAddNode> ));
            BOOST_STATIC_ASSERT(( boost::mpl::size<ExpectedBinaryAddIndices>::type::value == 2 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ExpectedBinaryAddIndices, 0>::type::value == 0 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ExpectedBinaryAddIndices, 1>::type::value == 1 ));

            // A*B -> A*B
            typedef expt::Node<ConstantNode, expt::MultiplyOp, ConstantNode> BinaryMultiply;
            typedef expt::RemoveUnecessaryTemporaries<BinaryMultiply, BinaryMultiply::Indices, 0>::TransformedNodeType ExpectedBinaryMultiplyNode;
            typedef expt::RemoveUnecessaryTemporaries<BinaryMultiply, BinaryMultiply::Indices, 0>::TransformedIndicesType ExpectedBinaryMultiplyIndices;

            BOOST_MPL_ASSERT(( boost::is_same<BinaryMultiply, ExpectedBinaryMultiplyNode> ));
            BOOST_STATIC_ASSERT(( boost::mpl::size<ExpectedBinaryMultiplyIndices>::type::value == 2 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ExpectedBinaryMultiplyIndices, 0>::type::value == 0 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ExpectedBinaryMultiplyIndices, 1>::type::value == 1 ));
        }

        ////////////////////////////////////////////////////
        // Tests for Expression 1: A+(B+C)
        ////////////////////////////////////////////////////
        BOOST_AUTO_TEST_CASE(Exp1NonAssociativeNonCommutative)
        {
            // A + (B+C) -> A+ (B+C)
            typedef ExpressionTemplateTestObject00 Obj;
            typedef expt::Node<Obj> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> BinaryNode;

            typedef expt::Node<ConstantNode, expt::AddOp, BinaryNode> Exp;
            typedef expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::TransformedNodeType ResultExpNode;
            typedef expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::TransformedIndicesType ResultExpIndices;

            BOOST_MPL_ASSERT(( boost::is_same<Exp, ResultExpNode> ));
            BOOST_STATIC_ASSERT(( boost::mpl::size<ResultExpIndices>::type::value == 3 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 0>::type::value == 0 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 1>::type::value == 1 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 2>::type::value == 2 ));
        }

        BOOST_AUTO_TEST_CASE(Exp1NonAssociativeCommutative)
        {
            // A + (B+C) -> (B+C) + A
            typedef ExpressionTemplateTestObject01 Obj;
            typedef expt::Node<Obj> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> BinaryNode;

            typedef expt::Node<ConstantNode, expt::AddOp, BinaryNode> Exp;
            typedef expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::TransformedNodeType ResultExpNode;
            typedef expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::TransformedIndicesType ResultExpIndices;
            typedef expt::Node<BinaryNode, expt::AddOp, ConstantNode> ExpectedResultNode;

            // Test the associative cluster sorting associated with this expression.
            BOOST_STATIC_ASSERT(( expt::SortAssociativeCommutativeClusters<ExpectedResultNode, Exp::Indices, 0>::SpecializationId == 0 ));
            BOOST_MPL_ASSERT(( boost::is_same<expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::Tree0Type, ExpectedResultNode> ));
            BOOST_MPL_ASSERT(( boost::is_same<expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::TransformedNodeType, ExpectedResultNode> ));
            BOOST_MPL_ASSERT(( boost::is_same<expt::SortAssociativeCommutativeClusters<ExpectedResultNode, Exp::Indices, 0>::TransformedNodeType, ExpectedResultNode> ));

            BOOST_MPL_ASSERT(( boost::is_same<ExpectedResultNode, ResultExpNode> ));
            BOOST_STATIC_ASSERT(( boost::mpl::size<ResultExpIndices>::type::value == 3 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 0>::type::value == 1 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 1>::type::value == 2 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 2>::type::value == 0 ));
        }


        BOOST_AUTO_TEST_CASE(Exp1AssociativeNonCommutative)
        {
            // A + (B+C) -> (A+B) + C
            typedef ExpressionTemplateTestObject10 Obj;
            typedef expt::Node<Obj> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> BinaryNode;

            typedef expt::Node<ConstantNode, expt::AddOp, BinaryNode> Exp;
            typedef expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::TransformedNodeType ResultExpNode;
            typedef expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::TransformedIndicesType ResultExpIndices;
            typedef expt::Node<BinaryNode, expt::AddOp, ConstantNode> ExpectedResultNode;

            BOOST_MPL_ASSERT(( boost::is_same<ExpectedResultNode, ResultExpNode> ));
            BOOST_STATIC_ASSERT(( boost::mpl::size<ResultExpIndices>::type::value == 3 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 0>::type::value == 0 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 1>::type::value == 1 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 2>::type::value == 2 ));
        }

        BOOST_AUTO_TEST_CASE(Exp1AssociativeCommutative)
        {
            // A + (B+C) -> (A+B) + C
            typedef ExpressionTemplateTestObject11 Obj;
            typedef expt::Node<Obj> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> BinaryNode;

            typedef expt::Node<ConstantNode, expt::AddOp, BinaryNode> Exp;
            typedef expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::TransformedNodeType ResultExpNode;
            typedef expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::TransformedIndicesType ResultExpIndices;
            typedef expt::Node<BinaryNode, expt::AddOp, ConstantNode> ExpectedResultNode;

            BOOST_MPL_ASSERT(( boost::is_same<ExpectedResultNode, ResultExpNode> ));
            BOOST_STATIC_ASSERT(( boost::mpl::size<ResultExpIndices>::type::value == 3 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 0>::type::value == 0 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 1>::type::value == 1 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 2>::type::value == 2 ));
        }


        ////////////////////////////////////////////////////
        // Expression 2: (A+B) + (C+D)
        ///////////////////////////////////////////////////
        BOOST_AUTO_TEST_CASE(Exp2NonAssociativeNonCommutative)
        {
            // (A+B) + (C+D) -> (A+B) + (C+D)
            typedef ExpressionTemplateTestObject00 Obj;
            typedef expt::Node<Obj> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> BinaryNode;

            typedef expt::Node<BinaryNode, expt::AddOp, BinaryNode> Exp;
            typedef expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::TransformedNodeType ResultExpNode;
            typedef expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::TransformedIndicesType ResultExpIndices;
            typedef Exp ExpectedResultNode;

            BOOST_MPL_ASSERT(( boost::is_same<ExpectedResultNode, ResultExpNode> ));
            BOOST_STATIC_ASSERT(( boost::mpl::size<ResultExpIndices>::type::value == 4 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 0>::type::value == 0 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 1>::type::value == 1 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 2>::type::value == 2 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 3>::type::value == 3 ));
        }

        BOOST_AUTO_TEST_CASE(Exp2NonAssociativeCommutative)
        {
            // (A+B) + (C+D) -> (A+B) + (C+D)
            typedef ExpressionTemplateTestObject01 Obj;
            typedef expt::Node<Obj> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> BinaryNode;

            typedef expt::Node<BinaryNode, expt::AddOp, BinaryNode> Exp;
            typedef expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::TransformedNodeType ResultExpNode;
            typedef expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::TransformedIndicesType ResultExpIndices;
            typedef Exp ExpectedResultNode;

            BOOST_MPL_ASSERT(( boost::is_same<ExpectedResultNode, ResultExpNode> ));
            BOOST_STATIC_ASSERT(( boost::mpl::size<ResultExpIndices>::type::value == 4 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 0>::type::value == 0 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 1>::type::value == 1 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 2>::type::value == 2 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 3>::type::value == 3 ));
        }

        BOOST_AUTO_TEST_CASE(Exp2AssociativeNonCommutative)
        {
            // (A+B) + (C+D) -> ((A+B)+C)+D
            typedef ExpressionTemplateTestObject10 Obj;
            typedef expt::Node<Obj> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> BinaryNode;

            typedef expt::Node<BinaryNode, expt::AddOp, BinaryNode> Exp;
            typedef expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::TransformedNodeType ResultExpNode;
            typedef expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::TransformedIndicesType ResultExpIndices;
            typedef expt::Node<BinaryNode, expt::AddOp, ConstantNode> T0;
            typedef expt::Node<T0, expt::AddOp, ConstantNode> ExpectedResultNode;

            BOOST_MPL_ASSERT(( boost::is_same<ExpectedResultNode, ResultExpNode> ));
            BOOST_STATIC_ASSERT(( boost::mpl::size<ResultExpIndices>::type::value == 4 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 0>::type::value == 0 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 1>::type::value == 1 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 2>::type::value == 2 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 3>::type::value == 3 ));
        }

        BOOST_AUTO_TEST_CASE(Exp2AssociativeCommutative)
        {
            // (A+B) + (C+D) -> ((A+B)+C)+D
            typedef ExpressionTemplateTestObject11 Obj;
            typedef expt::Node<Obj> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> BinaryNode;

            typedef expt::Node<BinaryNode, expt::AddOp, BinaryNode> Exp;
            typedef expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::TransformedNodeType ResultExpNode;
            typedef expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::TransformedIndicesType ResultExpIndices;
            typedef expt::Node<BinaryNode, expt::AddOp, ConstantNode> T0;
            typedef expt::Node<T0, expt::AddOp, ConstantNode> ExpectedResultNode;

            BOOST_MPL_ASSERT(( boost::is_same<ExpectedResultNode, ResultExpNode> ));
            BOOST_STATIC_ASSERT(( boost::mpl::size<ResultExpIndices>::type::value == 4 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 0>::type::value == 0 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 1>::type::value == 1 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 2>::type::value == 2 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 3>::type::value == 3 ));
        }

        BOOST_AUTO_TEST_CASE(NonCommutativeNonAssociativeOperators)
        {
            // This test case concerns itself with operators that are not commutative nor associative.
            // No optimizations are possible in this case and tree should not undergo any transformations.

            typedef ExpressionTemplateTestObject00 Obj;
            typedef expt::Node<Obj> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> BinaryNode;

            // A + (B+C) -> A+ (B+C)
            typedef expt::Node<ConstantNode, expt::AddOp, BinaryNode> Exp1;
            typedef expt::RemoveUnecessaryTemporaries<Exp1, Exp1::Indices, 0>::TransformedNodeType ExpectedExp1Node;
            typedef expt::RemoveUnecessaryTemporaries<Exp1, Exp1::Indices, 0>::TransformedIndicesType ExpectedExp1Indices;

            BOOST_MPL_ASSERT(( boost::is_same<Exp1, ExpectedExp1Node> ));
            BOOST_STATIC_ASSERT(( boost::mpl::size<ExpectedExp1Indices>::type::value == 3 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ExpectedExp1Indices, 0>::type::value == 0 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ExpectedExp1Indices, 1>::type::value == 1 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ExpectedExp1Indices, 2>::type::value == 2 ));

            // (A+B) + (C+D) -> (A+B) + (C+D)
            typedef expt::Node<BinaryNode, expt::AddOp, BinaryNode> Exp2;
            typedef expt::RemoveUnecessaryTemporaries<Exp2, Exp2::Indices, 0>::TransformedNodeType ExpectedExp2Node;
            typedef expt::RemoveUnecessaryTemporaries<Exp2, Exp2::Indices, 0>::TransformedIndicesType ExpectedExp2Indices;

            BOOST_MPL_ASSERT(( boost::is_same<Exp2, ExpectedExp2Node> ));
            BOOST_STATIC_ASSERT(( boost::mpl::size<ExpectedExp2Indices>::type::value == 4 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ExpectedExp2Indices, 0>::type::value == 0 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ExpectedExp2Indices, 1>::type::value == 1 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ExpectedExp2Indices, 2>::type::value == 2 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ExpectedExp2Indices, 3>::type::value == 3 ));

            
        }




        ////////////////////////////////////////////////////
        // Expression 3: (A + (B+C)) + (D+(E+F))
        ///////////////////////////////////////////////////
        BOOST_AUTO_TEST_CASE(Exp3NonAssociativeNonCommutative)
        {
            // (A + (B+C)) + (D+(E+F)) -> (A + (B+C)) + (D+(E+F))
            typedef ExpressionTemplateTestObject00 Obj;
            typedef expt::Node<Obj> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> BinaryNode;
            typedef expt::Node<ConstantNode, expt::AddOp, BinaryNode> T0;
            typedef expt::Node<T0, expt::AddOp, T0> Exp;
            typedef expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::TransformedNodeType ResultExpNode;
            typedef expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::TransformedIndicesType ResultExpIndices;

            typedef Exp ExpectedResultNode;
            
            BOOST_MPL_ASSERT(( boost::is_same<ExpectedResultNode, ResultExpNode> ));
            BOOST_STATIC_ASSERT(( boost::mpl::size<ResultExpIndices>::type::value == 6 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 0>::type::value == 0 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 1>::type::value == 1 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 2>::type::value == 2 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 3>::type::value == 3 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 4>::type::value == 4 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 5>::type::value == 5 ));
        }

        BOOST_AUTO_TEST_CASE(Exp3NonAssociativeCommutative)
        {
            // (A + (B+C)) + (D+(E+F)) -> ((B+C)+A) + ((E+F)+D)
            typedef ExpressionTemplateTestObject01 Obj;
            typedef expt::Node<Obj> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> BinaryNode;
            typedef expt::Node<ConstantNode, expt::AddOp, BinaryNode> T0;
            typedef expt::Node<T0, expt::AddOp, T0> Exp;
            typedef expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::TransformedNodeType ResultExpNode;
            typedef expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::TransformedIndicesType ResultExpIndices;

            typedef expt::Node<BinaryNode, expt::AddOp, ConstantNode> T1;
            typedef expt::Node<T1, expt::AddOp, T1> ExpectedResultNode;

            BOOST_MPL_ASSERT(( boost::is_same<ExpectedResultNode, ResultExpNode> ));
            BOOST_STATIC_ASSERT(( boost::mpl::size<ResultExpIndices>::type::value == 6 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 0>::type::value == 1 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 1>::type::value == 2 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 2>::type::value == 0 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 3>::type::value == 4 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 4>::type::value == 5 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 5>::type::value == 3 ));
        }

        BOOST_AUTO_TEST_CASE(Exp3AssociativeNonCommutative)
        {
            // (A + (B+C)) + (D+(E+F)) -> ((((A+B)+C)+D)+E)+F
            typedef ExpressionTemplateTestObject10 Obj;
            typedef expt::Node<Obj> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> BinaryNode;
            typedef expt::Node<ConstantNode, expt::AddOp, BinaryNode> T0;
            typedef expt::Node<T0, expt::AddOp, T0> Exp;
            typedef expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::TransformedNodeType ResultExpNode;
            typedef expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::TransformedIndicesType ResultExpIndices;


            typedef expt::Node<BinaryNode, expt::AddOp, ConstantNode> ABC;
            typedef expt::Node<ABC, expt::AddOp, ConstantNode> ABCD;
            typedef expt::Node<ABCD, expt::AddOp, ConstantNode> ABCDE;
            typedef expt::Node<ABCDE, expt::AddOp, ConstantNode> ExpectedResultNode;
            
            BOOST_MPL_ASSERT(( boost::is_same<ExpectedResultNode, ResultExpNode> ));
            BOOST_STATIC_ASSERT(( boost::mpl::size<ResultExpIndices>::type::value == 6 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 0>::type::value == 0 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 1>::type::value == 1 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 2>::type::value == 2 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 3>::type::value == 3 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 4>::type::value == 4 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 5>::type::value == 5 ));
        }

        BOOST_AUTO_TEST_CASE(Exp3AssociativeCommutative)
        {
            // (A + (B+C)) + (D+(E+F)) -> ((((A+B)+C)+D)+E)+F
            typedef ExpressionTemplateTestObject11 Obj;
            typedef expt::Node<Obj> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> BinaryNode;
            typedef expt::Node<ConstantNode, expt::AddOp, BinaryNode> T0;
            typedef expt::Node<T0, expt::AddOp, T0> Exp;
            typedef expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::TransformedNodeType ResultExpNode;
            typedef expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::TransformedIndicesType ResultExpIndices;


            typedef expt::Node<BinaryNode, expt::AddOp, ConstantNode> ABC;
            typedef expt::Node<ABC, expt::AddOp, ConstantNode> ABCD;
            typedef expt::Node<ABCD, expt::AddOp, ConstantNode> ABCDE;
            typedef expt::Node<ABCDE, expt::AddOp, ConstantNode> ExpectedResultNode;
            
            BOOST_MPL_ASSERT(( boost::is_same<ExpectedResultNode, ResultExpNode> ));
            BOOST_STATIC_ASSERT(( boost::mpl::size<ResultExpIndices>::type::value == 6 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 0>::type::value == 0 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 1>::type::value == 1 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 2>::type::value == 2 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 3>::type::value == 3 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 4>::type::value == 4 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 5>::type::value == 5 ));
        }

        ////////////////////////////////////////////////////
        // Expression 4: (A + (B+C)) + (D+E))
        ///////////////////////////////////////////////////
        BOOST_AUTO_TEST_CASE(Exp4AssociativeNonCommutative)
        {
            // (A + (B+C)) + (D+E)) -> (((A+B)+C)+D)+E
            typedef ExpressionTemplateTestObject10 Obj;
            typedef expt::Node<Obj> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> BinaryNode;
            typedef expt::Node<ConstantNode, expt::AddOp, BinaryNode> T0;
            typedef expt::Node<T0, expt::AddOp, BinaryNode> Exp;
            typedef expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::TransformedNodeType ResultExpNode;
            typedef expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::TransformedIndicesType ResultExpIndices;

            typedef expt::Node<BinaryNode, expt::AddOp, ConstantNode> ABC;
            typedef expt::Node<ABC, expt::AddOp, ConstantNode> ABCD;
            typedef expt::Node<ABCD, expt::AddOp, ConstantNode> ExpectedResultNode;

            BOOST_MPL_ASSERT(( boost::is_same<ExpectedResultNode, ResultExpNode> ));
            BOOST_STATIC_ASSERT(( boost::mpl::size<ResultExpIndices>::type::value == 5 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 0>::type::value == 0 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 1>::type::value == 1 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 2>::type::value == 2 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 3>::type::value == 3 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 4>::type::value == 4 ));
        }

        ////////////////////////////////////////////////////
        // Expression 5: (A + (B+C)) + D
        ///////////////////////////////////////////////////
        BOOST_AUTO_TEST_CASE(Exp5AssociativeNonCommutative)
        {
            // (A + (B+C)) + D -> ((A+B)+C)+D
            typedef ExpressionTemplateTestObject10 Obj;
            typedef expt::Node<Obj> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> BinaryNode;
            typedef expt::Node<ConstantNode, expt::AddOp, BinaryNode> T0;
            typedef expt::Node<T0, expt::AddOp, ConstantNode> Exp;
            typedef expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::TransformedNodeType ResultExpNode;
            typedef expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::TransformedIndicesType ResultExpIndices;

            typedef expt::Node<BinaryNode, expt::AddOp, ConstantNode> ABC;
            typedef expt::Node<ABC, expt::AddOp, ConstantNode> ExpectedResultNode;

            BOOST_MPL_ASSERT(( boost::is_same<ExpectedResultNode, ResultExpNode> ));
            BOOST_STATIC_ASSERT(( boost::mpl::size<ResultExpIndices>::type::value == 4 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 0>::type::value == 0 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 1>::type::value == 1 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 2>::type::value == 2 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 3>::type::value == 3 ));
        }
      
        ////////////////////////////////////////////////////
        // Expression 6: (A+B) + CD
        ///////////////////////////////////////////////////
        BOOST_AUTO_TEST_CASE(Exp6NonAssociativeNonCommutative)
        {
            // (A+B) + CD -> (A+B) + CD
            typedef ExpressionTemplateTestObject00 Obj;
            typedef expt::Node<Obj> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> AddNode;
            typedef expt::Node<ConstantNode, expt::MultiplyOp, ConstantNode> MultiplyNode;
            typedef expt::Node<AddNode, expt::AddOp, MultiplyNode> Exp;
            typedef expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::TransformedNodeType ResultExpNode;
            typedef expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::TransformedIndicesType ResultExpIndices;


            typedef Exp ExpectedResultNode;
            
            BOOST_MPL_ASSERT(( boost::is_same<ExpectedResultNode, ResultExpNode> ));
            BOOST_STATIC_ASSERT(( boost::mpl::size<ResultExpIndices>::type::value == 4 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 0>::type::value == 0 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 1>::type::value == 1 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 2>::type::value == 2 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 3>::type::value == 3 ));
        }

        BOOST_AUTO_TEST_CASE(Exp6NonAssociativeCommutative)
        {
            // (A+B) + CD -> (A+B) + CD
            typedef ExpressionTemplateTestObject01 Obj;
            typedef expt::Node<Obj> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> AddNode;
            typedef expt::Node<ConstantNode, expt::MultiplyOp, ConstantNode> MultiplyNode;
            typedef expt::Node<AddNode, expt::AddOp, MultiplyNode> Exp;
            typedef expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::TransformedNodeType ResultExpNode;
            typedef expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::TransformedIndicesType ResultExpIndices;

            typedef Exp ExpectedResultNode;
            
            BOOST_MPL_ASSERT(( boost::is_same<ExpectedResultNode, ResultExpNode> ));
            BOOST_STATIC_ASSERT(( boost::mpl::size<ResultExpIndices>::type::value == 4 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 0>::type::value == 0 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 1>::type::value == 1 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 2>::type::value == 2 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 3>::type::value == 3 ));
        }

        BOOST_AUTO_TEST_CASE(Exp6AssociativeNonCommutative)
        {
            // (A+B) + CD -> (A+B) + CD
            typedef ExpressionTemplateTestObject10 Obj;
            typedef expt::Node<Obj> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> AddNode;
            typedef expt::Node<ConstantNode, expt::MultiplyOp, ConstantNode> MultiplyNode;
            typedef expt::Node<AddNode, expt::AddOp, MultiplyNode> Exp;
            typedef expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::TransformedNodeType ResultExpNode;
            typedef expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::TransformedIndicesType ResultExpIndices;

            typedef Exp ExpectedResultNode;
            
            BOOST_MPL_ASSERT(( boost::is_same<ExpectedResultNode, ResultExpNode> ));
            BOOST_STATIC_ASSERT(( boost::mpl::size<ResultExpIndices>::type::value == 4 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 0>::type::value == 0 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 1>::type::value == 1 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 2>::type::value == 2 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 3>::type::value == 3 ));
        }

        BOOST_AUTO_TEST_CASE(Exp6AssociativeCommutative)
        {
            // (A+B) + CD -> CD + (A+B) 
//            typedef ExpressionTemplateTestObject11 Obj;
//            typedef expt::Node<Obj> ConstantNode;
//            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> AddNode;
//            typedef expt::Node<ConstantNode, expt::MultiplyOp, ConstantNode> MultiplyNode;
//            typedef expt::Node<AddNode, expt::AddOp, MultiplyNode> Exp;
//            typedef expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::TransformedNodeType ResultExpNode;
//            typedef expt::RemoveUnecessaryTemporaries<Exp, Exp::Indices, 0>::TransformedIndicesType ResultExpIndices;
//
//            typedef expt::Node<MultiplyNode, expt::AddOp, ConstantNode> T0;
//            typedef expt::Node<T0, expt::AddOp, ConstantNode> ExpectedResultNode;
//
//            // Test the sorting of the ac cluster.
//            BOOST_MPL_ASSERT(( boost::is_same<RemoveUnecessaryTemporariesInternal<Exp, Exp::Indices, 0>::Tree0Type, Exp> ));
//            BOOST_MPL_ASSERT(( boost::is_same<RemoveUnecessaryTemporariesInternal<Exp, Exp::Indices, 0>::Tree1Type, Exp> ));
//            BOOST_MPL_ASSERT(( boost::is_same<RemoveUnecessaryTemporariesInternal<Exp, Exp::Indices, 0>::Tree2Type, Exp> ));
//            BOOST_MPL_ASSERT(( boost::is_same<RemoveUnecessaryTemporariesInternal<Exp, Exp::Indices, 0>::Tree3Type, Exp> ));
//            BOOST_STATIC_ASSERT(( SortAssociativeCommutativeClusters<Exp, Exp::Indices, 0>::SpecializationId == 1 ));
//            
//            BOOST_MPL_ASSERT(( boost::is_same<ExpectedResultNode, ResultExpNode> ));
//            BOOST_STATIC_ASSERT(( boost::mpl::size<ResultExpIndices>::type::value == 4 ));
//            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 0>::type::value == 2 ));
//            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 1>::type::value == 3 ));
//            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 2>::type::value == 0 ));
//            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 3>::type::value == 1 ));
        }


        //////////////////////////////////////////////////////////
        // Expression 7: (A + BC) + EF -> (BC + A) + EF
        //////////////////////////////////////////////////////////
    }
}
