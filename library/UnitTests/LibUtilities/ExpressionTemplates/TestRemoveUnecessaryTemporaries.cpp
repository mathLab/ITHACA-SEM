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

using namespace expt;

namespace Nektar
{
    namespace UnitTests
    {
        BOOST_AUTO_TEST_CASE(BaseCases)
        {
            typedef NekMatrix<double> Matrix;

            // A -> A
            typedef expt::Node<Matrix> ConstantNode;
            typedef expt::RemoveUnecessaryTemporaries<ConstantNode>::TransformedNodeType ExpectedConstantNode;

            BOOST_MPL_ASSERT(( boost::is_same<ConstantNode, ExpectedConstantNode> ));

            // -A -> -A
            typedef expt::Node<ConstantNode, expt::NegateOp> UnaryNode;
            typedef expt::RemoveUnecessaryTemporaries<UnaryNode>::TransformedNodeType ExpectedUnaryNode;
            BOOST_MPL_ASSERT(( boost::is_same<UnaryNode, ExpectedUnaryNode> ));

            // A+B -> A+B
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> BinaryAdd;
            typedef expt::RemoveUnecessaryTemporaries<BinaryAdd>::TransformedNodeType ExpectedBinaryAddNode;
            typedef expt::RemoveUnecessaryTemporaries<BinaryAdd>::TransformedIndicesType ExpectedBinaryAddIndices;

            BOOST_MPL_ASSERT(( boost::is_same<BinaryAdd, ExpectedBinaryAddNode> ));
            BOOST_STATIC_ASSERT(( boost::mpl::size<ExpectedBinaryAddIndices>::type::value == 2 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ExpectedBinaryAddIndices, 0>::type::value == 0 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ExpectedBinaryAddIndices, 1>::type::value == 1 ));

            // A*B -> A*B
            typedef expt::Node<ConstantNode, expt::MultiplyOp, ConstantNode> BinaryMultiply;
            typedef expt::RemoveUnecessaryTemporaries<BinaryMultiply>::TransformedNodeType ExpectedBinaryMultiplyNode;
            typedef expt::RemoveUnecessaryTemporaries<BinaryMultiply>::TransformedIndicesType ExpectedBinaryMultiplyIndices;

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
            typedef TestObject Obj;
            typedef expt::Node<Obj> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> BinaryNode;

            typedef expt::Node<ConstantNode, expt::AddOp, BinaryNode> Exp;
            typedef expt::RemoveUnecessaryTemporaries<Exp>::TransformedNodeType ResultExpNode;
            typedef expt::RemoveUnecessaryTemporaries<Exp>::TransformedIndicesType ResultExpIndices;

            BOOST_MPL_ASSERT(( boost::is_same<Exp, ResultExpNode> ));
            BOOST_STATIC_ASSERT(( boost::mpl::size<ResultExpIndices>::type::value == 3 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 0>::type::value == 0 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 1>::type::value == 1 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 2>::type::value == 2 ));
        }

        BOOST_AUTO_TEST_CASE(Exp1NonAssociativeCommutative)
        {
            // A + (B+C) -> (B+C) + A
            typedef TestObjectC Obj;
            typedef expt::Node<Obj> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> BinaryNode;

            typedef expt::Node<ConstantNode, expt::AddOp, BinaryNode> Exp;
            typedef expt::RemoveUnecessaryTemporaries<Exp>::TransformedNodeType ResultExpNode;
            typedef expt::RemoveUnecessaryTemporaries<Exp>::TransformedIndicesType ResultExpIndices;
            typedef expt::Node<BinaryNode, expt::AddOp, ConstantNode> ExpectedResultNode;

            // Test the associative cluster sorting associated with this expression.
            BOOST_MPL_ASSERT(( boost::is_same<expt::RemoveUnecessaryTemporaries<Exp>::Tree0Type, ExpectedResultNode> ));
            BOOST_MPL_ASSERT(( boost::is_same<expt::RemoveUnecessaryTemporaries<Exp>::TransformedNodeType, ExpectedResultNode> ));

            BOOST_MPL_ASSERT(( boost::is_same<ExpectedResultNode, ResultExpNode> ));
            BOOST_STATIC_ASSERT(( boost::mpl::size<ResultExpIndices>::type::value == 3 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 0>::type::value == 1 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 1>::type::value == 2 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 2>::type::value == 0 ));
        }


        BOOST_AUTO_TEST_CASE(Exp1AssociativeNonCommutative)
        {
            // A + (B+C) -> (A+B) + C
            typedef TestObjectA Obj;
            typedef expt::Node<Obj> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> BinaryNode;

            typedef expt::Node<ConstantNode, expt::AddOp, BinaryNode> Exp;
            typedef expt::RemoveUnecessaryTemporaries<Exp>::TransformedNodeType ResultExpNode;
            typedef expt::RemoveUnecessaryTemporaries<Exp>::TransformedIndicesType ResultExpIndices;
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
            typedef TestObjectAC Obj;
            typedef expt::Node<Obj> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> BinaryNode;

            typedef expt::Node<ConstantNode, expt::AddOp, BinaryNode> Exp;
            typedef expt::RemoveUnecessaryTemporaries<Exp>::TransformedNodeType ResultExpNode;
            typedef expt::RemoveUnecessaryTemporaries<Exp>::TransformedIndicesType ResultExpIndices;
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
            typedef TestObject Obj;
            typedef expt::Node<Obj> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> BinaryNode;

            typedef expt::Node<BinaryNode, expt::AddOp, BinaryNode> Exp;
            typedef expt::RemoveUnecessaryTemporaries<Exp>::TransformedNodeType ResultExpNode;
            typedef expt::RemoveUnecessaryTemporaries<Exp>::TransformedIndicesType ResultExpIndices;
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
            typedef TestObjectC Obj;
            typedef expt::Node<Obj> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> BinaryNode;

            typedef expt::Node<BinaryNode, expt::AddOp, BinaryNode> Exp;
            typedef expt::RemoveUnecessaryTemporaries<Exp>::TransformedNodeType ResultExpNode;
            typedef expt::RemoveUnecessaryTemporaries<Exp>::TransformedIndicesType ResultExpIndices;
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
            typedef TestObjectA Obj;
            typedef expt::Node<Obj> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> BinaryNode;

            typedef expt::Node<BinaryNode, expt::AddOp, BinaryNode> Exp;
            typedef expt::RemoveUnecessaryTemporaries<Exp>::TransformedNodeType ResultExpNode;
            typedef expt::RemoveUnecessaryTemporaries<Exp>::TransformedIndicesType ResultExpIndices;
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
            typedef TestObjectAC Obj;
            typedef expt::Node<Obj> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> BinaryNode;

            typedef expt::Node<BinaryNode, expt::AddOp, BinaryNode> Exp;
            typedef expt::RemoveUnecessaryTemporaries<Exp>::TransformedNodeType ResultExpNode;
            typedef expt::RemoveUnecessaryTemporaries<Exp>::TransformedIndicesType ResultExpIndices;
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

            typedef TestObject Obj;
            typedef expt::Node<Obj> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> BinaryNode;

            // A + (B+C) -> A+ (B+C)
            typedef expt::Node<ConstantNode, expt::AddOp, BinaryNode> Exp1;
            typedef expt::RemoveUnecessaryTemporaries<Exp1>::TransformedNodeType ExpectedExp1Node;
            typedef expt::RemoveUnecessaryTemporaries<Exp1>::TransformedIndicesType ExpectedExp1Indices;

            BOOST_MPL_ASSERT(( boost::is_same<Exp1, ExpectedExp1Node> ));
            BOOST_STATIC_ASSERT(( boost::mpl::size<ExpectedExp1Indices>::type::value == 3 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ExpectedExp1Indices, 0>::type::value == 0 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ExpectedExp1Indices, 1>::type::value == 1 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ExpectedExp1Indices, 2>::type::value == 2 ));

            // (A+B) + (C+D) -> (A+B) + (C+D)
            typedef expt::Node<BinaryNode, expt::AddOp, BinaryNode> Exp2;
            typedef expt::RemoveUnecessaryTemporaries<Exp2>::TransformedNodeType ExpectedExp2Node;
            typedef expt::RemoveUnecessaryTemporaries<Exp2>::TransformedIndicesType ExpectedExp2Indices;

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
            typedef TestObject Obj;
            typedef expt::Node<Obj> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> BinaryNode;
            typedef expt::Node<ConstantNode, expt::AddOp, BinaryNode> T0;
            typedef expt::Node<T0, expt::AddOp, T0> Exp;
            typedef expt::RemoveUnecessaryTemporaries<Exp>::TransformedNodeType ResultExpNode;
            typedef expt::RemoveUnecessaryTemporaries<Exp>::TransformedIndicesType ResultExpIndices;

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
            typedef TestObjectC Obj;
            typedef expt::Node<Obj> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> BinaryNode;
            typedef expt::Node<ConstantNode, expt::AddOp, BinaryNode> T0;
            typedef expt::Node<T0, expt::AddOp, T0> Exp;
            typedef expt::RemoveUnecessaryTemporaries<Exp>::TransformedNodeType ResultExpNode;
            typedef expt::RemoveUnecessaryTemporaries<Exp>::TransformedIndicesType ResultExpIndices;

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
            typedef TestObjectA Obj;
            typedef expt::Node<Obj> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> BinaryNode;
            typedef expt::Node<ConstantNode, expt::AddOp, BinaryNode> T0;
            typedef expt::Node<T0, expt::AddOp, T0> Exp;
            typedef expt::RemoveUnecessaryTemporaries<Exp>::TransformedNodeType ResultExpNode;
            typedef expt::RemoveUnecessaryTemporaries<Exp>::TransformedIndicesType ResultExpIndices;


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
            typedef TestObjectAC Obj;
            typedef expt::Node<Obj> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> BinaryNode;
            typedef expt::Node<ConstantNode, expt::AddOp, BinaryNode> T0;
            typedef expt::Node<T0, expt::AddOp, T0> Exp;
            typedef expt::RemoveUnecessaryTemporaries<Exp>::TransformedNodeType ResultExpNode;
            typedef expt::RemoveUnecessaryTemporaries<Exp>::TransformedIndicesType ResultExpIndices;


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
            typedef TestObjectA Obj;
            typedef expt::Node<Obj> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> BinaryNode;
            typedef expt::Node<ConstantNode, expt::AddOp, BinaryNode> T0;
            typedef expt::Node<T0, expt::AddOp, BinaryNode> Exp;
            typedef expt::RemoveUnecessaryTemporaries<Exp>::TransformedNodeType ResultExpNode;
            typedef expt::RemoveUnecessaryTemporaries<Exp>::TransformedIndicesType ResultExpIndices;

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
            typedef TestObjectA Obj;
            typedef expt::Node<Obj> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> BinaryNode;
            typedef expt::Node<ConstantNode, expt::AddOp, BinaryNode> T0;
            typedef expt::Node<T0, expt::AddOp, ConstantNode> Exp;
            typedef expt::RemoveUnecessaryTemporaries<Exp>::TransformedNodeType ResultExpNode;
            typedef expt::RemoveUnecessaryTemporaries<Exp>::TransformedIndicesType ResultExpIndices;

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
            typedef TestObject Obj;
            typedef expt::Node<Obj> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> AddNode;
            typedef expt::Node<ConstantNode, expt::MultiplyOp, ConstantNode> MultiplyNode;
            typedef expt::Node<AddNode, expt::AddOp, MultiplyNode> Exp;
            typedef expt::RemoveUnecessaryTemporaries<Exp>::TransformedNodeType ResultExpNode;
            typedef expt::RemoveUnecessaryTemporaries<Exp>::TransformedIndicesType ResultExpIndices;


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
            typedef TestObjectC Obj;
            typedef expt::Node<Obj> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> AddNode;
            typedef expt::Node<ConstantNode, expt::MultiplyOp, ConstantNode> MultiplyNode;
            typedef expt::Node<AddNode, expt::AddOp, MultiplyNode> Exp;
            typedef expt::RemoveUnecessaryTemporaries<Exp>::TransformedNodeType ResultExpNode;
            typedef expt::RemoveUnecessaryTemporaries<Exp>::TransformedIndicesType ResultExpIndices;

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
            typedef TestObjectA Obj;
            typedef expt::Node<Obj> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> AddNode;
            typedef expt::Node<ConstantNode, expt::MultiplyOp, ConstantNode> MultiplyNode;
            typedef expt::Node<AddNode, expt::AddOp, MultiplyNode> Exp;
            typedef expt::RemoveUnecessaryTemporaries<Exp>::TransformedNodeType ResultExpNode;
            typedef expt::RemoveUnecessaryTemporaries<Exp>::TransformedIndicesType ResultExpIndices;

            typedef Exp ExpectedResultNode;
            
            BOOST_MPL_ASSERT(( boost::is_same<ExpectedResultNode, ResultExpNode> ));
            BOOST_STATIC_ASSERT(( boost::mpl::size<ResultExpIndices>::type::value == 4 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 0>::type::value == 0 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 1>::type::value == 1 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 2>::type::value == 2 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 3>::type::value == 3 ));
        }

        BOOST_AUTO_TEST_CASE(TestFullWithAssociativeCluster)
        {
            // ( A + (BC+D) ) + E
            typedef Node<NekVector<double> > V;
            typedef Node<V, MultiplyOp, V> M;
            typedef Node<M, AddOp, V> BCDNode;
            typedef Node<V, AddOp, BCDNode> ABCDNode;
            typedef Node<ABCDNode, AddOp, V> Test1;

            typedef RemoveUnecessaryTemporaries<Test1>::TransformedNodeType Test1NodeType;
            typedef RemoveUnecessaryTemporaries<Test1>::TransformedIndicesType Test1IndicesType;

            // Expected is  ((BC + A) + D) + E
            typedef Node<Node<Node<M, AddOp, V>, AddOp, V>, AddOp, V> ExpectedTest1NodeType;

            BOOST_MPL_ASSERT(( boost::is_same<Test1NodeType, ExpectedTest1NodeType> ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test1IndicesType, 0>::type::value == 1));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test1IndicesType, 1>::type::value == 2));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test1IndicesType, 2>::type::value == 0));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test1IndicesType, 3>::type::value == 3));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test1IndicesType, 4>::type::value == 4));
            
        }

        BOOST_AUTO_TEST_CASE(TestWithTwoFullWithAssociativeClusters)
        {
            // A + ( ( (BC) * (D+E) ) + (F+G) )
            typedef Node<NekVector<double> > V;
            typedef Node<V, MultiplyOp, V> M;
            typedef Node<V, AddOp, V> A;
            typedef Node<M, MultiplyOp, A> BCDENode;
            typedef Node<BCDENode, AddOp, A> RhsNode;
            typedef Node<V, AddOp, RhsNode> Test1;

            BOOST_STATIC_ASSERT(( TemporaryCount<Test1>::Value == 3 ));
            
            typedef RemoveUnecessaryTemporaries<Test1>::Tree0Type Tree0Type;
            typedef RemoveUnecessaryTemporaries<Test1>::T0 T0;
            typedef RemoveUnecessaryTemporaries<Test1>::T1 T1;
            typedef RemoveUnecessaryTemporaries<Test1>::T2 T2;
            typedef impl::RemoveUnecessaryTemporariesInternal<T2, T2::Indices, 0>::RightNode0Type RightNode0Type;
            typedef impl::RemoveUnecessaryTemporariesInternal<T2, T2::Indices, 0>::Tree0Type Tree0InternalType;
            typedef impl::RemoveUnecessaryTemporariesInternal<T2, T2::Indices, 0>::Tree1Type Tree1InternalType;
            typedef impl::RemoveUnecessaryTemporariesInternal<T2, T2::Indices, 0>::LeftNode1Type LeftNode1InternalType;
            typedef impl::RemoveUnecessaryTemporariesInternal<T2, T2::Indices, 0>::LeftNode2Type LeftNode2InternalType;
            typedef impl::RemoveUnecessaryTemporariesInternal<T2, T2::Indices, 0>::Tree3Type Tree3InternalType;

            BOOST_MPL_ASSERT(( boost::is_same<T0, Test1> ));
            BOOST_MPL_ASSERT(( boost::is_same<T1, Test1> ));
            BOOST_MPL_ASSERT(( boost::is_same<T2, Test1> ));

            typedef RemoveUnecessaryTemporaries<Test1>::TransformedNodeType Test1NodeType;
            typedef RemoveUnecessaryTemporaries<Test1>::TransformedIndicesType Test1IndicesType;

            BOOST_STATIC_ASSERT(( TemporaryCount<Test1NodeType>::Value == 0 ));
            BOOST_STATIC_ASSERT(( TemporaryCount<Tree0Type>::Value == 1 ));

            typedef Node<Node<BCDENode, AddOp, V>, AddOp, V> ExpectedRightNode0Type;
            typedef Node<V, AddOp, ExpectedRightNode0Type> ExpectedTree0InternalType;
            typedef Node<BCDENode, AddOp, V> Internal1;
            typedef Node<Node<V, AddOp, Internal1>, AddOp, V> IncorrectExpectedTree1InternalType;
            typedef Node<V, AddOp, Internal1> ExpectedLeftNode1InternalType;
            typedef Node<Node<BCDENode, AddOp, V>, AddOp, V> ExpectedLeftNode2InternalType;

            typedef Node<ExpectedLeftNode2InternalType, AddOp, V> ExectedTree3InternalType;

            typedef Node<A, MultiplyOp, V> A0;
            typedef Node<A0, MultiplyOp, V> A1;
            typedef Node<A1, AddOp, V> A2;
            typedef Node<A2, AddOp, V> A3;
            typedef Node<A3, AddOp, V> ExpectedTest1NodeType;


            BOOST_MPL_ASSERT(( boost::is_same<RightNode0Type, ExpectedRightNode0Type> ));
            BOOST_MPL_ASSERT(( boost::is_same<Tree0InternalType, ExpectedTree0InternalType> ));
            BOOST_MPL_ASSERT(( boost::is_same<Tree1InternalType, IncorrectExpectedTree1InternalType> ));
            BOOST_MPL_ASSERT(( boost::is_same<LeftNode1InternalType, ExpectedLeftNode1InternalType> ));
            BOOST_MPL_ASSERT(( boost::is_same<LeftNode2InternalType, ExpectedLeftNode2InternalType> ));
            BOOST_MPL_ASSERT(( boost::is_same<Tree3InternalType, ExectedTree3InternalType> ));
            BOOST_MPL_ASSERT((boost::is_same<Test1NodeType, ExpectedTest1NodeType> ));

            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test1IndicesType, 0>::type::value == 3));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test1IndicesType, 1>::type::value == 4));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test1IndicesType, 2>::type::value == 2));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test1IndicesType, 3>::type::value == 1));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test1IndicesType, 4>::type::value == 0));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test1IndicesType, 5>::type::value == 5));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test1IndicesType, 6>::type::value == 6));
            
        }

        BOOST_AUTO_TEST_CASE(TestWithTwoFullWithAssociativeClusters1)
        {
            // A + ( (BC*(D+E)) + F )
            typedef Node<NekVector<double> > V;
            typedef Node<V, MultiplyOp, V> M;
            typedef Node<V, AddOp, V> A;
            typedef Node<M, MultiplyOp, A> BCDENode;
            typedef Node<BCDENode, AddOp, V> RhsNode;
            typedef Node<V, AddOp, RhsNode> Test1;

            typedef RemoveUnecessaryTemporaries<Test1>::Tree0Type Tree0Type;
            typedef RemoveUnecessaryTemporaries<Test1>::T0 T0;
            typedef RemoveUnecessaryTemporaries<Test1>::T1 T1;
            typedef RemoveUnecessaryTemporaries<Test1>::T2 T2;
            typedef impl::RemoveUnecessaryTemporariesInternal<T2, T2::Indices, 0>::RightNode0Type RightNode0Type;
            typedef impl::RemoveUnecessaryTemporariesInternal<T2, T2::Indices, 0>::Tree0Type Tree0InternalType;
            typedef impl::RemoveUnecessaryTemporariesInternal<T2, T2::Indices, 0>::Tree1Type Tree1InternalType;
            typedef impl::RemoveUnecessaryTemporariesInternal<T2, T2::Indices, 0>::LeftNode1Type LeftNode1InternalType;
            typedef impl::RemoveUnecessaryTemporariesInternal<T2, T2::Indices, 0>::LeftNode2Type LeftNode2InternalType;
            typedef impl::RemoveUnecessaryTemporariesInternal<T2, T2::Indices, 0>::Tree3Type Tree3InternalType;
            typedef RemoveUnecessaryTemporaries<Test1>::TransformedNodeType Test1NodeType;
            typedef RemoveUnecessaryTemporaries<Test1>::TransformedIndicesType Test1IndicesType;


            BOOST_MPL_ASSERT(( boost::is_same<T0, Test1> ));
            BOOST_MPL_ASSERT(( boost::is_same<T1, Test1> ));
            BOOST_MPL_ASSERT(( boost::is_same<T2, Test1> ));
            BOOST_MPL_ASSERT(( boost::is_same<RhsNode, RightNode0Type> ));
            BOOST_MPL_ASSERT(( boost::is_same<Test1, Tree0InternalType> ));

            typedef Node<Node<V, AddOp, BCDENode>, AddOp, V> ExpectedTree1InternalType;
            BOOST_MPL_ASSERT(( boost::is_same<ExpectedTree1InternalType, Tree1InternalType> ));

            typedef Node<BCDENode, AddOp, V> ExpectedLeftNode2InternalType;
            BOOST_MPL_ASSERT(( boost::is_same<ExpectedLeftNode2InternalType, LeftNode2InternalType> ));

            typedef Node<A, MultiplyOp, V> A0;
            typedef Node<A0, MultiplyOp, V> A1;
            typedef Node<A1, AddOp, V> A2;
            typedef Node<A2, AddOp, V> ExpectedTest1NodeType;
            BOOST_MPL_ASSERT(( boost::is_same<ExpectedTest1NodeType, Test1NodeType> ));


            
        }

        BOOST_AUTO_TEST_CASE(Exp6AssociativeCommutative)
        {
            // (A+B) + CD -> CD + (A+B) 
            typedef TestObjectAC Obj;
            typedef expt::Node<Obj> ConstantNode;
            typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> AddNode;
            typedef expt::Node<ConstantNode, expt::MultiplyOp, ConstantNode> MultiplyNode;
            typedef expt::Node<AddNode, expt::AddOp, MultiplyNode> Exp;
            typedef expt::RemoveUnecessaryTemporaries<Exp>::TransformedNodeType ResultExpNode;
            typedef expt::RemoveUnecessaryTemporaries<Exp>::TransformedIndicesType ResultExpIndices;

            typedef expt::Node<MultiplyNode, expt::AddOp, ConstantNode> T0;
            typedef expt::Node<T0, expt::AddOp, ConstantNode> ExpectedResultNode;

            // Test the sorting of the ac cluster.
            //BOOST_MPL_ASSERT(( boost::is_same<RemoveUnecessaryTemporariesInternal<Exp>::Tree0Type, Exp> ));
            //BOOST_MPL_ASSERT(( boost::is_same<RemoveUnecessaryTemporariesInternal<Exp>::Tree1Type, Exp> ));
            //BOOST_MPL_ASSERT(( boost::is_same<RemoveUnecessaryTemporariesInternal<Exp>::Tree2Type, Exp> ));
            //BOOST_MPL_ASSERT(( boost::is_same<RemoveUnecessaryTemporariesInternal<Exp>::Tree3Type, Exp> ));
            
            BOOST_MPL_ASSERT(( boost::is_same<ExpectedResultNode, ResultExpNode> ));
            BOOST_STATIC_ASSERT(( boost::mpl::size<ResultExpIndices>::type::value == 4 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 0>::type::value == 2 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 1>::type::value == 3 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 2>::type::value == 1 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<ResultExpIndices, 3>::type::value == 0 ));
        }


        //////////////////////////////////////////////////////////
        // Expression 7: (A + BC) + EF -> (BC + A) + EF
        //////////////////////////////////////////////////////////

        // For the expression A-B, the forward inverse 
        // transform creates A + (-B)
        BOOST_AUTO_TEST_CASE(TestForwardAndBackwardTransforms)
        {
            typedef Node<NekVector<double> > V;
            typedef Node<V, SubtractOp, V> TestNode;

            typedef ForwardInverseTransform<TestNode>::Type Test1;
            typedef Test1::Indices Test1IndicesType;
            typedef Node<V, AddOp, Node<V, NegateOp> > ExpectedTest1;

            BOOST_MPL_ASSERT(( boost::is_same<Test1, ExpectedTest1> ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test1IndicesType, 0>::type::value == 0));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<Test1IndicesType, 1>::type::value == 1));

            typedef BackwardInverseTransform<Test1, Test1::Indices, 0>::TransformedNodeType Test2;
            typedef Node<Node<V, NegateOp>, AddOp, V> ExpectedTest2;

            BOOST_MPL_ASSERT(( boost::is_same<Test2, TestNode> ));

            typedef RemoveUnecessaryTemporaries<TestNode>::TransformedNodeType Test3;
            BOOST_MPL_ASSERT(( boost::is_same<Test3, TestNode> ));

            // Motivation was incorrect temporary in v-v
            double v0_buf[] = {1, 2};
            double v1_buf[] = {3, 4};
            NekVector<double> v0(2, v0_buf);
            NekVector<double> v1(2, v1_buf);
            NekVector<double> result = v0-v1;

            double r_buf[] = {1-3, 2-4};
            NekVector<double> expectedResult(2, r_buf);
            BOOST_CHECK_EQUAL(result, expectedResult);
        }
    }
}
