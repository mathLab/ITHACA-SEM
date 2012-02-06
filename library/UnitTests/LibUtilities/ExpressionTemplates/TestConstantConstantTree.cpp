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
#include <ExpressionTemplates/Node.hpp>
#include <ExpressionTemplates/RemoveAllUnecessaryTemporaries.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/type_traits.hpp>
#include <LibUtilities/LinearAlgebra/MatrixSize.hpp>

namespace Nektar
{
    namespace UnitTests
    {
        BOOST_AUTO_TEST_CASE(TestMatrixSizeConstantNode)
        {
            typedef NekMatrix<NekDouble> Matrix;
            typedef NekVector<NekDouble> Vector;

            typedef expt::Node<Matrix> LhsNode;
            typedef expt::Node<Vector> RhsNode;
            typedef expt::Node<LhsNode, expt::MultiplyOp, RhsNode> Expression;
            typedef Expression::Indices Indices;

            double m_buf[] = {1, 4, 2, 5, 3, 6};
            double v_buf[] = {1, 2, 3};
            Matrix m(2, 3, m_buf);
            Vector v(3, v_buf);

            Expression e = m*v;
            
            boost::tuple<unsigned int, unsigned int, unsigned int> sizes = 
                MatrixSize<LhsNode, Indices, 0>::GetRequiredSize(e.GetData());

            unsigned int rows = sizes.get<0>();
            unsigned int cols = sizes.get<1>();
            unsigned int bufferSize = sizes.get<2>();
            
            BOOST_CHECK_EQUAL(rows, 2u);
            BOOST_CHECK_EQUAL(cols, 3u);
            BOOST_CHECK_EQUAL(bufferSize, rows*cols);
    
            rows = 0;
            cols = 0;
            bufferSize = 0;
            sizes = MatrixSize<RhsNode, Indices, 1>::GetRequiredSize(e.GetData());
            rows = sizes.get<0>();
            cols = sizes.get<1>();
            bufferSize = sizes.get<2>();
            
            BOOST_CHECK_EQUAL(rows, 3u);
            BOOST_CHECK_EQUAL(cols, 1u);
            BOOST_CHECK_EQUAL(bufferSize, rows*cols);

            sizes = MatrixSize<Expression, Indices, 0>::GetRequiredSize(e.GetData());
            rows = sizes.get<0>();
            cols = sizes.get<1>();
            bufferSize = sizes.get<2>();
            BOOST_CHECK_EQUAL(rows, 2u);
            BOOST_CHECK_EQUAL(cols, 1u);
            BOOST_CHECK_EQUAL(bufferSize, rows*cols);
        }

        
        BOOST_AUTO_TEST_CASE(TestConstantConstantTree)
        {
            typedef expt::Node<CountedObject<double> > LhsNode;
            BOOST_MPL_ASSERT(( boost::is_same<LhsNode::ResultType, CountedObject<double> > ));
            BOOST_MPL_ASSERT(( boost::is_same<LhsNode::DataType, const CountedObject<double>&> ));
            
            typedef expt::Node<CountedObject<double> > RhsNode;
            typedef expt::Node<LhsNode, expt::AddOp, RhsNode> Expression;
            typedef Expression::Indices Indices;

            typedef expt::RemoveUnecessaryTemporaries<Expression>::TransformedNodeType TransformedNodeType;
            BOOST_MPL_ASSERT(( boost::is_same<Expression, TransformedNodeType> ));
            
            typedef expt::RemoveUnecessaryTemporaries<Expression>::TransformedIndicesType TransformedIndicesType;
            BOOST_STATIC_ASSERT(( boost::mpl::size<TransformedIndicesType>::type::value == 2 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<TransformedIndicesType, 0>::type::value == 0 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<TransformedIndicesType, 1>::type::value == 1 ));

#if BOOST_VERSION < 104200
            typedef boost::fusion::vector<const CountedObject<double>&, const CountedObject<double>&> ExpectedListType;
#else
            typedef boost::fusion::vector2<const CountedObject<double>&, const CountedObject<double>&> ExpectedListType;
#endif
            BOOST_MPL_ASSERT(( boost::is_same<Expression::VectorType, ExpectedListType> ));
            
            CountedObject<double> lhs(3);
            CountedObject<double> rhs(17);
            CountedObject<double>::ClearCounters();

            Expression e = lhs + rhs;
            CountedObject<double>::Check(0, 0, 0, 0, 0, 0);
            CountedObject<double> result = lhs + rhs;
            CountedObject<double>::Check(0, 0, 0, 0, 0, 0);
            BOOST_CHECK_EQUAL( CountedObject<double>::numberOfExpressionConstructions, 1 );
            BOOST_CHECK_EQUAL( CountedObject<double>::numberOfExpressionAssignments, 0 );

            BOOST_CHECK_EQUAL(result.value, lhs.value + rhs.value);
        }

        BOOST_AUTO_TEST_CASE(TestConstantConstantTreeVaryDataType)
        {
            typedef NekMatrix<NekDouble> Matrix;
            typedef NekVector<NekDouble> Vector;

            typedef expt::Node<Matrix> LhsNode;
            BOOST_MPL_ASSERT(( boost::is_same<LhsNode::ResultType, Matrix > ));
            BOOST_MPL_ASSERT(( boost::is_same<LhsNode::DataType, const Matrix&> ));
            
            typedef expt::Node<Vector> RhsNode;
            typedef expt::Node<LhsNode, expt::MultiplyOp, RhsNode> Expression;
            typedef Expression::Indices Indices;

            typedef expt::RemoveUnecessaryTemporaries<Expression>::TransformedNodeType TransformedNodeType;
            BOOST_MPL_ASSERT(( boost::is_same<Expression, TransformedNodeType> ));
            
            typedef expt::RemoveUnecessaryTemporaries<Expression>::TransformedIndicesType TransformedIndicesType;
            BOOST_STATIC_ASSERT(( boost::mpl::size<TransformedIndicesType>::type::value == 2 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<TransformedIndicesType, 0>::type::value == 0 ));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<TransformedIndicesType, 1>::type::value == 1 ));

#if BOOST_VERSION < 104200
            typedef boost::fusion::vector<const Matrix&, const Vector&> ExpectedListType;
#else
            typedef boost::fusion::vector2<const Matrix&, const Vector&> ExpectedListType;
#endif
            BOOST_MPL_ASSERT(( boost::is_same<Expression::VectorType, ExpectedListType> ));
            
            double m_buf[] = {1, 4, 2, 5, 3, 6};
            double v_buf[] = {1, 2, 3};
            Matrix m(2, 3, m_buf);
            Vector v(3, v_buf);

            Expression e = m*v;
            Vector result = m*v;

            double expected_result_buf[] = {1+4+9, 4+10+18};
            Vector expected_result(2, expected_result_buf); 
            BOOST_CHECK_EQUAL(result, expected_result);
        }

    }

    BOOST_AUTO_TEST_CASE(TestBinaryBinaryMatrixTree)
    {
        typedef NekMatrix<NekDouble> Matrix;
        double m_buf[] = {1, 4, 2, 5};
        Matrix m1(2, 2, m_buf);
        Matrix m2(2, 2, m_buf);
        Matrix m3(2, 2, m_buf);
        Matrix m4(2, 2, m_buf);

        Array<TwoD, const NekDouble> gmat;

        Matrix result = m1*m2 + m3*m4;
        //Matrix result1 = (gmat[0][0]*gmat[0][0]+gmat[2][0]*gmat[2][0]+gmat[4][0]*gmat[4][0])* (m1) + 
        //    (gmat[0][0]*gmat[1][0] + gmat[2][0]*gmat[3][0] + gmat[4][0]*gmat[5][0])*(m2 + Transpose(m2)) +
        //    (gmat[1][0]*gmat[1][0] + gmat[3][0]*gmat[3][0] + gmat[5][0]*gmat[5][0])* (m3);
        //Matrix result2 = (1.0)*(m1) + (2.0)*(m2) + (3.0)*(m3);
        
        //Matrix result3 = m1*m1 + m2*m2 + m3*m3;
        typedef expt::Node<expt::Node<Matrix>, expt::MultiplyOp, expt::Node<Matrix> > T;
        typedef expt::Node<T, expt::AddOp, T> Lhs;
        typedef T Rhs;
        typedef expt::Node<Lhs, expt::AddOp, Rhs> ExpressionType;
        ExpressionType e = m1*m1 + m2*m2 + m3*m3;
        typedef ExpressionType::Indices Indices;

        //typedef expt::RemoveUnecessaryTemporaries<ExpressionType, Indices>::TransformedNodeType TransformedNodeType;
        typedef expt::RemoveUnecessaryTemporaries<ExpressionType>::TransformedIndicesType TransformedIndicesType;
    }
}
