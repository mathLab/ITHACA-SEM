///////////////////////////////////////////////////////////////////////////////
//
// File: TestMatrixBufferSize.cpp
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
// Description: 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_USE_EXPRESSION_TEMPLATES
#define NEKTAR_USE_EXPRESSION_TEMPLATES
#endif //NEKTAR_USE_EXPRESSION_TEMPLATES

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <LibUtilities/BasicUtils/OperatorGenerators.hpp>
#include <UnitTests/CountedObject.h>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>

namespace Nektar
{

    BOOST_AUTO_TEST_CASE(TestLargerIntermediateMatrix)
    {
        double buf1[] = {1, 4,
                         2, 5,
                        3, 6};
        double buf2[] = {7, 10, 13,
                        8, 11, 14,
                        9, 12, 15,
                        16, 17, 18,
                        19, 20, 21};
        double buf3[] = {1, 2, 3, 4, 5};
        
        typedef NekMatrix<double> Matrix;
        typedef expt::Node<Matrix> MatrixNode;
        typedef expt::Node<MatrixNode, expt::MultiplyOp, MatrixNode> LhsTree;
        typedef MatrixNode RhsTree;
        typedef expt::Node<LhsTree, expt::MultiplyOp, RhsTree> Expression;
        typedef Expression::Indices Indices;
        
        NekMatrix<double> m1(2, 3, buf1);
        NekMatrix<double> m2(3, 5, buf2);
        NekMatrix<double> tempResult = m1*m2;
        
        double temp_result_buf[] = {66, 156,
                                    72, 171,
                                    78, 186,
                                    104, 257,
                                    122, 302};
        NekMatrix<double> temp_expected_result(2, 5, temp_result_buf);
        BOOST_CHECK_EQUAL(temp_expected_result, tempResult);
        
        NekMatrix<double> m3(5, 1, buf3);
        
        Expression exp = m1*m2*m3;
        
        boost::tuple<unsigned int, unsigned int, unsigned int> sizes = 
                MatrixSize<Expression, Indices, 0>::GetRequiredSize(exp.GetData());
                
        BOOST_CHECK_EQUAL(sizes.get<0>(), 2u);
        BOOST_CHECK_EQUAL(sizes.get<1>(), 1u);
        BOOST_CHECK_EQUAL(sizes.get<2>(), 10u);
        
        NekMatrix<double> result(exp);
        
        double result_buf[] = {1470, 3594};
        NekMatrix<double> expected_result(2, 1, result_buf);
        
        BOOST_CHECK_EQUAL(expected_result, result);
    }
    
}













