///////////////////////////////////////////////////////////////////////////////
//
// File: TestMatrixBlasOptimizations.cpp
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
    BOOST_AUTO_TEST_CASE(DgemmTest)
    {
        double a_buf[] = {1, 2, 3, 4};
        double b_buf[] = {4, 5, 6, 7};
        double c_buf[] = {8, 9, 10, 11};
        
        NekMatrix<double> a(2, 2, a_buf);
        NekMatrix<double> b(2, 2, b_buf);
        NekMatrix<double> c(2, 2, c_buf);
        
        NekMatrix<double> result = a*b + c;
        
        
        const Expression<BinaryExpressionPolicy
                    <
                        ConstantExpressionPolicy<NekMatrix<double> >, 
                        MultiplyOp, 
                        ConstantExpressionPolicy<NekMatrix<double> >
                    > >& t1 = a*b;
        
        const Expression
            <
                BinaryExpressionPolicy
                <
                    BinaryExpressionPolicy
                    <
                        ConstantExpressionPolicy<NekMatrix<double> >,
                        MultiplyOp,
                        ConstantExpressionPolicy<NekMatrix<double> >
                    >,
                    AddOp,
                    ConstantExpressionPolicy<NekMatrix<double> >
                >
            >& t2 = a*b + c;
        double expected_result_buf[] = {27, 37, 37, 51};
        NekMatrix<double> expected_result(2, 2, expected_result_buf);
        BOOST_CHECK_EQUAL(expected_result, result);
    }
    
}

