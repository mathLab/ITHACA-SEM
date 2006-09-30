///////////////////////////////////////////////////////////////////////////////
//
// File: testLinearSystem.cpp
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
// Description: Test code for NekVector
//
///////////////////////////////////////////////////////////////////////////////

#include <UnitTests/testLinearSystem.h>
#include <LibUtilities/LinearAlgebra/NekLinSys.hpp>

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

namespace Nektar
{
    namespace LinearSystemUnitTests
    {
        void testDiagonalSystem()
        {
            unsigned int matrix_buf[] = { 10, 5, 2 };

            unsigned int result_buf[] = { 20, 50, 10 };

            boost::shared_ptr<NekMatrix<unsigned int, eDiagonal> >  A(new NekMatrix<unsigned int, eDiagonal>(3,3,matrix_buf));
            boost::shared_ptr<NekVector<unsigned int, 3> > b(new NekVector<unsigned int, 3>(result_buf));

            LinearSystem<NekMatrix<unsigned int, eDiagonal>, NekVector<unsigned int, 3> > linsys(A, b);

            NekVector<unsigned int, 3> result = linsys.Solve();
            
            unsigned int expected_result_buf[] = { 2, 10, 5 };
            NekVector<unsigned int, 3> expectedResult(expected_result_buf);

            BOOST_CHECK_EQUAL(result, expectedResult);
        }
    }
}
