///////////////////////////////////////////////////////////////////////////////
//
// File: TestBandedMatrix.cpp
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

#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <boost/test/auto_unit_test.hpp>

namespace Nektar
{
    namespace BandedMatrixPolicyUnitTests
    {
        typedef MatrixStoragePolicy<NekDouble, BandedMatrixTag> Policy;
        typedef Policy::PolicySpecificDataHolderType DataHolderType;

        BOOST_AUTO_TEST_CASE(TestUninitializedDataConstructionWithDefaultDataHolder)
        {
            Array<OneD, NekDouble> result = Policy::Initialize(10, 10, DataHolderType());
            BOOST_CHECK_EQUAL(result.num_elements(), 100);

            BOOST_CHECK_THROW(Policy::Initialize(10, 9, DataHolderType()), ErrorUtil::NekError);
            BOOST_CHECK_THROW(Policy::Initialize(9, 10, DataHolderType()), ErrorUtil::NekError);
        }

        BOOST_AUTO_TEST_CASE(TestUninitializedDataConstructionWithUserDefinedNumberOfDiagonals)
        {
            DataHolderType d1(1,0);
            Array<OneD, NekDouble> result1 = Policy::Initialize(10, 10, d1);
            BOOST_CHECK_EQUAL(result1.num_elements(), 30);

            DataHolderType d2(0,1);
            Array<OneD, NekDouble> result2 = Policy::Initialize(10, 10, d2);
            BOOST_CHECK_EQUAL(result2.num_elements(), 20);

            DataHolderType d3(1,1);
            Array<OneD, NekDouble> result3 = Policy::Initialize(10, 10, d3);
            BOOST_CHECK_EQUAL(result3.num_elements(), 40);

            DataHolderType d4(3,1);
            Array<OneD, NekDouble> result4 = Policy::Initialize(10, 10, d4);
            BOOST_CHECK_EQUAL(result4.num_elements(), 80);

            DataHolderType d5(5,1);
            Array<OneD, NekDouble> result5 = Policy::Initialize(10, 10, d5);
            BOOST_CHECK_EQUAL(result5.num_elements(), 100);

        }

    }

    namespace BandedMatrixMultiplicationUnitTests
    {
        BOOST_AUTO_TEST_CASE(TestDirectBlasCall)
        {
           
            NekDouble x[] = {1, 2, 3, 4};
            NekDouble y[] = {0,0,0,0};

            NekDouble a[] = {0.0, 2.0, 5.0, 8.0,
                               1.0, 4.0, 7.0, 10.0,
                               3.0, 6.0, 9.0, 0.0 };
            Blas::Dgbmv('T', 4, 4, 1, 1, 1.0, a, 4, x, 1, 0.0, y, 1);

            // This works.
            //NekDouble a[] = {0, 1, 3, 2, 4, 6, 5, 7, 9, 8, 10, 0};
            //Blas::Dgbmv('N', 4, 4, 1, 1, 1.0, a, 3, x, 1, 0.0, y, 1);

            NekDouble expected_result_buf[] = { 5, 26, 65, 67 };
            NekVector<NekDouble, 4> expected_result(expected_result_buf);
            NekVector<NekDouble, 4> result(y);
            BOOST_CHECK_EQUAL(expected_result, result);
        }
        
        BOOST_AUTO_TEST_CASE(TestSquareDirectBlasCall)
        {
            NekDouble x[] = {1, 2, 3, 4};
            NekDouble y[] = {0,0,0,0};

            NekDouble a[] = {0, 0, 0, 4,
                             0, 0, 3, 8,
                             0, 2, 7, 12,
                             1, 6, 11, 16,
                             5, 10, 15, 0,
                             9, 14, 0, 0,
                             13, 0, 0, 0};
            int result = Blas::Dgbmv('T', 4, 4, 3, 3, 1.0, a, 3, x, 1, 0.0, y, 1);

            // This works.
            //NekDouble a[] = {0, 0, 0, 1, 5, 9, 13,
            //                 0, 0, 2, 6, 10, 14, 0,
            //                 0, 3, 7, 11, 15, 0, 0,
            //                 4, 8, 12, 16, 0, 0, 0};
            //Blas::Dgbmv('N', 4, 4, 3, 3, 1.0, a, 7, x, 1, 0.0, y, 1);


            NekDouble expected_result_buf[] = { 30, 70, 110, 150 };
            NekVector<NekDouble, 4> expected_result(expected_result_buf);
            NekVector<NekDouble, 4> result(y);
            BOOST_CHECK_EQUAL(expected_result, result);
        }

        BOOST_AUTO_TEST_CASE(TestStandardMatrixVectorMultiply)
        {
            typedef NekMatrix<NekDouble, BandedMatrixTag> MatrixType;

            {
                // This is the matrix
                // [ 1 2 0 0 ]
                // [ 3 4 5 0 ]
                // [ 0 6 7 8 ]
                // [ 0 0 9 10]
                NekDouble buf[] = { 0.0, 0.0, 0.0, 0.0,
                                    0.0, 2.0, 5.0, 8.0,
                                    1.0, 4.0, 7.0, 10.0,
                                    3.0, 6.0, 9.0, 0.0 };
    
                MatrixType m(4, 4, buf, MatrixType::PolicySpecificDataHolderType(1, 1));

                NekDouble vector_buf[] = {1.0, 2.0, 3.0, 4.0};
                NekVector<NekDouble, 4> v(vector_buf);

                NekVector<NekDouble, 4> result = m*v;

                NekDouble expected_result_buf[] = { 5, 26, 65, 67 };
                NekVector<NekDouble, 4> expected_result(expected_result_buf);
                BOOST_CHECK_EQUAL(expected_result, result);
            }
        }
    }
}


