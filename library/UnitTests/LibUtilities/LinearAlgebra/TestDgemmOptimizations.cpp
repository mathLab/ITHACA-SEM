///////////////////////////////////////////////////////////////////////////////
//
// File: TestDgemmOptimizations.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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


#include <boost/test/unit_test.hpp>

#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/auto_unit_test.hpp>

namespace Nektar
{
    
    BOOST_AUTO_TEST_CASE(TestABPlusbCScaled)
    {
        typedef NekMatrix<NekMatrix<double>, ScaledMatrixTag> ScaledMatrix;
        
        double abuf[] = {1, 2, 3, 4};
        double bbuf[] = {5, 6, 7, 8};
        double cbuf[] = {9, 10, 11, 12};
        double beta = 7.0;
        double alpha = 2.0;
        
        std::shared_ptr<NekMatrix<double> > innerA(new NekMatrix<double>(2, 2, abuf));
        std::shared_ptr<NekMatrix<double> > innerB(new NekMatrix<double>(2, 2, bbuf));
        std::shared_ptr<NekMatrix<double> > innerC(new NekMatrix<double>(2, 2, cbuf));
        
        ScaledMatrix A(2.0, innerA);
        ScaledMatrix B(3.0, innerB);
        ScaledMatrix C(4.0, innerC);
        double epsilon = 1e-11;
        
        NekMatrix<double> result1 = A*B + C;
        BOOST_CHECK_CLOSE(174.0, *result1(0,0), epsilon);
        BOOST_CHECK_CLOSE(244.0, *result1(1,0), epsilon);
        BOOST_CHECK_CLOSE(230.0, *result1(0,1), epsilon);
        BOOST_CHECK_CLOSE(324.0, *result1(1,1), epsilon);
        
        NekMatrix<double> result2 = A*B + beta*C;
        BOOST_CHECK_CLOSE(390.0, *result2(0,0), epsilon);
        BOOST_CHECK_CLOSE(484.0, *result2(1,0), epsilon);
        BOOST_CHECK_CLOSE(494.0, *result2(0,1), epsilon);
        BOOST_CHECK_CLOSE(612.0, *result2(1,1), epsilon);
        
        NekMatrix<double> result3 = alpha*A*B + C;
        BOOST_CHECK_CLOSE(312.0, *result3(0,0), epsilon);
        BOOST_CHECK_CLOSE(448.0, *result3(1,0), epsilon);
        BOOST_CHECK_CLOSE(416.0, *result3(0,1), epsilon);
        BOOST_CHECK_CLOSE(600.0, *result3(1,1), epsilon);        
        
        NekMatrix<double> result4 = alpha*A*B + beta*C;
        BOOST_CHECK_CLOSE(528.0, *result4(0,0), epsilon);
        BOOST_CHECK_CLOSE(688.0, *result4(1,0), epsilon);
        BOOST_CHECK_CLOSE(680.0, *result4(0,1), epsilon);
        BOOST_CHECK_CLOSE(888.0, *result4(1,1), epsilon);      
        
        NekMatrix<double> result5 = alpha*(A*B) + C;
        BOOST_CHECK_CLOSE(312.0, *result5(0,0), epsilon);
        BOOST_CHECK_CLOSE(448.0, *result5(1,0), epsilon);
        BOOST_CHECK_CLOSE(416.0, *result5(0,1), epsilon);
        BOOST_CHECK_CLOSE(600.0, *result5(1,1), epsilon); 
        
        NekMatrix<double> result6 = alpha*(A*B) + beta*C;
        BOOST_CHECK_CLOSE(528.0, *result6(0,0), epsilon);
        BOOST_CHECK_CLOSE(688.0, *result6(1,0), epsilon);
        BOOST_CHECK_CLOSE(680.0, *result6(0,1), epsilon);
        BOOST_CHECK_CLOSE(888.0, *result6(1,1), epsilon); 

    }
        
    BOOST_AUTO_TEST_CASE(TestABPlusbCBlock)
    {
        typedef NekMatrix<NekMatrix<double>, BlockMatrixTag> BlockMatrix;
        
        double abuf[] = {1, 2, 3, 4};
        double bbuf[] = {5, 6, 7, 8};
        double cbuf[] = {9, 10, 11, 12};
        double dbuf[] = {13, 14, 15, 16};
        double beta = 7.0;
        double alpha = 2.0;
        
        std::shared_ptr<NekMatrix<double> > innerA(new NekMatrix<double>(2, 2, abuf));
        std::shared_ptr<NekMatrix<double> > innerB(new NekMatrix<double>(2, 2, bbuf));
        std::shared_ptr<NekMatrix<double> > innerC(new NekMatrix<double>(2, 2, cbuf));
        std::shared_ptr<NekMatrix<double> > innerD(new NekMatrix<double>(2, 2, dbuf));
        
        BlockMatrix A(2, 2, 2, 2);
        BlockMatrix B(2, 2, 2, 2);
        BlockMatrix C(2, 2, 2, 2);
        
        A.SetBlock(0,0, innerA);
        A.SetBlock(0,1, innerB);
        A.SetBlock(1,0, innerC);
        A.SetBlock(1,1, innerD);
        
        B.SetBlock(0,0, innerA);
        B.SetBlock(0,1, innerB);
        B.SetBlock(1,0, innerC);
        B.SetBlock(1,1, innerD);
        
        C.SetBlock(0,0, innerA);
        C.SetBlock(0,1, innerB);
        C.SetBlock(1,0, innerC);
        C.SetBlock(1,1, innerD);
        double epsilon = 1e-11;
        
        NekMatrix<double> result1 = A*B + C;
        double expected_result_buf1[] = {123, 146, 307, 330, 
                                         157, 188, 405, 436, 
                                         191, 230, 503, 542, 
                                         225, 272, 601, 648};
        NekMatrix<double> expected_result1(4, 4, expected_result_buf1);
        BOOST_CHECK_EQUAL(expected_result1, result1);
        
        NekMatrix<double> result2 = A*B + beta*C;
        BOOST_CHECK_CLOSE(129.0, *result2(0,0), epsilon);
        BOOST_CHECK_CLOSE(158.0, *result2(1,0), epsilon);
        BOOST_CHECK_CLOSE(361.0, *result2(2,0), epsilon);
        BOOST_CHECK_CLOSE(390.0, *result2(3,0), epsilon);
        BOOST_CHECK_CLOSE(175.0, *result2(0,1), epsilon);
        BOOST_CHECK_CLOSE(212.0, *result2(1,1), epsilon);
        BOOST_CHECK_CLOSE(471.0, *result2(2,1), epsilon);
        BOOST_CHECK_CLOSE(508.0, *result2(3,1), epsilon);
        BOOST_CHECK_CLOSE(221.0, *result2(0,2), epsilon);
        BOOST_CHECK_CLOSE(266.0, *result2(1,2), epsilon);
        BOOST_CHECK_CLOSE(581.0, *result2(2,2), epsilon);
        BOOST_CHECK_CLOSE(626.0, *result2(3,2), epsilon);
        BOOST_CHECK_CLOSE(267.0, *result2(0,3), epsilon);
        BOOST_CHECK_CLOSE(320.0, *result2(1,3), epsilon);
        BOOST_CHECK_CLOSE(691.0, *result2(2,3), epsilon);
        BOOST_CHECK_CLOSE(744.0, *result2(3,3), epsilon);
        
        NekMatrix<double> result3 = alpha*A*B + C;
        BOOST_CHECK_CLOSE(245.0, *result3(0,0), epsilon);
        BOOST_CHECK_CLOSE(290.0, *result3(1,0), epsilon);
        BOOST_CHECK_CLOSE(605.0, *result3(2,0), epsilon);
        BOOST_CHECK_CLOSE(650.0, *result3(3,0), epsilon);
        BOOST_CHECK_CLOSE(311.0, *result3(0,1), epsilon);
        BOOST_CHECK_CLOSE(372.0, *result3(1,1), epsilon);
        BOOST_CHECK_CLOSE(799.0, *result3(2,1), epsilon);
        BOOST_CHECK_CLOSE(860.0, *result3(3,1), epsilon);
        BOOST_CHECK_CLOSE(377.0, *result3(0,2), epsilon);
        BOOST_CHECK_CLOSE(454.0, *result3(1,2), epsilon);
        BOOST_CHECK_CLOSE(993.0, *result3(2,2), epsilon);
        BOOST_CHECK_CLOSE(1070.0, *result3(3,2), epsilon);
        BOOST_CHECK_CLOSE(443.0, *result3(0,3), epsilon);
        BOOST_CHECK_CLOSE(536.0, *result3(1,3), epsilon);
        BOOST_CHECK_CLOSE(1187.0, *result3(2,3), epsilon);
        BOOST_CHECK_CLOSE(1280.0, *result3(3,3), epsilon);
        
        NekMatrix<double> result4 = alpha*A*B + beta*C;
        BOOST_CHECK_CLOSE(251.0, *result4(0,0), epsilon);
        BOOST_CHECK_CLOSE(302.0, *result4(1,0), epsilon);
        BOOST_CHECK_CLOSE(659.0, *result4(2,0), epsilon);
        BOOST_CHECK_CLOSE(710.0, *result4(3,0), epsilon);
        BOOST_CHECK_CLOSE(329.0, *result4(0,1), epsilon);
        BOOST_CHECK_CLOSE(396.0, *result4(1,1), epsilon);
        BOOST_CHECK_CLOSE(865.0, *result4(2,1), epsilon);
        BOOST_CHECK_CLOSE(932.0, *result4(3,1), epsilon);
        BOOST_CHECK_CLOSE(407.0, *result4(0,2), epsilon);
        BOOST_CHECK_CLOSE(490.0, *result4(1,2), epsilon);
        BOOST_CHECK_CLOSE(1071.0, *result4(2,2), epsilon);
        BOOST_CHECK_CLOSE(1154.0, *result4(3,2), epsilon);
        BOOST_CHECK_CLOSE(485.0, *result4(0,3), epsilon);
        BOOST_CHECK_CLOSE(584.0, *result4(1,3), epsilon);
        BOOST_CHECK_CLOSE(1277.0, *result4(2,3), epsilon);
        BOOST_CHECK_CLOSE(1376.0, *result4(3,3), epsilon);
        
        NekMatrix<double> result5 = alpha*(A*B) + C;
        BOOST_CHECK_CLOSE(245.0, *result5(0,0), epsilon);
        BOOST_CHECK_CLOSE(290.0, *result5(1,0), epsilon);
        BOOST_CHECK_CLOSE(605.0, *result5(2,0), epsilon);
        BOOST_CHECK_CLOSE(650.0, *result5(3,0), epsilon);
        BOOST_CHECK_CLOSE(311.0, *result5(0,1), epsilon);
        BOOST_CHECK_CLOSE(372.0, *result5(1,1), epsilon);
        BOOST_CHECK_CLOSE(799.0, *result5(2,1), epsilon);
        BOOST_CHECK_CLOSE(860.0, *result5(3,1), epsilon);
        BOOST_CHECK_CLOSE(377.0, *result5(0,2), epsilon);
        BOOST_CHECK_CLOSE(454.0, *result5(1,2), epsilon);
        BOOST_CHECK_CLOSE(993.0, *result5(2,2), epsilon);
        BOOST_CHECK_CLOSE(1070.0, *result5(3,2), epsilon);
        BOOST_CHECK_CLOSE(443.0, *result5(0,3), epsilon);
        BOOST_CHECK_CLOSE(536.0, *result5(1,3), epsilon);
        BOOST_CHECK_CLOSE(1187.0, *result5(2,3), epsilon);
        BOOST_CHECK_CLOSE(1280.0, *result5(3,3), epsilon);
        
        NekMatrix<double> result6 = alpha*(A*B) + beta*C;
        BOOST_CHECK_CLOSE(251.0, *result6(0,0), epsilon);
        BOOST_CHECK_CLOSE(302.0, *result6(1,0), epsilon);
        BOOST_CHECK_CLOSE(659.0, *result6(2,0), epsilon);
        BOOST_CHECK_CLOSE(710.0, *result6(3,0), epsilon);
        BOOST_CHECK_CLOSE(329.0, *result6(0,1), epsilon);
        BOOST_CHECK_CLOSE(396.0, *result6(1,1), epsilon);
        BOOST_CHECK_CLOSE(865.0, *result6(2,1), epsilon);
        BOOST_CHECK_CLOSE(932.0, *result6(3,1), epsilon);
        BOOST_CHECK_CLOSE(407.0, *result6(0,2), epsilon);
        BOOST_CHECK_CLOSE(490.0, *result6(1,2), epsilon);
        BOOST_CHECK_CLOSE(1071.0, *result6(2,2), epsilon);
        BOOST_CHECK_CLOSE(1154.0, *result6(3,2), epsilon);
        BOOST_CHECK_CLOSE(485.0, *result6(0,3), epsilon);
        BOOST_CHECK_CLOSE(584.0, *result6(1,3), epsilon);
        BOOST_CHECK_CLOSE(1277.0, *result6(2,3), epsilon);
        BOOST_CHECK_CLOSE(1376.0, *result6(3,3), epsilon);

    }
}
