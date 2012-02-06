///////////////////////////////////////////////////////////////////////////////
//
// File: testNekMatrixOperations.h
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
// Description: Tests NekMatrix functionality.
//
///////////////////////////////////////////////////////////////////////////////

#include "LibUtilitiesUnitTestsPrecompiledHeader.h"

#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <LibUtilities/BasicUtils/BoostUtil.hpp>
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/progress.hpp>
#include <iostream>
#include <boost/bind.hpp>
#include <functional>

namespace Nektar
{
    class DoAddition
    {
        public:
            template<typename LhsType, typename RhsType>
            typename expt::AddOp::ResultType<LhsType, RhsType>::type 
            operator()(const LhsType& lhs, const RhsType& rhs) const
            {
                return lhs + rhs;
            }
    };

    class DoSubtraction
    {
        public:
            template<typename LhsType, typename RhsType>
            typename expt::SubtractOp::ResultType<LhsType, RhsType>::type 
            operator()(const LhsType& lhs, const RhsType& rhs) const
            {
                return lhs - rhs;
            }
    };

    template<typename LhsScaledInnerMatrixType,
             typename LhsBlockInnerMatrixType,
             typename RhsScaledInnerMatrixType,
             typename RhsBlockInnerMatrixType,
             typename OpType>
    void RunAllTestCombinations(const NekMatrix<NekDouble, StandardMatrixTag>& l1,
                                const NekMatrix<NekMatrix<NekDouble, LhsScaledInnerMatrixType>, ScaledMatrixTag>& l2,
                                const NekMatrix<NekMatrix<NekDouble, LhsBlockInnerMatrixType>, BlockMatrixTag>& l3,
                                const NekMatrix<NekDouble, StandardMatrixTag>& r1,
                                const NekMatrix<NekMatrix<NekDouble, RhsScaledInnerMatrixType>, ScaledMatrixTag>& r2,
                                const NekMatrix<NekMatrix<NekDouble, RhsBlockInnerMatrixType>, BlockMatrixTag>& r3,
                                const NekMatrix<NekDouble, StandardMatrixTag>& result,
                                const OpType& f)
    {
        BOOST_CHECK_EQUAL(f(l1, r1), result);
        BOOST_CHECK_EQUAL(f(l1, r2), result);
        BOOST_CHECK_EQUAL(f(l1, r3), result);
        BOOST_CHECK_EQUAL(f(l2, r1), result);
        BOOST_CHECK_EQUAL(f(l2, r2), result);
        BOOST_CHECK_EQUAL(f(l2, r3), result);
        BOOST_CHECK_EQUAL(f(l3, r1), result);
        BOOST_CHECK_EQUAL(f(l3, r2), result);
        BOOST_CHECK_EQUAL(f(l3, r3), result);
    }

    
    void GenerateFullMatrices(double values[], double scale,
        boost::shared_ptr<NekMatrix<NekDouble, StandardMatrixTag> >& m1,
        boost::shared_ptr<NekMatrix<NekMatrix<NekDouble>, ScaledMatrixTag> >& m2,
        boost::shared_ptr<NekMatrix<NekMatrix<NekDouble>, BlockMatrixTag> >& m3)
    {
        m1 = MakePtr(new NekMatrix<NekDouble, StandardMatrixTag>(4, 4, values));
        
        double inner_values[16];
        std::transform(values, values+16, inner_values, boost::bind(std::divides<NekDouble>(), _1, scale));

        boost::shared_ptr<NekMatrix<NekDouble> > inner(
            new NekMatrix<NekDouble>(4, 4, inner_values)); 
        m2 = MakePtr(new NekMatrix<NekMatrix<NekDouble>, ScaledMatrixTag>(scale, inner));
        
        double block_1_values[] = {values[0], values[1], 
                            values[4], values[5]};
        double block_2_values[] = {values[2], values[3],
                            values[6], values[7]};
        double block_3_values[] = {values[8], values[9], 
                            values[12], values[13]};
        double block_4_values[] = {values[10], values[11],
                            values[14], values[15]};
        boost::shared_ptr<NekMatrix<NekDouble> > block1(new NekMatrix<NekDouble>(2, 2, block_1_values));
        boost::shared_ptr<NekMatrix<NekDouble> > block2(new NekMatrix<NekDouble>(2, 2, block_2_values));
        boost::shared_ptr<NekMatrix<NekDouble> > block3(new NekMatrix<NekDouble>(2, 2, block_3_values));
        boost::shared_ptr<NekMatrix<NekDouble> > block4(new NekMatrix<NekDouble>(2, 2, block_4_values));
        
        m3 = MakePtr(new NekMatrix<NekMatrix<NekDouble>, BlockMatrixTag>(2, 2, 2, 2));
        m3->SetBlock(0,0, block1);
        m3->SetBlock(1,0, block2);
        m3->SetBlock(0,1, block3);
        m3->SetBlock(1,1, block4);
    }

    void GenerateUpperTriangularMatrices(NekDouble values[], NekDouble scale,
        boost::shared_ptr<NekMatrix<NekDouble, StandardMatrixTag> >& m1,
        boost::shared_ptr<NekMatrix<NekMatrix<NekDouble>, ScaledMatrixTag> >& m2,
        boost::shared_ptr<NekMatrix<NekMatrix<NekDouble>, BlockMatrixTag> >& m3)
    {
        m1 = MakePtr(new NekMatrix<NekDouble, StandardMatrixTag>(4, 4, values, eUPPER_TRIANGULAR));
        
        double inner_values[10];
        std::transform(values, values+10, inner_values, boost::bind(std::divides<NekDouble>(), _1, scale));

        boost::shared_ptr<NekMatrix<NekDouble> > inner(
            new NekMatrix<NekDouble>(4, 4, inner_values, eUPPER_TRIANGULAR)); 
        m2 = MakePtr(new NekMatrix<NekMatrix<NekDouble>, ScaledMatrixTag>(scale, inner));
        
        double block_1_values[] = {values[0], 0.0,
                                   values[1], values[2]};
        double block_2_values[] = {values[3], values[4],
                                   values[6], values[7]};
        double block_4_values[] = {values[5], 0.0,
                                   values[8], values[9]};
        boost::shared_ptr<NekMatrix<NekDouble> > block1(new NekMatrix<NekDouble>(2, 2, block_1_values));
        boost::shared_ptr<NekMatrix<NekDouble> > block2(new NekMatrix<NekDouble>(2, 2, block_2_values));
        boost::shared_ptr<NekMatrix<NekDouble> > block4(new NekMatrix<NekDouble>(2, 2, block_4_values));
        
        m3 = MakePtr(new NekMatrix<NekMatrix<NekDouble>, BlockMatrixTag>(2, 2, 2, 2));
        m3->SetBlock(0,0, block1);
        m3->SetBlock(0,1, block2);
        m3->SetBlock(1,1, block4);
    }


    namespace MatrixOperationTests
    {
        
//         void GenerateDiagonalMatrices(double values[], double scale,
//             boost::shared_ptr<NekMatrix<NekDouble, FullMatrixTag, StandardMatrixTag> >& m1,
//             boost::shared_ptr<NekMatrix<NekMatrix<NekDouble>, FullMatrixTag, ScaledMatrixTag> >& m2,
//             boost::shared_ptr<NekMatrix<NekMatrix<NekDouble>, FullMatrixTag, BlockMatrixTag> >& m3)
//         {
//             m1 = MakePtr(new NekMatrix<NekDouble, DiagonalMatrixTag, StandardMatrixTag>(4, 4, values));
//             
//             double inner_values[4];
//             std::transform(values, values+4, inner_values, boost::bind(std::divides<NekDouble>(), _1, scale));
// 
//             boost::shared_ptr<NekMatrix<NekDouble, DiagonalMatrixTag> > inner(
//                 new NekMatrix<NekDouble, DiagonalMatrixTag>(4, 4, inner_values)); 
//             m2 = MakePtr(new NekMatrix<NekMatrix<NekDouble, DiagonalMatrixTag>, DiagonalMatrixTag, ScaledMatrixTag>(scale, inner));
//             
//             double block_1_values[] = {values[0], values[1], 
//                                 values[4], values[5]};
//             double block_2_values[] = {values[2], values[3],
//                                 values[6], values[7]};
//             double block_3_values[] = {values[8], values[9], 
//                                 values[12], values[13]};
//             double block_4_values[] = {values[10], values[11],
//                                 values[14], values[15]};
//             boost::shared_ptr<NekMatrix<NekDouble> > block1(new NekMatrix<NekDouble>(2, 2, block_1_values));
//             boost::shared_ptr<NekMatrix<NekDouble> > block2(new NekMatrix<NekDouble>(2, 2, block_2_values));
//             boost::shared_ptr<NekMatrix<NekDouble> > block3(new NekMatrix<NekDouble>(2, 2, block_3_values));
//             boost::shared_ptr<NekMatrix<NekDouble> > block4(new NekMatrix<NekDouble>(2, 2, block_4_values));
//             
//             m3 = MakePtr(new NekMatrix<NekMatrix<NekDouble>, FullMatrixTag, BlockMatrixTag>(2, 2, 2, 2));
//             m3->SetBlock(0,0, block1);
//             m3->SetBlock(0,1, block2);
//             m3->SetBlock(1,0, block3);
//             m3->SetBlock(1,1, block4);
//        }

        BOOST_AUTO_TEST_CASE(TestLhsFullRhsFull)
        {
            //double lhs_values[] = {2, 4, 6, 8,
            //                       10, 12, 14, 16,
            //                       18, 20, 22, 24,
            //                       26, 28, 30, 32};
            double lhs_values[] = {2, 10, 18, 26,
                                 4, 12, 20, 28,
                                 6, 14, 22, 30,
                                 8, 16, 24, 32};
                                   
            boost::shared_ptr<NekMatrix<NekDouble, StandardMatrixTag> > lhs1;
            boost::shared_ptr<NekMatrix<NekMatrix<NekDouble>, ScaledMatrixTag> > lhs2;
            boost::shared_ptr<NekMatrix<NekMatrix<NekDouble>, BlockMatrixTag> > lhs3;
            
            GenerateFullMatrices(lhs_values, 2.0, lhs1, lhs2, lhs3);
            //double rhs_values[] = {4, 8, 12, 16,
            //                       20, 24, 28, 32,
            //                       36, 40, 44, 48,
            //                       52, 56, 60, 64};
            double rhs_values[] = {4, 20, 36, 52,
                                   8, 24, 40, 56,
                                   12, 28, 44, 60,
                                   16, 32, 48, 64};
            boost::shared_ptr<NekMatrix<NekDouble, StandardMatrixTag> > rhs1;
            boost::shared_ptr<NekMatrix<NekMatrix<NekDouble>, ScaledMatrixTag> > rhs2;
            boost::shared_ptr<NekMatrix<NekMatrix<NekDouble>, BlockMatrixTag> > rhs3; 
            GenerateFullMatrices(rhs_values, 2.0, rhs1, rhs2, rhs3);
            
            double result_values[16];
            std::transform(lhs_values, lhs_values+16, rhs_values, result_values, std::plus<NekDouble>());
            NekMatrix<NekDouble> result(4, 4, result_values);
            
            RunAllTestCombinations(*lhs1, *lhs2, *lhs3, *rhs1, *rhs2, *rhs3, result, DoAddition());
        }
        
        BOOST_AUTO_TEST_CASE(TestLhsFullRhsDiagonal)
        {
//             double lhs_values[] = {2, 4, 6, 8,
//                                    10, 12, 14, 16,
//                                    18, 20, 22, 24,
//                                    26, 28, 30, 32};
//                                    
//             boost::shared_ptr<NekMatrix<NekDouble, FullMatrixTag, StandardMatrixTag> > lhs1;
//             boost::shared_ptr<NekMatrix<NekMatrix<NekDouble>, FullMatrixTag, ScaledMatrixTag> > lhs2;
//             boost::shared_ptr<NekMatrix<NekMatrix<NekDouble>, FullMatrixTag, BlockMatrixTag> > lhs3;
//             GenerateMatrices(lhs_values, 2.0, lhs1, lhs2, lhs3);
            
//             double rhs_values[] = {10, 20, 30, 40};
//             boost::shared_ptr<NekMatrix<NekDouble, DiagonalMatrixTag, StandardMatrixTag> > rhs1;
//             boost::shared_ptr<NekMatrix<NekMatrix<NekDouble>, DiagonalMatrixTag, ScaledMatrixTag> > rhs2;
//             boost::shared_ptr<NekMatrix<NekMatrix<NekDouble>, DiagonalMatrixTag, BlockMatrixTag> > rhs3;
//             GenerateDiagonalMatrices(rhs_values, 2.0, rhs1, rhs2, rhs3);
            
        }
        
        BOOST_AUTO_TEST_CASE(TestComboExpression)
        {
            {
                //double buf[] = {1.0, 2.0, 3.0, 4.0};
                double buf[] = {1.0, 3.0,
                                2.0, 4.0};
                SharedNekMatrixPtr inner1(new NekMatrix<NekDouble>(2,2,buf));
                SharedNekMatrixPtr inner2(new NekMatrix<NekDouble>(2,2,buf));
                
                DNekScalMat m1(2.0, inner1);
                DNekScalMat m2(3.0, inner2);
                
                (*inner1)*2.0;
                2.0*(*inner1);
                //DNekScalMat m3 = m1*2.0;
                NekMatrix<NekDouble> result = m1 + 1.0/6.0*m2;
                
                BOOST_CHECK_EQUAL(result(0,0), 2.5);
                BOOST_CHECK_EQUAL(result(0,1), 5.0);
                BOOST_CHECK_EQUAL(result(1,0), 7.5);
                BOOST_CHECK_EQUAL(result(1,1), 10.0);
            }

            {
                double buf[] = {1.0, 2.0, 3.0, 4.0};
                SharedNekMatrixPtr inner1(new NekMatrix<NekDouble>(2,2,buf));
                SharedNekMatrixPtr inner2(new NekMatrix<NekDouble>(2,2,buf));
                
                boost::shared_ptr<DNekScalMat> m1(new DNekScalMat(2.0, inner1));
                boost::shared_ptr<DNekScalMat> m2(new DNekScalMat(3.0, inner2));
                
                (*m1)*2.0;

                //BOOST_CHECK_EQUAL(result->GetValue(0,0), 2.0);
                //BOOST_CHECK_EQUAL(result->GetValue(0,1), 4.0);
                //BOOST_CHECK_EQUAL(result->GetValue(1,0), 6.0);
                //BOOST_CHECK_EQUAL(result->GetValue(1,1), 8.0);
            }

        }

        BOOST_AUTO_TEST_CASE(TestLhsUpperTriangularRhsUpperTriangular)
        {
            //double lhs_values[] = {2, 4,  6,  8,
            //                          12, 14, 16,
            //                              22, 24,
            //                                  32};
            double lhs_values[] = {2,
                                   4, 12,
                                   6, 14, 22,
                                   8, 16, 24, 32};

            boost::shared_ptr<NekMatrix<NekDouble, StandardMatrixTag> > lhs1;
            boost::shared_ptr<NekMatrix<NekMatrix<NekDouble>, ScaledMatrixTag> > lhs2;
            boost::shared_ptr<NekMatrix<NekMatrix<NekDouble>, BlockMatrixTag> > lhs3;
            GenerateUpperTriangularMatrices(lhs_values, 2.0, lhs1, lhs2, lhs3);

            //double rhs_values[] = {4,   8, 12, 16,
            //                           24, 28, 32,
            //                               44, 48,
            //                                   64};
            double rhs_values[] = {4,
                                   8, 24,
                                   12, 28, 44,
                                   16, 32, 48, 64};
            boost::shared_ptr<NekMatrix<NekDouble, StandardMatrixTag> > rhs1;
            boost::shared_ptr<NekMatrix<NekMatrix<NekDouble>, ScaledMatrixTag> > rhs2;
            boost::shared_ptr<NekMatrix<NekMatrix<NekDouble>, BlockMatrixTag> > rhs3; 
            GenerateUpperTriangularMatrices(rhs_values, 2.0, rhs1, rhs2, rhs3);
            
            double result_values[10];
            std::transform(lhs_values, lhs_values+10, rhs_values, result_values, std::plus<NekDouble>());
            NekMatrix<NekDouble, StandardMatrixTag> result(4, 4, result_values, eUPPER_TRIANGULAR);
            
            //RunAllTestCombinations(*lhs1, *lhs2, *lhs3, *rhs1, *rhs2, *rhs3, result, DoAddition());
        }

        BOOST_AUTO_TEST_CASE(TestLhsLowerTriangularRhsLowerTriangular)
        {
        }
    }

    namespace MatrixSubtractionTests
    {
        BOOST_AUTO_TEST_CASE(TestLhsFullRhsFull)
        {
            //double lhs_values[] = {2, 4, 6, 8,
            //                        10, 12, 14, 16,
            //                        18, 20, 22, 24,
            //                        26, 28, 30, 32};
                                    
            double lhs_values[] = {2, 10, 18, 26,
                                   4, 12, 20, 28,
                                   6, 14, 22, 30,
                                   8, 16, 24, 32};

            boost::shared_ptr<NekMatrix<NekDouble, StandardMatrixTag> > lhs1;
            boost::shared_ptr<NekMatrix<NekMatrix<NekDouble>, ScaledMatrixTag> > lhs2;
            boost::shared_ptr<NekMatrix<NekMatrix<NekDouble>, BlockMatrixTag> > lhs3;
            
            GenerateFullMatrices(lhs_values, 2.0, lhs1, lhs2, lhs3);
            //double rhs_values[] = {4, 8, 12, 16,
            //                        20, 24, 28, 32,
            //                        36, 40, 44, 48,
            //                        52, 56, 60, 64};
            double rhs_values[] = {4, 20, 36, 52,
                                   8, 24, 40, 56,
                                   12, 28, 44, 60,
                                   16, 32, 48, 64};
            boost::shared_ptr<NekMatrix<NekDouble, StandardMatrixTag> > rhs1;
            boost::shared_ptr<NekMatrix<NekMatrix<NekDouble>, ScaledMatrixTag> > rhs2;
            boost::shared_ptr<NekMatrix<NekMatrix<NekDouble>, BlockMatrixTag> > rhs3; 
            GenerateFullMatrices(rhs_values, 2.0, rhs1, rhs2, rhs3);
            
            double result_values[16];
            std::transform(lhs_values, lhs_values+16, rhs_values, result_values, std::minus<NekDouble>());
            NekMatrix<NekDouble> result(4, 4, result_values);
            RunAllTestCombinations(*lhs1, *lhs2, *lhs3, *rhs1, *rhs2, *rhs3, result, DoSubtraction());
        }
    }
}


