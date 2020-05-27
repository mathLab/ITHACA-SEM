///////////////////////////////////////////////////////////////////////////////
//
// File: testNekMatrix.cpp
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
// Description: Tests NekMatrix functionality.
//
///////////////////////////////////////////////////////////////////////////////

#include <UnitTests/util.h>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/progress.hpp>
#include <iostream>

namespace Nektar
{
    namespace MatrixUnitTests
    {
        BOOST_AUTO_TEST_CASE(TestFullMatrixInversion)
        {
            {
                double buf[] = {1.0, 3.0, 2.0, 4.0};
                NekMatrix<double> m(2,2,buf);
                m.Invert();

                BOOST_CHECK_CLOSE(*m(0,0), -2.0, .00001);
                BOOST_CHECK_CLOSE(*m(0,1), 1.0, .00001);
                BOOST_CHECK_CLOSE(*m(1,0), 3.0/2.0, .00001);
                BOOST_CHECK_CLOSE(*m(1,1), -1.0/2.0, .00001);
            }
            
            {
                double buf[] = {1.7, 5.6, 9.1, -3.4, -2.0,
                                4.5, 7.8, 8.2, -2.5, 5.0,
                                9.0, 3.6, 7.3, .98, 1.0,
                                -12.6, 2.9, 6.4, .02, -4.0,
                                -1, 1.7, 5.6, 2.98, 5.0};
                NekMatrix<double> m(5,5,buf);
                m.Invert();

                BOOST_CHECK_CLOSE(*m(0,0), 0.0005010995978, .000001);
                BOOST_CHECK_CLOSE(*m(0,1), -0.4704403712, .000001);
                BOOST_CHECK_CLOSE(*m(0,2), 0.2719063614 , .000001);
                BOOST_CHECK_CLOSE(*m(0,3), -0.3941557805 , .000001);
                BOOST_CHECK_CLOSE(*m(0,4), 0.09043166650, .000001);

                BOOST_CHECK_CLOSE(*m(1,0), -0.01087166322 , .000001);
                BOOST_CHECK_CLOSE(*m(1,1), 0.3242048735 , .000001);
                BOOST_CHECK_CLOSE(*m(1,2), -0.1605116333, .000001);
                BOOST_CHECK_CLOSE(*m(1,3), 0.09133673974 , .000001);
                BOOST_CHECK_CLOSE(*m(1,4), 0.01293234281, .000001);

                BOOST_CHECK_CLOSE(*m(2,0), 0.06465598241, .000001);
                BOOST_CHECK_CLOSE(*m(2,1), 0.2593895742 , .000001);
                BOOST_CHECK_CLOSE(*m(2,2), -0.09057233795, .000001);
                BOOST_CHECK_CLOSE(*m(2,3), 0.3106562459, .000001);
                BOOST_CHECK_CLOSE(*m(2,4), -0.1589713628, .000001);

                BOOST_CHECK_CLOSE(*m(3,0), -0.03464982687, .000001);
                BOOST_CHECK_CLOSE(*m(3,1), 0.2655177914, .000001);
                BOOST_CHECK_CLOSE(*m(3,2), -0.1016866584, .000001);
                BOOST_CHECK_CLOSE(*m(3,3), 0.2125363445, .000001);
                BOOST_CHECK_CLOSE(*m(3,4), -0.1099886183, .000001);

                BOOST_CHECK_CLOSE(*m(4,0), -0.02957895492, .000001);
                BOOST_CHECK_CLOSE(*m(4,1), -0.3518447037, .000001);
                BOOST_CHECK_CLOSE(*m(4,2), 0.2060393188, .000001);
                BOOST_CHECK_CLOSE(*m(4,3), -0.1411012255, .000001);
                BOOST_CHECK_CLOSE(*m(4,4), 0.1670437017, .000001);
            }
            
            {
            }
        }
        
        BOOST_AUTO_TEST_CASE(TestDiagonalMatrixInversion)
        {
            double buf[] = {1.0, 2.0, 3.0, 4.0};
            NekMatrix<double> m(4, 4, buf, eDIAGONAL);
            m.Invert();
            
            BOOST_CHECK_EQUAL(m(0,0), 1.0/1.0);
            BOOST_CHECK_EQUAL(m(1,1), 1.0/2.0);
            BOOST_CHECK_EQUAL(m(2,2), 1.0/3.0);
            BOOST_CHECK_EQUAL(m(3,3), 1.0/4.0);
        }
        
       
        BOOST_AUTO_TEST_CASE(TestNekMatrixConstruction)
        {
            {
                std::shared_ptr<Nektar::Matrix<double> > a(new Nektar::NekMatrix<double>(3,4));
                std::shared_ptr<Nektar::NekMatrix<double> > b(new Nektar::NekMatrix<double>(5,6));
                
                BOOST_CHECK_EQUAL(a->GetRows(), 3u);
                BOOST_CHECK_EQUAL(a->GetColumns(), 4u);
                BOOST_CHECK_EQUAL(a->GetSize()[0], 3u);
                BOOST_CHECK_EQUAL(a->GetSize()[1], 4u);
                
                BOOST_CHECK_EQUAL(b->GetRows(), 5u);
                BOOST_CHECK_EQUAL(b->GetColumns(), 6u);
                BOOST_CHECK_EQUAL(b->GetSize()[0], 5u);
                BOOST_CHECK_EQUAL(b->GetSize()[1], 6u);
                
                BOOST_CHECK_EQUAL(a->GetStorageType(), eFULL);
                BOOST_CHECK_EQUAL(b->GetStorageType(), eFULL);
                
                BOOST_CHECK_EQUAL(a->GetStorageSize(), 3*4);
                BOOST_CHECK_EQUAL(b->GetStorageSize(), 5*6);
            }
            
            {
                std::shared_ptr<Nektar::Matrix<double> > a(new Nektar::NekMatrix<double>(3,3, eDIAGONAL));
                std::shared_ptr<Nektar::NekMatrix<double> > b(new Nektar::NekMatrix<double>(5,5,eDIAGONAL));
                
                BOOST_CHECK_EQUAL(a->GetRows(), 3);
                BOOST_CHECK_EQUAL(a->GetColumns(), 3);
                BOOST_CHECK_EQUAL(a->GetSize()[0], 3);
                BOOST_CHECK_EQUAL(a->GetSize()[1], 3);
                
                BOOST_CHECK_EQUAL(b->GetRows(), 5);
                BOOST_CHECK_EQUAL(b->GetColumns(), 5);
                BOOST_CHECK_EQUAL(b->GetSize()[0], 5);
                BOOST_CHECK_EQUAL(b->GetSize()[1], 5);
                
                BOOST_CHECK_EQUAL(a->GetStorageType(), eDIAGONAL);
                BOOST_CHECK_EQUAL(b->GetStorageType(), eDIAGONAL);
                
                BOOST_CHECK_EQUAL(a->GetStorageSize(), 3);
                BOOST_CHECK_EQUAL(b->GetStorageSize(), 5);
            }
        }
        
        BOOST_AUTO_TEST_CASE(TestFullNekMatrixGetValue)
        {
            UnitTests::RedirectCerrIfNeeded();
            double data[] = {1.0, 89.0, 0.0, 45.12,
                             2.0, -12.3, 892.2532, 76.12,
                             3.0, -56.7, 211.22, 45.23};

            Nektar::NekMatrix<double> m(4, 3, data);
            
            BOOST_CHECK_EQUAL( m(0,0), 1.0 );
            BOOST_CHECK_EQUAL( m(0,1), 2.0 );
            BOOST_CHECK_EQUAL( m(0,2), 3.0 );
            
            BOOST_CHECK_EQUAL( m(1,0), 89.0 );
            BOOST_CHECK_EQUAL( m(1,1), -12.3 );
            BOOST_CHECK_EQUAL( m(1,2), -56.7 );
            
            BOOST_CHECK_EQUAL( m(2,0), 0.0 );
            BOOST_CHECK_EQUAL( m(2,1), 892.2532 );
            BOOST_CHECK_EQUAL( m(2,2), 211.22 );
            
            BOOST_CHECK_EQUAL( m(3,0), 45.12 );
            BOOST_CHECK_EQUAL( m(3,1), 76.12 );
            BOOST_CHECK_EQUAL( m(3,2), 45.23 );
            
            #ifdef NEKTAR_FULLDEBUG
            BOOST_CHECK_THROW( m(0,3), ErrorUtil::NekError);
            BOOST_CHECK_THROW( m(4,1), ErrorUtil::NekError);
            BOOST_CHECK_THROW( m(4,3), ErrorUtil::NekError);
            #endif //NEKTAR_FULLDEBUG
            
            std::shared_ptr<Nektar::Matrix<double> > m1(new Nektar::NekMatrix<double>(4, 3, data));
            std::shared_ptr<Nektar::NekMatrix<double> > m2(new Nektar::NekMatrix<double>(4, 3, data));
            
            BOOST_CHECK_EQUAL( (*m1)(0,0), 1.0 );
            BOOST_CHECK_EQUAL( (*m1)(0,1), 2.0 );
            BOOST_CHECK_EQUAL( (*m1)(0,2), 3.0 );
            
            BOOST_CHECK_EQUAL( (*m1)(1,0), 89.0 );
            BOOST_CHECK_EQUAL( (*m1)(1,1), -12.3 );
            BOOST_CHECK_EQUAL( (*m1)(1,2), -56.7 );
            
            BOOST_CHECK_EQUAL( (*m1)(2,0), 0.0 );
            BOOST_CHECK_EQUAL( (*m1)(2,1), 892.2532 );
            BOOST_CHECK_EQUAL( (*m1)(2,2), 211.22 );
            
            BOOST_CHECK_EQUAL( (*m1)(3,0), 45.12 );
            BOOST_CHECK_EQUAL( (*m1)(3,1), 76.12 );
            BOOST_CHECK_EQUAL( (*m1)(3,2), 45.23 );
            
            #ifdef NEKTAR_FULLDEBUG
            BOOST_CHECK_THROW( (*m1)(0,3), ErrorUtil::NekError);
            BOOST_CHECK_THROW( (*m1)(4,1), ErrorUtil::NekError);
            BOOST_CHECK_THROW( (*m1)(4,3), ErrorUtil::NekError);
            #endif //NEKTAR_FULLDEBUG
            
            BOOST_CHECK_EQUAL( (*m2)(0,0), 1.0 );
            BOOST_CHECK_EQUAL( (*m2)(0,1), 2.0 );
            BOOST_CHECK_EQUAL( (*m2)(0,2), 3.0 );
            
            BOOST_CHECK_EQUAL( (*m2)(1,0), 89.0 );
            BOOST_CHECK_EQUAL( (*m2)(1,1), -12.3 );
            BOOST_CHECK_EQUAL( (*m2)(1,2), -56.7 );
            
            BOOST_CHECK_EQUAL( (*m2)(2,0), 0.0 );
            BOOST_CHECK_EQUAL( (*m2)(2,1), 892.2532 );
            BOOST_CHECK_EQUAL( (*m2)(2,2), 211.22 );
            
            BOOST_CHECK_EQUAL( (*m2)(3,0), 45.12 );
            BOOST_CHECK_EQUAL( (*m2)(3,1), 76.12 );
            BOOST_CHECK_EQUAL( (*m2)(3,2), 45.23 );
            
            #ifdef NEKTAR_FULLDEBUG
            BOOST_CHECK_THROW( (*m2)(0,3), ErrorUtil::NekError);
            BOOST_CHECK_THROW( (*m2)(4,1), ErrorUtil::NekError);
            BOOST_CHECK_THROW( (*m2)(4,3), ErrorUtil::NekError);
            #endif //NEKTAR_FULLDEBUG
            
        }
        
        BOOST_AUTO_TEST_CASE(TestDiagonalMatrixGetValue)
        {
            UnitTests::RedirectCerrIfNeeded();
            double data[] = {8.9, 3.4, 5.7};
            Nektar::NekMatrix<double> m1(3, 3, data, eDIAGONAL);
            std::shared_ptr<Nektar::Matrix<double> > m2(new Nektar::NekMatrix<double>(3, 3, data, eDIAGONAL));
            std::shared_ptr<Nektar::NekMatrix<double > > m3(new Nektar::NekMatrix<double>(3, 3, data, eDIAGONAL));
            
            BOOST_CHECK_EQUAL(m1(0,0), 8.9);
            BOOST_CHECK_EQUAL((*m2)(0,0), 8.9);
            BOOST_CHECK_EQUAL((*m3)(0,0), 8.9);
            
            BOOST_CHECK_EQUAL(m1(0,1), 0.0);
            BOOST_CHECK_EQUAL((*m2)(0,1), 0.0);
            BOOST_CHECK_EQUAL((*m3)(0,1), 0.0);
            
            BOOST_CHECK_EQUAL(m1(0,2), 0.0);
            BOOST_CHECK_EQUAL((*m2)(0,2), 0.0);
            BOOST_CHECK_EQUAL((*m3)(0,2), 0.0);
            
            BOOST_CHECK_EQUAL(m1(1,0), 0.0);
            BOOST_CHECK_EQUAL((*m2)(1,0), 0.0);
            BOOST_CHECK_EQUAL((*m3)(1,0), 0.0);
            
            BOOST_CHECK_EQUAL(m1(1,1), 3.4);
            BOOST_CHECK_EQUAL((*m2)(1,1), 3.4);
            BOOST_CHECK_EQUAL((*m3)(1,1), 3.4);
            
            BOOST_CHECK_EQUAL(m1(1,2), 0.0);
            BOOST_CHECK_EQUAL((*m2)(1,2), 0.0);
            BOOST_CHECK_EQUAL((*m3)(1,2), 0.0);
            
            BOOST_CHECK_EQUAL(m1(2,0), 0.0);
            BOOST_CHECK_EQUAL((*m2)(2,0), 0.0);
            BOOST_CHECK_EQUAL((*m3)(2,0), 0.0);
            
            BOOST_CHECK_EQUAL(m1(2,1), 0.0);
            BOOST_CHECK_EQUAL((*m2)(2,1), 0.0);
            BOOST_CHECK_EQUAL((*m3)(2,1), 0.0);
            
            BOOST_CHECK_EQUAL(m1(2,2), 5.7);
            BOOST_CHECK_EQUAL((*m2)(2,2), 5.7);
            BOOST_CHECK_EQUAL((*m3)(2,2), 5.7);
            
            #ifdef NEKTAR_FULLDEBUG
            BOOST_CHECK_THROW(m1(0,4), ErrorUtil::NekError);
            BOOST_CHECK_THROW((*m2)(0,4), ErrorUtil::NekError);
            BOOST_CHECK_THROW((*m3)(0,4), ErrorUtil::NekError);
            
            BOOST_CHECK_THROW(m1(4,0), ErrorUtil::NekError);
            BOOST_CHECK_THROW((*m2)(4,0), ErrorUtil::NekError);
            BOOST_CHECK_THROW((*m3)(4,0), ErrorUtil::NekError);
            
            BOOST_CHECK_THROW(m1(4,4), ErrorUtil::NekError);
            BOOST_CHECK_THROW((*m2)(4,4), ErrorUtil::NekError);
            BOOST_CHECK_THROW((*m3)(4,4), ErrorUtil::NekError);
            #endif //NEKTAR_FULLDEBUG
            
            
        }
        
        BOOST_AUTO_TEST_CASE(TestFullNekMatrixSetValue)
        {
            double data[] = {1.0, 3.0, 5.0,
                             2.0, 4.0, 6.0};
            Nektar::NekMatrix<double> m1(3, 2, data);
            std::shared_ptr<Nektar::Matrix<double> > m2(new Nektar::NekMatrix<double>(3, 2, data));
            std::shared_ptr<Nektar::NekMatrix<double> > m3(new Nektar::NekMatrix<double>(3, 2, data));
            
            m1.SetValue(0,0,-1.0);
            m2->SetValue(1,1,-2.0);
            m3->SetValue(2,1, -3.0);
            
            BOOST_CHECK_EQUAL(m1(0,0), -1.0);
            BOOST_CHECK_EQUAL(m1(0,1), 2.0);
            BOOST_CHECK_EQUAL(m1(1,0), 3.0);
            BOOST_CHECK_EQUAL(m1(1,1), 4.0);
            BOOST_CHECK_EQUAL(m1(2,0), 5.0);
            BOOST_CHECK_EQUAL(m1(2,1), 6.0);
            
            BOOST_CHECK_EQUAL((*m2)(0,0), 1.0);
            BOOST_CHECK_EQUAL((*m2)(0,1), 2.0);
            BOOST_CHECK_EQUAL((*m2)(1,0), 3.0);
            BOOST_CHECK_EQUAL((*m2)(1,1), -2.0);
            BOOST_CHECK_EQUAL((*m2)(2,0), 5.0);
            BOOST_CHECK_EQUAL((*m2)(2,1), 6.0);
            
            BOOST_CHECK_EQUAL((*m3)(0,0), 1.0);
            BOOST_CHECK_EQUAL((*m3)(0,1), 2.0);
            BOOST_CHECK_EQUAL((*m3)(1,0), 3.0);
            BOOST_CHECK_EQUAL((*m3)(1,1), 4.0);
            BOOST_CHECK_EQUAL((*m3)(2,0), 5.0);
            BOOST_CHECK_EQUAL((*m3)(2,1), -3.0);
        }
        
        BOOST_AUTO_TEST_CASE(TestDiagonalNekMatrixSetValue)
        {
            UnitTests::RedirectCerrIfNeeded();
            double data[] = {8.9, 3.4, 5.7};
            Nektar::NekMatrix<double> m1(3, 3, data, eDIAGONAL);
            std::shared_ptr<Nektar::Matrix<double> > m2(new Nektar::NekMatrix<double>(3, 3, data, eDIAGONAL));
            std::shared_ptr<Nektar::NekMatrix<double > > m3(new Nektar::NekMatrix<double>(3, 3, data, eDIAGONAL));
            
            m1.SetValue(0,0,1.0);
            m2->SetValue(1,1,2.0);
            m3->SetValue(2,2,3.0);
            
            BOOST_CHECK_EQUAL(m1.GetStorageType(), eDIAGONAL);
            BOOST_CHECK_EQUAL(m2->GetStorageType(), eDIAGONAL);
            BOOST_CHECK_EQUAL(m3->GetStorageType(), eDIAGONAL);
            
            BOOST_CHECK_EQUAL(m1(0,0), 1.0);
            BOOST_CHECK_EQUAL(m1(1,1), 3.4);
            BOOST_CHECK_EQUAL(m1(2,2), 5.7);
            
            BOOST_CHECK_EQUAL((*m2)(0,0), 8.9);
            BOOST_CHECK_EQUAL((*m2)(1,1), 2.0);
            BOOST_CHECK_EQUAL((*m2)(2,2), 5.7);
            
            BOOST_CHECK_EQUAL((*m3)(0,0), 8.9);
            BOOST_CHECK_EQUAL((*m3)(1,1), 3.4);
            BOOST_CHECK_EQUAL((*m3)(2,2), 3.0);
            
            #ifdef NEKTAR_FULLDEBUG
            BOOST_CHECK_THROW(m1.SetValue(0,3,2.0), ErrorUtil::NekError);
            BOOST_CHECK_THROW(m1.SetValue(3,0,2.0), ErrorUtil::NekError);
            BOOST_CHECK_THROW(m1.SetValue(3,3,2.0), ErrorUtil::NekError);
            
            BOOST_CHECK_THROW(m2->SetValue(0,3,2.0), ErrorUtil::NekError);
            BOOST_CHECK_THROW(m2->SetValue(3,0,2.0), ErrorUtil::NekError);
            BOOST_CHECK_THROW(m2->SetValue(3,3,2.0), ErrorUtil::NekError);
            
            BOOST_CHECK_THROW(m3->SetValue(0,3,2.0), ErrorUtil::NekError);
            BOOST_CHECK_THROW(m3->SetValue(3,0,2.0), ErrorUtil::NekError);
            BOOST_CHECK_THROW(m3->SetValue(3,3,2.0), ErrorUtil::NekError);
            #endif //NEKTAR_FULLDEBUG
        }
        
        
        
        BOOST_AUTO_TEST_CASE(TestFullFullMatrixAddition)
        {
            double lhs_data[] = {1.0, 4.0,
                                 2.0, 5.0,
                                 3.0, 6.0};
            double rhs_data[] = {10.0, 13.0,
                                 11.0, 14.0,
                                 12.0, 15.0};
            
            Nektar::NekMatrix<double> lhs(2,3,lhs_data);
            Nektar::NekMatrix<double> rhs(2,3,rhs_data);
            Nektar::NekMatrix<double> result = lhs + rhs;
            
            BOOST_CHECK_EQUAL(result.GetRows(), 2);
            BOOST_CHECK_EQUAL(result.GetColumns(), 3);
            BOOST_CHECK_EQUAL(result.GetStorageType(), eFULL);
            BOOST_CHECK_EQUAL(result(0,0), 11.0);
            BOOST_CHECK_EQUAL(result(0,1), 13.0);
            BOOST_CHECK_EQUAL(result(0,2), 15.0);
            BOOST_CHECK_EQUAL(result(1,0), 17.0);
            BOOST_CHECK_EQUAL(result(1,1), 19.0);
            BOOST_CHECK_EQUAL(result(1,2), 21.0);
        }
    }
    
    namespace UnitTests
    {
        
        using namespace Nektar;
                
        BOOST_AUTO_TEST_CASE(testNekMatrixConstruction)
        {
            // Basic, dense matrix construction.
            {
                double buf[] = {1.0, 4.0, 7.0, 10.0,
                                2.0, 5.0, 8.0, 11.0,
                                3.0, 6.0, 9.0, 12.0};
                NekMatrix<double> dynamic_matrix(4, 3, buf);

                BOOST_CHECK(dynamic_matrix.GetRows() == 4);
                BOOST_CHECK(dynamic_matrix.GetColumns() == 3);

                for(unsigned int i = 0; i < 4; ++i)
                {
                    for(unsigned int j = 0; j < 3; ++j)
                    {
                        BOOST_CHECK(dynamic_matrix(i,j) == buf[4*j + i]);
                    }
                }
            }
        }



        BOOST_AUTO_TEST_CASE(testNekMatrixBasicMath)
        {
             // Addition tests.
             {
                  double buf[] = {1.0, 4.0, 7.0,
                                 2.0, 5.0, 8.0,
                                 3.0, 6.0, 9.0};
                 NekMatrix<double> m1(3, 3, buf);
                 NekMatrix<double> m2(3, 3, buf);
                 NekMatrix<double> m3 = m1 + m2;
 
                 for(unsigned int i = 0; i < 3; ++i)
                 {
                     for(unsigned int j = 0; j < 3; ++j)
                     {
                         BOOST_CHECK(m3(i,j) == buf[3*j+i] + buf[3*j+i]);
                     }
                 }
 
                 NekMatrix<double> m4(3, 3, buf);
                 NekMatrix<double> m5(3, 3, buf);
                 NekMatrix<double> m6 = m4+m5;
 
                 for(unsigned int i = 0; i < 3; ++i)
                 {
                     for(unsigned int j = 0; j < 3; ++j)
                     {
                         BOOST_CHECK(m6(i,j) == buf[3*j+i] + buf[3*j+i]);
                     }
                 }
             }
             // Multiply

             {
                 double buf1[] = {1, 4, 7,
                                  2, 5, 8,
                                  3, 6, 9};
                 double buf2[] = {10, 15, 19,
                                  11, 16, 20,
                                  12, 17, 21};
                 NekMatrix<double> lhs(3, 3, buf1);
                 NekMatrix<double> rhs(3, 3, buf2);
 
                 NekMatrix<double> result = lhs*rhs;
 
                 BOOST_CHECK(result.GetRows() == 3);
                 BOOST_CHECK(result.GetColumns() == 3);
 
                 double epsilon = 1e-12;
                 BOOST_CHECK_CLOSE(*result(0,0), 97.0, epsilon);
                 BOOST_CHECK_CLOSE(*result(0,1), 103.0, epsilon);
                 BOOST_CHECK_CLOSE(*result(0,2), 109.0, epsilon);
 
                 BOOST_CHECK_CLOSE(*result(1,0), 229.0, epsilon);
                 BOOST_CHECK_CLOSE(*result(1,1), 244.0, epsilon);
                 BOOST_CHECK_CLOSE(*result(1,2), 259.0, epsilon);
 
                 BOOST_CHECK_CLOSE(*result(2,0), 361.0, epsilon);
                 BOOST_CHECK_CLOSE(*result(2,1), 385.0, epsilon);
                 BOOST_CHECK_CLOSE(*result(2,2), 409.0, epsilon);
             }
         
             {
                 double buf1[] = {1, 4, 7,
                                  2, 5, 8,
                                  3, 6, 9};
                 double buf2[] = {10, 15, 19,
                                  11, 16, 20,
                                  12, 17, 21,
                                  14, 18, 22};
                 NekMatrix<double> lhs(3, 3, buf1);
                 NekMatrix<double> rhs(3, 4, buf2);
 
                 NekMatrix<double> result = lhs*rhs;
 
                 BOOST_CHECK(result.GetRows() == 3);
                 BOOST_CHECK(result.GetColumns() == 4);
 
                 double epsilon = 1e-12;
                 BOOST_CHECK_CLOSE(*result(0,0), 97.0, epsilon);
                 BOOST_CHECK_CLOSE(*result(0,1), 103.0, epsilon);
                 BOOST_CHECK_CLOSE(*result(0,2), 109.0, epsilon);
                 BOOST_CHECK_CLOSE(*result(0,3), 116.0, epsilon);
 
                 BOOST_CHECK_CLOSE(*result(1,0), 229.0, epsilon);
                 BOOST_CHECK_CLOSE(*result(1,1), 244.0, epsilon);
                 BOOST_CHECK_CLOSE(*result(1,2), 259.0, epsilon);
                 BOOST_CHECK_CLOSE(*result(1,3), 278.0, epsilon);
 
                 BOOST_CHECK_CLOSE(*result(2,0), 361.0, epsilon);
                 BOOST_CHECK_CLOSE(*result(2,1), 385.0, epsilon);
                 BOOST_CHECK_CLOSE(*result(2,2), 409.0, epsilon);
                 BOOST_CHECK_CLOSE(*result(2,3), 440.0, epsilon);
             }
            
 

// 
//             // Transpose
// 
//             // Determinant.
// 
//             // Invert/Check for singularity.
// 
//             // Eigenvalues/vectors
// 
//             // Condition number wrt various norms.
// 
//             // Various norm computations.
// 
//             // LU Decomposition?  More appropriate in LinAlg?
        }




    }
}
