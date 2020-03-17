///////////////////////////////////////////////////////////////////////////////
//
// File: BlockMatrixUnitTests.cpp
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

#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/progress.hpp>
#include <iostream>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>

namespace Nektar
{
    namespace BlockMatrixUnitTests
    {
        BOOST_AUTO_TEST_CASE(TestEqualSizedBlockConstruction)
        {
            {
                typedef NekMatrix<NekMatrix<NekDouble>, BlockMatrixTag> BlockMatrixType;
                BlockMatrixType m1(3, 2, 2, 2);
                std::shared_ptr<ConstMatrix<NekDouble> > m2(new BlockMatrixType(2, 3, 2, 2));
                
                BOOST_CHECK_EQUAL(m1.GetRows(), 6u);
                BOOST_CHECK_EQUAL(m1.GetColumns(), 4u);
                BOOST_CHECK_EQUAL(m2->GetRows(), 4u);
                BOOST_CHECK_EQUAL(m2->GetColumns(), 6u);               
            }
            
            // Bad Input
            {
            }
            
        }
                
        BOOST_AUTO_TEST_CASE(TestElementAccessBlockMatrix)
        {
            {
                typedef NekMatrix<NekMatrix<NekDouble>, BlockMatrixTag> BlockMatrixType;
                BlockMatrixType m1(3, 2, 2, 2);
                std::shared_ptr<ConstMatrix<NekDouble> > m2(new BlockMatrixType(3, 2, 2, 2));
                
                double vals1[] = {1.0, 3.0, 2.0, 4.0};
                double vals2[] = {5.0, 7.0, 6.0, 8.0};
                double vals3[] = {9.0, 11.0, 10.0, 12.0};
                double vals4[] = {13.0, 15.0, 14.0, 16.0};
                double vals5[] = {17.0, 19.0, 18.0, 20.0};
                double vals6[] = {21.0, 23.0, 22.0, 24.0};
                
                std::shared_ptr<NekMatrix<NekDouble> > sub1(new NekMatrix<NekDouble>(2, 2, vals1));
                std::shared_ptr<NekMatrix<NekDouble> > sub2(new NekMatrix<NekDouble>(2, 2, vals2));
                std::shared_ptr<NekMatrix<NekDouble> > sub3(new NekMatrix<NekDouble>(2, 2, vals3));
                std::shared_ptr<NekMatrix<NekDouble> > sub4(new NekMatrix<NekDouble>(2, 2, vals4));
                std::shared_ptr<NekMatrix<NekDouble> > sub5(new NekMatrix<NekDouble>(2, 2, vals5));
                std::shared_ptr<NekMatrix<NekDouble> > sub6(new NekMatrix<NekDouble>(2, 2, vals6));
                
                m1.SetBlock(0,0, sub1);
                m1.SetBlock(0,1, sub2);
                m1.SetBlock(1,0, sub3);
                m1.SetBlock(1,1, sub4);
                m1.SetBlock(2,0, sub5);
                m1.SetBlock(2,1, sub6);
                
                std::shared_ptr<BlockMatrixType> m2_cast = std::dynamic_pointer_cast<BlockMatrixType>(m2);
                
                m2_cast->SetBlock(0,0, sub1);
                m2_cast->SetBlock(0,1, sub2);
                m2_cast->SetBlock(1,0, sub3);
                m2_cast->SetBlock(1,1, sub4);
                m2_cast->SetBlock(2,0, sub5);
                m2_cast->SetBlock(2,1, sub6);
                               
                BOOST_CHECK_EQUAL(m1(0,0), 1.0);               
                BOOST_CHECK_EQUAL(m1(0,1), 2.0);               
                BOOST_CHECK_EQUAL(m1(0,2), 5.0);               
                BOOST_CHECK_EQUAL(m1(0,3), 6.0);               
                BOOST_CHECK_EQUAL(m1(1,0), 3.0);               
                BOOST_CHECK_EQUAL(m1(1,1), 4.0);               
                BOOST_CHECK_EQUAL(m1(1,2), 7.0);               
                BOOST_CHECK_EQUAL(m1(1,3), 8.0);               
                BOOST_CHECK_EQUAL(m1(2,0), 9.0);               
                BOOST_CHECK_EQUAL(m1(2,1), 10.0);               
                BOOST_CHECK_EQUAL(m1(2,2), 13.0);               
                BOOST_CHECK_EQUAL(m1(2,3), 14.0);               
                BOOST_CHECK_EQUAL(m1(3,0), 11.0);               
                BOOST_CHECK_EQUAL(m1(3,1), 12.0);               
                BOOST_CHECK_EQUAL(m1(3,2), 15.0);               
                BOOST_CHECK_EQUAL(m1(3,3), 16.0);               
                BOOST_CHECK_EQUAL(m1(4,0), 17.0);               
                BOOST_CHECK_EQUAL(m1(4,1), 18.0);               
                BOOST_CHECK_EQUAL(m1(4,2), 21.0);               
                BOOST_CHECK_EQUAL(m1(4,3), 22.0);               
                BOOST_CHECK_EQUAL(m1(5,0), 19.0);               
                BOOST_CHECK_EQUAL(m1(5,1), 20.0);               
                BOOST_CHECK_EQUAL(m1(5,2), 23.0);               
                BOOST_CHECK_EQUAL(m1(5,3), 24.0);               
                
                BOOST_CHECK_EQUAL((*m2)(0,0), 1.0);               
                BOOST_CHECK_EQUAL((*m2)(0,1), 2.0);               
                BOOST_CHECK_EQUAL((*m2)(0,2), 5.0);               
                BOOST_CHECK_EQUAL((*m2)(0,3), 6.0);               
                BOOST_CHECK_EQUAL((*m2)(1,0), 3.0);               
                BOOST_CHECK_EQUAL((*m2)(1,1), 4.0);               
                BOOST_CHECK_EQUAL((*m2)(1,2), 7.0);               
                BOOST_CHECK_EQUAL((*m2)(1,3), 8.0);               
                BOOST_CHECK_EQUAL((*m2)(2,0), 9.0);               
                BOOST_CHECK_EQUAL((*m2)(2,1), 10.0);               
                BOOST_CHECK_EQUAL((*m2)(2,2), 13.0);               
                BOOST_CHECK_EQUAL((*m2)(2,3), 14.0);               
                BOOST_CHECK_EQUAL((*m2)(3,0), 11.0);               
                BOOST_CHECK_EQUAL((*m2)(3,1), 12.0);               
                BOOST_CHECK_EQUAL((*m2)(3,2), 15.0);               
                BOOST_CHECK_EQUAL((*m2)(3,3), 16.0);               
                BOOST_CHECK_EQUAL((*m2)(4,0), 17.0);               
                BOOST_CHECK_EQUAL((*m2)(4,1), 18.0);               
                BOOST_CHECK_EQUAL((*m2)(4,2), 21.0);               
                BOOST_CHECK_EQUAL((*m2)(4,3), 22.0);               
                BOOST_CHECK_EQUAL((*m2)(5,0), 19.0);               
                BOOST_CHECK_EQUAL((*m2)(5,1), 20.0);               
                BOOST_CHECK_EQUAL((*m2)(5,2), 23.0);               
                BOOST_CHECK_EQUAL((*m2)(5,3), 24.0);     
                         
            }
            
            
            {
                typedef NekMatrix<NekMatrix<NekDouble>, BlockMatrixTag> BlockMatrixType;
                unsigned int rowSizes[] = {2, 2, 2};
                unsigned int columnSizes[] = {2, 2};
                BlockMatrixType m1(3, 2, rowSizes, columnSizes);
                std::shared_ptr<ConstMatrix<NekDouble> > m2(new BlockMatrixType(3, 2, rowSizes, columnSizes));
                
                double vals1[] = {1.0, 3.0, 2.0, 4.0};
                double vals2[] = {5.0, 7.0, 6.0, 8.0};
                double vals3[] = {9.0, 11.0, 10.0, 12.0};
                double vals4[] = {13.0, 15.0, 14.0, 16.0};
                double vals5[] = {17.0, 19.0, 18.0, 20.0};
                double vals6[] = {21.0, 23.0, 22.0, 24.0};
                
                std::shared_ptr<NekMatrix<NekDouble> > sub1(new NekMatrix<NekDouble>(2, 2, vals1));
                std::shared_ptr<NekMatrix<NekDouble> > sub2(new NekMatrix<NekDouble>(2, 2, vals2));
                std::shared_ptr<NekMatrix<NekDouble> > sub3(new NekMatrix<NekDouble>(2, 2, vals3));
                std::shared_ptr<NekMatrix<NekDouble> > sub4(new NekMatrix<NekDouble>(2, 2, vals4));
                std::shared_ptr<NekMatrix<NekDouble> > sub5(new NekMatrix<NekDouble>(2, 2, vals5));
                std::shared_ptr<NekMatrix<NekDouble> > sub6(new NekMatrix<NekDouble>(2, 2, vals6));
                
                m1.SetBlock(0,0, sub1);
                m1.SetBlock(0,1, sub2);
                m1.SetBlock(1,0, sub3);
                m1.SetBlock(1,1, sub4);
                m1.SetBlock(2,0, sub5);
                m1.SetBlock(2,1, sub6);
                
                std::shared_ptr<BlockMatrixType> m2_cast = std::dynamic_pointer_cast<BlockMatrixType>(m2);
                
                m2_cast->SetBlock(0,0, sub1);
                m2_cast->SetBlock(0,1, sub2);
                m2_cast->SetBlock(1,0, sub3);
                m2_cast->SetBlock(1,1, sub4);
                m2_cast->SetBlock(2,0, sub5);
                m2_cast->SetBlock(2,1, sub6);
                               
                BOOST_CHECK_EQUAL(m1(0,0), 1.0);               
                BOOST_CHECK_EQUAL(m1(0,1), 2.0);               
                BOOST_CHECK_EQUAL(m1(0,2), 5.0);               
                BOOST_CHECK_EQUAL(m1(0,3), 6.0);               
                BOOST_CHECK_EQUAL(m1(1,0), 3.0);               
                BOOST_CHECK_EQUAL(m1(1,1), 4.0);               
                BOOST_CHECK_EQUAL(m1(1,2), 7.0);               
                BOOST_CHECK_EQUAL(m1(1,3), 8.0);               
                BOOST_CHECK_EQUAL(m1(2,0), 9.0);               
                BOOST_CHECK_EQUAL(m1(2,1), 10.0);               
                BOOST_CHECK_EQUAL(m1(2,2), 13.0);               
                BOOST_CHECK_EQUAL(m1(2,3), 14.0);               
                BOOST_CHECK_EQUAL(m1(3,0), 11.0);               
                BOOST_CHECK_EQUAL(m1(3,1), 12.0);               
                BOOST_CHECK_EQUAL(m1(3,2), 15.0);               
                BOOST_CHECK_EQUAL(m1(3,3), 16.0);               
                BOOST_CHECK_EQUAL(m1(4,0), 17.0);               
                BOOST_CHECK_EQUAL(m1(4,1), 18.0);               
                BOOST_CHECK_EQUAL(m1(4,2), 21.0);               
                BOOST_CHECK_EQUAL(m1(4,3), 22.0);               
                BOOST_CHECK_EQUAL(m1(5,0), 19.0);               
                BOOST_CHECK_EQUAL(m1(5,1), 20.0);               
                BOOST_CHECK_EQUAL(m1(5,2), 23.0);               
                BOOST_CHECK_EQUAL(m1(5,3), 24.0);               
                
                BOOST_CHECK_EQUAL((*m2)(0,0), 1.0);               
                BOOST_CHECK_EQUAL((*m2)(0,1), 2.0);               
                BOOST_CHECK_EQUAL((*m2)(0,2), 5.0);               
                BOOST_CHECK_EQUAL((*m2)(0,3), 6.0);               
                BOOST_CHECK_EQUAL((*m2)(1,0), 3.0);               
                BOOST_CHECK_EQUAL((*m2)(1,1), 4.0);               
                BOOST_CHECK_EQUAL((*m2)(1,2), 7.0);               
                BOOST_CHECK_EQUAL((*m2)(1,3), 8.0);               
                BOOST_CHECK_EQUAL((*m2)(2,0), 9.0);               
                BOOST_CHECK_EQUAL((*m2)(2,1), 10.0);               
                BOOST_CHECK_EQUAL((*m2)(2,2), 13.0);               
                BOOST_CHECK_EQUAL((*m2)(2,3), 14.0);               
                BOOST_CHECK_EQUAL((*m2)(3,0), 11.0);               
                BOOST_CHECK_EQUAL((*m2)(3,1), 12.0);               
                BOOST_CHECK_EQUAL((*m2)(3,2), 15.0);               
                BOOST_CHECK_EQUAL((*m2)(3,3), 16.0);               
                BOOST_CHECK_EQUAL((*m2)(4,0), 17.0);               
                BOOST_CHECK_EQUAL((*m2)(4,1), 18.0);               
                BOOST_CHECK_EQUAL((*m2)(4,2), 21.0);               
                BOOST_CHECK_EQUAL((*m2)(4,3), 22.0);               
                BOOST_CHECK_EQUAL((*m2)(5,0), 19.0);               
                BOOST_CHECK_EQUAL((*m2)(5,1), 20.0);               
                BOOST_CHECK_EQUAL((*m2)(5,2), 23.0);               
                BOOST_CHECK_EQUAL((*m2)(5,3), 24.0);     
                          
                
            }
            
            {
                typedef NekMatrix<NekMatrix<NekDouble>, BlockMatrixTag> BlockMatrixType;
                unsigned int rowSizesBuf[] = {2, 2, 2};
                unsigned int columnSizesBuf[] = {2, 2};
                Array<OneD, unsigned int> rowSizes(3, rowSizesBuf);
                Array<OneD, unsigned int> columnSizes(2, columnSizesBuf);
                
                BlockMatrixType m1(3, 2, rowSizes, columnSizes);
                std::shared_ptr<ConstMatrix<NekDouble> > m2(new BlockMatrixType(3, 2, rowSizes, columnSizes));
                
                double vals1[] = {1.0, 3.0, 2.0, 4.0};
                double vals2[] = {5.0, 7.0, 6.0, 8.0};
                double vals3[] = {9.0, 11.0, 10.0, 12.0};
                double vals4[] = {13.0, 15.0, 14.0, 16.0};
                double vals5[] = {17.0, 19.0, 18.0, 20.0};
                double vals6[] = {21.0, 23.0, 22.0, 24.0};
                
                std::shared_ptr<NekMatrix<NekDouble> > sub1(new NekMatrix<NekDouble>(2, 2, vals1));
                std::shared_ptr<NekMatrix<NekDouble> > sub2(new NekMatrix<NekDouble>(2, 2, vals2));
                std::shared_ptr<NekMatrix<NekDouble> > sub3(new NekMatrix<NekDouble>(2, 2, vals3));
                std::shared_ptr<NekMatrix<NekDouble> > sub4(new NekMatrix<NekDouble>(2, 2, vals4));
                std::shared_ptr<NekMatrix<NekDouble> > sub5(new NekMatrix<NekDouble>(2, 2, vals5));
                std::shared_ptr<NekMatrix<NekDouble> > sub6(new NekMatrix<NekDouble>(2, 2, vals6));
                
                m1.SetBlock(0,0, sub1);
                m1.SetBlock(0,1, sub2);
                m1.SetBlock(1,0, sub3);
                m1.SetBlock(1,1, sub4);
                m1.SetBlock(2,0, sub5);
                m1.SetBlock(2,1, sub6);
                
                std::shared_ptr<BlockMatrixType> m2_cast = std::dynamic_pointer_cast<BlockMatrixType>(m2);
                
                m2_cast->SetBlock(0,0, sub1);
                m2_cast->SetBlock(0,1, sub2);
                m2_cast->SetBlock(1,0, sub3);
                m2_cast->SetBlock(1,1, sub4);
                m2_cast->SetBlock(2,0, sub5);
                m2_cast->SetBlock(2,1, sub6);
                               
                BOOST_CHECK_EQUAL(m1(0,0), 1.0);               
                BOOST_CHECK_EQUAL(m1(0,1), 2.0);               
                BOOST_CHECK_EQUAL(m1(0,2), 5.0);               
                BOOST_CHECK_EQUAL(m1(0,3), 6.0);               
                BOOST_CHECK_EQUAL(m1(1,0), 3.0);               
                BOOST_CHECK_EQUAL(m1(1,1), 4.0);               
                BOOST_CHECK_EQUAL(m1(1,2), 7.0);               
                BOOST_CHECK_EQUAL(m1(1,3), 8.0);               
                BOOST_CHECK_EQUAL(m1(2,0), 9.0);               
                BOOST_CHECK_EQUAL(m1(2,1), 10.0);               
                BOOST_CHECK_EQUAL(m1(2,2), 13.0);               
                BOOST_CHECK_EQUAL(m1(2,3), 14.0);               
                BOOST_CHECK_EQUAL(m1(3,0), 11.0);               
                BOOST_CHECK_EQUAL(m1(3,1), 12.0);               
                BOOST_CHECK_EQUAL(m1(3,2), 15.0);               
                BOOST_CHECK_EQUAL(m1(3,3), 16.0);               
                BOOST_CHECK_EQUAL(m1(4,0), 17.0);               
                BOOST_CHECK_EQUAL(m1(4,1), 18.0);               
                BOOST_CHECK_EQUAL(m1(4,2), 21.0);               
                BOOST_CHECK_EQUAL(m1(4,3), 22.0);               
                BOOST_CHECK_EQUAL(m1(5,0), 19.0);               
                BOOST_CHECK_EQUAL(m1(5,1), 20.0);               
                BOOST_CHECK_EQUAL(m1(5,2), 23.0);               
                BOOST_CHECK_EQUAL(m1(5,3), 24.0);               
                
                BOOST_CHECK_EQUAL((*m2)(0,0), 1.0);               
                BOOST_CHECK_EQUAL((*m2)(0,1), 2.0);               
                BOOST_CHECK_EQUAL((*m2)(0,2), 5.0);               
                BOOST_CHECK_EQUAL((*m2)(0,3), 6.0);               
                BOOST_CHECK_EQUAL((*m2)(1,0), 3.0);               
                BOOST_CHECK_EQUAL((*m2)(1,1), 4.0);               
                BOOST_CHECK_EQUAL((*m2)(1,2), 7.0);               
                BOOST_CHECK_EQUAL((*m2)(1,3), 8.0);               
                BOOST_CHECK_EQUAL((*m2)(2,0), 9.0);               
                BOOST_CHECK_EQUAL((*m2)(2,1), 10.0);               
                BOOST_CHECK_EQUAL((*m2)(2,2), 13.0);               
                BOOST_CHECK_EQUAL((*m2)(2,3), 14.0);               
                BOOST_CHECK_EQUAL((*m2)(3,0), 11.0);               
                BOOST_CHECK_EQUAL((*m2)(3,1), 12.0);               
                BOOST_CHECK_EQUAL((*m2)(3,2), 15.0);               
                BOOST_CHECK_EQUAL((*m2)(3,3), 16.0);               
                BOOST_CHECK_EQUAL((*m2)(4,0), 17.0);               
                BOOST_CHECK_EQUAL((*m2)(4,1), 18.0);               
                BOOST_CHECK_EQUAL((*m2)(4,2), 21.0);               
                BOOST_CHECK_EQUAL((*m2)(4,3), 22.0);               
                BOOST_CHECK_EQUAL((*m2)(5,0), 19.0);               
                BOOST_CHECK_EQUAL((*m2)(5,1), 20.0);               
                BOOST_CHECK_EQUAL((*m2)(5,2), 23.0);               
                BOOST_CHECK_EQUAL((*m2)(5,3), 24.0);     
                          
                
            }
        }

        BOOST_AUTO_TEST_CASE(TestEmptyRowElementAccess)
        {
            typedef NekMatrix<NekMatrix<NekDouble>, BlockMatrixTag> BlockMatrixType;
            unsigned int rowSizes[] = {2, 0, 2, 0, 0, 2};
            unsigned int columnSizes[] = {2, 2};
            BlockMatrixType m1(6, 2, rowSizes, columnSizes);
            std::shared_ptr<ConstMatrix<NekDouble> > m2(new BlockMatrixType(6, 2, rowSizes, columnSizes));
            
            double vals1[] = {1.0, 3.0, 2.0, 4.0};
            double vals2[] = {5.0, 7.0, 6.0, 8.0};
            double vals3[] = {9.0, 11.0, 10.0, 12.0};
            double vals4[] = {13.0, 15.0, 14.0, 16.0};
            double vals5[] = {17.0, 19.0, 18.0, 20.0};
            double vals6[] = {21.0, 23.0, 22.0, 24.0};
            
            std::shared_ptr<NekMatrix<NekDouble> > sub1(new NekMatrix<NekDouble>(2, 2, vals1));
            std::shared_ptr<NekMatrix<NekDouble> > sub2(new NekMatrix<NekDouble>(2, 2, vals2));
            std::shared_ptr<NekMatrix<NekDouble> > sub3(new NekMatrix<NekDouble>(2, 2, vals3));
            std::shared_ptr<NekMatrix<NekDouble> > sub4(new NekMatrix<NekDouble>(2, 2, vals4));
            std::shared_ptr<NekMatrix<NekDouble> > sub5(new NekMatrix<NekDouble>(2, 2, vals5));
            std::shared_ptr<NekMatrix<NekDouble> > sub6(new NekMatrix<NekDouble>(2, 2, vals6));
            
            m1.SetBlock(0,0, sub1);
            m1.SetBlock(0,1, sub2);
            m1.SetBlock(2,0, sub3);
            m1.SetBlock(2,1, sub4);
            m1.SetBlock(5,0, sub5);
            m1.SetBlock(5,1, sub6);
            
            std::shared_ptr<BlockMatrixType> m2_cast = std::dynamic_pointer_cast<BlockMatrixType>(m2);
            
            m2_cast->SetBlock(0,0, sub1);
            m2_cast->SetBlock(0,1, sub2);
            m2_cast->SetBlock(2,0, sub3);
            m2_cast->SetBlock(2,1, sub4);
            m2_cast->SetBlock(5,0, sub5);
            m2_cast->SetBlock(5,1, sub6);
                           
            BOOST_CHECK_EQUAL(m1(0,0), 1.0);               
            BOOST_CHECK_EQUAL(m1(0,1), 2.0);               
            BOOST_CHECK_EQUAL(m1(0,2), 5.0);               
            BOOST_CHECK_EQUAL(m1(0,3), 6.0);               
            BOOST_CHECK_EQUAL(m1(1,0), 3.0);               
            BOOST_CHECK_EQUAL(m1(1,1), 4.0);               
            BOOST_CHECK_EQUAL(m1(1,2), 7.0);               
            BOOST_CHECK_EQUAL(m1(1,3), 8.0);               
            BOOST_CHECK_EQUAL(m1(2,0), 9.0);               
            BOOST_CHECK_EQUAL(m1(2,1), 10.0);               
            BOOST_CHECK_EQUAL(m1(2,2), 13.0);               
            BOOST_CHECK_EQUAL(m1(2,3), 14.0);               
            BOOST_CHECK_EQUAL(m1(3,0), 11.0);               
            BOOST_CHECK_EQUAL(m1(3,1), 12.0);               
            BOOST_CHECK_EQUAL(m1(3,2), 15.0);               
            BOOST_CHECK_EQUAL(m1(3,3), 16.0);               
            BOOST_CHECK_EQUAL(m1(4,0), 17.0);               
            BOOST_CHECK_EQUAL(m1(4,1), 18.0);               
            BOOST_CHECK_EQUAL(m1(4,2), 21.0);               
            BOOST_CHECK_EQUAL(m1(4,3), 22.0);               
            BOOST_CHECK_EQUAL(m1(5,0), 19.0);               
            BOOST_CHECK_EQUAL(m1(5,1), 20.0);               
            BOOST_CHECK_EQUAL(m1(5,2), 23.0);               
            BOOST_CHECK_EQUAL(m1(5,3), 24.0);               
            
            BOOST_CHECK_EQUAL((*m2)(0,0), 1.0);               
            BOOST_CHECK_EQUAL((*m2)(0,1), 2.0);               
            BOOST_CHECK_EQUAL((*m2)(0,2), 5.0);               
            BOOST_CHECK_EQUAL((*m2)(0,3), 6.0);               
            BOOST_CHECK_EQUAL((*m2)(1,0), 3.0);               
            BOOST_CHECK_EQUAL((*m2)(1,1), 4.0);               
            BOOST_CHECK_EQUAL((*m2)(1,2), 7.0);               
            BOOST_CHECK_EQUAL((*m2)(1,3), 8.0);               
            BOOST_CHECK_EQUAL((*m2)(2,0), 9.0);               
            BOOST_CHECK_EQUAL((*m2)(2,1), 10.0);               
            BOOST_CHECK_EQUAL((*m2)(2,2), 13.0);               
            BOOST_CHECK_EQUAL((*m2)(2,3), 14.0);               
            BOOST_CHECK_EQUAL((*m2)(3,0), 11.0);               
            BOOST_CHECK_EQUAL((*m2)(3,1), 12.0);               
            BOOST_CHECK_EQUAL((*m2)(3,2), 15.0);               
            BOOST_CHECK_EQUAL((*m2)(3,3), 16.0);               
            BOOST_CHECK_EQUAL((*m2)(4,0), 17.0);               
            BOOST_CHECK_EQUAL((*m2)(4,1), 18.0);               
            BOOST_CHECK_EQUAL((*m2)(4,2), 21.0);               
            BOOST_CHECK_EQUAL((*m2)(4,3), 22.0);               
            BOOST_CHECK_EQUAL((*m2)(5,0), 19.0);               
            BOOST_CHECK_EQUAL((*m2)(5,1), 20.0);               
            BOOST_CHECK_EQUAL((*m2)(5,2), 23.0);               
            BOOST_CHECK_EQUAL((*m2)(5,3), 24.0);     
        }

        BOOST_AUTO_TEST_CASE(TestEmptyRowElementMultiplication)
        {
            typedef NekMatrix<NekMatrix<NekDouble>, BlockMatrixTag> BlockMatrixType;
            unsigned int rowSizes[] = {0, 2, 0, 2, 0, 0, 2};
            unsigned int columnSizes[] = {0, 2, 2};
            BlockMatrixType m1(7, 3, rowSizes, columnSizes);
            
            double vals1[] = {1.0, 3.0, 2.0, 4.0};
            double vals2[] = {5.0, 7.0, 6.0, 8.0};
            double vals3[] = {9.0, 11.0, 10.0, 12.0};
            double vals4[] = {13.0, 15.0, 14.0, 16.0};
            double vals5[] = {17.0, 19.0, 18.0, 20.0};
            double vals6[] = {21.0, 23.0, 22.0, 24.0};
            
            std::shared_ptr<NekMatrix<NekDouble> > sub1(new NekMatrix<NekDouble>(2, 2, vals1));
            std::shared_ptr<NekMatrix<NekDouble> > sub2(new NekMatrix<NekDouble>(2, 2, vals2));
            std::shared_ptr<NekMatrix<NekDouble> > sub3(new NekMatrix<NekDouble>(2, 2, vals3));
            std::shared_ptr<NekMatrix<NekDouble> > sub4(new NekMatrix<NekDouble>(2, 2, vals4));
            std::shared_ptr<NekMatrix<NekDouble> > sub5(new NekMatrix<NekDouble>(2, 2, vals5));
            std::shared_ptr<NekMatrix<NekDouble> > sub6(new NekMatrix<NekDouble>(2, 2, vals6));
            
            std::shared_ptr<NekMatrix<NekDouble> > empty(new NekMatrix<NekDouble>(0, 0));
            m1.SetBlock(0,0, empty);
            m1.SetBlock(0,1, empty);
            m1.SetBlock(0,2, empty);
            m1.SetBlock(1,0, empty);
            m1.SetBlock(1,1, sub1);
            m1.SetBlock(1,2, sub2);
            m1.SetBlock(2,0, empty);
            m1.SetBlock(3,1, sub3);
            m1.SetBlock(3,2, sub4);
            m1.SetBlock(6,1, sub5);
            m1.SetBlock(6,2, sub6);
                                       
            double x_buf[] = {1, 2, 3, 4};
            NekVector<double> x(4, x_buf);

            double expected_result_buf[] = { 44, 64, 124, 144, 204, 224};
            NekVector<double> expected_result(6, expected_result_buf);

            NekVector<double> result1 = m1*x;
            BOOST_CHECK_EQUAL(result1, expected_result);

        }

    }
}
