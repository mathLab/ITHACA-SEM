///////////////////////////////////////////////////////////////////////////////
//
// File: ScaledMatrixUnitTests.cpp
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

#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/progress.hpp>
#include <iostream>

namespace Nektar
{
    namespace ScaledMatrixUnitTests
    {
        BOOST_AUTO_TEST_CASE(TestConstruction)
        {
            typedef NekMatrix<double> OwnedType;
            double d[] = {1.0, 2.0, 3.0, 4.0,
                          5.0, 6.0, 7.0, 8.0,
                          9.0, 10.0, 11.0, 12.0};
            std::shared_ptr<OwnedType> o(new OwnedType(3, 4, d));
            NekMatrix<OwnedType, ScaledMatrixTag> m1(2.0, o);
            std::shared_ptr<ConstMatrix<double> > m2(new NekMatrix<OwnedType, ScaledMatrixTag>(3.0, o));
            std::shared_ptr<NekMatrix<OwnedType, ScaledMatrixTag> > m3(new NekMatrix<OwnedType, ScaledMatrixTag>(0.0, o));
            
            BOOST_CHECK_EQUAL(m1.Scale(), 2.0);
            //BOOST_CHECK_EQUAL(m2->Scale(), 3.0);
            BOOST_CHECK_EQUAL(m3->Scale(), 0.0);
                        
            BOOST_CHECK_EQUAL(m1.GetOwnedMatrix(), o);            
            //BOOST_CHECK_EQUAL(m2->GetOwnedMatrix(), o);            
            BOOST_CHECK_EQUAL(m3->GetOwnedMatrix(), o);            
            
        }
        
        BOOST_AUTO_TEST_CASE(TestElementAccessScaledMatrix)
        {
            typedef NekMatrix<double> OwnedType;
            double d[] = {1.0, 5.0, 9.0,
                          2.0, 6.0, 10.0,
                          3.0, 7.0, 11.0,
                          4.0, 8.0, 12.0};
            std::shared_ptr<OwnedType> o(new OwnedType(3, 4, d));
            NekMatrix<OwnedType, ScaledMatrixTag> m1(2.0, o);
            std::shared_ptr<ConstMatrix<double> > m2(new NekMatrix<OwnedType, ScaledMatrixTag>(3.0, o));
            std::shared_ptr<NekMatrix<OwnedType, ScaledMatrixTag> > m3(new NekMatrix<OwnedType, ScaledMatrixTag>(0.0, o));
            
            BOOST_CHECK_EQUAL(m1(0,0), 2.0);
            BOOST_CHECK_EQUAL(m1(0,1), 4.0);
            BOOST_CHECK_EQUAL(m1(0,2), 6.0);
            BOOST_CHECK_EQUAL(m1(0,3), 8.0);
            BOOST_CHECK_EQUAL(m1(1,0), 10.0);
            BOOST_CHECK_EQUAL(m1(1,1), 12.0);
            BOOST_CHECK_EQUAL(m1(1,2), 14.0);
            BOOST_CHECK_EQUAL(m1(1,3), 16.0);
            BOOST_CHECK_EQUAL(m1(2,0), 18.0);
            BOOST_CHECK_EQUAL(m1(2,1), 20.0);
            BOOST_CHECK_EQUAL(m1(2,2), 22.0);
            BOOST_CHECK_EQUAL(m1(2,3), 24.0);
            
            BOOST_CHECK_EQUAL((*m2)(0,0), 3.0);
            BOOST_CHECK_EQUAL((*m2)(0,1), 6.0);
            BOOST_CHECK_EQUAL((*m2)(0,2), 9.0);
            BOOST_CHECK_EQUAL((*m2)(0,3), 12.0);
            BOOST_CHECK_EQUAL((*m2)(1,0), 15.0);
            BOOST_CHECK_EQUAL((*m2)(1,1), 18.0);
            BOOST_CHECK_EQUAL((*m2)(1,2), 21.0);
            BOOST_CHECK_EQUAL((*m2)(1,3), 24.0);
            BOOST_CHECK_EQUAL((*m2)(2,0), 27.0);
            BOOST_CHECK_EQUAL((*m2)(2,1), 30.0);
            BOOST_CHECK_EQUAL((*m2)(2,2), 33.0);
            BOOST_CHECK_EQUAL((*m2)(2,3), 36.0);
            
            BOOST_CHECK_EQUAL((*m3)(0,0), 0.0);
            BOOST_CHECK_EQUAL((*m3)(0,1), 0.0);
            BOOST_CHECK_EQUAL((*m3)(0,2), 0.0);
            BOOST_CHECK_EQUAL((*m3)(0,3), 0.0);
            BOOST_CHECK_EQUAL((*m3)(1,0), 0.0);
            BOOST_CHECK_EQUAL((*m3)(1,1), 0.0);
            BOOST_CHECK_EQUAL((*m3)(1,2), 0.0);
            BOOST_CHECK_EQUAL((*m3)(1,3), 0.0);
            BOOST_CHECK_EQUAL((*m3)(2,0), 0.0);
            BOOST_CHECK_EQUAL((*m3)(2,1), 0.0);
            BOOST_CHECK_EQUAL((*m3)(2,2), 0.0);
            BOOST_CHECK_EQUAL((*m3)(2,3), 0.0);
            
        }
        
        
        BOOST_AUTO_TEST_CASE(TestGetNumElements)
        {
            typedef NekMatrix<double> OwnedType;
            double d[] = {1.0, 5.0, 9.0,
                          2.0, 6.0, 10.0,
                          3.0, 7.0, 11.0,
                          4.0, 8.0, 12.0};
            std::shared_ptr<OwnedType> o(new OwnedType(3, 4, d));
            NekMatrix<OwnedType, ScaledMatrixTag> m1(2.0, o);
            std::shared_ptr<ConstMatrix<double> > m2(new NekMatrix<OwnedType, ScaledMatrixTag>(3.0, o));
            std::shared_ptr<NekMatrix<OwnedType, ScaledMatrixTag> > m3(new NekMatrix<OwnedType, ScaledMatrixTag>(0.0, o));
            
            BOOST_CHECK_EQUAL(m1.GetStorageSize(), 12);
            BOOST_CHECK_EQUAL(m2->GetStorageSize(), 12);
            BOOST_CHECK_EQUAL(m3->GetStorageSize(), 12);
        }
        
        BOOST_AUTO_TEST_CASE(TestGetStorageType)
        {
            typedef NekMatrix<double> OwnedType;
            double d[] = {1.0, 5.0, 9.0,
                          2.0, 6.0, 10.0,
                          3.0, 7.0, 11.0,
                          4.0, 8.0, 12.0};
            std::shared_ptr<OwnedType> o(new OwnedType(3, 4, d));
            NekMatrix<OwnedType, ScaledMatrixTag> m1(2.0, o);
            std::shared_ptr<ConstMatrix<double> > m2(new NekMatrix<OwnedType, ScaledMatrixTag>(3.0, o));
            std::shared_ptr<NekMatrix<OwnedType, ScaledMatrixTag> > m3(new NekMatrix<OwnedType, ScaledMatrixTag>(0.0, o));
            
            BOOST_CHECK_EQUAL(m1.GetStorageType(), eFULL);
            BOOST_CHECK_EQUAL(m2->GetStorageType(), eFULL);
            BOOST_CHECK_EQUAL(m3->GetStorageType(), eFULL);
        }
    }
}
