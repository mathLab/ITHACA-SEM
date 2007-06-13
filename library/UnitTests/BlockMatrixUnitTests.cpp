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

#include <UnitTests/testNekMatrix.h>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/progress.hpp>
#include <iostream>

namespace Nektar
{
    namespace BlockMatrixUnitTests
    {
        void TestEqualSizedBlockConstruction()
        {
            {
                typedef NekMatrix<NekMatrix<NekDouble>, FullMatrixTag, BlockMatrixTag> BlockMatrixType;
                BlockMatrixType m1(3, 2, 2, 2);
                boost::shared_ptr<ConstMatrix<NekDouble> > m2(new BlockMatrixType(2, 3, 2, 2));
                boost::shared_ptr<Matrix<NekDouble> > m3(new BlockMatrixType(3, 3, 2, 2));
                
                BOOST_CHECK_EQUAL(m1.GetRows(), 6);
                BOOST_CHECK_EQUAL(m1.GetColumns(), 4);
                BOOST_CHECK_EQUAL(m2->GetRows(), 4);
                BOOST_CHECK_EQUAL(m2->GetColumns(), 6);
                BOOST_CHECK_EQUAL(m3->GetRows(), 6);
                BOOST_CHECK_EQUAL(m3->GetColumns(), 6);
               
            }
            
            // Bad Input
            {
            }
            
        }
        
        void TestVariableSizedBlockConstruction()
        {
        }
        
        void TestElementAccess()
        {
            
            
        }
        
        
        void TestGetNumElements()
        {
            
        }
        
        void TestGetStorageType()
        {
            
        }
    }
}
