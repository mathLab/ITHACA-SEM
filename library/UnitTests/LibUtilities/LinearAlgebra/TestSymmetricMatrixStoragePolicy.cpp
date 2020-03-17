///////////////////////////////////////////////////////////////////////////////
//
// File: TestSymmetricMatrixStoragePolicy.cpp
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
#include <UnitTests/CountedObject.h>
#include <UnitTests/util.h>
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <boost/test/auto_unit_test.hpp>

namespace Nektar
{
    namespace SymmetricMatrixStoragePolicyUnitTests
    {
        typedef SymmetricMatrixFuncs Policy;

                    
        BOOST_AUTO_TEST_CASE(TestAdvanceSymmetric)
        {
            UnitTests::RedirectCerrIfNeeded();
            {

                unsigned int curRow = 0; 
                unsigned int curColumn = 0;
                std::tie(curRow, curColumn) = Policy::Advance(3, 3, curRow, curColumn);
                BOOST_CHECK_EQUAL(1, curRow);
                BOOST_CHECK_EQUAL(0, curColumn);

                std::tie(curRow, curColumn) = Policy::Advance(3, 3, curRow, curColumn);
                BOOST_CHECK_EQUAL(2, curRow);
                BOOST_CHECK_EQUAL(0, curColumn);

                std::tie(curRow, curColumn) = Policy::Advance(3, 3, curRow, curColumn);
                BOOST_CHECK_EQUAL(0, curRow);
                BOOST_CHECK_EQUAL(1, curColumn);

                std::tie(curRow, curColumn) = Policy::Advance(3, 3, curRow, curColumn);
                BOOST_CHECK_EQUAL(1, curRow);
                BOOST_CHECK_EQUAL(1, curColumn);

                std::tie(curRow, curColumn) = Policy::Advance(3, 3, curRow, curColumn);
                BOOST_CHECK_EQUAL(2, curRow);
                BOOST_CHECK_EQUAL(1, curColumn);

                std::tie(curRow, curColumn) = Policy::Advance(3, 3, curRow, curColumn);
                BOOST_CHECK_EQUAL(0, curRow);
                BOOST_CHECK_EQUAL(2, curColumn);

                std::tie(curRow, curColumn) = Policy::Advance(3, 3, curRow, curColumn);
                BOOST_CHECK_EQUAL(1, curRow);
                BOOST_CHECK_EQUAL(2, curColumn);

                std::tie(curRow, curColumn) = Policy::Advance(3, 3, curRow, curColumn);
                BOOST_CHECK_EQUAL(2, curRow);
                BOOST_CHECK_EQUAL(2, curColumn);

                std::tie(curRow, curColumn) = Policy::Advance(3, 3, curRow, curColumn);
                BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), curRow);
                BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), curColumn);
            }

            {

                unsigned int curRow = 0; 
                unsigned int curColumn = 0;
                std::tie(curRow, curColumn) = Policy::Advance(1, 1, curRow, curColumn);
                BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), curRow);
                BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), curColumn);
            }

            {

                unsigned int curRow = 0; 
                unsigned int curColumn = 0;
                std::tie(curRow, curColumn) = Policy::Advance(2, 2, curRow, curColumn);
                BOOST_CHECK_EQUAL(1, curRow);
                BOOST_CHECK_EQUAL(0, curColumn);

                std::tie(curRow, curColumn) = Policy::Advance(2, 2, curRow, curColumn);
                BOOST_CHECK_EQUAL(0, curRow);
                BOOST_CHECK_EQUAL(1, curColumn);

                std::tie(curRow, curColumn) = Policy::Advance(2, 2, curRow, curColumn);
                BOOST_CHECK_EQUAL(1, curRow);
                BOOST_CHECK_EQUAL(1, curColumn);

                std::tie(curRow, curColumn) = Policy::Advance(2, 2, curRow, curColumn);
                BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), curRow);
                BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), curColumn);
            }
        }
    }
}


