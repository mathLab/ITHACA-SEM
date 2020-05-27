///////////////////////////////////////////////////////////////////////////////
//
// File: testNekLinAlgAlgorithms.cpp
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

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <LibUtilities/LinearAlgebra/NekLinAlgAlgorithms.hpp>
#include <LibUtilities/LinearAlgebra/NekVector.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <vector>

namespace Nektar
{
    namespace NekLinAlgTests
    {
        BOOST_AUTO_TEST_CASE(TestGramSchmidtOrthogonalizationBookExample)
        {
            double buf_0[] = {1.0, 0.0, 2.0};
            double buf_1[] = {2.0, 3.0, 0.0};
            
            std::vector<NekVector<double> > x;
            x.push_back(NekVector<double>(3, buf_0));
            x.push_back(NekVector<double>(3, buf_1));
            
            std::vector<NekVector<double> > q = GramSchmidtOrthogonalization(x);
            
            BOOST_CHECK_EQUAL(q.size(), 2);
            
            double epsilon = 1e-1;
            BOOST_CHECK_CLOSE(q[0][0], .4472, epsilon);
            BOOST_CHECK_CLOSE(q[0][1], .0, epsilon);
            BOOST_CHECK_CLOSE(q[0][2], .8942, epsilon);
            
            BOOST_CHECK_CLOSE(q[1][0], .45811, epsilon);
            BOOST_CHECK_CLOSE(q[1][1], .85890, epsilon);
            BOOST_CHECK_CLOSE(q[1][2], -.22890, epsilon);
        }
    }
}
