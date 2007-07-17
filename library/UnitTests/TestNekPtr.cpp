///////////////////////////////////////////////////////////////////////////////
//
// File: TestNekPtr.cpp
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

#include <UnitTests/TestNekPtr.h>
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/progress.hpp>
#include <iostream>

#include <LibUtilities/BasicUtils/SharedPtr.hpp>
#include <UnitTests/CountedObject.h>

namespace Nektar
{
    namespace PtrUnitTests
    {
        void TestScope()
        {
            CountedObject<double>::ClearCounters();

            {
                ptr<CountedObject<double> > p(new CountedObject<double>());
                ptr<CountedObject<double> > p1(p);

                {
                    ptr<CountedObject<double> > p2(new CountedObject<double>());
                    CountedObject<double>::check(2, 0, 0, 0, 0, 0);
                }
                CountedObject<double>::check(2, 0, 1, 0, 0, 0);
            }
            CountedObject<double>::check(2, 0, 2, 0, 0, 0);
        }

        void TestConstConversions()
        {
            CountedObject<double>::ClearCounters();
            {
                ptr<CountedObject<double> > p(new CountedObject<double>());
                const_ptr<CountedObject<double> > p1(p);

                CountedObject<double>::check(1, 0, 0, 0, 0, 0);
            }
            CountedObject<double>::check(1, 0, 1, 0, 0, 0);
        }
    }
}
