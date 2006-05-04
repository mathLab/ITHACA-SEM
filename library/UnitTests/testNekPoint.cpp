///////////////////////////////////////////////////////////////////////////////
//
// File: testNekPoint.cpp
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
// Description: Test code for NekPoint
//
///////////////////////////////////////////////////////////////////////////////

#include <UnitTests/testNekPoint.h>
#include <LibUtilities/NekPoint.hpp>

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/progress.hpp>

namespace Nektar
{
    namespace UnitTests
    {
        class PointTestClass
        {
            public:
                PointTestClass() : m_dataValue(0) {}
                explicit PointTestClass(int data) : m_dataValue(data) {}
                ~PointTestClass(){}

                PointTestClass operator=(const PointTestClass& rhs)
                {
                    m_dataValue = rhs.m_dataValue;
                    return *this;
                }

                int value() const { return m_dataValue; }

            private:
                int m_dataValue;
        };

        bool operator==(const PointTestClass& lhs, const PointTestClass& rhs)
        {
            return lhs.value() == rhs.value();
        }

        bool operator!=(const PointTestClass& lhs, const PointTestClass& rhs)
        {
            return !(lhs == rhs);
        }


        void testNekPointConstruction()
        {
            using namespace Nektar::LibUtilities;

            // Test the constructors on a numeric type.
            {
                NekPoint<double, 3> p1;
                BOOST_CHECK(p1.x() == 0.0);
                BOOST_CHECK(p1.y() == 0.0);
                BOOST_CHECK(p1.z() == 0.0);

                NekPoint<int, 10> p2(-2);
                for(int i = 0; i < 10; ++i)
                {
                    BOOST_CHECK(p2[i] == -2);
                }

                NekPoint<double, 3> p3(p1);
                BOOST_CHECK(p1==p3);
            }

            // Now test it on an arbitrary type.
            {
                NekPoint<PointTestClass, 3> p1;
                BOOST_CHECK(p1.x() == PointTestClass(0));
                BOOST_CHECK(p1.y() == PointTestClass(0));
                BOOST_CHECK(p1.z() == PointTestClass(0));

                NekPoint<PointTestClass, 10> p2(PointTestClass(-2));
                for(int i = 0; i < 10; ++i)
                {
                    BOOST_CHECK(p2[i] == PointTestClass(-2));
                }

                NekPoint<PointTestClass, 3> p3(p1);
                BOOST_CHECK(p1==p3);
            }
        }

        void testNekPointArithmetic()
        {
        }
    }
}

/**
    $Log: testNekPoint.cpp,v $
    Revision 1.3  2006/04/11 02:02:13  bnelson
    Added more tests.

    Revision 1.2  2006/04/06 03:45:50  bnelson
    Added some more point tests and fixed compiler errors.

    Revision 1.1  2006/01/31 14:06:23  bnelson
    Added the new UnitTest project.

**/
