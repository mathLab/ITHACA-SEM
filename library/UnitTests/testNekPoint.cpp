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

#include <UnitTests/util.h>
#include <LibUtilities/LinearAlgebra/NekPoint.hpp>
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
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
                PointTestClass(const PointTestClass&) = default;
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

        class TestPoint : public Nektar::NekPoint<double>
        {
            public:
                TestPoint() {}
                TestPoint(const TestPoint& rhs) :
                Nektar::NekPoint<double>(rhs)
                {
                }
        };

        void test()
        {
            TestPoint p;
            TestPoint p1(p);
        }

        struct TenD 
        {
            static const unsigned int Value = 10;
        };
        
        BOOST_AUTO_TEST_CASE(testNekPointConstruction)
        {
            using namespace Nektar;

            // Test the constructors on a numeric type.
            {
                NekPoint<double> p1;
                BOOST_CHECK(p1.x() == 0.0);
                BOOST_CHECK(p1.y() == 0.0);
                BOOST_CHECK(p1.z() == 0.0);


                NekPoint<double> p3(p1);
                BOOST_CHECK(p1==p3);
            }

            // Now test it on an arbitrary type.
            {
                NekPoint<PointTestClass> p1;
                BOOST_CHECK(p1.x() == PointTestClass(0));
                BOOST_CHECK(p1.y() == PointTestClass(0));
                BOOST_CHECK(p1.z() == PointTestClass(0));

                BOOST_CHECK(p1.x() == PointTestClass());
                BOOST_CHECK(p1.y() == PointTestClass());
                BOOST_CHECK(p1.z() == PointTestClass());

                NekPoint<PointTestClass> p3(p1);
                BOOST_CHECK(p1==p3);
            }

        }

        

    }
}
