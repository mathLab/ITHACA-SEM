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

        class TestPoint : public Nektar::NekPoint<double, 3, 0>
        {
            public:
                TestPoint() {}
                TestPoint(const TestPoint& rhs) :
                Nektar::NekPoint<double, 3, 0>(rhs)
                {
                }
        };

        void test()
        {
            TestPoint p;
            TestPoint p1(p);
        }

        BOOST_AUTO_TEST_CASE(testNekPointConstruction)
        {
            using namespace Nektar;

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

                BOOST_CHECK(p1.x() == PointTestClass());
                BOOST_CHECK(p1.y() == PointTestClass());
                BOOST_CHECK(p1.z() == PointTestClass());

                NekPoint<PointTestClass, 10> p2(PointTestClass(-2));
                for(int i = 0; i < 10; ++i)
                {
                    BOOST_CHECK(p2[i] == PointTestClass(-2));
                }

                NekPoint<PointTestClass, 3> p3(p1);
                BOOST_CHECK(p1==p3);
            }

            // Test string construction.
            {
                std::string testValues("<1,2>");

                NekPoint<int, 2> p2;
                NekPoint<int, 3> p3;

                BOOST_CHECK(fromString(testValues, p2));
                BOOST_CHECK(p2.x() == 1);
                BOOST_CHECK(p2.y() == 2);

                BOOST_CHECK(!fromString(testValues, p3));

                try
                {
                    NekPoint<int, 3> p4(testValues);
                    BOOST_CHECK(false);
                }
                catch(std::runtime_error&)
                {
                    BOOST_CHECK(true);
                }

            }

            {
                NekPoint<int, 3> p(1, 2, 3);
                BOOST_CHECK(p.x() == 1);
                BOOST_CHECK(p.y() == 2);
                BOOST_CHECK(p.z() == 3);

                NekPoint<int, 3> p1(p);
                BOOST_CHECK(p == p1);

                NekPoint<int, 3> p2;
                p2 = p;
                BOOST_CHECK(p2 == p);
            }
        }

        BOOST_AUTO_TEST_CASE(testNekPointDataAccess)
        {
            using namespace Nektar;

            NekPoint<int, 3> p(1,2,3);
            BOOST_CHECK((p.x() == p.a()) && (p.x() == p.r()));
            BOOST_CHECK((p.y() == p.b()) && (p.y() == p.s()));
            BOOST_CHECK((p.z() == p.c()) && (p.z() == p.t()));

            p.SetX(10);
            p.SetY(11);
            p.SetZ(12);

            BOOST_CHECK((p.x() == p.a()) && (p.x() == p.r()) && (p.x() == 10));
            BOOST_CHECK((p.y() == p.b()) && (p.y() == p.s()) && (p.y() == 11));
            BOOST_CHECK((p.z() == p.c()) && (p.z() == p.t()) && (p.z() == 12));

            BOOST_CHECK(p.x() == p(0));
            BOOST_CHECK(p.y() == p(1));
            BOOST_CHECK(p.z() == p(2));

            BOOST_CHECK(p.x() == p[0]);
            BOOST_CHECK(p.y() == p[1]);
            BOOST_CHECK(p.z() == p[2]);

            BOOST_CHECK_THROW(p(3), std::runtime_error);
            BOOST_CHECK_NO_THROW(p[3]);

            // Test constant versions.
            const NekPoint<int, 3> p2(p);

            BOOST_CHECK((p2.x() == p2.a()) && (p2.x() == p2.r()) && (p2.x() == 10));
            BOOST_CHECK((p2.y() == p2.b()) && (p2.y() == p2.s()) && (p2.y() == 11));
            BOOST_CHECK((p2.z() == p2.c()) && (p2.z() == p2.t()) && (p2.z() == 12));

            BOOST_CHECK(p2.x() == p2(0));
            BOOST_CHECK(p2.y() == p2(1));
            BOOST_CHECK(p2.z() == p2(2));

            BOOST_CHECK(p2.x() == p2[0]);
            BOOST_CHECK(p2.y() == p2[1]);
            BOOST_CHECK(p2.z() == p2[2]);

            BOOST_CHECK_THROW(p2(3), std::runtime_error);
            BOOST_CHECK_NO_THROW(p2[3]);


        }

        BOOST_AUTO_TEST_CASE(testNekPointPointerManipulation)
        {
            using namespace Nektar;
            NekPoint<int, 3> p(1,2,3);
            const int* ptr = p.GetPtr();
            int expected[] = {1, 2, 3};
            BOOST_CHECK(memcmp(expected, ptr, sizeof(int)*3) == 0);
        }

        BOOST_AUTO_TEST_CASE(testNekPointComparison)
        {
            using namespace Nektar;
            NekPoint<int, 3> lhs(1, 2, 3);
            NekPoint<int, 3> rhs(lhs);
            NekPoint<int, 3> ne(3, 2, 1);

            BOOST_CHECK(lhs == rhs);
            BOOST_CHECK(lhs != ne);
        }

        BOOST_AUTO_TEST_CASE(testNekPointOperators)
        {
            using Nektar::NekPoint;

            {
                NekPoint<int, 3> p(1,2,3);
                NekPoint<int, 3> negated = -p;

                BOOST_CHECK(p.x() == -negated.x());
                BOOST_CHECK(p.y() == -negated.y());
                BOOST_CHECK(p.z() == -negated.z());
            }

            {
               NekPoint<int, 3> p(1,2,3);
               p += p;
               BOOST_CHECK(p.x() == 2);
               BOOST_CHECK(p.y() == 4);
               BOOST_CHECK(p.z() == 6);

               p += 8;
               BOOST_CHECK(p.x() == 10);
               BOOST_CHECK(p.y() == 12);
               BOOST_CHECK(p.z() == 14);

               NekPoint<int, 3> sub(9, 3, -5);
               p -= sub;
               BOOST_CHECK(p.x() == 1);
               BOOST_CHECK(p.y() == 9);
               BOOST_CHECK(p.z() == 19);

               p -= 1;
               BOOST_CHECK(p.x() == 0);
               BOOST_CHECK(p.y() == 8);
               BOOST_CHECK(p.z() == 18);

               p *= 3;
               BOOST_CHECK(p.x() == 0);
               BOOST_CHECK(p.y() == 24);
               BOOST_CHECK(p.z() == 54);

               p[0] = 2;
               p /= 2;
               BOOST_CHECK(p.x() == 1);
               BOOST_CHECK(p.y() == 12);
               BOOST_CHECK(p.z() == 27);


               std::string s = p.AsString();
               BOOST_CHECK(s == "(1, 12, 27)");
            }

            {
                Nektar::NekPoint<int, 3, 0> p1(1, 2, 3);
                Nektar::NekPoint<int, 3, 0> p2(10, 20, 30);

                Nektar::NekPoint<int, 3> p3 = p1 + p2;
                BOOST_CHECK(p3.x() == 11);
                BOOST_CHECK(p3.y() == 22);
                BOOST_CHECK(p3.z() == 33);

                Nektar::NekPoint<int, 3> p4 = p1 + 2;
                BOOST_CHECK(p4.x() == 3);
                BOOST_CHECK(p4.y() == 4);
                BOOST_CHECK(p4.z() == 5);

                Nektar::NekPoint<int, 3> p5 = 2 + p1;
                BOOST_CHECK(p5 == p4);

                Nektar::NekPoint<int, 3> p6 = p2 - p1;
                BOOST_CHECK(p6.x() == 9);
                BOOST_CHECK(p6.y() == 18);
                BOOST_CHECK(p6.z() == 27);

                Nektar::NekPoint<int, 3> p7 = p1 - 2;
                BOOST_CHECK(p7.x() == -1);
                BOOST_CHECK(p7.y() == 0);
                BOOST_CHECK(p7.z() == 1);

                Nektar::NekPoint<int, 3> p8 = 2 - p1;
                BOOST_CHECK(p8.x() == 1);
                BOOST_CHECK(p8.y() == 0);
                BOOST_CHECK(p8.z() == -1);

                Nektar::operator*<int, 3, 0>(2, p1);

                p1*(int)2;
                NekPoint<int, 3> p9 = p1*(int)2;
                BOOST_CHECK(p9.x() == 2);
                BOOST_CHECK(p9.y() == 4);
                BOOST_CHECK(p9.z() == 6);

                NekPoint<int, 3> p10 = 2*p1;
                BOOST_CHECK(p9 == p10);

                NekPoint<int, 3> p11 = p2/2;
                BOOST_CHECK(p11.x() == 5);
                BOOST_CHECK(p11.y() == 10);
                BOOST_CHECK(p11.z() == 15);
            }

        }

        BOOST_AUTO_TEST_CASE(testNekPointMisc)
        {
            using namespace Nektar;

            NekPoint<int, 3> p1;
            NekPoint<int, 100> p2;

            BOOST_CHECK(p1.dimension() == 3);
            BOOST_CHECK(p2.dimension() == 100);

            NekPoint<double, 2> source(0.0, 0.0);
            NekPoint<double, 2> dest(1.0, 1.0);

            BOOST_CHECK(distanceBetween(source, dest) == sqrt(2.0));
            BOOST_CHECK(distanceBetween(dest, source) == sqrt(2.0));
        }

        BOOST_AUTO_TEST_CASE(testNekPointAssignment)
        {
            using namespace Nektar;
            NekPoint<int, 3> lhs(1,2,3);
            NekPoint<int, 3> rhs(10, 11, 12);

            lhs = rhs;

            BOOST_CHECK(lhs == rhs);
            BOOST_CHECK(lhs.x() == rhs.x() && lhs.x() == 10);
            BOOST_CHECK(lhs.y() == rhs.y() && lhs.y() == 11);
            BOOST_CHECK(lhs.z() == rhs.z() && lhs.z() == 12);
        }

    }
}
