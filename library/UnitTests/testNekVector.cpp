///////////////////////////////////////////////////////////////////////////////
//
// File: testNekVector.cpp
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
// Description: Test code for NekVector
//
///////////////////////////////////////////////////////////////////////////////

#include <UnitTests/testNekVector.h>
#include <LibUtilities/LinearAlgebra/NekVector.hpp>

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/progress.hpp>

#include <iostream>

using std::cout;
using std::endl;

namespace Nektar
{
    namespace UnitTests
    {
        using namespace Nektar;

        class VectorTestClass
        {
            public:
                VectorTestClass() : m_dataValue(0) {}
                explicit VectorTestClass(int data) : m_dataValue(data) {}
                ~VectorTestClass(){}

                VectorTestClass operator=(const VectorTestClass& rhs)
                {
                    m_dataValue = rhs.m_dataValue;
                    return *this;
                }

                int value() const { return m_dataValue; }

            private:
                int m_dataValue;
        };

        bool operator==(const VectorTestClass& lhs, const VectorTestClass& rhs)
        {
            return lhs.value() == rhs.value();
        }

        bool operator!=(const VectorTestClass& lhs, const VectorTestClass& rhs)
        {
            return !(lhs == rhs);
        }


        void testNekVectorConstruction()
        {
            // Test the constructors on a numeric type.
            {
                NekVector<double> p1(10, 2.7);
                BOOST_CHECK(p1.GetDimension() == 10);
                for(unsigned int i = 0; i < p1.GetDimension(); ++i)
                {
                    BOOST_CHECK(p1(i) == 2.7);
                    BOOST_CHECK(p1[i] == 2.7);
                }

                NekVector<double> p2(1,0.0);
                BOOST_CHECK(p2.GetDimension() == 1);
                BOOST_CHECK(p2(0) == 0.0);

                NekVector<double> p3(p1);
                BOOST_CHECK(p3.GetDimension() == 10);
                for(unsigned int i = 0; i < p3.GetDimension(); ++i)
                {
                    BOOST_CHECK(p3(i) == 2.7);
                    BOOST_CHECK(p3[i] == 2.7);
                }

                p2 = p3;
                BOOST_CHECK(p2.GetDimension() == 10);
                for(unsigned int i = 0; i < p2.GetDimension(); ++i)
                {
                    BOOST_CHECK_EQUAL(p2(i), 2.7);
                    BOOST_CHECK_EQUAL(p2[i], 2.7);
                }

                NekVector<int> v10(10, -2);
                for(int i = 0; i < 10; ++i)
                {
                    BOOST_CHECK(v10[i] == -2);
                }
            }

            {
                // Now do the same tests with constant sized vectors.
                NekVector<double, 10> p1(2.7);
                BOOST_CHECK(p1.GetDimension() == 10);
                for(unsigned int i = 0; i < p1.GetDimension(); ++i)
                {
                    BOOST_CHECK(p1(i) == 2.7);
                    BOOST_CHECK(p1[i] == 2.7);
                }

                NekVector<double, 1> p2;
                BOOST_CHECK(p2.GetDimension() == 1);
                BOOST_CHECK(p2(0) == 0.0);

                NekVector<double, 10> p3(p1);
                BOOST_CHECK(p3.GetDimension() == 10);
                for(unsigned int i = 0; i < p3.GetDimension(); ++i)
                {
                    BOOST_CHECK(p3(i) == 2.7);
                    BOOST_CHECK(p3[i] == 2.7);
                }

                NekVector<double, 10> p4;
                p4 = p3;
                BOOST_CHECK(p4.GetDimension() == 10);
                for(unsigned int i = 0; i < p4.GetDimension(); ++i)
                {
                    BOOST_CHECK(p4(i) == 2.7);
                    BOOST_CHECK(p4[i] == 2.7);
                }

                NekVector<int> v10(10, -2);
                for(int i = 0; i < 10; ++i)
                {
                    BOOST_CHECK(v10[i] == -2);
                }
            }


            // Now test it on an arbitrary type.
            {
                NekVector<VectorTestClass, 3> v1;
                BOOST_CHECK(v1.x() == VectorTestClass(0));
                BOOST_CHECK(v1.y() == VectorTestClass(0));
                BOOST_CHECK(v1.z() == VectorTestClass(0));

                NekVector<VectorTestClass, 10> p2(VectorTestClass(-2));
                for(int i = 0; i < 10; ++i)
                {
                    BOOST_CHECK(p2[i] == VectorTestClass(-2));
                }

                NekVector<VectorTestClass, 3> p3(v1);
                BOOST_CHECK(v1==p3);
            }
        }

        void testNekVectorOperators()
        {
            NekVector<double> v1(3, 1.0);
            v1(0) = 1.1;
            v1(1) = 1.2;
            v1(2) = 1.3;

            BOOST_CHECK(v1(0) == 1.1);
            BOOST_CHECK(v1(1) == 1.2);
            BOOST_CHECK(v1(2) == 1.3);

            v1[0] = 1.4;
            v1[1] = 1.5;
            v1[2] = 1.6;

            BOOST_CHECK(v1[0] == 1.4);
            BOOST_CHECK(v1[1] == 1.5);
            BOOST_CHECK(v1[2] == 1.6);

            v1.x() = 1.7;
            v1.y() = 1.8;
            v1.z() = 1.9;

            BOOST_CHECK(v1[0] == 1.7);
            BOOST_CHECK(v1[1] == 1.8);
            BOOST_CHECK(v1[2] == 1.9);

            BOOST_CHECK(v1.x() == 1.7);
            BOOST_CHECK(v1.y() == 1.8);
            BOOST_CHECK(v1.z() == 1.9);

            NekVector<double> v2(3, 1.0);
            v2.x() = 1.7;
            v2.y() = 1.8;
            v2.z() = 1.9;

            BOOST_CHECK(v1 == v2);
            BOOST_CHECK(!(v1 != v2));

            NekVector<double> v3(4, 2.0);
            BOOST_CHECK(v3 != v1);

            NekVector<double> v4(3, 0.0);
            BOOST_CHECK(v4 != v1);

        }

        void testNekVectorArithmetic()
        {
            {
                NekVector<double, 3> v1(1.0);
                v1(0) = 1.0;
                v1(1) = 2.0;
                v1(2) = 3.0;

                NekVector<double, 3> v2 = -v1;
                BOOST_CHECK(v2(0) == -v1(0));
                BOOST_CHECK(v2(1) == -v1(1));
                BOOST_CHECK(v2(2) == -v1(2));
            }

            {
                NekVector<double, 3> lhs;
                NekVector<double, 3> rhs;

                lhs(0) = 1.0;
                lhs(1) = 2.0;
                lhs(2) = 3.0;
                rhs(0) = 4.0;
                rhs(1) = 5.0;
                rhs(2) = 6.0;

                rhs += lhs;

                BOOST_CHECK(rhs(0) == 5.0);
                BOOST_CHECK(rhs(1) == 7.0);
                BOOST_CHECK(rhs(2) == 9.0);

                rhs -= lhs;

                BOOST_CHECK(rhs(0) == 4.0);
                BOOST_CHECK(rhs(1) == 5.0);
                BOOST_CHECK(rhs(2) == 6.0);

                NekVector<double, 3> add_result = rhs+lhs;
                NekVector<double, 3> sub_result = rhs-lhs;

                BOOST_CHECK(add_result(0) == 5.0);
                BOOST_CHECK(add_result(1) == 7.0);
                BOOST_CHECK(add_result(2) == 9.0);

                BOOST_CHECK(sub_result(0) == 3.0);
                BOOST_CHECK(sub_result(1) == 3.0);
                BOOST_CHECK(sub_result(2) == 3.0);
            }

            {
                NekVector<double, 3> v1;
                v1(0) = 1.0;
                v1(1) = 2.0;
                v1(2) = 3.0;

                NekVector<double, 3> v2(v1);

                v1 *= 1.2;
                BOOST_CHECK(v1(0) == 1.0*1.2);
                BOOST_CHECK(v1(1) == 2.0*1.2);
                BOOST_CHECK(v1(2) == 3.0*1.2);

                NekVector<double, 3> v3 = v2*1.2;
                NekVector<double, 3> v4 = 1.2*v2;
                BOOST_CHECK(v3 == v1);
                BOOST_CHECK(v4 == v1);
            }

            {
                NekVector<double, 3> v1;
                v1(0) = 2.0;
                v1(1) = 6.0;
                v1(2) = 9.0;

                NekVector<double, 3> v2 = v1;
                Normalize(v2);
                BOOST_CHECK(v1.Magnitude() == 11.0);
                v1.Normalize();
                BOOST_CHECK(v1(0) == 2.0/11.0);
                BOOST_CHECK(v1(1) == 6.0/11.0);
                BOOST_CHECK(v1(2) == 9.0/11.0);
                BOOST_CHECK(v2 == v1);
            }

            {
                NekVector<double, 3> v1;
                v1(0) = 1.0;
                v1(1) = 10.0;
                v1(2) = -3.4;

                BOOST_CHECK(v1.Magnitude() == sqrt((1.0*1.0) + (10.0*10.0) + (-3.4*-3.4)));

                //NekVector<double> v2;
                //BOOST_CHECK(v2.magnitude() == 0.0);
            }

            {
                NekVector<double, 3> v1;
                v1(0) = 1.0;
                v1(1) = 10.0;
                v1(2) = -3.4;

                NekVector<double, 3> v2;
                v2(0) = 8.9;
                v2(1) = -4.5;
                v2(2) = 9.12;

                BOOST_CHECK(v1.Dot(v2) == 1.0*8.9 + 10.0*-4.5 + 9.12*-3.4);
            }

        }
        
        void testNorms()
        {
            double vals[] = {1,-2,3};
            NekVector<double, 3> v(vals);
                
            double epsilon = 1e-11;
            BOOST_CHECK_EQUAL(v.L1Norm(), 1.0 + 2.0 + 3.0);
            BOOST_CHECK_CLOSE(v.L2Norm(), sqrt(1.0+4.0+9.0), epsilon);
            BOOST_CHECK_EQUAL(v.InfinityNorm(), 3);
        }
    }
}

/**
    $Log: testNekVector.cpp,v $
    Revision 1.5  2006/11/18 17:18:46  bnelson
    Added L1, L2, and Infinity norm tests.

    Revision 1.4  2006/09/30 15:38:29  bnelson
    no message

    Revision 1.3  2006/09/15 02:01:16  bnelson
    no message

    Revision 1.2  2006/06/05 02:23:17  bnelson
    Updates for the reorganization of LibUtilities.

    Revision 1.1  2006/05/04 18:59:56  kirby
    *** empty log message ***

    Revision 1.1  2006/04/11 02:02:13  bnelson
    Added more tests.


**/
