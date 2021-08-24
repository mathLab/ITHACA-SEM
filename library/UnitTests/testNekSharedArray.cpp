///////////////////////////////////////////////////////////////////////////////
//
// File: testNekSharedArray.h
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
// Description: Test code for Nektar::shared_array
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <UnitTests/CountedObject.h>
#include <UnitTests/util.h>

#include <LibUtilities/LinearAlgebra/NekVector.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>
#include <iostream>

namespace Nektar
{
    namespace SharedArrayUnitTests
    {


        class ParameterTestClass
        {
            public:
                ParameterTestClass() :
                    a(5),
                    b(a)
                {
                    a[0] = 0.0;
                    a[1] = 1.0;
                    a[2] = 2.0;
                    a[3] = 3.0;
                    a[4] = 4.0;
                }

                // out should refer to the same array as a, so any
                // changes in out should be reflected in a.
                void getNonConstByReference(Array<OneD, NekDouble>& out)
                {
                    out = a;
                }

                // Does nothing - the assignment goes into a temporary and is
                // immediately lost.
                void getNonConstByValue(Array<OneD, NekDouble> out)
                {
                    out = a;
                }

                // The following should not compile because the assignment it going into a
                // C++ const variable.
//                 void getNonConstByConstReference(const Array<OneD, NekDouble>& out)
//                 {
//                     out = a;
//                 }
//
//                 void getNonConstByConstValue(const Array<OneD, NekDouble> out)
//                 {
//                     out = a;
//                 }

                // out should refer to the same array as a, so any
                // changes in out should be reflected in a.
                void getConstByReference(Array<OneD, const NekDouble>& out)
                {
                    out = a;
                }

                // Does nothing - the assignment goes into a temporary and is
                // immediately lost.
                void getConstByValue(Array<OneD, const NekDouble> out)
                {
                    out = a;
                }

                Array<OneD, NekDouble> a;
                Array<OneD, const NekDouble> b;
        };

        BOOST_AUTO_TEST_CASE(TestParameterPopulation)
        {
            ParameterTestClass obj;

            Array<OneD, NekDouble> zero(5, 0.0);
            Array<OneD, NekDouble> temp(zero);
            Array<OneD, const NekDouble> const_zero(5, 0.0);
            Array<OneD, const NekDouble> const_temp(const_zero);

            obj.getNonConstByReference(temp);
            BOOST_CHECK(temp == obj.a);

            temp = zero;
            obj.getNonConstByValue(temp);
            BOOST_CHECK(temp == zero);
            BOOST_CHECK(temp != obj.a);

            // should not compile.
            //obj.getNonConstByReference(const_temp);
            //obj.getNonConstByValue(const_temp);

            // Now the ConstArray versions.
            obj.getConstByReference(const_temp);
            BOOST_CHECK(const_temp == obj.a);
            BOOST_CHECK(const_temp != const_zero);
            const_temp = const_zero;
            BOOST_CHECK(const_temp == const_zero);

            obj.getConstByValue(const_temp);
            BOOST_CHECK(const_temp == const_zero);
            BOOST_CHECK(const_temp != obj.a );

        }

        BOOST_AUTO_TEST_CASE(TestEmptyConstructor)
        {
            {
                CountedObject<double>::ClearCounters();
                Array<OneD, const CountedObject<double> > a;
                Array<OneD, CountedObject<double> > b;
                CountedObject<double>::Check(0, 0, 0, 0, 0, 0);

                BOOST_CHECK(a.begin() == a.end());
                BOOST_CHECK(b.begin() == b.end());

                BOOST_CHECK_EQUAL(a.size(), 0);
                BOOST_CHECK_EQUAL(b.size(), 0);

                BOOST_CHECK_EQUAL(a.num_dimensions(), 1);
                BOOST_CHECK_EQUAL(b.num_dimensions(), 1);
            }
            CountedObject<double>::Check(0, 0, 0, 0, 0, 0);

            {
                CountedObject<double>::ClearCounters();
                Array<TwoD, const CountedObject<double> > a;
                Array<TwoD, CountedObject<double> > b;
                CountedObject<double>::Check(0, 0, 0, 0, 0, 0);

                BOOST_CHECK(a.begin() == a.end());
                BOOST_CHECK(b.begin() == b.end());

                BOOST_CHECK_EQUAL(a.size(), 0);
                BOOST_CHECK_EQUAL(b.size(), 0);

                BOOST_CHECK_EQUAL(a.num_dimensions(), 2);
                BOOST_CHECK_EQUAL(b.num_dimensions(), 2);
            }
            CountedObject<double>::Check(0, 0, 0, 0, 0, 0);
        }

        BOOST_AUTO_TEST_CASE(TestUninitializedConstructor)
        {
            {
                CountedObject<double>::ClearCounters();

                Array<OneD, const CountedObject<double> > a(5);
                Array<OneD, CountedObject<double> > b(10);
                CountedObject<double>::Check(15, 0, 0, 0, 0, 0);

                BOOST_CHECK(a.begin() != a.end());
                BOOST_CHECK(b.begin() != b.end());

                BOOST_CHECK_EQUAL(a.size(), 5);
                BOOST_CHECK_EQUAL(b.size(), 10);

                BOOST_CHECK_EQUAL(a.num_dimensions(), 1);
                BOOST_CHECK_EQUAL(b.num_dimensions(), 1);
            }
            CountedObject<double>::Check(15, 0, 15, 0, 0, 0);

            {
                CountedObject<double>::ClearCounters();

                Array<TwoD, const CountedObject<double> > a(5, 10);
                Array<TwoD, CountedObject<double> > b(10, 10);
                CountedObject<double>::Check(150, 0, 0, 0, 0, 0);

                BOOST_CHECK(a.begin() != a.end());
                BOOST_CHECK(b.begin() != b.end());

                BOOST_CHECK_EQUAL(a.size(), 50);
                BOOST_CHECK_EQUAL(b.size(), 100);

                BOOST_CHECK_EQUAL(a.shape()[0], 5);
                BOOST_CHECK_EQUAL(a.shape()[1], 10);
                BOOST_CHECK_EQUAL(b.shape()[0], 10);
                BOOST_CHECK_EQUAL(b.shape()[1], 10);

                BOOST_CHECK_EQUAL(a.num_dimensions(), 2);
                BOOST_CHECK_EQUAL(b.num_dimensions(), 2);
            }
            CountedObject<double>::Check(150, 0, 150, 0, 0, 0);
        }

        BOOST_AUTO_TEST_CASE(TestSingleValueInitialization)
        {
            {
                CountedObject<double> initValue(98);
                CountedObject<double>::ClearCounters();

                Array<OneD, const CountedObject<double> > a(5, initValue);
                Array<OneD, CountedObject<double> > b(10, initValue);
                CountedObject<double>::Check(0, 0, 0, 15, 0, 0);

                BOOST_CHECK(a.begin() != a.end());
                BOOST_CHECK(b.begin() != b.end());

                BOOST_CHECK_EQUAL(a.size(), 5);
                BOOST_CHECK_EQUAL(b.size(), 10);

                BOOST_CHECK_EQUAL(a.num_dimensions(), 1);
                BOOST_CHECK_EQUAL(b.num_dimensions(), 1);

                for(Array<OneD,
                    const CountedObject<double> >::const_iterator iter = a.begin();
                    iter != a.end(); ++iter)
                {
                    BOOST_CHECK(*iter == initValue);
                }

                for(Array<OneD, CountedObject<double> >::iterator iter = b.begin();
                    iter != b.end(); ++iter)
                {
                    BOOST_CHECK(*iter == initValue);
                }
            }
            CountedObject<double>::Check(0, 0, 16, 15, 0, 0);

            {
                CountedObject<double> initValue(98);
                CountedObject<double>::ClearCounters();

                Array<TwoD, const CountedObject<double> > a(5, 10, initValue);
                Array<TwoD, CountedObject<double> > b(10, 10, initValue);
                CountedObject<double>::Check(0, 0, 0, 150, 0, 0);

                BOOST_CHECK(a.begin() != a.end());
                BOOST_CHECK(b.begin() != b.end());

                BOOST_CHECK_EQUAL(a.size(), 50);
                BOOST_CHECK_EQUAL(b.size(), 100);

                BOOST_CHECK_EQUAL(a.shape()[0], 5);
                BOOST_CHECK_EQUAL(a.shape()[1], 10);
                BOOST_CHECK_EQUAL(b.shape()[0], 10);
                BOOST_CHECK_EQUAL(b.shape()[1], 10);

                BOOST_CHECK_EQUAL(a.num_dimensions(), 2);
                BOOST_CHECK_EQUAL(b.num_dimensions(), 2);

                using size_type = Array<TwoD, double >::size_type;
                for(size_type i = 0; i < a.shape()[0]; ++i)
                {
                    for(size_type j = 0; j < a.shape()[1]; ++j)
                    {
                        BOOST_CHECK(a[i][j] == initValue);
                    }
                }

                for(size_type i = 0; i < b.shape()[0]; ++i)
                {
                    for(size_type j = 0; j < b.shape()[1]; ++j)
                    {
                        BOOST_CHECK(b[i][j] == initValue);
                    }
                }
            }
            CountedObject<double>::Check(0, 0, 151, 150, 0, 0);

        }

        BOOST_AUTO_TEST_CASE(TestPopulationFromCArray)
        {
            {
                CountedObject<double> a_array[] = { CountedObject<double>(1),
                                                    CountedObject<double>(2),
                                                    CountedObject<double>(3),
                                                    CountedObject<double>(4) };
                CountedObject<double> b_array[] = { CountedObject<double>(1),
                                                    CountedObject<double>(2),
                                                    CountedObject<double>(3),
                                                    CountedObject<double>(4),
                                                    CountedObject<double>(5)};
                CountedObject<double>::ClearCounters();

                Array<OneD, const CountedObject<double> > a(4, a_array);
                Array<OneD, CountedObject<double> > b(5, b_array);
                CountedObject<double>::Check(0, 0, 0, 9, 0, 0);

                BOOST_CHECK(a.begin() != a.end());
                BOOST_CHECK(b.begin() != b.end());

                BOOST_CHECK_EQUAL(a.size(), 4);
                BOOST_CHECK_EQUAL(b.size(), 5);

                BOOST_CHECK_EQUAL(a.num_dimensions(), 1);
                BOOST_CHECK_EQUAL(b.num_dimensions(), 1);

                using size_type = Array<OneD, double >::size_type;
                for(size_type i = 0; i < a.size(); ++i)
                {
                    BOOST_CHECK(a[i] == a_array[i]);
                }

                for(size_type i = 0; i < b.size(); ++i)
                {
                    BOOST_CHECK(b[i] == b_array[i]);
                }
            }
            CountedObject<double>::Check(0, 0, 18, 9, 0, 0);

            {
                double a_array[] = {1.0, 2.0, 3.0, 4.0};
                double b_array[] = {5.0, 6.0, 7.0, 8.0, 9.0};
                Array<OneD, double> a(4, a_array);
                Array<OneD, double> b(5, b_array);

                BOOST_CHECK(a.size() == 4);
                BOOST_CHECK(b.size() == 5);

                using size_type = Array<OneD, double >::size_type;
                for(size_type i = 0; i < a.size(); ++i)
                {
                    BOOST_CHECK(a[i] == a_array[i]);
                }

                for(size_type i = 0; i < b.size(); ++i)
                {
                    BOOST_CHECK(b[i] == b_array[i]);
                }
            }
        }

        BOOST_AUTO_TEST_CASE(TestCopyConstruction)
        {
            {
                CountedObject<double> a_array[] = { CountedObject<double>(1),
                                                    CountedObject<double>(2),
                                                    CountedObject<double>(3),
                                                    CountedObject<double>(4) };
                CountedObject<double> b_array[] = { CountedObject<double>(1),
                                                    CountedObject<double>(2),
                                                    CountedObject<double>(3),
                                                    CountedObject<double>(4),
                                                    CountedObject<double>(5)};
                CountedObject<double>::ClearCounters();

                Array<OneD, const CountedObject<double> > a(4, a_array);
                Array<OneD, CountedObject<double> > b(5, b_array);
                CountedObject<double>::Check(0, 0, 0, 9, 0, 0);

                {
                    Array<OneD, const CountedObject<double> > c(a);
                    Array<OneD, CountedObject<double> > d(b);
                    Array<OneD, const CountedObject<double> > e(a);

                    BOOST_CHECK_EQUAL(c.size(), a.size());
                    BOOST_CHECK_EQUAL(e.size(), a.size());
                    BOOST_CHECK_EQUAL(d.size(), b.size());
                    BOOST_CHECK(c.data() == a.data());
                    BOOST_CHECK(e.data() == a.data());
                    BOOST_CHECK(d.data() == b.data());

                    CountedObject<double>::Check(0, 0, 0, 9, 0, 0);
                }
                CountedObject<double>::Check(0, 0, 0, 9, 0, 0);
            }
            CountedObject<double>::Check(0, 0, 18, 9, 0, 0);

            {
                // Test offset copy constructor.
                CountedObject<double> a_array[] = { CountedObject<double>(1),
                                                    CountedObject<double>(2),
                                                    CountedObject<double>(3),
                                                    CountedObject<double>(4) };
                CountedObject<double> b_array[] = { CountedObject<double>(1),
                                                    CountedObject<double>(2),
                                                    CountedObject<double>(3),
                                                    CountedObject<double>(4),
                                                    CountedObject<double>(5)};
                CountedObject<double>::ClearCounters();
                Array<OneD, const CountedObject<double> > a(4, a_array);
                Array<OneD, CountedObject<double> > b(5, b_array);
                CountedObject<double>::Check(0, 0, 0, 9, 0, 0);

                Array<OneD, const CountedObject<double> > a_off = a + 1;
                Array<OneD, const CountedObject<double> > bb_off = b + 2;
                Array<OneD, CountedObject<double> > b_off = b + 2;

                BOOST_CHECK_EQUAL(a_off[0].value, a[1].value);
                BOOST_CHECK_EQUAL(a_off[1].value, a[2].value);
                BOOST_CHECK_EQUAL(a_off[2].value, a[3].value);
                BOOST_CHECK_EQUAL(a_off.size(), 3);

                BOOST_CHECK_EQUAL(b_off[0].value, b[2].value);
                BOOST_CHECK_EQUAL(b_off[1].value, b[3].value);
                BOOST_CHECK_EQUAL(b_off[2].value, b[4].value);
                BOOST_CHECK_EQUAL(b_off.size(), 3);


                BOOST_CHECK_EQUAL(bb_off[0].value, b[2].value);
                BOOST_CHECK_EQUAL(bb_off[1].value, b[3].value);
                BOOST_CHECK_EQUAL(bb_off[2].value, b[4].value);
                BOOST_CHECK_EQUAL(bb_off.size(), 3);

            }
        }

        BOOST_AUTO_TEST_CASE(Test1DAssignmentOperator)
        {
            {
                // 1D
                CountedObject<double> a_array[] = { CountedObject<double>(1),
                                                    CountedObject<double>(2),
                                                    CountedObject<double>(3),
                                                    CountedObject<double>(4) };
                CountedObject<double> b_array[] = { CountedObject<double>(1),
                                                    CountedObject<double>(2),
                                                    CountedObject<double>(3),
                                                    CountedObject<double>(4),
                                                    CountedObject<double>(5)};
                CountedObject<double>::ClearCounters();
                Array<OneD, const CountedObject<double> > a(4, a_array);
                Array<OneD, CountedObject<double> > b(5, b_array);
                CountedObject<double>::Check(0, 0, 0, 9, 0, 0);

                {
                    Array<OneD, const CountedObject<double> > lhs_a;
                    Array<OneD, CountedObject<double> > lhs_b1;
                    Array<OneD, const CountedObject<double> > lhs_b2;
                    CountedObject<double>::Check(0, 0, 0, 9, 0, 0);

                    BOOST_CHECK_EQUAL(lhs_a.size(), 0);
                    BOOST_CHECK_EQUAL(lhs_b1.size(), 0);
                    BOOST_CHECK_EQUAL(lhs_b2.size(), 0);

                    {
                        lhs_a = a;
                        lhs_b1 = b;
                        lhs_b2 = b;
                        CountedObject<double>::Check(0, 0, 0, 9, 0, 0);

                        BOOST_CHECK_EQUAL(lhs_a.size(), a.size());
                        BOOST_CHECK_EQUAL(lhs_b1.size(), b.size());
                        BOOST_CHECK_EQUAL(lhs_b2.size(), b.size());
                        BOOST_CHECK_EQUAL(lhs_a.data(), a.data());
                        BOOST_CHECK_EQUAL(lhs_b1.data(), b.data());
                        BOOST_CHECK_EQUAL(lhs_b2.data(), b.data());

                        using size_type = Array<OneD, double >::size_type;
                        for(size_type i = 0; i < a.size(); ++i)
                        {
                            BOOST_CHECK(lhs_a[i] == a[i]);
                        }

                            for(size_type i = 0; i < b.size(); ++i)
                            {
                                BOOST_CHECK(lhs_b1[i] == b[i]);
                                BOOST_CHECK(lhs_b2[i] == b[i]);
                            }
                            CountedObject<double>::Check(0, 0, 0, 9, 0, 0);
                        }
                        CountedObject<double>::Check(0, 0, 0, 9, 0, 0);
                    }
                    CountedObject<double>::Check(0, 0, 0, 9, 0, 0);
            }
        }

        BOOST_AUTO_TEST_CASE(Test2DAssignmentOperator)
        {
            {
                // 2D
                CountedObject<double> initValue(98);
                CountedObject<double>::ClearCounters();

                Array<TwoD, const CountedObject<double> > a(5, 10, initValue);
                Array<TwoD, CountedObject<double> > b(10, 10, initValue);
                CountedObject<double>::Check(0, 0, 0, 150, 0, 0);

                BOOST_CHECK(a.begin() != a.end());
                BOOST_CHECK(b.begin() != b.end());

                BOOST_CHECK_EQUAL(a.size(), 50);
                BOOST_CHECK_EQUAL(b.size(), 100);
                BOOST_CHECK_EQUAL(a.shape()[0], 5);
                BOOST_CHECK_EQUAL(a.shape()[1], 10);
                BOOST_CHECK_EQUAL(b.shape()[0], 10);
                BOOST_CHECK_EQUAL(b.shape()[1], 10);

                BOOST_CHECK_EQUAL(a.num_dimensions(), 2);
                BOOST_CHECK_EQUAL(b.num_dimensions(), 2);

                {
                    Array<TwoD, const CountedObject<double> > lhs_a;
                    Array<TwoD, const CountedObject<double> > lhs_b1;
                    Array<TwoD, CountedObject<double> > lhs_b2;

                    BOOST_CHECK_EQUAL(lhs_a.size(), 0);
                    BOOST_CHECK_EQUAL(lhs_b1.size(), 0);
                    BOOST_CHECK_EQUAL(lhs_b2.size(), 0);

                    lhs_a = a;
                    lhs_b1 = b;
                    lhs_b2 = b;
                    CountedObject<double>::Check(0, 0, 0, 150, 0, 0);

                    BOOST_CHECK_EQUAL(lhs_a.size(), a.size());
                    BOOST_CHECK_EQUAL(lhs_b1.size(), b.size());
                    BOOST_CHECK_EQUAL(lhs_b2.size(), b.size());
                    BOOST_CHECK_EQUAL(lhs_a.shape()[0], a.shape()[0]);
                    BOOST_CHECK_EQUAL(lhs_a.shape()[1], a.shape()[1]);

                    BOOST_CHECK_EQUAL(lhs_b1.shape()[0], b.shape()[0]);
                    BOOST_CHECK_EQUAL(lhs_b1.shape()[1], b.shape()[1]);

                    BOOST_CHECK_EQUAL(lhs_b2.shape()[0], b.shape()[0]);
                    BOOST_CHECK_EQUAL(lhs_b2.shape()[1], b.shape()[1]);

                    using size_type = Array<TwoD, double >::size_type;
                    for(size_type i = 0; i < a.shape()[0]; ++i)
                    {
                        for(size_type j = 0; j < a.shape()[1]; ++j)
                        {
                            BOOST_CHECK(a[i][j] == initValue);
                            BOOST_CHECK(lhs_a[i][j] == a[i][j]);
                        }
                    }

                    for(size_type i = 0; i < b.shape()[0]; ++i)
                    {
                        for(size_type j = 0; j < b.shape()[1]; ++j)
                        {
                            BOOST_CHECK(b[i][j] == initValue);
                            BOOST_CHECK(lhs_b1[i][j] == b[i][j]);
                            BOOST_CHECK(lhs_b2[i][j] == b[i][j]);
                        }
                    }
                }
                CountedObject<double>::Check(0, 0, 0, 150, 0, 0);
            }

            CountedObject<double>::Check(0, 0, 151, 150, 0, 0);
        }

        BOOST_AUTO_TEST_CASE(TestOffsetAssignmentOperator)
        {
            double a[] = {1.0, 2.0, 3.0, 4.0, 5.0};
            double b[] = {10.0, 20.0, 30.0, 40.0 };

            Array<OneD, NekDouble> rhs_a(5, a);
            Array<OneD, const NekDouble> rhs_b(4, b);

            Array<OneD, NekDouble> offset_a = rhs_a + 2;;
            Array<OneD, const NekDouble> offset_b = rhs_a + 2;

            Array<OneD, NekDouble> lhs_a;
            Array<OneD, const NekDouble> lhs_b;

            lhs_a = offset_a;
            lhs_b = offset_b;

            BOOST_CHECK_EQUAL(lhs_a.size(), offset_a.size());
            BOOST_CHECK_EQUAL(lhs_b.size(), offset_b.size());

            using size_type = Array<TwoD, double >::size_type;
            for (size_type i = 0; i < lhs_a.size(); ++i)
            {
                BOOST_CHECK(lhs_a[i] == offset_a[i]);
            }

            for (size_type i = 0; i < lhs_b.size(); ++i)
            {
                BOOST_CHECK(lhs_b[i] == offset_b[i]);
            }
        }

        BOOST_AUTO_TEST_CASE(Test1DAccessOperator)
        {
            UnitTests::RedirectCerrIfNeeded();
            // Normal access.
            double a[] = {1.0, 2.0, 3.0, 4.0, 5.0};
            double b[] = {10.0, 20.0, 30.0, 40.0 };

            Array<OneD, NekDouble> rhs_a(5, a);
            Array<OneD, const NekDouble> rhs_b(4, b);

            BOOST_CHECK_EQUAL(rhs_a[0], 1.0);
            BOOST_CHECK_EQUAL(rhs_a[1], 2.0);
            BOOST_CHECK_EQUAL(rhs_a[2], 3.0);
            BOOST_CHECK_EQUAL(rhs_a[3], 4.0);
            BOOST_CHECK_EQUAL(rhs_a[4], 5.0);

            BOOST_CHECK_EQUAL(rhs_b[0], 10.0);
            BOOST_CHECK_EQUAL(rhs_b[1], 20.0);
            BOOST_CHECK_EQUAL(rhs_b[2], 30.0);
            BOOST_CHECK_EQUAL(rhs_b[3], 40.0);
            #if defined(NEKTAR_DEBUG) || defined(NEKTAR_FULL_DEBUG)
                BOOST_CHECK_NO_THROW(rhs_a[4]);
                BOOST_CHECK_THROW(rhs_a[5], ErrorUtil::NekError);
                BOOST_CHECK_NO_THROW(rhs_b[3]);
                BOOST_CHECK_THROW(rhs_b[4], ErrorUtil::NekError);
            #endif
        }

        BOOST_AUTO_TEST_CASE(Test2DAccessOperator)
        {
            {
                NekDouble a_vals[] = {1.0, 2.0,
                                    3.0, 4.0,
                                    5.0, 6.0,
                                    7.0, 8.0,
                                    9.0, 10.0};
                Array<TwoD, NekDouble> a(5, 2, a_vals);

                BOOST_CHECK_EQUAL(a.GetRows(), 5);
                BOOST_CHECK_EQUAL(a.GetColumns(), 2);

                using size_type = Array<TwoD, double >::size_type;
                for(size_type i = 0; i < a.GetRows(); ++i)
                {
                    for(size_type j = 0; j < a.GetColumns(); ++j)
                    {
                        BOOST_CHECK_EQUAL(a[i][j], a_vals[i*a.GetColumns() + j]);
                    }
                }

                // Now test assignment.
                NekDouble rhs_vals[] = {18.0, -76.2,
                                        45.2, 1352.246,
                                        -46.346, -2463.346,
                                        26.347, 1.0,
                                        0.0, 23465.3};
                for(size_type i = 0; i < a.GetRows(); ++i)
                {
                    for(size_type j = 0; j < a.GetColumns(); ++j)
                    {
                        a[i][j] = rhs_vals[i*a.GetColumns() + j];
                    }
                }

                for(size_type i = 0; i < a.GetRows(); ++i)
                {
                    for(size_type j = 0; j < a.GetColumns(); ++j)
                    {
                        BOOST_CHECK_EQUAL(a[i][j], rhs_vals[i*a.GetColumns() + j]);
                    }
                }
            }

            {
                // Test 1 Column.
                NekDouble a_vals[] = {1.0,
                                    3.0,
                                    5.0,
                                    7.0,
                                    9.0};
                Array<TwoD, NekDouble> a(5, 1, a_vals);
                BOOST_CHECK_EQUAL(a.GetRows(), 5);
                BOOST_CHECK_EQUAL(a.GetColumns(), 1);

                using size_type = Array<TwoD, double >::size_type;
                for(size_type i = 0; i < a.GetRows(); ++i)
                {
                    BOOST_CHECK_EQUAL(a[i][0], a_vals[i]);
                }
            }

            {
                Array<TwoD, NekDouble> a(1,1);
                a[0][0] = 17.23;
                BOOST_CHECK_EQUAL(a[0][0], 17.23);
            }
        }

        BOOST_AUTO_TEST_CASE(TestSharedPtr)
        {
            std::shared_ptr<double> a(new double[10]);
            std::shared_ptr<const double> b(a);
            std::shared_ptr<const double> c;
            c = a;

        }

        BOOST_AUTO_TEST_CASE(Test1DEqualOperator)
        {
            {
                // integer types
                using int_type = int;

                int_type a[] = {1, 2, 3, 4, 5};

                Array<OneD, int_type> nek_a(5, a);
                Array<OneD, int_type> nek_b(5, a);
                Array<OneD, const int_type> nek_c(5, a);

                // still equal
                BOOST_CHECK(nek_a == nek_b);
                // equal and only one const
                BOOST_CHECK(nek_a == nek_c);

                // not equal
                nek_b[3] += 1;
                BOOST_CHECK(!(nek_a == nek_b));
            }

            {
                // float types
                using float_type = double;

                float_type a[] = {1.0, 2.0, 3.0, 4.0, 5.0};

                Array<OneD, float_type> nek_a(5, a);
                Array<OneD, float_type> nek_b(5, a);
                Array<OneD, const float_type> nek_c(5, a);

                // still equal
                nek_b[3] *= (1.0 + std::numeric_limits<float_type>::epsilon());
                BOOST_CHECK(nek_a == nek_b);
                 // equal and only one const
                BOOST_CHECK(nek_a == nek_c);

                // not equal
                nek_b[4] *=
                    (1.0 + 100 * std::numeric_limits<float_type>::epsilon());
                BOOST_CHECK(!(nek_a == nek_b));

            }

            // This should not compile
            // {
            //     // mixed types should not compare
            //     using float_type = double;
            //     using int_type = int;

            //     float_type a[] = {1.0, 2.0, 3.0, 4.0, 5.0};
            //     int_type b[] = {1, 2, 3, 4, 5};

            //     Array<OneD, float_type> nek_a(5, a);
            //     Array<OneD, int_type> nek_b(5, b);

            //     // still equal
            //     BOOST_CHECK(nek_a == nek_b);
            // }
        }

        BOOST_AUTO_TEST_CASE(Test2DEqualOperator)
        {
            {
                // integer types
                using int_type = int;

                int_type a[] = {1, 2, 3, 4, 5};
                Array<TwoD, int_type> nek_a(5, 1, a);
                Array<TwoD, int_type> nek_b(5, 1, a);
                Array<TwoD, const int_type> nek_c(5, 1, a);

                // still equal
                BOOST_CHECK(nek_a == nek_b);
                // equal and only one const
                BOOST_CHECK(nek_a == nek_c);

                // not equal
                nek_b[3][0] += 1;
                BOOST_CHECK(!(nek_a == nek_b));
            }

            {
                // floating types
                using float_type = double;

                float_type a[] = {1., 2., 3., 4., 5.};
                Array<TwoD, float_type> nek_a(5, 1, a);
                Array<TwoD, float_type> nek_b(5, 1, a);
                Array<TwoD, const float_type> nek_c(5, 1, a);

                // still equal
                nek_b[3][0] *=
                    (1.0 + std::numeric_limits<float_type>::epsilon());
                BOOST_CHECK(nek_a == nek_b);
                 // equal and only one const
                BOOST_CHECK(nek_a == nek_c);

                // not equal
                nek_b[4][0] *=
                    (1.0 + 100 * std::numeric_limits<float_type>::epsilon());
                BOOST_CHECK(!(nek_a == nek_b));
            }
        }



    } // End SharedArrayUnitTests
} // End Nektar


