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
// Description: Test code for Nektar::shared_array
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_UNIT_TESTS_SHARED_ARRAY_HPP
#define NEKTAR_UNIT_TESTS_SHARED_ARRAY_HPP

#include <UnitTests/testNekSharedArray.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <UnitTests/CountedObject.h>

#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>

namespace Nektar
{
    namespace SharedArrayUnitTests
    {
        SharedArray<double> testFunc1()
        {
            SharedArray<double> result(new double[10], 10);
            return result;
        }

        class TestClass
        {
            public:
                TestClass() : result(new double[10], 10) {}

                SharedArray<double> foo() { return result; }
                const SharedArray<double>& bar() const { return result; }

            private:
                SharedArray<double> result;
        };

        // Nothing really runs here - it just checks valid compilation.  The commented out 
        // lines should not compile.
        void testGet()
        {
            Nektar::SharedArray<double> d1(new double[10], 10);
            const Nektar::SharedArray<double> d2(new double[10], 10);
            const Nektar::SharedArray<double>& d3 = d1;

            double* a1 = d1.get();
            const double* a2 = d1.get();

            //double* a3 = d2.get();
            const double* a4 = d2.get();

            //double* a5 = d3.get();
            const double* a6 = d3.get();

            //Nektar::SharedArray<double> d4(d2);
            //d4[0] = 2.7;

            const Nektar::SharedArray<double> d5(d1);
            //d5 = d3;
            //d5[0] = 3.4;

            Nektar::SharedArray<double> d6;
            d6 = d1;

            Nektar::SharedArray<double> d7 = testFunc1();

            TestClass t;
            Nektar::SharedArray<double> d8 = t.foo();
            //Nektar::SharedArray<double> d9 = t.bar();
            const Nektar::SharedArray<double>& d9 = t.bar();
        }

        void testAccessOperator()
        {
            Nektar::SharedArray<double> d1(new double[10], 10);
            const Nektar::SharedArray<double> d2(new double[10], 10);
            const Nektar::SharedArray<double>& d3 = d1;

            double& a1 = d1[0];
            const double& a2 = d1[0]; 

            //double& a3 = d2[0];
            const double& a4 = d2[0];

            //double& a5 = d3[0];
            const double& a6 = d3[0];
        }

        void testOffset()
        {
            Nektar::SharedArray<double> d1(new double[10], 10);
            for(unsigned int i = 0; i < 10; ++i)
            {
                d1[i] = (double)i + .7;
            }

            Nektar::SharedArray<double> d2 = d1 + 5;
            BOOST_CHECK_EQUAL(d1[5], d2[0]);


            const Nektar::SharedArray<double>& d3 = d2;
            //Nektar::SharedArray<double> d4 = d3 + 2;
        }

        class TestReference
        {
            public:
                TestReference() : obj(new double[10], 10) {}
                const SharedArray<double>& foo() const { return obj; }

                SharedArray<double> obj;
        };
        
        void testConstruction()
        {
            {
                CountedObject<double>::ClearCounters();

                SharedArray<CountedObject<double> > a1 = 
                    MemoryManager::AllocateSharedArray<CountedObject<double> >(10);
                SharedArray<const CountedObject<double> > a2 = 
                    MemoryManager::AllocateSharedArray<const CountedObject<double> >(10);
                const SharedArray<CountedObject<double> > a3 = 
                    MemoryManager::AllocateSharedArray<CountedObject<double> >(10);
                const SharedArray<const CountedObject<double> > a4 = 
                    MemoryManager::AllocateSharedArray<const CountedObject<double> >(10);

                {
                    // All should compile.
                    SharedArray<const CountedObject<double> > b1(a1);
                    SharedArray<const CountedObject<double> > b2(a2);
                    SharedArray<const CountedObject<double> > b3(a3);
                    SharedArray<const CountedObject<double> > b4(a4);
                }

                {
                    // All should compile.
                    const SharedArray<const CountedObject<double> > b1(a1);
                    const SharedArray<const CountedObject<double> > b2(a2);
                    const SharedArray<const CountedObject<double> > b3(a3);
                    const SharedArray<const CountedObject<double> > b4(a4);
                }

                {
                    // Commented out should not compile
                    SharedArray<CountedObject<double> > b1(a1);
                    //SharedArray<CountedObject<double> > b2(a2);
                    //SharedArray<CountedObject<double> > b3(a3);
                    //SharedArray<CountedObject<double> > b4(a4);
                }

                {
                    // Commented out should not compile
                    const SharedArray<CountedObject<double> > b1(a1);
                    //const SharedArray<CountedObject<double> > b2(a2);
                    //const SharedArray<CountedObject<double> > b3(a3);
                    //const SharedArray<CountedObject<double> > b4(a4);
                }
                BOOST_CHECK_EQUAL(CountedObject<double>::numberDestroyed, 0);
            }
            BOOST_CHECK_EQUAL(CountedObject<double>::numberDestroyed, 40);

            {
                // TODO - Find a way to prohibit data not on the heap.
                double* data1 = new double[10];
                double* data2 = new double[10];
                for(int i = 0; i < 10; ++i)
                {
                    data1[i] = i+1;
                    data2[i] = i+1;
                }
                
                SharedArray<double> a1(data1, 10);
                SharedArray<const double> a2(data2, 10);
                
                {
                    SharedArray<const double> b1(a1);
                    a1[0] = 17.0;
                    BOOST_CHECK_EQUAL(b1[0], 17.0);
                    BOOST_CHECK_EQUAL(b1[0], a1[0]);

                    // Should not compile.
                    //b1[0] = 15;
                    //b1.get()[1] = 15;
                }

                {
                    SharedArray<double> b1(a1);
                    b1[0] = 23.0;
                    BOOST_CHECK_EQUAL(b1[0], 23.0);
                    BOOST_CHECK_EQUAL(b1[0], a1[0]);
                }
            }
        }

        void testAssignment()
        {
            {
                CountedObject<double>::ClearCounters();

                SharedArray<CountedObject<double> > a1 = 
                    MemoryManager::AllocateSharedArray<CountedObject<double> >(10);
                SharedArray<const CountedObject<double> > a2 = 
                    MemoryManager::AllocateSharedArray<const CountedObject<double> >(10);
                const SharedArray<CountedObject<double> > a3 = 
                    MemoryManager::AllocateSharedArray<CountedObject<double> >(10);
                const SharedArray<const CountedObject<double> > a4 = 
                    MemoryManager::AllocateSharedArray<const CountedObject<double> >(10);

                {
                    // All should compile.
                    SharedArray<const CountedObject<double> > b1;
                    SharedArray<const CountedObject<double> > b2;
                    SharedArray<const CountedObject<double> > b3;
                    SharedArray<const CountedObject<double> > b4;

                    b1 = a1;
                    b2 = a2;
                    b3 = a3;
                    b4 = a4;
                }

                {
                    // None should compile, we can't assign into a constant 
                    // Shared Array
                    const SharedArray<const CountedObject<double> > b1;
                    const SharedArray<const CountedObject<double> > b2;
                    const SharedArray<const CountedObject<double> > b3;
                    const SharedArray<const CountedObject<double> > b4;

                    //b1 = a1;
                    //b2 = a2;
                    //b3 = a3;
                    //b4 = a4;
                }

                {
                    // Commented out should not compile
                    SharedArray<CountedObject<double> > b1;
                    SharedArray<CountedObject<double> > b2;
                    SharedArray<CountedObject<double> > b3;
                    SharedArray<CountedObject<double> > b4;

                    b1 = a1;
                    //b2 = a2;
                    //b3 = a3;
                    //b4 = a4;
                }

                {
                    // None should compile, we can't assign into a constant 
                    // Shared Array
                    const SharedArray<CountedObject<double> > b1;
                    const SharedArray<CountedObject<double> > b2;
                    const SharedArray<CountedObject<double> > b3;
                    const SharedArray<CountedObject<double> > b4;

                    //b1 = a1;
                    //b2 = a2;
                    //b3 = a3;
                    //b4 = a4;
                }
                BOOST_CHECK_EQUAL(CountedObject<double>::numberDestroyed, 0);
            }
            BOOST_CHECK_EQUAL(CountedObject<double>::numberDestroyed, 40);
        }

        void foo(SharedArray<double>& f)
        {
            SharedArray<double> d = MemoryManager::AllocateSharedArray<double>(10);
            d[0] = 7.1;
            d[1] = 8.2;

            f = d;
        }

        void testFunctionCall()
        {
            SharedArray<double> a;
            foo(a);

            BOOST_CHECK_EQUAL(a[0], 7.1);
            BOOST_CHECK_EQUAL(a[1], 8.2);
        }
    }
}

#endif //NEKTAR_UNIT_TESTS_SHARED_ARRAY_HPP
