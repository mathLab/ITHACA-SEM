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
#include <LibUtilities/LinearAlgebra/NekVector.hpp>
#include <LibUtilities/LinearAlgebra/NekVectorVariableSized.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>
#include <iostream>

namespace Nektar
{
    namespace SharedArrayUnitTests
    {
        

        void testNewOffset()
        {
        }
        
        void testConstantResultType()
        {
        }
        
        
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
                void getConstByReference(ConstArray<OneD, NekDouble>& out)
                {
                    out = a;
                }
                
                // Does nothing - the assignment goes into a temporary and is 
                // immediately lost.
                void getConstByValue(ConstArray<OneD, NekDouble> out)
                {
                    out = a;
                }
                
                Array<OneD, NekDouble> a;
                ConstArray<OneD, NekDouble> b;
        };
        
        void testParameterPopulation()
        {
            ParameterTestClass obj;
            
            Array<OneD, NekDouble> zero(5, 0.0);
            Array<OneD, NekDouble> temp(zero);
            ConstArray<OneD, NekDouble> const_zero(5, 0.0);
            ConstArray<OneD, NekDouble> const_temp(const_zero);
            
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
        
        void testEmptyConstructor()
        {
            {
                CountedObject<double>::ClearCounters();
                ConstArray<OneD, CountedObject<double> > a;
                Array<OneD, CountedObject<double> > b;
                CountedObject<double>::check(0, 0, 0, 0, 0, 0);
                
                BOOST_CHECK(a.begin() == a.end());
                BOOST_CHECK(b.begin() == b.end());
                
                BOOST_CHECK_EQUAL(a.size(), 0);
                BOOST_CHECK_EQUAL(b.size(), 0);
                
                BOOST_CHECK_EQUAL(a.num_dimensions(), 1);
                BOOST_CHECK_EQUAL(b.num_dimensions(), 1);
            }
            CountedObject<double>::check(0, 0, 0, 0, 0, 0);
            
            {
                CountedObject<double>::ClearCounters();
                ConstArray<TwoD, CountedObject<double> > a;
                Array<TwoD, CountedObject<double> > b;
                CountedObject<double>::check(0, 0, 0, 0, 0, 0);
                
                BOOST_CHECK(a.begin() == a.end());
                BOOST_CHECK(b.begin() == b.end());
                
                BOOST_CHECK_EQUAL(a.size(), 0);
                BOOST_CHECK_EQUAL(b.size(), 0);
                
                BOOST_CHECK_EQUAL(a.num_dimensions(), 2);
                BOOST_CHECK_EQUAL(b.num_dimensions(), 2);
            }
            CountedObject<double>::check(0, 0, 0, 0, 0, 0);
        }
        
        void testUninitializedConstructor()
        {
            {
                CountedObject<double>::ClearCounters();
                
                ConstArray<OneD, CountedObject<double> > a(5);
                Array<OneD, CountedObject<double> > b(10);
                CountedObject<double>::check(15, 0, 0, 0, 0, 0);
                
                BOOST_CHECK(a.begin() != a.end());
                BOOST_CHECK(b.begin() != b.end());
                
                BOOST_CHECK_EQUAL(a.size(), 5);
                BOOST_CHECK_EQUAL(b.size(), 10);
                
                BOOST_CHECK_EQUAL(a.num_dimensions(), 1);
                BOOST_CHECK_EQUAL(b.num_dimensions(), 1);
            }
            CountedObject<double>::check(15, 0, 15, 0, 0, 0);
            
            {
                CountedObject<double>::ClearCounters();
                
                ConstArray<TwoD, CountedObject<double> > a(5, 10);
                Array<TwoD, CountedObject<double> > b(10, 10);
                CountedObject<double>::check(150, 0, 0, 0, 0, 0);
                
                BOOST_CHECK(a.begin() != a.end());
                BOOST_CHECK(b.begin() != b.end());
                
                BOOST_CHECK_EQUAL(a.size(), 5);
                BOOST_CHECK_EQUAL(b.size(), 10);
                BOOST_CHECK_EQUAL(a.shape()[0], 5);
                BOOST_CHECK_EQUAL(a.shape()[1], 10);
                BOOST_CHECK_EQUAL(b.shape()[0], 10);
                BOOST_CHECK_EQUAL(b.shape()[1], 10);
                
                BOOST_CHECK_EQUAL(a.num_dimensions(), 2);
                BOOST_CHECK_EQUAL(b.num_dimensions(), 2);
            }
            CountedObject<double>::check(150, 0, 150, 0, 0, 0);
            
            {
                CountedObject<double>::ClearCounters();
                
                ConstArray<ThreeD, CountedObject<double> > a(1, 2, 3);
                Array<ThreeD, CountedObject<double> > b(4, 5, 6);
                CountedObject<double>::check(126, 0, 0, 0, 0, 0);
                
                BOOST_CHECK(a.begin() != a.end());
                BOOST_CHECK(b.begin() != b.end());
                
                BOOST_CHECK_EQUAL(a.size(), 1);
                BOOST_CHECK_EQUAL(b.size(), 4);
                BOOST_CHECK_EQUAL(a.shape()[0], 1);
                BOOST_CHECK_EQUAL(a.shape()[1], 2);
                BOOST_CHECK_EQUAL(a.shape()[2], 3);
                BOOST_CHECK_EQUAL(b.shape()[0], 4);
                BOOST_CHECK_EQUAL(b.shape()[1], 5);
                BOOST_CHECK_EQUAL(b.shape()[2], 6);
               
                BOOST_CHECK_EQUAL(a.num_dimensions(), 3);
                BOOST_CHECK_EQUAL(b.num_dimensions(), 3);
            }
            CountedObject<double>::check(126, 0, 126, 0, 0, 0);
        }
        
        void testSingleValueInitialization()
        {
            {
                CountedObject<double> initValue(98);
                CountedObject<double>::ClearCounters();
                
                ConstArray<OneD, CountedObject<double> > a(5, initValue);
                Array<OneD, CountedObject<double> > b(10, initValue);
                CountedObject<double>::check(0, 0, 0, 15, 0, 0);
                
                BOOST_CHECK(a.begin() != a.end());
                BOOST_CHECK(b.begin() != b.end());
                
                BOOST_CHECK_EQUAL(a.size(), 5);
                BOOST_CHECK_EQUAL(b.size(), 10);
                
                BOOST_CHECK_EQUAL(a.num_dimensions(), 1);
                BOOST_CHECK_EQUAL(b.num_dimensions(), 1);
                
                for(ConstArray<OneD, CountedObject<double> >::const_iterator iter = a.begin(); iter != a.end(); ++iter)
                {
                    BOOST_CHECK(*iter == initValue);
                }
                    
                for(Array<OneD, CountedObject<double> >::iterator iter = b.begin(); iter != b.end(); ++iter)
                {
                    BOOST_CHECK(*iter == initValue);
                }
            }
            CountedObject<double>::check(0, 0, 16, 15, 0, 0);
            
            {
                CountedObject<double> initValue(98);
                CountedObject<double>::ClearCounters();
                
                ConstArray<TwoD, CountedObject<double> > a(5, 10, initValue);
                Array<TwoD, CountedObject<double> > b(10, 10, initValue);
                CountedObject<double>::check(0, 0, 0, 150, 0, 0);
                
                BOOST_CHECK(a.begin() != a.end());
                BOOST_CHECK(b.begin() != b.end());
                
                BOOST_CHECK_EQUAL(a.size(), 5);
                BOOST_CHECK_EQUAL(b.size(), 10);
                BOOST_CHECK_EQUAL(a.shape()[0], 5);
                BOOST_CHECK_EQUAL(a.shape()[1], 10);
                BOOST_CHECK_EQUAL(b.shape()[0], 10);
                BOOST_CHECK_EQUAL(b.shape()[1], 10);
                
                BOOST_CHECK_EQUAL(a.num_dimensions(), 2);
                BOOST_CHECK_EQUAL(b.num_dimensions(), 2);
                    
                for(unsigned int i = 0; i < a.shape()[0]; ++i)
                {
                    for(unsigned int j = 0; j < a.shape()[1]; ++j)
                    {
                        BOOST_CHECK(a[i][j] == initValue);
                    }
                }
                
                for(unsigned int i = 0; i < b.shape()[0]; ++i)
                {
                    for(unsigned int j = 0; j < b.shape()[1]; ++j)
                    {
                        BOOST_CHECK(b[i][j] == initValue);
                    }
                }
            }
            CountedObject<double>::check(0, 0, 151, 150, 0, 0);
            
            {
                CountedObject<double> initValue(98);
                CountedObject<double>::ClearCounters();
                
                ConstArray<ThreeD, CountedObject<double> > a(1, 2, 3, initValue);
                Array<ThreeD, CountedObject<double> > b(4, 5, 6, initValue);
                CountedObject<double>::check(0, 0, 0, 126, 0, 0);
                
                BOOST_CHECK(a.begin() != a.end());
                BOOST_CHECK(b.begin() != b.end());
                
                BOOST_CHECK_EQUAL(a.size(), 1);
                BOOST_CHECK_EQUAL(b.size(), 4);
                BOOST_CHECK_EQUAL(a.shape()[0], 1);
                BOOST_CHECK_EQUAL(a.shape()[1], 2);
                BOOST_CHECK_EQUAL(a.shape()[2], 3);
                BOOST_CHECK_EQUAL(b.shape()[0], 4);
                BOOST_CHECK_EQUAL(b.shape()[1], 5);
                BOOST_CHECK_EQUAL(b.shape()[2], 6);
               
                BOOST_CHECK_EQUAL(a.num_dimensions(), 3);
                BOOST_CHECK_EQUAL(b.num_dimensions(), 3);
                    
                for(unsigned int i = 0; i < a.shape()[0]; ++i)
                {
                    for(unsigned int j = 0; j < a.shape()[1]; ++j)
                    {
                        for(unsigned int k = 0; k < a.shape()[2]; ++k)
                        {
                            BOOST_CHECK(a[i][j][k] == initValue);
                        }
                    }
                }
                
                for(unsigned int i = 0; i < b.shape()[0]; ++i)
                {
                    for(unsigned int j = 0; j < b.shape()[1]; ++j)
                    {
                        for(unsigned int k = 0; k < b.shape()[2]; ++k)
                        {
                            BOOST_CHECK(b[i][j][k] == initValue);
                        }
                    }
                }
            }
            CountedObject<double>::check(0, 0, 127, 126, 0, 0);
        }
        
        void testPopulationFromCArray()
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
            
                ConstArray<OneD, CountedObject<double> > a(4, a_array);
                Array<OneD, CountedObject<double> > b(5, b_array);
                CountedObject<double>::check(0, 0, 0, 9, 0, 0);
            
                BOOST_CHECK(a.begin() != a.end());
                BOOST_CHECK(b.begin() != b.end());
            
                BOOST_CHECK_EQUAL(a.size(), 4);
                BOOST_CHECK_EQUAL(b.size(), 5);
            
                BOOST_CHECK_EQUAL(a.num_dimensions(), 1);
                BOOST_CHECK_EQUAL(b.num_dimensions(), 1);
            
                for(unsigned int i = 0; i < a.size(); ++i)
                {
                    BOOST_CHECK(a[i] == a_array[i]);
                }
                
                for(unsigned int i = 0; i < b.size(); ++i)
                {
                    BOOST_CHECK(b[i] == b_array[i]);
                }
            }
            CountedObject<double>::check(0, 0, 18, 9, 0, 0);

            {
                double a_array[] = {1.0, 2.0, 3.0, 4.0};
                double b_array[] = {5.0, 6.0, 7.0, 8.0, 9.0};
                Array<OneD, double> a(4, a_array);
                Array<OneD, double> b(5, b_array);

                BOOST_CHECK(a.size() == 4);
                BOOST_CHECK(b.size() == 5);
                for(unsigned int i = 0; i < a.size(); ++i)
                {
                    BOOST_CHECK(a[i] == a_array[i]);
                }
                
                for(unsigned int i = 0; i < b.size(); ++i)
                {
                    BOOST_CHECK(b[i] == b_array[i]);
                }
            }
        }
    }
}

#endif //NEKTAR_UNIT_TESTS_SHARED_ARRAY_HPP
