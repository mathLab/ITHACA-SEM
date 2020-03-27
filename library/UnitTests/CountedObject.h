///////////////////////////////////////////////////////////////////////////////
//
// File: CountedObject.h
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

#ifndef NEKTAR_UNIT_TESTS_COUNTED_OBJECT_H
#define NEKTAR_UNIT_TESTS_COUNTED_OBJECT_H


#include <boost/test/unit_test.hpp>
#include <iostream>

namespace Nektar
{
    
    template<typename DerivedType>
    class CountedObject
    {
        public:
            CountedObject() :
                value(0)
            {
                ++numberDefaultConstructed;
            }

            explicit CountedObject(unsigned int v) :
                value(v)
            {
                ++numberConstructedFromInt;
                ++numberOf1ParameterConstructions;
            }

            CountedObject(unsigned int a1, unsigned int a2) :
                value(a1 + a2)
            {
                ++numberOf2ParameterConstructions;
            }

            CountedObject(unsigned int a1, unsigned int a2, unsigned int a3) :
                value(a1 + a2 + a3)
            {
                ++numberOf3ParameterConstructions;
            }

            CountedObject(const CountedObject<DerivedType>& rhs) :
                value(rhs.value)
            {
                ++numberCopied;
            }

            virtual ~CountedObject()
            {
                ++numberDestroyed;
            }

            CountedObject<DerivedType>& operator=(const CountedObject<DerivedType>& rhs)
            {
                ++numberAssigned;
                value = rhs.value;
                return *this;
            }

            bool operator==(const CountedObject<DerivedType>& rhs) const
            {
                return value == rhs.value;
            }

            bool operator!=(const CountedObject<DerivedType>& rhs) const
            {
                return !(*this == rhs);
            }

            CountedObject<DerivedType>* clone() const
            {
                ++numberCloned;
                return new CountedObject<DerivedType>(*this);
            }

            static void ClearCounters()
            {
                numberDefaultConstructed = 0;
                numberConstructedFromInt = 0;
                numberDestroyed = 0;
                numberCopied = 0;
                numberAssigned = 0;
                numberCloned = 0;
                numberOf1ParameterConstructions = 0;
                numberOf2ParameterConstructions = 0;
                numberOf3ParameterConstructions = 0;
                numberOfExpressionConstructions = 0;
                numberOfExpressionAssignments = 0;
            }

            static void Check(unsigned int expectedDefaultConstructed, unsigned int expectedConstructedFromInt,
                                unsigned int expectedDestroyed, unsigned int expectedCopied, unsigned int expectedCloned,
                                unsigned int expectedAssigned)
            {
                BOOST_CHECK_EQUAL(numberDefaultConstructed, expectedDefaultConstructed);
                BOOST_CHECK_EQUAL(numberConstructedFromInt, expectedConstructedFromInt);
                BOOST_CHECK_EQUAL(numberDestroyed, expectedDestroyed);
                BOOST_CHECK_EQUAL(numberCopied, expectedCopied);
                BOOST_CHECK_EQUAL(numberAssigned, expectedAssigned);
                BOOST_CHECK_EQUAL(numberCloned, expectedCloned);
            }

            unsigned int GetValue() const { return value; }
            operator unsigned int() const { return value; }
            unsigned int value;

            static unsigned int numberDefaultConstructed;
            static unsigned int numberConstructedFromInt;
            static unsigned int numberDestroyed;
            static unsigned int numberCopied;
            static unsigned int numberAssigned;
            static unsigned int numberCloned;
            static unsigned int numberOf1ParameterConstructions;
            static unsigned int numberOf2ParameterConstructions;
            static unsigned int numberOf3ParameterConstructions;
            static unsigned int numberOfExpressionConstructions;
            static unsigned int numberOfExpressionAssignments;
    };
    
    template<typename T>
    std::ostream& operator<<(std::ostream& os, const CountedObject<T>& lhs)
    {
        os << lhs.GetValue();
        return os;
    }

    template<typename DerivedType>
    unsigned int CountedObject<DerivedType>::numberDefaultConstructed = 0;
    
    template<typename DerivedType>
    unsigned int CountedObject<DerivedType>::numberConstructedFromInt = 0;
    
    template<typename DerivedType>
    unsigned int CountedObject<DerivedType>::numberDestroyed = 0;
    
    template<typename DerivedType>
    unsigned int CountedObject<DerivedType>::numberCopied = 0;
    
    template<typename DerivedType>
    unsigned int CountedObject<DerivedType>::numberAssigned = 0;
    
    template<typename DerivedType>
    unsigned int CountedObject<DerivedType>::numberCloned = 0;

    template<typename DerivedType>
    unsigned int CountedObject<DerivedType>::numberOf1ParameterConstructions = 0;

    template<typename DerivedType>
    unsigned int CountedObject<DerivedType>::numberOf2ParameterConstructions = 0;

    template<typename DerivedType>
    unsigned int CountedObject<DerivedType>::numberOf3ParameterConstructions = 0;

    template<typename DerivedType>
    unsigned int CountedObject<DerivedType>::numberOfExpressionConstructions = 0;

    template<typename DerivedType>
    unsigned int CountedObject<DerivedType>::numberOfExpressionAssignments = 0;
}

#endif //NEKTAR_UNIT_TESTS_COUNTED_OBJECT_H
