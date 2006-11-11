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
// Description: 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_UNIT_TESTS_COUNTED_OBJECT_H
#define NEKTAR_UNIT_TESTS_COUNTED_OBJECT_H

#include <boost/test/unit_test.hpp>

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

            CountedObject(unsigned int v) :
                value(v)
            {
                ++numberConstructedFromInt;
            }

            CountedObject(const FICounterObject& rhs) :
                value(rhs.value)
            {
                ++numberCopied;
            }

            ~CountedObject()
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

            static void clearCounters()
            {
                numberDefaultConstructed = 0;
                numberConstructedFromInt = 0;
                numberDestroyed = 0;
                numberCopied = 0;
                numberAssigned = 0;
                numberCloned = 0;
            }

            static void check(unsigned int expectedDefaultConstructed, unsigned int expectedConstructedFromInt,
                                unsigned int expectedDestroyed, unsigned int expectedCopied, unsigned int expectedCloned,
                                unsigned int expectedAssigned)
            {
                BOOST_CHECK(numberDefaultConstructed == expectedDefaultConstructed);
                BOOST_CHECK(numberConstructedFromInt == expectedConstructedFromInt);
                BOOST_CHECK(numberDestroyed == expectedDestroyed);
                BOOST_CHECK(numberCopied == expectedCopied);
                BOOST_CHECK(numberAssigned == expectedAssigned);
                BOOST_CHECK(numberCloned == expectedCloned);
            }

            unsigned int value;

            static unsigned int numberDefaultConstructed ;
            static unsigned int numberConstructedFromInt;
            static unsigned int numberDestroyed;
            static unsigned int numberCopied;
            static unsigned int numberAssigned;
            static unsigned int numberCloned;
};

#endif //NEKTAR_UNIT_TESTS_COUNTED_OBJECT_H

/**
    $Log: $
**/
