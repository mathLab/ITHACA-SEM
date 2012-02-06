///////////////////////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////////////////////

#include "ExpressionTemplateObjects.h"

namespace expt
{

    TestObject Add(const TestObject& lhs, const TestObject& rhs)
    {
        return TestObject(lhs.value + rhs.value);
    }

    void AddEqual(TestObject& accumulator, const TestObject& rhs)
    {
        accumulator.value += rhs.value;
    }

    void Add(TestObject& accumulator, const TestObject& lhs, const TestObject& rhs)
    {
        accumulator.value = lhs.value + rhs.value;
    }
    
    TestObject Multiply(const TestObject& lhs, const TestObject& rhs)
    {
        return TestObject(lhs.value * rhs.value);
    }

    void MultiplyEqual(TestObject& accumulator, const TestObject& rhs)
    {
        accumulator.value *= rhs.value;
    }

    void Multiply(TestObject& accumulator, const TestObject& lhs, const TestObject& rhs)
    {
        accumulator.value = lhs.value * rhs.value;
    }

    GENERATE_ADDITION_OPERATOR(TestObject, 0, TestObject, 0);
    GENERATE_MULTIPLICATION_OPERATOR(TestObject, 0, TestObject, 0);

    std::ostream& operator<<(std::ostream& os, const TestObject& obj)
    {
        os << obj.value;
        return os;
    }

    bool operator==(const TestObject& lhs, const TestObject& rhs)
    {
        return lhs.value == rhs.value;
    }


    TestObjectC Add(const TestObjectC& lhs, const TestObjectC& rhs)
    {
        return TestObjectC(lhs.value + rhs.value);
    }

    void AddEqual(TestObjectC& accumulator, const TestObjectC& rhs)
    {
        accumulator.value += rhs.value;
    }

    void Add(TestObjectC& accumulator, const TestObjectC& lhs, const TestObjectC& rhs)
    {
        accumulator.value = lhs.value + rhs.value;
    }

    TestObjectC Multiply(const TestObjectC& lhs, const TestObjectC& rhs)
    {
        return TestObjectC(lhs.value * rhs.value);
    }

    void MultiplyEqual(TestObjectC& accumulator, const TestObjectC& rhs)
    {
        accumulator.value *= rhs.value;
    }

    void Multiply(TestObjectC& accumulator, const TestObjectC& lhs, const TestObjectC& rhs)
    {
        accumulator.value = lhs.value * rhs.value;
    }

    GENERATE_ADDITION_OPERATOR(TestObjectC, 0, TestObjectC, 0);    
    GENERATE_MULTIPLICATION_OPERATOR(TestObjectC, 0, TestObjectC, 0);    

    std::ostream& operator<<(std::ostream& os, const TestObjectC& obj)
    {
        os << obj.value;
        return os;
    }

    bool operator==(const TestObjectC& lhs, const TestObjectC& rhs)
    {
        return lhs.value == rhs.value;
    }


    TestObjectA Add(const TestObjectA& lhs, const TestObjectA& rhs)
    {
        return TestObjectA(lhs.value + rhs.value);
    }

    void AddEqual(TestObjectA& accumulator, const TestObjectA& rhs)
    {
        accumulator.value += rhs.value;
    }

    void Add(TestObjectA& accumulator, const TestObjectA& lhs, const TestObjectA& rhs)
    {
        accumulator.value = lhs.value + rhs.value;
    }

    TestObjectA Multiply(const TestObjectA& lhs, const TestObjectA& rhs)
    {
        return TestObjectA(lhs.value * rhs.value);
    }

    void MultiplyEqual(TestObjectA& accumulator, const TestObjectA& rhs)
    {
        accumulator.value *= rhs.value;
    }

    void Multiply(TestObjectA& accumulator, const TestObjectA& lhs, const TestObjectA& rhs)
    {
        accumulator.value = lhs.value * rhs.value;
    }

    GENERATE_ADDITION_OPERATOR(TestObjectA, 0, TestObjectA, 0);
    GENERATE_MULTIPLICATION_OPERATOR(TestObjectA, 0, TestObjectA, 0);

    std::ostream& operator<<(std::ostream& os, const TestObjectA& obj)
    {
        os << obj.value;
        return os;
    }

    bool operator==(const TestObjectA& lhs, const TestObjectA& rhs)
    {
        return lhs.value == rhs.value;
    }


    TestObjectAC Add(const TestObjectAC& lhs, const TestObjectAC& rhs)
    {
        return TestObjectAC(lhs.value + rhs.value);
    }

    void AddEqual(TestObjectAC& accumulator, const TestObjectAC& rhs)
    {
        accumulator.value += rhs.value;
    }

    void Add(TestObjectAC& accumulator, const TestObjectAC& lhs, const TestObjectAC& rhs)
    {
        accumulator.value = lhs.value + rhs.value;
    }

    TestObjectAC Multiply(const TestObjectAC& lhs, const TestObjectAC& rhs)
    {
        return TestObjectAC(lhs.value * rhs.value);
    }

    void MultiplyEqual(TestObjectAC& accumulator, const TestObjectAC& rhs)
    {
        accumulator.value *= rhs.value;
    }

    void Multiply(TestObjectAC& accumulator, const TestObjectAC& lhs, const TestObjectAC& rhs)
    {
        accumulator.value = lhs.value * rhs.value;
    }

    GENERATE_ADDITION_OPERATOR(TestObjectAC, 0, TestObjectAC, 0);
    GENERATE_MULTIPLICATION_OPERATOR(TestObjectAC, 0, TestObjectAC, 0);

    std::ostream& operator<<(std::ostream& os, const TestObjectAC& obj)
    {
        os << obj.value;
        return os;
    }

    bool operator==(const TestObjectAC& lhs, const TestObjectAC& rhs)
    {
        return lhs.value == rhs.value;
    }
}
    

