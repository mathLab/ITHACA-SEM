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

#ifndef NEKTAR_UNIT_TESTS_EXPRESSION_TEMPLATE_UNIT_TESTS_EXPRESSION_TEMPLATE_OBJECTS_H
#define NEKTAR_UNIT_TESTS_EXPRESSION_TEMPLATE_UNIT_TESTS_EXPRESSION_TEMPLATE_OBJECTS_H

#include <ExpressionTemplates/Node.hpp>
#include <ExpressionTemplates/Operators.hpp>
#include <ExpressionTemplates/AssociativeTraits.hpp>
#include <ExpressionTemplates/CommutativeTraits.hpp>
#include <ExpressionTemplates/ExpressionTemplates.hpp>
#include <LibUtilities/BasicUtils/OperatorGenerators.hpp>

#include <boost/typeof/typeof.hpp>
#include  BOOST_TYPEOF_INCREMENT_REGISTRATION_GROUP()

namespace expt
{
    //The operators associated with this object are arbitrarily assigned to associative and/or
    //commutative categories for testing purposes.
    //
    // The suffix encodes the properties the object posseses
    // C - Commutative
    // A - Associative

    class TestObject
    {
        public:
            TestObject() : value(0.0) {}
            TestObject(double v) : value(v) {}

            template<typename T1, typename Op, typename T2>
            TestObject(const Node<T1, Op, T2>& node)
            {
                ExpressionEvaluator::Evaluate(node, *this);
            }

            double value;
    };

    TestObject Add(const TestObject& lhs, const TestObject& rhs);

    void AddEqual(TestObject& accumulator, const TestObject& rhs);

    void Add(TestObject& accumulator, const TestObject& lhs, const TestObject& rhs);
    
    TestObject Multiply(const TestObject& lhs, const TestObject& rhs);

    void MultiplyEqual(TestObject& accumulator, const TestObject& rhs);

    void Multiply(TestObject& accumulator, const TestObject& lhs, const TestObject& rhs);

    template<>
    struct AssociativeTraits<TestObject, expt::AddOp, TestObject, expt::AddOp> : public boost::false_type {};
    template<>
    struct CommutativeTraits<TestObject, expt::AddOp, TestObject> : public boost::false_type {};

    template<>
    struct AssociativeTraits<TestObject, expt::MultiplyOp, TestObject, expt::MultiplyOp> : public boost::false_type {};
    template<>
    struct CommutativeTraits<TestObject, expt::MultiplyOp, TestObject> : public boost::false_type {};

    BOOST_TYPEOF_REGISTER_TYPE(TestObject);
    DECLARE_ADDITION_OPERATOR(TestObject, 0, TestObject, 0);
    DECLARE_MULTIPLICATION_OPERATOR(TestObject, 0, TestObject, 0);

    std::ostream& operator<<(std::ostream& os, const TestObject& obj);

    bool operator==(const TestObject& lhs, const TestObject& rhs);

    class TestObjectC
    {
        public:
            TestObjectC() : value(0.0) {}
            TestObjectC(double v) : value(v) {}

            template<typename T1, typename Op, typename T2>
            TestObjectC(const Node<T1, Op, T2>& node)
            {
                ExpressionEvaluator::Evaluate(node, *this);
            }

            double value;
    };

    TestObjectC Add(const TestObjectC& lhs, const TestObjectC& rhs);

    void AddEqual(TestObjectC& accumulator, const TestObjectC& rhs);

    void Add(TestObjectC& accumulator, const TestObjectC& lhs, const TestObjectC& rhs);

    TestObjectC Multiply(const TestObjectC& lhs, const TestObjectC& rhs);

    void MultiplyEqual(TestObjectC& accumulator, const TestObjectC& rhs);

    void Multiply(TestObjectC& accumulator, const TestObjectC& lhs, const TestObjectC& rhs);

    template<>
    struct AssociativeTraits<TestObjectC, expt::AddOp, TestObjectC, expt::AddOp> : public boost::false_type {};
    template<>
    struct CommutativeTraits<TestObjectC, expt::AddOp, TestObjectC> : public boost::true_type {};
    template<>
    struct AssociativeTraits<TestObjectC, expt::MultiplyOp, TestObjectC, expt::MultiplyOp> : public boost::false_type {};
    template<>
    struct CommutativeTraits<TestObjectC, expt::MultiplyOp, TestObjectC> : public boost::false_type {};

    BOOST_TYPEOF_REGISTER_TYPE(TestObjectC);
    DECLARE_ADDITION_OPERATOR(TestObjectC, 0, TestObjectC, 0);    
    DECLARE_MULTIPLICATION_OPERATOR(TestObjectC, 0, TestObjectC, 0);    

    std::ostream& operator<<(std::ostream& os, const TestObjectC& obj);

    bool operator==(const TestObjectC& lhs, const TestObjectC& rhs);

    class TestObjectA
    {
        public:
            TestObjectA() : value(0.0) {}
            TestObjectA(double v) : value(v) {}

            template<typename T1, typename Op, typename T2>
            TestObjectA(const Node<T1, Op, T2>& node)
            {
                ExpressionEvaluator::Evaluate(node, *this);
            }

            double value;
    };

    TestObjectA Add(const TestObjectA& lhs, const TestObjectA& rhs);

    void AddEqual(TestObjectA& accumulator, const TestObjectA& rhs);

    void Add(TestObjectA& accumulator, const TestObjectA& lhs, const TestObjectA& rhs);

    TestObjectA Multiply(const TestObjectA& lhs, const TestObjectA& rhs);

    void MultiplyEqual(TestObjectA& accumulator, const TestObjectA& rhs);

    void Multiply(TestObjectA& accumulator, const TestObjectA& lhs, const TestObjectA& rhs);

    template<>
    struct AssociativeTraits<TestObjectA, expt::AddOp, TestObjectA, expt::AddOp> : public boost::true_type {};
    template<>
    struct CommutativeTraits<TestObjectA, expt::AddOp, TestObjectA> : public boost::false_type {};
    template<>
    struct AssociativeTraits<TestObjectA, expt::MultiplyOp, TestObjectA, expt::MultiplyOp> : public boost::false_type {};
    template<>
    struct CommutativeTraits<TestObjectA, expt::MultiplyOp, TestObjectA> : public boost::false_type {};

    BOOST_TYPEOF_REGISTER_TYPE(TestObjectA);
    DECLARE_ADDITION_OPERATOR(TestObjectA, 0, TestObjectA, 0);
    DECLARE_MULTIPLICATION_OPERATOR(TestObjectA, 0, TestObjectA, 0);

    std::ostream& operator<<(std::ostream& os, const TestObjectA& obj);

    bool operator==(const TestObjectA& lhs, const TestObjectA& rhs);

    class TestObjectAC
    {
        public:
            TestObjectAC() : value(0.0) {}
            TestObjectAC(double v) : value(v) {}
            
            template<typename T1, typename Op, typename T2>
            TestObjectAC(const Node<T1, Op, T2>& node)
            {
                ExpressionEvaluator::Evaluate(node, *this);
            }

            double value;
    };

    TestObjectAC Add(const TestObjectAC& lhs, const TestObjectAC& rhs);

    void AddEqual(TestObjectAC& accumulator, const TestObjectAC& rhs);

    void Add(TestObjectAC& accumulator, const TestObjectAC& lhs, const TestObjectAC& rhs);

    TestObjectAC Multiply(const TestObjectAC& lhs, const TestObjectAC& rhs);

    void MultiplyEqual(TestObjectAC& accumulator, const TestObjectAC& rhs);

    void Multiply(TestObjectAC& accumulator, const TestObjectAC& lhs, const TestObjectAC& rhs);

    template<>
    struct AssociativeTraits<TestObjectAC, expt::AddOp, TestObjectAC, expt::AddOp> : public boost::true_type {};
    template<>
    struct CommutativeTraits<TestObjectAC, expt::AddOp, TestObjectAC> : public boost::true_type {};
    template<>
    struct AssociativeTraits<TestObjectAC, expt::MultiplyOp, TestObjectAC, expt::MultiplyOp> : public boost::false_type {};
    template<>
    struct CommutativeTraits<TestObjectAC, expt::MultiplyOp, TestObjectAC> : public boost::false_type {};

    BOOST_TYPEOF_REGISTER_TYPE(TestObjectAC);
    DECLARE_ADDITION_OPERATOR(TestObjectAC, 0, TestObjectAC, 0);
    DECLARE_MULTIPLICATION_OPERATOR(TestObjectAC, 0, TestObjectAC, 0);

    std::ostream& operator<<(std::ostream& os, const TestObjectAC& obj);

    bool operator==(const TestObjectAC& lhs, const TestObjectAC& rhs);
}
    


#endif //NEKTAR_UNIT_TESTS_EXPRESSION_TEMPLATE_UNIT_TESTS_EXPRESSION_TEMPLATE_OBJECTS_H
