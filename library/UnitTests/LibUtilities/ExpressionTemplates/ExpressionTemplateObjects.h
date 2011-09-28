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
#include <LibUtilities/BasicUtils/OperatorGenerators.hpp>

#include <boost/typeof/typeof.hpp>
#include  BOOST_TYPEOF_INCREMENT_REGISTRATION_GROUP()

namespace Nektar
{
     //The operators associated with this object are arbitrarily assigned to associative and/or
     //commutative categories for testing purposes.

    class ExpressionTemplateTestObject00
    {
        public:
            ExpressionTemplateTestObject00() : value(0.0) {}
            ExpressionTemplateTestObject00(double v) : value(v) {}
            double value;
    };

    ExpressionTemplateTestObject00 Add(const ExpressionTemplateTestObject00& lhs, const ExpressionTemplateTestObject00& rhs)
    {
        return ExpressionTemplateTestObject00(lhs.value + rhs.value);
    }

    void AddEqual(ExpressionTemplateTestObject00& accumulator, const ExpressionTemplateTestObject00& rhs)
    {
        accumulator.value += rhs.value;
    }

    void Add(ExpressionTemplateTestObject00& accumulator, const ExpressionTemplateTestObject00& lhs, const ExpressionTemplateTestObject00& rhs)
    {
        accumulator.value = lhs.value + rhs.value;
    }
    
    ExpressionTemplateTestObject00 Multiply(const ExpressionTemplateTestObject00& lhs, const ExpressionTemplateTestObject00& rhs)
    {
        return ExpressionTemplateTestObject00(lhs.value * rhs.value);
    }

    void MultiplyEqual(ExpressionTemplateTestObject00& accumulator, const ExpressionTemplateTestObject00& rhs)
    {
        accumulator.value *= rhs.value;
    }

    void Multiply(ExpressionTemplateTestObject00& accumulator, const ExpressionTemplateTestObject00& lhs, const ExpressionTemplateTestObject00& rhs)
    {
        accumulator.value = lhs.value * rhs.value;
    }
}

namespace expt
{
    template<>
    struct AssociativeTraits<Nektar::ExpressionTemplateTestObject00, expt::AddOp, Nektar::ExpressionTemplateTestObject00, expt::AddOp> : public boost::false_type {};
    template<>
    struct CommutativeTraits<Nektar::ExpressionTemplateTestObject00, expt::AddOp, Nektar::ExpressionTemplateTestObject00> : public boost::false_type {};

    template<>
    struct AssociativeTraits<Nektar::ExpressionTemplateTestObject00, expt::MultiplyOp, Nektar::ExpressionTemplateTestObject00, expt::MultiplyOp> : public boost::false_type {};
    template<>
    struct CommutativeTraits<Nektar::ExpressionTemplateTestObject00, expt::MultiplyOp, Nektar::ExpressionTemplateTestObject00> : public boost::false_type {};
}

namespace Nektar
{
    BOOST_TYPEOF_REGISTER_TYPE(ExpressionTemplateTestObject00);
    GENERATE_ADDITION_OPERATOR(ExpressionTemplateTestObject00, 0, ExpressionTemplateTestObject00, 0);
    GENERATE_MULTIPLICATION_OPERATOR(ExpressionTemplateTestObject00, 0, ExpressionTemplateTestObject00, 0);

    class ExpressionTemplateTestObject01
    {
        public:
            ExpressionTemplateTestObject01() : value(0.0) {}
            ExpressionTemplateTestObject01(double v) : value(v) {}
            double value;
    };

    ExpressionTemplateTestObject01 Add(const ExpressionTemplateTestObject01& lhs, const ExpressionTemplateTestObject01& rhs)
    {
        return ExpressionTemplateTestObject01(lhs.value + rhs.value);
    }

    void AddEqual(ExpressionTemplateTestObject01& accumulator, const ExpressionTemplateTestObject01& rhs)
    {
        accumulator.value += rhs.value;
    }

    void Add(ExpressionTemplateTestObject01& accumulator, const ExpressionTemplateTestObject01& lhs, const ExpressionTemplateTestObject01& rhs)
    {
        accumulator.value = lhs.value + rhs.value;
    }

    ExpressionTemplateTestObject01 Multiply(const ExpressionTemplateTestObject01& lhs, const ExpressionTemplateTestObject01& rhs)
    {
        return ExpressionTemplateTestObject01(lhs.value * rhs.value);
    }

    void MultiplyEqual(ExpressionTemplateTestObject01& accumulator, const ExpressionTemplateTestObject01& rhs)
    {
        accumulator.value *= rhs.value;
    }

    void Multiply(ExpressionTemplateTestObject01& accumulator, const ExpressionTemplateTestObject01& lhs, const ExpressionTemplateTestObject01& rhs)
    {
        accumulator.value = lhs.value * rhs.value;
    }
}

namespace expt
{
    template<>
    struct AssociativeTraits<Nektar::ExpressionTemplateTestObject01, expt::AddOp, Nektar::ExpressionTemplateTestObject01, expt::AddOp> : public boost::false_type {};
    template<>
    struct CommutativeTraits<Nektar::ExpressionTemplateTestObject01, expt::AddOp, Nektar::ExpressionTemplateTestObject01> : public boost::true_type {};
    template<>
    struct AssociativeTraits<Nektar::ExpressionTemplateTestObject01, expt::MultiplyOp, Nektar::ExpressionTemplateTestObject01, expt::MultiplyOp> : public boost::false_type {};
    template<>
    struct CommutativeTraits<Nektar::ExpressionTemplateTestObject01, expt::MultiplyOp, Nektar::ExpressionTemplateTestObject01> : public boost::false_type {};
}

namespace Nektar
{
    BOOST_TYPEOF_REGISTER_TYPE(ExpressionTemplateTestObject01);
    GENERATE_ADDITION_OPERATOR(ExpressionTemplateTestObject01, 0, ExpressionTemplateTestObject01, 0);    
    GENERATE_MULTIPLICATION_OPERATOR(ExpressionTemplateTestObject01, 0, ExpressionTemplateTestObject01, 0);    

    class ExpressionTemplateTestObject10
    {
        public:
            ExpressionTemplateTestObject10() : value(0.0) {}
            ExpressionTemplateTestObject10(double v) : value(v) {}
            double value;
    };

    ExpressionTemplateTestObject10 Add(const ExpressionTemplateTestObject10& lhs, const ExpressionTemplateTestObject10& rhs)
    {
        return ExpressionTemplateTestObject10(lhs.value + rhs.value);
    }

    void AddEqual(ExpressionTemplateTestObject10& accumulator, const ExpressionTemplateTestObject10& rhs)
    {
        accumulator.value += rhs.value;
    }

    void Add(ExpressionTemplateTestObject10& accumulator, const ExpressionTemplateTestObject10& lhs, const ExpressionTemplateTestObject10& rhs)
    {
        accumulator.value = lhs.value + rhs.value;
    }

    ExpressionTemplateTestObject10 Multiply(const ExpressionTemplateTestObject10& lhs, const ExpressionTemplateTestObject10& rhs)
    {
        return ExpressionTemplateTestObject10(lhs.value * rhs.value);
    }

    void MultiplyEqual(ExpressionTemplateTestObject10& accumulator, const ExpressionTemplateTestObject10& rhs)
    {
        accumulator.value *= rhs.value;
    }

    void Multiply(ExpressionTemplateTestObject10& accumulator, const ExpressionTemplateTestObject10& lhs, const ExpressionTemplateTestObject10& rhs)
    {
        accumulator.value = lhs.value * rhs.value;
    }
}

namespace expt
{
    template<>
    struct AssociativeTraits<Nektar::ExpressionTemplateTestObject10, expt::AddOp, Nektar::ExpressionTemplateTestObject10, expt::AddOp> : public boost::true_type {};
    template<>
    struct CommutativeTraits<Nektar::ExpressionTemplateTestObject10, expt::AddOp, Nektar::ExpressionTemplateTestObject10> : public boost::false_type {};
    template<>
    struct AssociativeTraits<Nektar::ExpressionTemplateTestObject10, expt::MultiplyOp, Nektar::ExpressionTemplateTestObject10, expt::MultiplyOp> : public boost::false_type {};
    template<>
    struct CommutativeTraits<Nektar::ExpressionTemplateTestObject10, expt::MultiplyOp, Nektar::ExpressionTemplateTestObject10> : public boost::false_type {};
}

namespace Nektar
{
    BOOST_TYPEOF_REGISTER_TYPE(ExpressionTemplateTestObject10);
    GENERATE_ADDITION_OPERATOR(ExpressionTemplateTestObject10, 0, ExpressionTemplateTestObject10, 0);
    GENERATE_MULTIPLICATION_OPERATOR(ExpressionTemplateTestObject10, 0, ExpressionTemplateTestObject10, 0);

    class ExpressionTemplateTestObject11
    {
        public:
            ExpressionTemplateTestObject11() : value(0.0) {}
            ExpressionTemplateTestObject11(double v) : value(v) {}
            double value;
    };

    ExpressionTemplateTestObject11 Add(const ExpressionTemplateTestObject11& lhs, const ExpressionTemplateTestObject11& rhs)
    {
        return ExpressionTemplateTestObject11(lhs.value + rhs.value);
    }

    void AddEqual(ExpressionTemplateTestObject11& accumulator, const ExpressionTemplateTestObject11& rhs)
    {
        accumulator.value += rhs.value;
    }

    void Add(ExpressionTemplateTestObject11& accumulator, const ExpressionTemplateTestObject11& lhs, const ExpressionTemplateTestObject11& rhs)
    {
        accumulator.value = lhs.value + rhs.value;
    }

    ExpressionTemplateTestObject11 Multiply(const ExpressionTemplateTestObject11& lhs, const ExpressionTemplateTestObject11& rhs)
    {
        return ExpressionTemplateTestObject11(lhs.value * rhs.value);
    }

    void MultiplyEqual(ExpressionTemplateTestObject11& accumulator, const ExpressionTemplateTestObject11& rhs)
    {
        accumulator.value *= rhs.value;
    }

    void Multiply(ExpressionTemplateTestObject11& accumulator, const ExpressionTemplateTestObject11& lhs, const ExpressionTemplateTestObject11& rhs)
    {
        accumulator.value = lhs.value * rhs.value;
    }
}

namespace expt
{
    template<>
    struct AssociativeTraits<Nektar::ExpressionTemplateTestObject11, expt::AddOp, Nektar::ExpressionTemplateTestObject11, expt::AddOp> : public boost::true_type {};
    template<>
    struct CommutativeTraits<Nektar::ExpressionTemplateTestObject11, expt::AddOp, Nektar::ExpressionTemplateTestObject11> : public boost::true_type {};
    template<>
    struct AssociativeTraits<Nektar::ExpressionTemplateTestObject11, expt::MultiplyOp, Nektar::ExpressionTemplateTestObject11, expt::MultiplyOp> : public boost::false_type {};
    template<>
    struct CommutativeTraits<Nektar::ExpressionTemplateTestObject11, expt::MultiplyOp, Nektar::ExpressionTemplateTestObject11> : public boost::false_type {};

}

namespace Nektar
{
    BOOST_TYPEOF_REGISTER_TYPE(ExpressionTemplateTestObject11);
    GENERATE_ADDITION_OPERATOR(ExpressionTemplateTestObject11, 0, ExpressionTemplateTestObject11, 0);
    GENERATE_MULTIPLICATION_OPERATOR(ExpressionTemplateTestObject11, 0, ExpressionTemplateTestObject11, 0);

}
    


#endif //NEKTAR_UNIT_TESTS_EXPRESSION_TEMPLATE_UNIT_TESTS_EXPRESSION_TEMPLATE_OBJECTS_H
