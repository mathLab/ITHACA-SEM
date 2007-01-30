///////////////////////////////////////////////////////////////////////////////
//
// File: ExpressionTemplates.hpp
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
// Three steps to integrate your classes with the expression templates.
//
// 1.  Provide a copy constructor and an assignment operator for your class which
//     takes a binary expression as the argument.
//
// 2.  Define a BinaryExpressionTraits object for each possible combination of
//     operators using your type.
//
// 3.  Define op_add, op_subtract, op_divide, op_multiply.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_HPP
#define NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_HPP

#include <LibUtilities/ExpressionTemplates/Expression.hpp>
#include <LibUtilities/ExpressionTemplates/ExpressionReturnType.hpp>
#include <LibUtilities/ExpressionTemplates/ConstantExpression.hpp>
#include <LibUtilities/ExpressionTemplates/UnaryExpression.hpp>
#include <LibUtilities/ExpressionTemplates/UnaryExpressionTraits.hpp>
#include <LibUtilities/ExpressionTemplates/NegateOp.hpp>
#include <LibUtilities/ExpressionTemplates/NullOp.hpp>

#include <LibUtilities/ExpressionTemplates/ConstantExpressionTraits.hpp>
#include <LibUtilities/ExpressionTemplates/BinaryExpressionTraits.hpp>
#include <LibUtilities/ExpressionTemplates/BinaryExpression.hpp>
#include <LibUtilities/ExpressionTemplates/BinaryExpressionOperators.hpp>


#define ENABLE_EXPRESSION_TEMPLATE_OPERATOR_ADD(lhsClassName, rhsClassName) \
    expt::BinaryExpressionType<lhsClassName, expt::AddOp, rhsClassName >::Type \
    operator+(const lhsClassName& lhs, const rhsClassName& rhs) \
{ \
    return expt::CreateBinaryExpression<expt::AddOp>(lhs, rhs); \
}
//         typedef expt::ConstantExpressionPolicy<lhsClassName > LhsExpressionPolicy; \
//         typedef expt::ConstantExpressionPolicy<rhsClassName > RhsExpressionPolicy; \
//         typedef expt::BinaryExpressionPolicy<LhsExpressionPolicy, RhsExpressionPolicy, expt::AddOp> ResultExpressionPolicyType; \
//         typedef expt::Expression<ResultExpressionPolicyType> ResultExpressionType; \
//         \
//         ResultExpressionPolicyType::DataType d((expt::Expression<LhsExpressionPolicy>(lhs)), expt::Expression<RhsExpressionPolicy>(rhs)); \
//         return ResultExpressionType(d); \
// }

#define ENABLE_EXPRESSION_TEMPLATE_OPERATOR_SUBTRACT(lhsClassName, rhsClassName) \
    expt::BinaryExpressionType<lhsClassName, expt::SubtractOp, rhsClassName >::Type \
    operator-(const lhsClassName& lhs, const rhsClassName& rhs) \
{ \
    return expt::CreateBinaryExpression<expt::SubtractOp>(lhs, rhs); \
}
//         typedef expt::ConstantExpressionPolicy<lhsClassName > LhsExpressionPolicy; \
//         typedef expt::ConstantExpressionPolicy<rhsClassName > RhsExpressionPolicy; \
//         typedef expt::BinaryExpressionPolicy<LhsExpressionPolicy, RhsExpressionPolicy, expt::SubtractOp> ResultExpressionPolicyType; \
//         typedef expt::Expression<ResultExpressionPolicyType> ResultExpressionType; \
//         \
//         ResultExpressionPolicyType::DataType d((expt::Expression<LhsExpressionPolicy>(lhs)), expt::Expression<RhsExpressionPolicy>(rhs)); \
//         return ResultExpressionType(d); \
// }

#define ENABLE_EXPRESSION_TEMPLATE_OPERATOR_MULTIPLY(lhsClassName, rhsClassName) \
    expt::BinaryExpressionType<lhsClassName, expt::MultiplyOp, rhsClassName >::Type \
    operator*(const lhsClassName& lhs, const rhsClassName& rhs) \
{ \
    return expt::CreateBinaryExpression<expt::MultiplyOp>(lhs, rhs); \
}
//         typedef expt::ConstantExpressionPolicy<lhsClassName > LhsExpressionPolicy; \
//         typedef expt::ConstantExpressionPolicy<rhsClassName > RhsExpressionPolicy; \
//         typedef expt::BinaryExpressionPolicy<LhsExpressionPolicy, RhsExpressionPolicy, expt::MultiplyOp> ResultExpressionPolicyType; \
//         typedef expt::Expression<ResultExpressionPolicyType> ResultExpressionType; \
//         \
//         ResultExpressionPolicyType::DataType d((expt::Expression<LhsExpressionPolicy>(lhs)), expt::Expression<RhsExpressionPolicy>(rhs)); \
//         return ResultExpressionType(d); \
// }

#define ENABLE_EXPRESSION_TEMPLATE_OPERATOR_DIVIDE(lhsClassName, rhsClassName) \
    expt::BinaryExpressionType<lhsClassName, expt::DivideOp, rhsClassName >::Type \
    operator/(const lhsClassName& lhs, const rhsClassName& rhs) \
{ \
    return expt::CreateBinaryExpression<expt::DivideOp>(lhs, rhs); \
}
//         typedef expt::ConstantExpressionPolicy<lhsClassName > LhsExpressionPolicy; \
//         typedef expt::ConstantExpressionPolicy<rhsClassName > RhsExpressionPolicy; \
//         typedef expt::BinaryExpressionPolicy<LhsExpressionPolicy, RhsExpressionPolicy, expt::DivideOp> ResultExpressionPolicyType; \
//         typedef expt::Expression<ResultExpressionPolicyType> ResultExpressionType; \
//         \
//         ResultExpressionPolicyType::DataType d((expt::Expression<LhsExpressionPolicy>(lhs)), expt::Expression<RhsExpressionPolicy>(rhs)); \
//         return ResultExpressionType(d); \
// }

// Utility macros to generate the operators needed by a class supporting expression templates.
#define ENABLE_EXPRESSION_TEMPLATE_OPERATORS2(lhsClassName, rhsClassName) \
    ENABLE_EXPRESSION_TEMPLATE_OPERATOR_DIVIDE(lhsClassName, rhsClassName) \
    ENABLE_EXPRESSION_TEMPLATE_OPERATOR_ADD(lhsClassName, rhsClassName) \
    ENABLE_EXPRESSION_TEMPLATE_OPERATOR_SUBTRACT(lhsClassName, rhsClassName) \
    ENABLE_EXPRESSION_TEMPLATE_OPERATOR_MULTIPLY(lhsClassName, rhsClassName) \





#define ENABLE_EXPRESSION_TEMPLATE_OPERATORS(className) \
    ENABLE_EXPRESSION_TEMPLATE_OPERATORS2(className, className)




/// Requirements to enable your type to do expressions.
/// For each operation, define the appropriate AdditionTraits, SubtractionTraits, DivisionTraits, and MultiplicationTraits objects.
///
/// Override commutative traits if they do not follow standard practice.
/// Override associative traits if they do not follow standard practice.
namespace Nektar
{
    namespace expt
    {
        template<typename ExpressionPolicyType>
        std::ostream& operator<<(std::ostream& os, const Expression<ExpressionPolicyType>& exp)
        {
            exp.Print(os);
            return os;
        }
    }
}

// Is this the correct place to put high level expression template documentation?
// To create an object that can use expression template:
//
// If you want to be able to add two objects together, you must
// 1. Provide an add method that takes three parameters: lhs, rhs, result.
// 2. Provide a += method.
#endif // NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_HPP

/**
    $Log: ExpressionTemplates.hpp,v $
    Revision 1.8  2007/01/16 17:37:55  bnelson
    Wrapped everything with #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

    Revision 1.7  2007/01/16 05:29:50  bnelson
    Major improvements for expression templates.

    Revision 1.6  2006/11/12 17:58:52  bnelson
    *** empty log message ***

    Revision 1.5  2006/09/11 03:24:24  bnelson
    Updated expression templates so they are all specializations of an Expression object, using policies to differentiate.

    Revision 1.4  2006/08/28 02:39:53  bnelson
    *** empty log message ***

    Revision 1.3  2006/08/25 01:33:47  bnelson
    no message

    Revision 1.2  2006/06/01 13:44:28  kirby
    *** empty log message ***

    Revision 1.1  2006/06/01 09:20:56  kirby
    *** empty log message ***

    Revision 1.1  2006/05/04 18:57:42  kirby
    *** empty log message ***

    Revision 1.2  2006/01/31 13:51:12  bnelson
    Updated for new configure.

    Revision 1.1  2006/01/10 14:50:31  bnelson
    Initial Revision.

**/
