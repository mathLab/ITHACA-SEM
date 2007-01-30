///////////////////////////////////////////////////////////////////////////////
//
// File: UnaryExpression.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_UNARY_EXPRESSION_HPP
#define NEKTAR_LIB_UTILITIES_UNARY_EXPRESSION_HPP

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>

#include <LibUtilities/ExpressionTemplates/Expression.hpp>
#include <LibUtilities/ExpressionTemplates/UnaryExpressionTraits.hpp>

#include <LibUtilities/ExpressionTemplates/NegateOp.hpp>
#include <LibUtilities/ExpressionTemplates/Accumulator.hpp>


namespace Nektar
{
    namespace expt
    {

        template<typename InputExpressionPolicyType, template <typename> class OpType>
        class UnaryExpressionPolicy
        {
            public:
                typedef Expression<InputExpressionPolicyType> InputExpressionType;
                typedef typename InputExpressionType::ResultType InputExpressionResultType;
                typedef typename OpType<InputExpressionPolicyType>::ResultType ResultType;
                
                typedef InputExpressionType DataType;
        };

        // OpType - A class with a statis method called Apply that takes a single
        // parameter and returns a result of the same or different type.
        // A template parameter to allow a single OpType templated class to be
        // used for a variety of types.
        // ParameterType - A type which follows the expression interface.
        //template<template <typename> class OpType, typename InputExpressionType>
        //class UnaryExpression
        //template<unsigned int NumberOfArguments, typename Argument1Type, typename Argument2Type, template <typename> class ArgumentUnaryOpType, template <typename, typename> class BinaryOpType, typename OpType>
        //class Expression<1, Expression<NumberOfArguments, Argument1Type, Argument2Type, ArgumentUnaryOpType, BinaryOpType>, void, OpType>
//         template<typename InputExpressionPolicyType, template <typename> class OpType>
//         class Expression<UnaryExpressionPolicy<Expression<InputExpressionPolicyType>, OpType> >
//         {
//             public:
//                 // Defined by the user who codes the operation.  They need to tell us what
//                 // the result type of the operation will be.
//                 typedef UnaryExpressionPolicy<Expression<InputExpressionPolicyType>, OpType > PolicyType;
//                 typedef Expression<InputExpressionPolicyType> InputExpressionType;
//                 typedef typename InputExpressionType::ResultType InputExpressionResultType;
//                 typedef typename OpType<InputExpressionResultType>::ResultType ResultType;
// 
//                 typedef typename OpType<InputExpressionResultType>::TraitsType OpTraitsType;
//                 //typedef typename ExpressionMetadataChooser<OpTraitsType>::MetadataType InputExpressionMetadataType;
//                 typedef typename ExpressionMetadataChooser<OpTraitsType>::MetadataType MetadataType;
//                 
//                 typedef Expression<UnaryExpressionPolicy<Expression<InputExpressionPolicyType>, OpType> > ThisType;
// 
//             public:
//                 explicit Expression(const InputExpressionType& value) :
//                     m_value(value),
//                     m_metadata(value.GetMetadata())
//                 {
//                 }
// 
//                 Expression(const ThisType& rhs) :
//                     m_value(rhs.m_value),
//                     m_metadata(rhs.m_metadata)
//                 {
//                 }
// 
//                 ~Expression() {}
// 
//                 // Two cases for the apply method.
//                 // 1.  Result and Parameter types are the same.
//                 // 2.  Result and Parameter types are different.
//                 void Apply(typename boost::call_traits<ResultType>::reference result) const
//                 {
//                     Accumulator<ResultType> accum(result);
//                     Apply(accum);
//                 }
// 
//                 template<typename IncomingOpType>
//                 void ApplyEqual(typename boost::call_traits<ResultType>::reference result) const
//                 {
//                     m_value.ApplyEqual<IncomingOpType>(result);
//                     OpType<ResultType>::Apply(result);
//                 }
// 
//                 void Apply(Accumulator<ResultType>& accum) const
//                 {
//                     // Evaluate the expression up to this point.
//                     m_value.Apply(accum);
// 
//                     // Now apply the operator to the result.
//                     OpType<ResultType>::Apply(accum.GetData());
//                 }
// 
//                 template<typename IncomingOpType>
//                 void ApplyEqual(Accumulator<ResultType>& accum) const
//                 {
//                     m_value.ApplyEqual<IncomingOpType>(accum.GetData());
//                     OpType<ResultType>::Apply(accum.GetData());
//                 }
//                 
//                 const MetadataType& GetMetadata() const
//                 {
//                     return m_metadata;
//                 }
// 
//                 void Print(std::ostream& os) const
//                 {
//                     os << OpType<ResultType>::AsString() << m_value;
//                 }
//                 
//             private:
//                 ThisType& operator=(const ThisType& rhs);
// 
//                 InputExpressionType m_value;
//                 MetadataType m_metadata;
//         };

        template<typename InputExpressionPolicyType>
        Expression<UnaryExpressionPolicy<Expression<InputExpressionPolicyType>, NegateOp> > operator-(const Expression<InputExpressionPolicyType>& rhs)
        {
            return Expression<UnaryExpressionPolicy<Expression<InputExpressionPolicyType>,NegateOp> >(rhs);
        }
    }
}

#endif // NEKTAR_LIB_UTILITIES_UNARY_EXPRESSION_HPP

/**
    $Log: UnaryExpression.hpp,v $
    Revision 1.9  2007/01/16 17:37:56  bnelson
    Wrapped everything with #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

    Revision 1.8  2007/01/16 05:29:50  bnelson
    Major improvements for expression templates.

    Revision 1.7  2006/11/12 17:58:53  bnelson
    *** empty log message ***

    Revision 1.6  2006/09/16 23:52:56  bnelson
    Changed unary operations to rely on methods outside the class (to aid adding expression templates to classes that can't be modified)

    Revision 1.5  2006/09/14 02:08:59  bnelson
    Fixed gcc compile errors

    Revision 1.4  2006/09/11 03:24:24  bnelson
    Updated expression templates so they are all specializations of an Expression object, using policies to differentiate.

    Revision 1.3  2006/08/28 02:39:53  bnelson
    *** empty log message ***

    Revision 1.2  2006/08/27 02:11:30  bnelson
    Added support for negating an expression.

    Revision 1.1  2006/08/25 01:33:48  bnelson
    no message

**/
