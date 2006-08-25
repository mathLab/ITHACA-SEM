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

namespace Nektar
{
    // OpType - A class with a statis method called Apply that takes a single 
    // parameter and returns a result of the same or different type.
    // A template parameter to allow a single OpType templated class to be 
    // used for a variety of types.
    template<template <typename> class OpType, typename ParameterType>
    class UnaryExpression : public Expression<typename OpType<ParameterType>::ResultType>
    {
        public:
            // Defined by the user who codes the operation.  They need to tell us what
            // the result type of the operation will be.
            typedef typename OpType<ParameterType>::ResultType ResultType;

        public:
            explicit UnaryExpression(typename boost::call_traits<ParameterType>::const_reference value) :
                m_value(value)
            {
            }

            UnaryExpression(const UnaryExpression<OpType, ParameterType>& rhs) :
                m_value(rhs.m_value)
            {
            }

        private:
            virtual void DoApply(typename boost::call_traits<ResultType>::reference result) const
            {
                OpType<ParameterType>::Apply(m_value, result);
            }

            UnaryExpression<OpType, ParameterType>& operator=(const UnaryExpression<OpType, ParameterType>& rhs);

            // This should always be an expression, but it really doesn't matter
            // as long as the interface is kept.
            ParameterType m_value;
    };

    //template<typename ResultType>
    //UnaryExpression<NegateOp, ResultType> operator-(const Expression<ResultType>& rhs)
    //{
    //    
    //}


    //// To evaluate an expression, I need to know if the lhs or rhs of the binary tree
    //// matches the overall result type.  This can be done with enable_if and disable_if,
    //// but for some reason the work for global functions but not member functions.  
    //// Therefore, these global functions exist solely so the binary functions can forward 
    //// the evaluation.
    //
    //template<template <typename, typename> class OpType, typename LhsType, typename RhsType>
    //class BinaryExpression
    //{
    //    
    //};

    //template<template <typename, typename> class OpType, typename LhsType, typename RhsType>
    //class BinaryExpression<OpType, ConstantExpression<LhsType>, ConstantExpression<RhsType> >
    //{
    //    public:
    //        typedef typename GetReturnType<LhsType>::result_type lhs_result_type;
    //        typedef typename GetReturnType<RhsType>::result_type rhs_result_type;

    //        typedef typename OpType<lhs_result_type, rhs_result_type>::result_type result_type;

    //    public:
    //        BinaryExpression(const ConstantExpression<LhsType>& lhs, const ConstantExpression<RhsType>& rhs) :
    //            m_lhs(lhs),
    //            m_rhs(rhs)
    //        {
    //        }

    //        BinaryExpression(const BinaryExpression& rhs) :
    //            m_lhs(rhs.m_lhs),
    //            m_rhs(rhs.m_rhs)
    //        {
    //        }

    //        ~BinaryExpression() {}

    //        BinaryExpression& operator=(const BinaryExpression& rhs)
    //        {
    //            m_lhs = rhs.m_lhs;
    //            m_rhs = rhs.m_rhs;
    //        }

    //        void apply(typename boost::call_traits<result_type>::reference result) const
    //        {
    //            OpType<lhs_result_type, rhs_result_type>::apply(m_lhs.getValue(), m_rhs.getValue(), result);
    //        }

    //    private:
    //        ConstantExpression<LhsType> m_lhs;
    //        ConstantExpression<RhsType> m_rhs;
    //};

    //// Case 1 - The lhs does match the result type.
    //template<template <typename, typename> class OpType, typename LhsType, typename RhsType, typename ResultType>
    //void evaluate_expression(typename boost::call_traits<LhsType>::param_type lhs,
    //                         typename boost::call_traits<RhsType>::param_type rhs,
    //                         typename boost::call_traits<ResultType>::reference result,
    //                         typename boost::enable_if<boost::is_same<typename GetReturnType<LhsType>::result_type, ResultType> >::type* v = 0)
    //{
    //    typedef typename GetReturnType<LhsType>::result_type lhs_result_type;
    //    typedef typename GetReturnType<RhsType>::result_type rhs_result_type;

    //    lhs.apply(result);
    //    OpType<lhs_result_type, rhs_result_type>::apply(result, rhs.getValue(), result);
    //}

    //template<template <typename, typename> class OpType, typename LhsType, typename RhsType, typename ResultType>
    //void evaluate_expression(typename boost::call_traits<LhsType>::param_type lhs,
    //                         typename boost::call_traits<RhsType>::param_type rhs,
    //                         typename boost::call_traits<ResultType>::reference result,
    //                         typename boost::disable_if<boost::is_same<typename GetReturnType<LhsType>::result_type, ResultType> >::type* v = 0)
    //{
    //    typedef typename GetReturnType<LhsType>::result_type LhsResultType;
    //    typedef typename GetReturnType<RhsType>::result_type RhsResultType;

    //    LhsResultType temp;
    //    lhs.apply(temp);
    //    OpType<LhsResultType, RhsResultType>::apply(temp, rhs.getValue(), result);
    //}

    //template<template <typename, typename> class OpType, typename LhsType, typename RhsDataType>
    //class BinaryExpression<OpType, LhsType, ConstantExpression<RhsDataType> >
    //{
    //    public:
    //        typedef ConstantExpression<RhsDataType> RhsType;
    //        typedef typename GetReturnType<LhsType>::result_type LhsResultType;
    //        typedef typename GetReturnType<RhsType>::result_type RhsResultType;

    //        typedef typename OpType<LhsResultType, RhsResultType>::result_type ResultType;

    //    public:
    //        BinaryExpression(typename boost::call_traits<LhsType>::param_type lhs, 
    //                         typename boost::call_traits<RhsType>::param_type rhs) :
    //            m_lhs(lhs),
    //            m_rhs(rhs)
    //        {
    //        }

    //        BinaryExpression(const BinaryExpression& rhs) :
    //            m_lhs(rhs.m_lhs),
    //            m_rhs(rhs.m_rhs)
    //        {
    //        }

    //        ~BinaryExpression() {}

    //        BinaryExpression& operator=(const BinaryExpression& rhs)
    //        {
    //            m_lhs = rhs.m_lhs;
    //            m_rhs = rhs.m_rhs;
    //        }

    //        // Case 1 - The return type of the lhs is the same as the result, so we can accumulate.
    //        void apply(typename boost::call_traits<ResultType>::reference result) const
    //        {
    //            evaluate_expression<OpType, LhsType, RhsType, ResultType>(m_lhs, m_rhs, result);
    //        }

    //        //// Case 2 - rhs result type is the same as the final result type and the lhs is not.
    //        //void apply(typename boost::call_traits<result_type>::reference result, 
    //        //    typename boost::disable_if<boost::is_same<int, int> >::type* param1 = 0,
    //        //    double t = 0) const
    //        //{
    //        //    // This assumes the the result type of the lhs is the same as the result type of 
    //        //    // the rhs.
    //        //    //m_lhs.apply(result);
    //        //    //m_rhs.apply(result);
    //        //    //OpType<typename LhsType::result_type, RhsType>::apply(result, m_rhs.getValue(), result);
    //        //}

    //    private:
    //        LhsType m_lhs;
    //        RhsType m_rhs;
    //};

    //template<template <typename, typename> class OpType, typename LhsType, typename RhsType>
    //class GetReturnType<BinaryExpression<OpType, LhsType, RhsType> >
    //{
    //    public:
    //        typedef typename BinaryExpression<OpType, LhsType, RhsType>::result_type result_type;
    //};
}

#endif // NEKTAR_LIB_UTILITIES_UNARY_EXPRESSION_HPP

/**
    $Log: $
**/
