///////////////////////////////////////////////////////////////////////////////
//
// File: BinaryExpressionEvaluator.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_BINARY_EXPRESSION_EVALUATOR_HPP
#define NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_BINARY_EXPRESSION_EVALUATOR_HPP

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

#include <LibUtilities/ExpressionTemplates/BinaryExpressionEvaluatorFwd.hpp>
#include <LibUtilities/ExpressionTemplates/ConstantExpression.hpp>
#include <LibUtilities/ExpressionTemplates/BinaryExpressionPolicyFwd.hpp>
#include <LibUtilities/ExpressionTemplates/AssociativeTraits.hpp>
#include <LibUtilities/ExpressionTemplates/CommutativeTraits.hpp>
#include <LibUtilities/ExpressionTemplates/AssignableTraits.hpp>
#include <LibUtilities/ExpressionTemplates/BinaryOperators.hpp>
#include <LibUtilities/ExpressionTemplates/BinaryExpressionEvaluator_BinaryConstant.hpp>
#include <LibUtilities/ExpressionTemplates/BinaryExpressionEvaluator_BinaryBinary.hpp>

#include <boost/utility/enable_if.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/not.hpp>
#include <iostream>

namespace Nektar
{
    ///////////////////////////////////////////////////////////////
    // Constant-Constant Specializations
    //
    // These specializations deal with binary expressions in which 
    // both leaves are constant expressions.  
    ///////////////////////////////////////////////////////////////
    
    /// \brief Constant-Constant specialization 0
    /// 
    /// With no incoming parent type, this evaluation is simply A op B, such 
    /// as A+B or A-B.  We'll just call the appropriate NekOp function which 
    /// takes 3 parameters.
    template<typename LhsDataType, typename RhsDataType, typename ResultType, 
             template <typename, typename> class OpType>
    struct BinaryExpressionEvaluator<ConstantExpressionPolicy<LhsDataType>,
                                     ConstantExpressionPolicy<RhsDataType>,
                                     ResultType, OpType, BinaryNullOp>
    {
        static void Eval(const Expression<ConstantExpressionPolicy<LhsDataType> >& lhs, 
                         const Expression<ConstantExpressionPolicy<RhsDataType> >& rhs,
                         Accumulator<ResultType>& result)
        {
            OpType<LhsDataType, RhsDataType>::Apply(result, *lhs, *rhs);
        }
        
        #ifdef NEKTAR_UNIT_TESTS
        static const unsigned int ClassId = 100;
        #endif //NEKTAR_UNIT_TESTS

    };
        
    /// \brief Constant-Constant specialization 1
    /// 
    /// This case occurs as part of a larger expression, such as 
    /// A + (B+C).  By definition, since there is a parent op type, it can't happen in isolation.
    ///
    /// Examples:
    ///
    /// A + (B+C) will enter this class with ParentOpType = +, and OpType = +, 
    /// and evaluate to A += B, A += C.
    ///
    /// A + (B - C) will enter with ParentOpType = +, OpType = -, and evaluate to 
    /// A += B, A -= C.
    ///
    /// A - (B + C) will enter with ParentOpType = -, OPType = +, and evaluate to 
    /// A -= B, A -= C.
    template<typename LhsResultType, typename RhsResultType, typename ResultType, 
             template <typename, typename> class OpType,
             template <typename, typename> class ParentOpType>
    struct BinaryExpressionEvaluator<ConstantExpressionPolicy<LhsResultType>,
                                     ConstantExpressionPolicy<RhsResultType>,
                                     ResultType, OpType, ParentOpType>
    {   
        static void Eval(const Expression<ConstantExpressionPolicy<LhsResultType> >& lhs, 
                         const Expression<ConstantExpressionPolicy<RhsResultType> >& rhs,
                         Accumulator<ResultType>& accum)
        {
            typedef Expression<ConstantExpressionPolicy<LhsResultType> > LhsExpressionType;
            typedef Expression<ConstantExpressionPolicy<RhsResultType> > RhsExpressionType;

            ParentOpType<ResultType, LhsResultType>::ApplyEqual(accum, *lhs);
            typedef AssociativeTraits<ConstantExpressionPolicy<ResultType>, ParentOpType,
                                      BinaryExpressionPolicy<LhsExpressionType, OpType, RhsExpressionType> > TraitsType;

            typedef typename TraitsType::template RhsOpType<ResultType, RhsResultType>::type NewRhsOpType;
            NewRhsOpType::ApplyEqual(accum, *rhs);
        }
        
        #ifdef NEKTAR_UNIT_TESTS
        static const unsigned int ClassId = 101;
        #endif //NEKTAR_UNIT_TESTS
    };



//        ////////////////////////////////////////////////////////////////////////
//        // Constant-Binary Expressions
//        //
//        // The main difficulty with these expressions comes from needing to 
//        // evaluate the constant.
//        // first.  If we pass the accumulator to the left side first, then we 
//        // must create a temporary 
//        // on the right side unless the expression is associative.  If the expression 
//        // is associative 
//        // then there is no problem passing the accumulator with a ParentOpType.
//        //
//        // Another option occurs if the expression is commutative and the binary
//        // expression has the same 
//        // type as the accumulator.  In this case, pass the accumulator to the 
//        // right side and then 
//        // op= the left side.
//        //
//        // A last condition assumes the user has defined a left equal operator 
//        // if the expression is 
//        // not commutative nor associative.
//        //
//        // Therefore, the variables of interest are:
//        // a - Expression is Associative.
//        // b - Expression is Commutative.
//        // c - Left Equal Operator is defined.
//        // d - Constant Expression is result type.
//        // e - Binary Expression is result type.
//        // f - In A + (B-C), (A+B) is the result type.
//        // 
//        // Case 0 - Both sides need temporaries.  
//        // Case 1 - Rhs needs temporary.
//        // Case 2 - Evaluate Rhs, Left Op equal lhs.
//        // Case 3 - Swap operations.  Evaluate new lhs, op= new rhs.
//        // Case 4 - Evaluate lhs, evaluate rhs with parent op type.
//        // Case 5 - Evaluate new (A+B) then op= C
//        
//        // .i 7
//        // .o 6
//        // .ilb a b c d e f g
//        // .ob c0 c1 c2 c3 c4 c5
//        // .p 11
//        // -000-0- 100000
//        // ---00-- 100000
//        // 0000--- 100000
//        // 1--1--0 100000
//        
//        // 0--10-- 010000
//        // 0001--- 010000
//        
//        // -0101-- 001000
//        // 001-1-- 001000
//        
//        // -1-01-- 000100
//        // 01--1-- 000100
//        
//        // 1--1--1 000010
//         
//        // 100011- 000001
//        
//
//        ////////////////////////////////////////////////////////////////////////
//        
//        
        // Constant-Binary specialization 0 - Both sides need temporaries.
        // Therefore, the variables of interest are:
        // a - Expression is Associative.
        // b - Expression is Commutative.
        // c - Left Equal Operator is defined.
        // d - Constant Expression is result type.
        // e - Binary Expression is result type.
        // f - In A + (B-C), (A+B) is the result type.
        //
        // -000-0- 100000
        // ---00-- 100000
        // 0000--- 100000
        // 1--1--0 100000
        // This case is already handled by the general purpose specializations.
        template<typename LhsResultType, typename RhsLhsPolicyType, typename RhsRhsPolicyType,
                 template <typename, typename> class RhsOpType,
                 template <typename, typename> class OpType,
                 typename ResultType>
        struct BinaryExpressionEvaluator<ConstantExpressionPolicy<LhsResultType>,
                                        BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>,
                                        ResultType, OpType, BinaryNullOp,
                                        typename boost::enable_if
                                        <
                                            boost::mpl::or_
                                            <
                                                boost::mpl::and_
                                                <
                                                    boost::mpl::not_<IsCommutative<ConstantExpressionPolicy<LhsResultType>, OpType, BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> > >,
                                                    boost::mpl::not_<typename OpType<LhsResultType, typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>::ResultType>::HasOpLeftEqual>,
                                                    boost::mpl::not_<boost::is_same<LhsResultType, ResultType> >,
                                                    boost::mpl::not_<boost::is_same<ResultType, typename AssociativeTraits<ConstantExpressionPolicy<LhsResultType>, OpType, BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >::LhsAsBinaryExpressionResultType> >
                                                >,
                                                boost::mpl::and_
                                                <
                                                    boost::mpl::not_< typename AssociativeTraits<ConstantExpressionPolicy<LhsResultType>, OpType, BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >::IsAssociative>,
                                                    boost::mpl::not_<IsCommutative<ConstantExpressionPolicy<LhsResultType>, OpType, BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> > >,
                                                    boost::mpl::not_<typename OpType<LhsResultType, typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>::ResultType>::HasOpLeftEqual>,
                                                    boost::mpl::not_<boost::is_same<LhsResultType, ResultType> >
                                                >,
                                                boost::mpl::and_
                                                <
                                                    boost::mpl::not_<boost::is_same<ResultType, typename ConstantExpressionPolicy<LhsResultType>::ResultType> >,
                                                    boost::mpl::not_<boost::is_same<ResultType, typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>::ResultType> >
                                                >,
                                                boost::mpl::and_
                                                <
                                                    typename AssociativeTraits<ConstantExpressionPolicy<LhsResultType>, OpType, BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >::IsAssociative,
                                                    boost::is_same<ResultType, typename ConstantExpressionPolicy<LhsResultType>::ResultType>, 
                                                    boost::mpl::not_<typename AssociativeTraits<ConstantExpressionPolicy<LhsResultType>, OpType, BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >::OpEqualsAreDefined>
                                                >
                                            >
                                         >::type
                                         >
        {
            static void Eval(const Expression<ConstantExpressionPolicy<LhsResultType> >& lhs,
                             const Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >& rhs, 
                             Accumulator<ResultType>& accum)
            {
                typedef typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>::ResultType RhsResultType;
                
                LhsResultType lhsTemp = lhs.Evaluate();
                RhsResultType rhsTemp = rhs.Evaluate();
                OpType<LhsResultType, RhsResultType>::Apply(accum, lhsTemp, rhsTemp);
            }
            
            #ifdef NEKTAR_UNIT_TESTS
            static const unsigned int ClassNum = 1;
            #endif //NEKTAR_UNIT_TESTS
        };


        // Constant-Binary specialization 0 - Both sides need temporaries.
        // Therefore, the variables of interest are:
        // a - Expression is Associative.
        // b - Expression is Commutative.
        // c - Left Equal Operator is defined.
        // d - Constant Expression is result type.
        // e - Binary Expression is result type.
        // f - In A + (B-C), (A+B) is the result type.
        //
        // 0--10- 010000
        // 0001-- 010000
        // This case is already handled by the general purpose specializations.
        template<typename LhsResultType, typename RhsLhsPolicyType, typename RhsRhsPolicyType,
                 template <typename, typename> class RhsOpType,
                 template <typename, typename> class OpType,
                 typename ResultType>
        struct BinaryExpressionEvaluator<ConstantExpressionPolicy<LhsResultType>,
                                        BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>,
                                        ResultType, OpType, BinaryNullOp,
                                        typename boost::enable_if
                                        <
                                            boost::mpl::or_
                                            <
                                                boost::mpl::and_
                                                <
                                                    boost::mpl::not_< typename AssociativeTraits<ConstantExpressionPolicy<LhsResultType>, OpType, BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >::IsAssociative>,
                                                    boost::is_same<LhsResultType, ResultType>,
                                                    boost::mpl::not_<boost::is_same<ResultType, typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>::ResultType> >
                                                >,
                                                boost::mpl::and_
                                                <
                                                    boost::mpl::not_< typename AssociativeTraits<ConstantExpressionPolicy<LhsResultType>, OpType, BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >::IsAssociative>,
                                                    boost::mpl::not_<IsCommutative<ConstantExpressionPolicy<LhsResultType>, OpType, BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> > >,
                                                    boost::mpl::not_<typename OpType<LhsResultType, typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>::ResultType>::HasOpLeftEqual>,
                                                    boost::is_same<LhsResultType, ResultType>
                                                >
                                            >
                                         >::type
                                         >
        {
            static void Eval(const Expression<ConstantExpressionPolicy<LhsResultType> >& lhs,
                             const Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >& rhs, 
                             Accumulator<ResultType>& accum)
            {
                typedef typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>::ResultType RhsResultType;
                
                lhs.Evaluate(accum);
                RhsResultType rhsTemp = rhs.Evaluate();
                OpType<LhsResultType, RhsResultType>::ApplyEqual(accum, rhsTemp);
            }
            
            #ifdef NEKTAR_UNIT_TESTS
            static const unsigned int ClassNum = 2;
            #endif //NEKTAR_UNIT_TESTS
        };


        // a - Expression is Associative.
        // b - Expression is Commutative.
        // c - Left Equal Operator is defined.
        // d - Constant Expression is result type.
        // e - Binary Expression is result type.
        // f - In A + (B-C), (A+B) is the result type.
        // Case 2 - Evaluate Rhs, Left Op equal lhs.
        // -0101- 001000
        // 001-1- 001000
        template<typename LhsResultType, typename RhsLhsPolicyType, typename RhsRhsPolicyType,
                 template <typename, typename> class RhsOpType,
                 template <typename, typename> class OpType,
                 typename ResultType>
        struct BinaryExpressionEvaluator<ConstantExpressionPolicy<LhsResultType>,
                                        BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>,
                                        ResultType, OpType, BinaryNullOp,
                                        typename boost::enable_if
                                        <
                                            boost::mpl::or_
                                            <
                                                boost::mpl::and_
                                                <
                                                    boost::mpl::not_<IsCommutative<ConstantExpressionPolicy<LhsResultType>, OpType, BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> > >,
                                                    typename OpType<LhsResultType, typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>::ResultType>::HasOpLeftEqual,
                                                    boost::mpl::not_<boost::is_same<LhsResultType, ResultType> >,
                                                    boost::is_same<ResultType, typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>::ResultType>
                                                >,
                                                boost::mpl::and_
                                                <
                                                    boost::mpl::not_< typename AssociativeTraits<ConstantExpressionPolicy<LhsResultType>, OpType, BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >::IsAssociative>,
                                                    typename OpType<LhsResultType, typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>::ResultType>::HasOpLeftEqual,
                                                    boost::mpl::not_<boost::is_same<LhsResultType, ResultType> >,
                                                    boost::is_same<ResultType, typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>::ResultType>
                                                >                                                
                                            >
                                         >::type
                                         >
        {
            static void Eval(const Expression<ConstantExpressionPolicy<LhsResultType> >& lhs,
                             const Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >& rhs, 
                             Accumulator<ResultType>& accum)
            {
                typedef typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>::ResultType RhsResultType;
                
                rhs.Evaluate(accum);
                OpType<LhsResultType, RhsResultType>::ApplyLeftEqual(accum, *lhs);
            }
            
            #ifdef NEKTAR_UNIT_TESTS
            static const unsigned int ClassNum = 3;
            #endif //NEKTAR_UNIT_TESTS
        };

        // a - Expression is Associative.
        // b - Expression is Commutative.
        // c - Left Equal Operator is defined.
        // d - Constant Expression is result type.
        // e - Binary Expression is result type.
        // f - In A + (B-C), (A+B) is the result type.
        // 
        // Case 3 - Swap operations.  Evaluate new lhs, op= new rhs.
        // .i 6
        // .o 6
        // .ilb a b c d e f
        // .ob c0 c1 c2 c3 c4 c5
        // .p 11
        // -1-01- 000100
        // 01--1- 000100
        template<typename LhsResultType, typename RhsLhsPolicyType, typename RhsRhsPolicyType,
                 template <typename, typename> class RhsOpType,
                 template <typename, typename> class OpType,
                 typename ResultType>
        struct BinaryExpressionEvaluator<ConstantExpressionPolicy<LhsResultType>,
                                        BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>,
                                        ResultType, OpType, BinaryNullOp,
                                        typename boost::enable_if
                                        <
                                            boost::mpl::or_
                                            <
                                                boost::mpl::and_
                                                <
                                                    IsCommutative<ConstantExpressionPolicy<LhsResultType>, OpType, BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >,
                                                    boost::mpl::not_<boost::is_same<ResultType, typename ConstantExpressionPolicy<LhsResultType>::ResultType> >,
                                                    boost::is_same<ResultType, typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>::ResultType> 
                                                >,
                                                boost::mpl::and_
                                                <
                                                    boost::mpl::not_< typename AssociativeTraits<ConstantExpressionPolicy<LhsResultType>, OpType, BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >::IsAssociative>,
                                                    IsCommutative<ConstantExpressionPolicy<LhsResultType>, OpType, BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >,
                                                    boost::is_same<ResultType, typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>::ResultType> 
                                                >
                                            >
                                         >::type
                                         >
        {
            static void Eval(const Expression<ConstantExpressionPolicy<LhsResultType> >& lhs,
                             const Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >& rhs, 
                             Accumulator<ResultType>& accum)
            {
                typedef typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>::ResultType RhsResultType;
                
                rhs.Evaluate(accum);
                OpType<RhsResultType, LhsResultType>::ApplyEqual(accum, *lhs);
            }
            
            #ifdef NEKTAR_UNIT_TESTS
            static const unsigned int ClassNum = 4;
            #endif //NEKTAR_UNIT_TESTS
        };


        // a - Expression is Associative.
        // b - Expression is Commutative.
        // c - Left Equal Operator is defined.
        // d - Constant Expression is result type.
        // e - Binary Expression is result type.
        // f - In A + (B-C), (A+B) is the result type.
        // g - In R = A + (B-C), R += B and R -= C is defined.
        // Case 4 - Evaluate lhs, evaluate rhs with parent op type.
        // 1--1--1 000010     
        template<typename LhsResultType, typename RhsLhsPolicyType, typename RhsRhsPolicyType,
                 template <typename, typename> class RhsOpType,
                 template <typename, typename> class OpType,
                 typename ResultType>
        struct BinaryExpressionEvaluator<ConstantExpressionPolicy<LhsResultType>,
                                        BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>,
                                        ResultType, OpType, BinaryNullOp,
                                        typename boost::enable_if
                                        <
                                            boost::mpl::and_
                                            <
                                                typename AssociativeTraits<ConstantExpressionPolicy<LhsResultType>, OpType, BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >::IsAssociative,
                                                boost::is_same<ResultType, typename ConstantExpressionPolicy<LhsResultType>::ResultType>,
                                                typename AssociativeTraits<ConstantExpressionPolicy<LhsResultType>, OpType, BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >::OpEqualsAreDefined
                                            >
                                         >::type
                                         >
        {
            static void Eval(const Expression<ConstantExpressionPolicy<LhsResultType> >& lhs,
                             const Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >& rhs, 
                             Accumulator<ResultType>& accum)
            {
                typedef typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>::ResultType RhsResultType;
                
                lhs.Evaluate(accum);
                rhs.template Evaluate<OpType>(accum);
            }
            
            #ifdef NEKTAR_UNIT_TESTS
            static const unsigned int ClassNum = 5;
            #endif //NEKTAR_UNIT_TESTS
        };


        // a - Expression is Associative.
        // b - Expression is Commutative.
        // c - Left Equal Operator is defined.
        // d - Constant Expression is result type.
        // e - Binary Expression is result type.
        // f - In A + (B-C), (A+B) is the result type.
        // 
        // Case 5 - Evaluate new (A+B) then op= C
        // 100011 000001
        template<typename LhsResultType, typename RhsLhsPolicyType, typename RhsRhsPolicyType,
                 template <typename, typename> class RhsOpType,
                 template <typename, typename> class OpType,
                 typename ResultType>
        struct BinaryExpressionEvaluator<ConstantExpressionPolicy<LhsResultType>,
                                        BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>,
                                        ResultType, OpType, BinaryNullOp,
                                        typename boost::enable_if
                                        <
                                            boost::mpl::and_
                                            <
                                                boost::mpl::and_
                                                <
                                                    typename AssociativeTraits<ConstantExpressionPolicy<LhsResultType>, OpType, BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >::IsAssociative,
                                                    boost::mpl::not_<IsCommutative<ConstantExpressionPolicy<LhsResultType>, OpType, BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> > >,
                                                    boost::mpl::not_<typename OpType<LhsResultType, typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>::ResultType>::HasOpLeftEqual>,
                                                    boost::mpl::not_<boost::is_same<ResultType, typename ConstantExpressionPolicy<LhsResultType>::ResultType> >,
                                                    boost::is_same<ResultType, typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>::ResultType>
                                                >,
                                                boost::is_same<ResultType, typename AssociativeTraits<ConstantExpressionPolicy<LhsResultType>, OpType, BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >::LhsAsBinaryExpressionResultType>
                                            >                                            
                                         >::type
                                         >
        {
            static void Eval(const Expression<ConstantExpressionPolicy<LhsResultType> >& lhs,
                             const Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >& rhs, 
                             Accumulator<ResultType>& accum)
            {
                typedef typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>::ResultType RhsResultType;
                
                // TODO.
            }
            
            #ifdef NEKTAR_UNIT_TESTS
            static const unsigned int ClassNum = 6;
            #endif //NEKTAR_UNIT_TESTS
        };


        ///////////////////////////////////////////////////////////////////////
        // Constant-Binary Specializations With Parent Op Types
        ///////////////////////////////////////////////////////////////////////
        //
        // These are expressions of the form R + (A + (B-C)), where R is our 
        // incoming accumulator, and + is the ParentOpType.  These specializations
        // are only valid if the expression is associative and R += A and 
        // R += typeof(B-C) is defined.  We don't check for these conditions 
        // explicitly since it should have already been handled elsewhere.
        //
        // There are only two cases, whether or not a temporary is needed for B-C.
        // This will only be necessary if the subexpression (A + (B-C)) is not 
        // associative.
        
        
        // Binary-Constant With Parent Op Type Specialization 0 
        // No temporaries required.
        template<typename LhsResultType, typename RhsLhsPolicyType, typename RhsRhsPolicyType,
                 template <typename, typename> class RhsOpType,
                 template <typename, typename> class OpType,
                 template <typename, typename> class ParentOpType,
                 typename ResultType>
        struct BinaryExpressionEvaluator<ConstantExpressionPolicy<LhsResultType>,
                                         BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>,
                                         ResultType, OpType, ParentOpType,
                                         typename boost::enable_if
                                         <
                                            boost::mpl::and_
                                            <
                                                typename AssociativeTraits<ConstantExpressionPolicy<LhsResultType>, OpType, BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >::IsAssociative,
                                                boost::mpl::not_<boost::is_same<BinaryNullOp<int, int>, ParentOpType<int, int> > >
                                            >
                                         >::type >
        {
            static void Eval(const Expression<ConstantExpressionPolicy<LhsResultType> >& lhs,
                             const Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >& rhs, 
                             Accumulator<ResultType>& accum)
            {
                typedef typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>::ResultType RhsResultType;
                typedef Expression<ConstantExpressionPolicy<LhsResultType> > LhsExpressionType;
                typedef Expression<ConstantExpressionPolicy<RhsResultType> > RhsExpressionType;

                ParentOpType<ResultType, LhsResultType>::ApplyEqual(accum, *lhs);
                typedef AssociativeTraits<ConstantExpressionPolicy<ResultType>, ParentOpType,
                                          BinaryExpressionPolicy<LhsExpressionType, OpType, RhsExpressionType> > TraitsType;

                typedef typename TraitsType::template LhsOpType<ResultType, typename RhsLhsPolicyType::ResultType, ParentOpType>::type NewLhsOpType;
                typedef typename TraitsType::template RhsOpType<ResultType, typename RhsRhsPolicyType::ResultType, ParentOpType>::type NewRhsOpType;
                
                NewLhsOpType::ApplyEqual(accum, (*rhs).first);
                NewRhsOpType::ApplyEqual(accum, (*rhs).second);
            }
        };
        
        // Binary-Constant With Parent Op Type Specialization 1 
        // Rhs needs a temporary.
        template<typename LhsResultType, typename RhsLhsPolicyType, typename RhsRhsPolicyType,
                 template <typename, typename> class RhsOpType,
                 template <typename, typename> class OpType,
                 template <typename, typename> class ParentOpType,
                 typename ResultType>
        struct BinaryExpressionEvaluator<ConstantExpressionPolicy<LhsResultType>,
                                         BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>,
                                         ResultType, OpType, ParentOpType,
                                         typename boost::enable_if
                                         <
                                            boost::mpl::and_
                                            <
                                                boost::mpl::not_<typename AssociativeTraits<ConstantExpressionPolicy<LhsResultType>, OpType, BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >::IsAssociative>,
                                                boost::mpl::not_<boost::is_same<BinaryNullOp<int, int>, ParentOpType<int, int> > >
                                            >
                                         >::type >
        {
            static void Eval(const Expression<ConstantExpressionPolicy<LhsResultType> >& lhs,
                             const Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >& rhs, 
                             Accumulator<ResultType>& accum)
            {
                typedef typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>::ResultType RhsResultType;
                typedef Expression<ConstantExpressionPolicy<LhsResultType> > LhsExpressionType;
                typedef Expression<ConstantExpressionPolicy<RhsResultType> > RhsExpressionType;

                ParentOpType<ResultType, LhsResultType>::ApplyEqual(accum, *lhs);
                
                RhsResultType rhsTemp = rhs.Evaluate();
                
                typedef AssociativeTraits<ConstantExpressionPolicy<ResultType>, ParentOpType,
                                          BinaryExpressionPolicy<LhsExpressionType, OpType, RhsExpressionType> > TraitsType;

                typedef typename TraitsType::template LhsOpType<ResultType, RhsResultType, ParentOpType>::type NewLhsOpType;
                
                NewLhsOpType::ApplyEqual(accum, rhsTemp);
            }
        };


                      
                      
        //////////////////////////////////////////////////
        // Binary-Binary expressions
        // 
        // R = (A+B) - (C*D)
        //
        // a - typeof(A+B) == typeof(R)
        // b - typeof(A+B) - (C*D) is associative.
        //
        ///////////////////////////////////////////////////////////
//        template<typename LhsLhsPolicyType, template<typename,typename> class LhsOpType, typename LhsRhsPolicyType,
//                 typename RhsLhsPolicyType, template<typename,typename> class RhsOpType, typename RhsRhsPolicyType,
//                 typename ResultType,
//                 template<typename,typename> class OpType>
//        struct BinaryExpressionEvaluator<BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType>,
//                                        BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>,
//                                        ResultType, OpType, BinaryNullOp,
//                                        typename boost::enable_if
//                                        <
//                                            boost::mpl::and_
//                                            <
//                                                boost::is_same<ResultType, typename BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType>::ResultType>,
//                                                typename AssociativeTraits
//                                                <
//                                                    BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType>,
//                                                    OpType,
//                                                    BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>
//                                                >::IsAssociative
//                                            >
//                                        >::type >
//        {
//            static void Eval(const Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType> >& lhs,
//                             const Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >& rhs,
//                             Accumulator<ResultType>& accum)
//            {
//                lhs.Evaluate(accum);
//                rhs.template Evaluate<OpType>(accum);
//            }
//
//        };
//        
//        template<typename LhsLhsPolicyType, template<typename,typename> class LhsOpType, typename LhsRhsPolicyType,
//                 typename RhsLhsPolicyType, template<typename,typename> class RhsOpType, typename RhsRhsPolicyType,
//                 typename ResultType,
//                 template<typename,typename> class OpType>
//        struct BinaryExpressionEvaluator<BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType>,
//                                        BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>,
//                                        ResultType, OpType, BinaryNullOp,
//                                        typename boost::enable_if
//                                        <
//                                            boost::mpl::not_<boost::is_same<ResultType, typename BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType>::ResultType> >
//                                        >::type >
//        {
//            static void Eval(const Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType> >& lhs,
//                             const Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >& rhs,
//                             Accumulator<ResultType>& accum)
//            {
//                typedef typename BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType>::ResultType LhsResultType;
//                typedef typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>::ResultType RhsResultType;
//                
//                LhsResultType lhsTemp = lhs.Evaluate();
//                RhsResultType rhsTemp = rhs.Evaluate();
//                OpType<LhsResultType, RhsResultType>::Apply(accum, lhs, rhs);
//            }
//
//        };
//
//        template<typename LhsLhsPolicyType, template<typename,typename> class LhsOpType, typename LhsRhsPolicyType,
//                 typename RhsLhsPolicyType, template<typename,typename> class RhsOpType, typename RhsRhsPolicyType,
//                 typename ResultType,
//                 template<typename,typename> class OpType>
//        struct BinaryExpressionEvaluator<BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType>,
//                                        BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>,
//                                        ResultType, OpType, BinaryNullOp,
//                                        typename boost::enable_if
//                                        <
//                                            boost::mpl::and_
//                                            <
//                                                boost::is_same<ResultType, typename BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType>::ResultType>,
//                                                boost::mpl::not_
//                                                <
//                                                    typename AssociativeTraits
//                                                    <
//                                                        BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType>,
//                                                        OpType,
//                                                        BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>
//                                                    >::IsAssociative
//                                                >
//                                            >
//                                        >::type >
//        {
//            static void Eval(const Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType> >& lhs,
//                             const Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >& rhs,
//                             Accumulator<ResultType>& accum)
//            {
//                typedef typename BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType>::ResultType LhsResultType;
//                typedef typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>::ResultType RhsResultType;
//            
//                lhs.Evaluate(accum);
//                RhsResultType rhsTemp = rhs.Evaluate();
//                OpType<LhsResultType, RhsResultType>::ApplyEqual(accum, rhsTemp);
//
//            }
//
//        };

        






    // Four cases where there is an incoming ParentOpType.
//    .i 4
//    .o 4
//    .ilb LhsIsResultType RhsIsResultType LhsIsAssociative RhsIsAssociative
//    .ob TwoTemps LhsTemp RhsTemp NoTemps
//    .p 9
//    --00 1000
//    0--0 1000
//    -00- 1000
//    00-- 1000
//    01-1 0100
//    -101 0100
//    101- 0010
//    1-10 0010
//    1111 0001
//    .e

    //    .ilb LhsIsResultType RhsIsResultType LhsIsAssociative RhsIsAssociative
    //    .ob TwoTemps LhsTemp RhsTemp NoTemps
    //    --00 1000
    //    0--0 1000
    //    -00- 1000
    //    00-- 1000
    template<typename LhsPolicyType, typename RhsPolicyType, typename ResultType, 
                 template <typename, typename> class OpType,
                 template <typename, typename> class ParentOpType>
    struct BinaryExpressionEvaluator<LhsPolicyType,
                                    RhsPolicyType,
                                    ResultType, OpType, ParentOpType,
                                    typename boost::enable_if
                                    <
                                        boost::mpl::and_
                                        <
                                            boost::mpl::and_
                                            <
                                                IsBinaryExpressionPolicy<LhsPolicyType>,
                                                IsBinaryExpressionPolicy<RhsPolicyType>
                                            >,
                                            boost::mpl::not_<boost::is_same<BinaryNullOp<int, int>, ParentOpType<int, int> > >,
                                            boost::mpl::or_
                                            <
                                                boost::mpl::and_
                                                <
                                                    boost::mpl::not_<boost::is_same<ResultType, typename LhsPolicyType::ResultType> >,
                                                    boost::mpl::not_<boost::is_same<ResultType, typename RhsPolicyType::ResultType> >
                                                >,
                                                boost::mpl::and_
                                                <
                                                     boost::mpl::not_<boost::is_same<ResultType, typename LhsPolicyType::ResultType> >,
                                                     boost::mpl::not_<typename AssociativeTraits<ConstantExpressionPolicy<ResultType>, ParentOpType, RhsPolicyType>::IsAssociative>
                                                >,
                                                boost::mpl::and_
                                                <
                                                    boost::mpl::not_<boost::is_same<ResultType, typename RhsPolicyType::ResultType> >,
                                                    boost::mpl::not_<typename AssociativeTraits<ConstantExpressionPolicy<ResultType>, ParentOpType, LhsPolicyType>::IsAssociative>
                                                >,
                                                boost::mpl::and_
                                                <
                                                    boost::mpl::not_<typename AssociativeTraits<ConstantExpressionPolicy<ResultType>, ParentOpType, RhsPolicyType>::IsAssociative>,
                                                    boost::mpl::not_<typename AssociativeTraits<ConstantExpressionPolicy<ResultType>, ParentOpType, LhsPolicyType>::IsAssociative>
                                                >
                                            >
                                        >
                                    >::type>
    {
        static void Eval(const Expression<LhsPolicyType>& lhs, 
                         const Expression<RhsPolicyType>& rhs,
                         Accumulator<ResultType>& result)
        {
            typedef typename LhsPolicyType::ResultType LhsResultType;
            typedef typename RhsPolicyType::ResultType RhsResultType;
            
            LhsResultType lhsTemp = lhs.Evaluate();
            RhsResultType rhsTemp = rhs.Evaluate();
            
            ParentOpType<ResultType, LhsResultType>::ApplyEqual(result, lhsTemp);
            
            typedef AssociativeTraits<ConstantExpressionPolicy<ResultType>, ParentOpType,
                                      BinaryExpressionPolicy<LhsPolicyType, OpType, RhsPolicyType> > TraitsType;

            typedef typename TraitsType::template RhsOpType<ResultType, RhsResultType, ParentOpType>::type NewRhsOpType;
            
            NewRhsOpType::ApplyEqual(result, rhsTemp);

        }
        
        #ifdef NEKTAR_UNIT_TESTS
        static const unsigned int ClassNum = 200;
        #endif //NEKTAR_UNIT_TESTS

    };
    
    //    .ilb LhsIsResultType RhsIsResultType LhsIsAssociative RhsIsAssociative
    //    .ob TwoTemps LhsTemp RhsTemp NoTemps
    //    01-1 0100
    //    -101 0100
    template<typename LhsPolicyType, typename RhsPolicyType, typename ResultType, 
                 template <typename, typename> class OpType,
                 template <typename, typename> class ParentOpType>
    struct BinaryExpressionEvaluator<LhsPolicyType,
                                    RhsPolicyType,
                                    ResultType, OpType, ParentOpType,
                                    typename boost::enable_if
                                    <
                                        boost::mpl::and_
                                        <
                                            boost::mpl::and_
                                            <
                                                IsBinaryExpressionPolicy<LhsPolicyType>,
                                                IsBinaryExpressionPolicy<RhsPolicyType>
                                            >,
                                            boost::mpl::not_<boost::is_same<BinaryNullOp<int, int>, ParentOpType<int, int> > >,
                                            boost::mpl::or_
                                            <
                                                boost::mpl::and_
                                                <
                                                    boost::mpl::not_<boost::is_same<ResultType, typename LhsPolicyType::ResultType> >,
                                                    boost::is_same<ResultType, typename RhsPolicyType::ResultType>,
                                                    typename AssociativeTraits<ConstantExpressionPolicy<ResultType>, ParentOpType, RhsPolicyType>::IsAssociative
                                                >,
                                                boost::mpl::and_
                                                <
                                                     boost::is_same<ResultType, typename RhsPolicyType::ResultType>,
                                                     boost::mpl::not_<typename AssociativeTraits<ConstantExpressionPolicy<ResultType>, ParentOpType, LhsPolicyType>::IsAssociative>,
                                                     typename AssociativeTraits<ConstantExpressionPolicy<ResultType>, ParentOpType, RhsPolicyType>::IsAssociative
                                                >
                                            >
                                        >
                                    >::type>
    {
        static void Eval(const Expression<LhsPolicyType>& lhs, 
                         const Expression<RhsPolicyType>& rhs,
                         Accumulator<ResultType>& result)
        {
            typedef typename LhsPolicyType::ResultType LhsResultType;
            typedef typename RhsPolicyType::ResultType RhsResultType;
            LhsResultType lhsTemp = lhs.Evaluate();
            
            ParentOpType<ResultType, LhsResultType>::ApplyEqual(result, lhsTemp);
            
            
            typedef AssociativeTraits<ConstantExpressionPolicy<ResultType>, ParentOpType,
                                      BinaryExpressionPolicy<LhsPolicyType, OpType, RhsPolicyType> > TraitsType;

            TraitsType::EvaluateNewRhs(rhs, result);
        }
        
        #ifdef NEKTAR_UNIT_TESTS
        static const unsigned int ClassNum = 201;
        #endif //NEKTAR_UNIT_TESTS
    };

//    .ilb LhsIsResultType RhsIsResultType LhsIsAssociative RhsIsAssociative
    //    .ob TwoTemps LhsTemp RhsTemp NoTemps
//    101- 0010
//    1-10 0010
    template<typename LhsPolicyType, typename RhsPolicyType, typename ResultType, 
                 template <typename, typename> class OpType,
                 template <typename, typename> class ParentOpType>
    struct BinaryExpressionEvaluator<LhsPolicyType,
                                    RhsPolicyType,
                                    ResultType, OpType, ParentOpType,
                                    typename boost::enable_if
                                    <
                                        boost::mpl::and_
                                        <
                                            boost::mpl::and_
                                            <
                                                IsBinaryExpressionPolicy<LhsPolicyType>,
                                                IsBinaryExpressionPolicy<RhsPolicyType>
                                            >,
                                            boost::mpl::not_<boost::is_same<BinaryNullOp<int, int>, ParentOpType<int, int> > >,
                                            boost::mpl::or_
                                            <
                                                boost::mpl::and_
                                                <
                                                    boost::is_same<ResultType, typename LhsPolicyType::ResultType>,
                                                    boost::mpl::not_<boost::is_same<ResultType, typename RhsPolicyType::ResultType> >,
                                                    typename AssociativeTraits<ConstantExpressionPolicy<ResultType>, ParentOpType, LhsPolicyType>::IsAssociative
                                                >,
                                                boost::mpl::and_
                                                <
                                                     boost::is_same<ResultType, typename LhsPolicyType::ResultType>,
                                                     boost::mpl::not_<typename AssociativeTraits<ConstantExpressionPolicy<ResultType>, ParentOpType, RhsPolicyType>::IsAssociative>,
                                                     typename AssociativeTraits<ConstantExpressionPolicy<ResultType>, ParentOpType, LhsPolicyType>::IsAssociative
                                                >
                                            >
                                        >
                                    >::type>
    {
        static void Eval(const Expression<LhsPolicyType>& lhs, 
                         const Expression<RhsPolicyType>& rhs,
                         Accumulator<ResultType>& result)
        {
            typedef typename LhsPolicyType::ResultType LhsResultType;
            typedef typename RhsPolicyType::ResultType RhsResultType;
            lhs.template Evaluate<ParentOpType>(result);
            
            RhsResultType rhsTemp = rhs.Evaluate();
                      
            typedef AssociativeTraits<ConstantExpressionPolicy<ResultType>, ParentOpType,
                                      BinaryExpressionPolicy<LhsPolicyType, OpType, RhsPolicyType> > TraitsType;

            typedef typename TraitsType::template RhsOpType<ResultType, typename RhsPolicyType::ResultType, ParentOpType>::type NewRhsOpType;
            NewRhsOpType::ApplyEqual(result, rhsTemp);
        }
    };

//    .ilb LhsIsResultType RhsIsResultType LhsIsAssociative RhsIsAssociative
    //    .ob TwoTemps LhsTemp RhsTemp NoTemps
//    1111 0001
    template<typename LhsPolicyType, typename RhsPolicyType, typename ResultType, 
                 template <typename, typename> class OpType,
                 template <typename, typename> class ParentOpType>
    struct BinaryExpressionEvaluator<LhsPolicyType,
                                    RhsPolicyType,
                                    ResultType, OpType, ParentOpType,
                                    typename boost::enable_if
                                    <
                                        boost::mpl::and_
                                        <
                                            boost::mpl::and_
                                            <
                                                IsBinaryExpressionPolicy<LhsPolicyType>,
                                                IsBinaryExpressionPolicy<RhsPolicyType>
                                            >,
                                            boost::mpl::not_<boost::is_same<BinaryNullOp<int, int>, ParentOpType<int, int> > >,
                                            boost::mpl::and_
                                            <
                                                boost::is_same<ResultType, typename LhsPolicyType::ResultType>,
                                                boost::is_same<ResultType, typename RhsPolicyType::ResultType>,
                                                typename AssociativeTraits<ConstantExpressionPolicy<ResultType>, ParentOpType, LhsPolicyType>::IsAssociative,
                                                typename AssociativeTraits<ConstantExpressionPolicy<ResultType>, ParentOpType, RhsPolicyType>::IsAssociative
                                            >
                                        >
                                    >::type>
    {
        static void Eval(const Expression<LhsPolicyType>& lhs, 
                         const Expression<RhsPolicyType>& rhs,
                         Accumulator<ResultType>& result)
        {
            lhs.template Evaluate<ParentOpType>(result);
            
            typedef AssociativeTraits<ConstantExpressionPolicy<ResultType>, ParentOpType,
                                      BinaryExpressionPolicy<LhsPolicyType, OpType, RhsPolicyType> > TraitsType;

            TraitsType::EvaluateNewRhs(rhs, result);
        }
    };
}

#endif //NEKTAR_USE_EXPRESSION_TEMPLATES
#endif //NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_BINARY_EXPRESSION_EVALUATOR_HPP

