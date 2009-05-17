///////////////////////////////////////////////////////////////////////////////
//
// File: BinaryExpressionEvaluator_BinaryConstant.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_BINARY_EXPRESSION_EVALUATOR_BINARY_BINARY_HPP
#define NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_BINARY_EXPRESSION_EVALUATOR_BINARY_BINARY_HPP

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

#include <LibUtilities/ExpressionTemplates/BinaryExpressionEvaluatorFwd.hpp>
#include <LibUtilities/ExpressionTemplates/ConstantExpression.hpp>
#include <LibUtilities/ExpressionTemplates/BinaryExpressionPolicyFwd.hpp>
#include <LibUtilities/ExpressionTemplates/AssociativeTraits.hpp>
#include <LibUtilities/ExpressionTemplates/CommutativeTraits.hpp>
#include <LibUtilities/ExpressionTemplates/AssignableTraits.hpp>
#include <LibUtilities/ExpressionTemplates/BinaryOperators.hpp>

namespace Nektar
{
    namespace Impl
    {
        template<typename LhsExpressionPolicyType, typename RhsExpressionPolicyType, typename ResultType, 
            template <typename, typename> class OpType, 
            template <typename, typename> class ParentOpType = BinaryNullOp,
            typename enabled = void>
        struct BinaryBinaryExpressionEvaluator;
        
        // Binary Expression Evaluator with no parent op type.
        // L Op R
        //
        // Input 0 - Lhs has different type
        // Input 1 - Rhs has different type
        // Input 2 - L Op R is commutative.
        // Input 3 - L Op R is associative.
        // Input 4 - R Op L is associative.
        // 
        // Case 0 - Temps for both sides.
        // Case 1 - Evaluate lhs, then Op = rhs.
        // Case 2 - Evaluate rhs, then op = lhs.
        // Case 3 - Eval lhs, temp rhs, then op= temp
        // Case 4 - Eval rhs, temp lhs, then op= temp
        
        //.i 5
        //.o 5
        //.ilb LhsIsResultType RhsIsResultType IsCommutative IsLeftAssociative IsRightAss
        //.ob TwoTemps EvalLhsTempRhs EvalRhsTempLhs EvalLhsParentOpRhs EvalRhsParentOpLhs
        //00--- 10000
        //0100- 10000
        //0101- 10000
        //0110- 00100
        //0111- 00001
        //1000- 01000
        //1001- 00010
        //1010- 01000
        //1011- 00010
        //1100- 01000
        //1101- 00010
        //1110- 01000
        //11110 00001
        //11111 00010
        //.e
        
        // Output is
    //    .i 4
    //    .o 5
    //    .ilb LhsIsResultType RhsIsResultType IsCommutative IsAssociative
    //    .ob TwoTemps EvalLhsTempRhs EvalRhsTempLhs EvalLhsParentOpRhs EvalRhsParentOpLhs
        //0-0-- 10000
        //00--- 10000
        //10-0- 01000
        //1--00 01000
        //1-00- 01000
        //011-0 00100
        //1--1- 00010
        //011-1 00001
        //-1101 00001
        //    .e


        // Specialization 0 - Need to create two temporaries
        //
        //    .ilb LhsIsResultType RhsIsResultType IsCommutative IsAssociative
        //    0-0-- 10000
        //    00--- 10000
        template<typename LhsPolicyType, typename RhsPolicyType, typename ResultType, 
                     template <typename, typename> class OpType>
        struct BinaryBinaryExpressionEvaluator<LhsPolicyType,
                                        RhsPolicyType,
                                        ResultType, OpType, BinaryNullOp,
                                        typename boost::enable_if
                                        <
                                            boost::mpl::and_
                                            <
                                                boost::mpl::and_
                                                <
                                                    IsBinaryExpressionPolicy<LhsPolicyType>,
                                                    IsBinaryExpressionPolicy<RhsPolicyType>
                                                    //boost::mpl::not_<IsConstantExpressionPolicy<LhsPolicyType> >,
                                                    //boost::mpl::not_<IsConstantExpressionPolicy<RhsPolicyType> >
                                                >,
                                                boost::mpl::or_
                                                <
                                                    boost::mpl::and_
                                                    <
                                                        boost::mpl::not_<boost::is_same<ResultType, typename LhsPolicyType::ResultType> >,
                                                        boost::mpl::not_<IsCommutative<LhsPolicyType, OpType, RhsPolicyType> >
                                                    >,
                                                    boost::mpl::and_
                                                    <
                                                        boost::mpl::not_<boost::is_same<ResultType, typename LhsPolicyType::ResultType> >,
                                                        boost::mpl::not_<boost::is_same<ResultType, typename RhsPolicyType::ResultType> >
                                                    >
                                                >
                                            >
                                        >::type>
        {
            static void Eval(const Expression<LhsPolicyType>& lhs, 
                             const Expression<RhsPolicyType>& rhs,
                             Accumulator<ResultType>& result)
            {
                typedef typename LhsPolicyType::ResultType LhsType;
                typedef typename RhsPolicyType::ResultType RhsType;
                
                LhsType lhsTemp;
                lhs.Evaluate(lhsTemp);
                
                RhsType rhsTemp;
                rhs.Evaluate(rhsTemp);

                OpType<LhsType, RhsType>::Apply(result, lhsTemp, rhsTemp);
            }
            
            #ifdef NEKTAR_UNIT_TESTS
            static const unsigned int ClassNum = 301;
            #endif //NEKTAR_UNIT_TESTS
        };


        // Specialization 1 - Evaluate the lhs with the accumulator, but create a 
        // temporary for the rhs.
        //
        //    .ilb LhsIsResultType RhsIsResultType IsCommutative IsAssociative
        //10-0- 01000
        //1--00 01000
        //1-00- 01000
        template<typename LhsPolicyType, typename RhsPolicyType, typename ResultType, 
                     template <typename, typename> class OpType>
        struct BinaryBinaryExpressionEvaluator<LhsPolicyType,
                                        RhsPolicyType,
                                        ResultType, OpType, BinaryNullOp,
                                        typename boost::enable_if
                                        <
                                            boost::mpl::and_
                                            <
                                                boost::mpl::and_
                                                <
                                                    IsBinaryExpressionPolicy<LhsPolicyType>,
                                                    IsBinaryExpressionPolicy<RhsPolicyType>
                                                    //boost::mpl::not_<IsConstantExpressionPolicy<LhsPolicyType> >,
                                                    //boost::mpl::not_<IsConstantExpressionPolicy<RhsPolicyType> >
                                                >,
                                                boost::mpl::or_
                                                <
                                                    boost::mpl::and_
                                                    <
                                                        boost::is_same<ResultType, typename LhsPolicyType::ResultType>,
                                                        boost::mpl::not_<boost::is_same<ResultType, typename RhsPolicyType::ResultType> >,
                                                        boost::mpl::not_<typename AssociativeTraits<LhsPolicyType, OpType, RhsPolicyType>::IsLeftAssociative>
                                                    >, 
                                                    boost::mpl::and_
                                                    <
                                                        boost::is_same<ResultType, typename LhsPolicyType::ResultType>, 
                                                        boost::mpl::not_<typename AssociativeTraits<LhsPolicyType, OpType, RhsPolicyType>::IsLeftAssociative>,
                                                        boost::mpl::not_<typename AssociativeTraits<LhsPolicyType, OpType, RhsPolicyType>::IsRightAssociative>
                                                    >,
                                                    boost::mpl::and_
                                                    <
                                                        boost::is_same<ResultType, typename LhsPolicyType::ResultType>, 
                                                        boost::mpl::not_<IsCommutative<LhsPolicyType, OpType, RhsPolicyType> >,
                                                        boost::mpl::not_<typename AssociativeTraits<LhsPolicyType, OpType, RhsPolicyType>::IsLeftAssociative>
                                                    >
                                                >
                                            >
                                        >::type>
        {
            static void Eval(const Expression<LhsPolicyType>& lhs, 
                             const Expression<RhsPolicyType>& rhs,
                             Accumulator<ResultType>& result)
            {
                typedef typename RhsPolicyType::ResultType RhsType;
                
                lhs.Evaluate(result);
                RhsType rhsTemp = rhs.Evaluate();

                OpType<ResultType, RhsType>::ApplyEqual(result, rhsTemp);
            }
        };
        
        // Specialization 2 - Evaluate the rhs, create a temporary for the left, 
        // then op= the result.  This is possible because the expression is 
        // commutative.
        //
        //    .ilb LhsIsResultType RhsIsResultType IsCommutative IsAssociative
        //011-0 00100
        template<typename LhsPolicyType, typename RhsPolicyType, typename ResultType, 
                     template <typename, typename> class OpType>
        struct BinaryBinaryExpressionEvaluator<LhsPolicyType,
                                        RhsPolicyType,
                                        ResultType, OpType, BinaryNullOp,
                                        typename boost::enable_if
                                        <
                                            boost::mpl::and_
                                            <
                                                boost::mpl::and_
                                                <
                                                    IsBinaryExpressionPolicy<LhsPolicyType>,
                                                    IsBinaryExpressionPolicy<RhsPolicyType>
                                                >,
                                                boost::mpl::and_
                                                <
                                                    boost::mpl::not_<boost::is_same<ResultType, typename LhsPolicyType::ResultType> >,
                                                    boost::is_same<ResultType, typename RhsPolicyType::ResultType>,
                                                    IsCommutative<LhsPolicyType, OpType, RhsPolicyType>,
                                                    boost::mpl::not_<typename AssociativeTraits<LhsPolicyType, OpType, RhsPolicyType>::IsRightAssociative>
                                                >
                                            >
                                        >::type>
        {
            static void Eval(const Expression<LhsPolicyType>& lhs, 
                             const Expression<RhsPolicyType>& rhs,
                             Accumulator<ResultType>& result)
            {
                typedef typename LhsPolicyType::ResultType LhsType;
                
                rhs.Evaluate(result);
                LhsType lhsTemp = lhs.Evaluate();

                OpType<ResultType, LhsType>::ApplyEqual(result, lhsTemp);
            }
        };
        



        // Specialization 3 - Evaluate lhs, then op= rhs
        //
        //    .ilb LhsIsResultType RhsIsResultType IsCommutative IsAssociative
        //1--1- 00010
        template<typename LhsPolicyType, typename RhsPolicyType, typename ResultType, 
                     template <typename, typename> class OpType>
        struct BinaryBinaryExpressionEvaluator<LhsPolicyType,
                                        RhsPolicyType,
                                        ResultType, OpType, BinaryNullOp,
                                        typename boost::enable_if
                                        <
                                            boost::mpl::and_
                                            <
                                                boost::mpl::and_
                                                <
                                                    IsBinaryExpressionPolicy<LhsPolicyType>,
                                                    IsBinaryExpressionPolicy<RhsPolicyType>
                                                >,
                                                boost::mpl::and_
                                                <
                                                    boost::is_same<ResultType, typename LhsPolicyType::ResultType>,
                                                    typename AssociativeTraits<LhsPolicyType, OpType, RhsPolicyType>::IsLeftAssociative
                                                >
                                            >
                                        >::type>
        {
            static void Eval(const Expression<LhsPolicyType>& lhs, 
                             const Expression<RhsPolicyType>& rhs,
                             Accumulator<ResultType>& result)
            {
                lhs.Evaluate(result);
                rhs.template Evaluate<OpType>(result);
                
            }
        };

        // Specialization 4 - Evaluate lhs, then op= rhs
        //
        //    .ilb LhsIsResultType RhsIsResultType IsCommutative IsAssociative
        //-1101 00001
        template<typename LhsPolicyType, typename RhsPolicyType, typename ResultType, 
                     template <typename, typename> class OpType>
        struct BinaryBinaryExpressionEvaluator<LhsPolicyType,
                                        RhsPolicyType,
                                        ResultType, OpType, BinaryNullOp,
                                        typename boost::enable_if
                                        <
                                            boost::mpl::and_
                                            <   
                                                boost::mpl::and_
                                                <
                                                    IsBinaryExpressionPolicy<LhsPolicyType>,
                                                    IsBinaryExpressionPolicy<RhsPolicyType>
                                                >,
                                                boost::mpl::and_
                                                <
                                                    boost::is_same<ResultType, typename RhsPolicyType::ResultType>,
                                                    IsCommutative<LhsPolicyType, OpType, RhsPolicyType>,
                                                    boost::mpl::not_<typename AssociativeTraits<LhsPolicyType, OpType, RhsPolicyType>::IsLeftAssociative>,
                                                    typename AssociativeTraits<LhsPolicyType, OpType, RhsPolicyType>::IsRightAssociative
                                                >
                                            >
                                        >::type>
        {
            static void Eval(const Expression<LhsPolicyType>& lhs, 
                             const Expression<RhsPolicyType>& rhs,
                             Accumulator<ResultType>& result)
            {
                rhs.Evaluate(result);
                lhs.template Evaluate<OpType>(result);
                
            }
        };
    }

    template<typename LhsLhsPolicyType, typename LhsRhsPolicyType, template<typename, typename> class LhsOpType,
         typename RhsLhsPolicyType, typename RhsRhsPolicyType, template<typename, typename> class RhsOpType,
         typename ResultType, template<typename, typename> class OpType>
    struct BinaryExpressionEvaluator<BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType>,
                                 BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>,
                                 ResultType, OpType, BinaryNullOp>
    {
        typedef BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType> LhsType;
        typedef BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> RhsType;
        
        static void Eval(const Expression<LhsType>& lhs, const Expression<RhsType>& rhs, 
                         Accumulator<ResultType>& result)
        {
            return Impl::BinaryBinaryExpressionEvaluator<LhsType, RhsType, ResultType, OpType, BinaryNullOp>::Eval(lhs, rhs, result);
        }
    };
    
}

#endif //NEKTAR_USE_EXPRESSION_TEMPLATES
#endif //NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_BINARY_EXPRESSION_EVALUATOR_BINARY_BINARY_HPP
