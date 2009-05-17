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

#ifndef NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_BINARY_EXPRESSION_EVALUATOR_BINARY_CONSTANT_HPP
#define NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_BINARY_EXPRESSION_EVALUATOR_BINARY_CONSTANT_HPP

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
            struct BinaryConstantBinaryExpressionEvaluator;
            
            ///////////////////////////////////////////////////////////////////////
            // Binary-Constant Specializations : (A+B) - C
            ///////////////////////////////////////////////////////////////////////
            // a - expression is associative.
            // b - (A+B) is the result type
            //
            // Case 1 - Evaluate LHS, the -= rhs.
            // Case 2 - Create a temporary for the lhs.
            // Case 3 - Restructure to A + (B-C) and pass to the Constant-Binary specialization.
            
            // .i 2
            // .o 3
            // .ilb a b 
            // .ob c0 c1 c2 
            // 00 010
            // -1 100
            // 10 001
            
            // Case 1 - Evaluate LHS, then -= rhs
            // a - expression is associative.
            // b - (A+B) is the result type
            // -1 100
            template<typename LhsLhsPolicyType, template<typename, typename> class LhsOpType, typename LhsRhsPolicyType,
                     typename RhsResultType, typename ResultType,
                     template<typename, typename> class OpType>
            struct BinaryConstantBinaryExpressionEvaluator<BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType>,
                                            ConstantExpressionPolicy<RhsResultType>,
                                            ResultType, OpType, BinaryNullOp,
                                            typename boost::enable_if
                                            <
                                                boost::is_same<typename BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType>::ResultType, ResultType>
                                            >::type >
            {
                static void Eval(const Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType> >& lhs,
                                 const Expression<ConstantExpressionPolicy<RhsResultType> >& rhs,
                                 Accumulator<ResultType>& accum)
                {
                    typedef typename BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType>::ResultType LhsResultType;
                    lhs.Evaluate(accum);
                    OpType<LhsResultType, RhsResultType>::ApplyEqual(accum, *rhs);
                }
                
                #ifdef NEKTAR_UNIT_TESTS
                static const unsigned int ClassNum = 6;
                #endif //NEKTAR_UNIT_TESTS
            };

            // Case 2 - Create a temporary for the lhs.
            // a - expression is associative.
            // b - (A+B) is the result type
            // 00 010
            template<typename LhsLhsPolicyType, template<typename, typename> class LhsOpType, typename LhsRhsPolicyType,
                     typename RhsResultType, typename ResultType,
                     template<typename, typename> class OpType>
            struct BinaryConstantBinaryExpressionEvaluator<BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType>,
                                            ConstantExpressionPolicy<RhsResultType>,
                                            ResultType, OpType, BinaryNullOp,
                                            typename boost::enable_if
                                            <
                                                boost::mpl::and_
                                                <
                                                    boost::mpl::not_<typename AssociativeTraits<BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType>, OpType, ConstantExpressionPolicy<RhsResultType> >::IsAssociative>,
                                                    boost::mpl::not_<boost::is_same<typename BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType>::ResultType, ResultType> >
                                                >
                                            >::type >
            {
                static void Eval(const Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType> >& lhs,
                                 const Expression<ConstantExpressionPolicy<RhsResultType> >& rhs,
                                 Accumulator<ResultType>& accum)
                {
                    typedef typename BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType>::ResultType LhsResultType;
                    
                    LhsResultType lhsTemp = lhs.Evaluate();
                    OpType<LhsResultType, RhsResultType>::Apply(accum, lhsTemp, *rhs);
                }
                
                #ifdef NEKTAR_UNIT_TESTS
                static const unsigned int ClassNum = 7;
                #endif //NEKTAR_UNIT_TESTS
            };
                      
            // Case 3 - Restructure to A + (B-C) and pass to the Constant-Binary specialization.
            // a - expression is associative.
            // b - (A+B) is the result type
            // .ilb a b 
            // 10 001
            template<typename LhsLhsPolicyType, template<typename, typename> class LhsOpType, typename LhsRhsPolicyType,
                     typename RhsResultType, typename ResultType,
                     template<typename, typename> class OpType>
            struct BinaryConstantBinaryExpressionEvaluator<BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType>,
                                            ConstantExpressionPolicy<RhsResultType>,
                                            ResultType, OpType, BinaryNullOp,
                                            typename boost::enable_if
                                            <
                                                boost::mpl::and_
                                                <
                                                    typename AssociativeTraits<BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType>, OpType, ConstantExpressionPolicy<RhsResultType> >::IsAssociative,
                                                    boost::mpl::not_<boost::is_same<typename BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType>::ResultType, ResultType> >
                                                >
                                            >::type >
            {
                static void Eval(const Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType> >& lhs,
                                 const Expression<ConstantExpressionPolicy<RhsResultType> >& rhs,
                                 Accumulator<ResultType>& accum)
                {
                    typedef AssociativeTraits<BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType>, OpType, ConstantExpressionPolicy<RhsResultType> > Traits;
                    typedef typename Traits::AdjustedExpressionEvaluatorType AdjustedExpressionEvaluatorType;
                    typedef typename Traits::AdjustedLhsType AdjustedLhsType;
                    typedef typename Traits::AdjustedRhsType AdjustedRhsType;
                    
                    AdjustedLhsType newLhs = Traits::CreateNewLhs(lhs, rhs);
                    AdjustedRhsType newRhs = Traits::CreateNewRhs(lhs, rhs);
                    AdjustedExpressionEvaluatorType::Eval(newLhs, newRhs, accum);
                }
                
                #ifdef NEKTAR_UNIT_TESTS
                static const unsigned int ClassNum = 8;
                #endif //NEKTAR_UNIT_TESTS
            };
        }
        
        template<typename LhsLhsPolicyType, template<typename, typename> class LhsOpType, typename LhsRhsPolicyType,
                     typename RhsResultType, typename ResultType,
                     template<typename, typename> class OpType>
        struct BinaryExpressionEvaluator<BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType>,
                                        ConstantExpressionPolicy<RhsResultType>,
                                        ResultType, OpType, BinaryNullOp >
        {
            static void Eval(const Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType> >& lhs,
                             const Expression<ConstantExpressionPolicy<RhsResultType> >& rhs,
                             Accumulator<ResultType>& accum)
            {
                return Impl::BinaryConstantBinaryExpressionEvaluator<BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType>,
                                        ConstantExpressionPolicy<RhsResultType>,
                                        ResultType, OpType, BinaryNullOp>::Eval(lhs, rhs, accum);
            }
            
            #ifdef NEKTAR_UNIT_TESTS
            static const unsigned int ClassNum = 8;
            #endif //NEKTAR_UNIT_TESTS
        };

}

#endif //NEKTAR_USE_EXPRESSION_TEMPLATES
#endif //NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_BINARY_EXPRESSION_EVALUATOR_BINARY_CONSTANT_HPP
