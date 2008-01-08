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

#include <boost/utility/enable_if.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/not.hpp>

namespace Nektar
{

        
        ////////////////////////////////////////////////
        // 1 op Specializations
        ////////////////////////////////////////////////
        
        /// \brief Constant-Constant specialization 1
        /// 
        /// With no incoming parent type, this evaluation is simply A op B, such 
        /// as A+B or A-B.  We'll just call the appropriate NekOp function which 
        /// takes 3 parameters.
        template<typename LhsDataType, typename RhsDataType, typename ResultType, 
                 template <typename, typename> class OpType>
        class BinaryExpressionEvaluator<ConstantExpressionPolicy<LhsDataType>,
                                        ConstantExpressionPolicy<RhsDataType>,
                                        ResultType, OpType, BinaryNullOp>
        {
            public:
            
                static void Eval(const Expression<ConstantExpressionPolicy<LhsDataType> >& lhs, 
                                 const Expression<ConstantExpressionPolicy<RhsDataType> >& rhs,
                                 Accumulator<ResultType>& result)
                {
                    OpType<LhsDataType, RhsDataType>::Apply(result, *lhs, *rhs);
                }
        };
           
        // R = -A + B
        //
        // Case 1 - typeof(A) == typeof(R)
        //          Pass accumulator to left side, then += rhs.
        //
        // Case 2 - typeof(A) != typeof(R)  
        //          typeof(A) + typeof(B) is commutative.
        //          - is convertible to a binary operation.
        //          Create new expression (B-A) and pass to constant expression evaluator.
        //
        // Case 3 - Create a temporary for the left, and then eval total expression.
        
        /// \brief Unary-Constant Specialization 1
        ///
    
//        class BinaryExpressionEvaluator<UnaryExpressionPolicy<LhsPolicyType, LhsOpType>,
//                                        ConstantExpressionPolicy<RhdDataType>,
//                                        ResultType, OpType, BinaryNulOp>
//        {
//            public:
//                typedef Expression<UnaryExpressionPolicy<LhsPolicyType, LhsOpType> > LhsExpressionType;
//                typedef typename LhsExpressionType::ResultType LhsResultType;
//                
//                static void Eval(const LhsExpressionType& lhs,
//                                 const Expression<ConstantExpressionPolicy<RhsDataType> >& rhs,
//                                 Accumulator<ResultTpye>& result)
//                {
//                    LhsResultType t = lhs;
//                    OpType<LhsResultType, RhsDataType
//                }
//        };
        
        // Temporary
        //template<typename C, typename LhsLhsPolicyType, typename LhsRhsPolicyType, typename R,
        //         template <typename, typename> class LhsOpType, 
        //         template <typename, typename> class OpType>
        //class BinaryExpressionEvaluator<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType>,
        //                                ConstantExpressionPolicy<C>,
        //                                R, OpType, BinaryNullOp, void>
        //{
        //    public:
        //        static void Eval(const Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >& lhs,
        //                         const Expression<ConstantExpressionPolicy<C> >& rhs, 
        //                         Accumulator<R>& accum)
        //        {
        //            typedef typename BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType>::ResultType LhsResultType;
        //            lhs.template Evaluate<BinaryNullOp>(accum);
        //            OpType<LhsResultType, C>::ApplyEqual(accum, *rhs);
        //        }
        //};


        /// \brief Constant-Constant specialization 2
        /// 
        /// This case occurs as part of a larger expression, such as 
        /// A + (B+C).  By definition, since there is a parent op type, it can't happen in isolation.
        ///
        /// For accumulator A, the following actions are performed:
        /// A op= lhs
        /// A op= rhs
        ///
        /// Change of operators due to associativity is detected.
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
//         template<typename LhsDataType, typename RhsDataType, typename ResultType, 
//                  typename <typename, typename> class OpType,
//                  template <typename, typename> class ParentOpType>
//         class BinaryExpressionEvaluator<ConstantExpressionPolicy<LhsDataType>,
//                                         ConstantExpressionPolicy<RhsDataType>,
//                                         ResultType, OpType, ParentOpType,
//                                         typename boost::enable_if_c
//                                         <
//                                                 ParentOpType<ResultType, LhsDataType>::HasOpEqual &&
//                                                 OpType<ResultType, RhsDataType>::HasOpEqual &&
//                                                 AssociativeTraits<ResultType, ParentOpType, LhsDataType, OpType, RhsDataType>::IsAssociative &&
//                                                 !boost::is_same<BinaryNullOp<void, void>, typename ParentOpType<LhsDataType, RhsDataType>::template Rebind<void, void>::type >::value,
//                                                 void
//                                         >::type >
//         {
//             public:
//                 
//                 static void Eval(const Expression<ConstantExpressionPolicy<LhsDataType> >& lhs, 
//                                 const Expression<ConstantExpressionPolicy<RhsDataType> >& rhs,
//                                 Accumulator<ResultType>& accum)
//                 {
//                     ParentOpType<ResultType, LhsDataType>::ApplyEqual(accum, *lhs);
// 
//                     typedef AssociativeTraits<ResultType, ParentOpType, LhsDataType, OpType, RhsDataType> AssociativeTraits;
//                     typedef typename ChooseOpChangeType<AssociativeTraits, OpType<LhsDataType, RhsDataType> >::Type OpChangeType;
//                     //typedef typename AssociativeTraits<ResultType, ParentOpType, LhsDataType, OpType, RhsDataType>::OpChangeType OpChangeType;
//                     OpChangeType::ApplyEqual(accum, *rhs);
//                 }
//         };
        
        // The only other case is where the OpEquals are not defined, but this is such a bizarre case that
        // I imagine if they aren't defined they are more likely user errors than a true part of the math.
        // So I will leave it undefined to help catch compiler errors.  
        // TODO - It may be worth defining it just so I can put a boost static assert to aid debugging.


        ////////////////////////////////////////////////////////////////////////
        // Constant-Binary Expressions
        //
        // The main difficulty with these expressions comes from needing to 
        // evaluate the constant.
        // first.  If we pass the accumulator to the left side first, then we 
        // must create a temporary 
        // on the right side unless the expression is associative.  If the expression 
        // is associative 
        // then there is no problem passing the accumulator with a ParentOpType.
        //
        // Another option occurs if the expression is commutative and the binary
        // expression has the same 
        // type as the accumulator.  In this case, pass the accumulator to the 
        // right side and then 
        // op= the left side.
        //
        // A last condition assumes the user has defined a left equal operator 
        // if the expression is 
        // not commutative nor associative.
        //
        // Therefore, the variables of interest are:
        // a - Expression is Associative.
        // b - Expression is Commutative.
        // c - Left Equal Operator is defined.
        // d - Constant Expression is result type.
        // e - Binary Expression is result type.
        // f - In A + (B-C), (A+B) is the result type.
        // 
        // Case 0 - Both sides need temporaries.  
        // Case 1 - Rhs needs temporary.
        // Case 2 - Evaluate Rhs, Left Op equal lhs.
        // Case 3 - Swap operations.  Evaluate new lhs, op= new rhs.
        // Case 4 - Evaluate lhs, evaluate rhs with parent op type.
        // Case 5 - Evaluate new (A+B) then op= C
        
        // .i 6
        // .o 6
        // .ilb a b c d e f
        // .ob c0 c1 c2 c3 c4 c5
        // .p 11
        // -000-0 100000
        // 0000-- 100000
        // ---00- 100000
        // 0--10- 010000
        // 0001-- 010000
        // -0101- 001000
        // 001-1- 001000
        // -1-01- 000100
        // 01--1- 000100
        // 1--1-- 000010       
        // 100011 000001
        // .e
        ////////////////////////////////////////////////////////////////////////
        
        
        // Constant-Binary specialization 0 - Both sides need temporaries.
        // Therefore, the variables of interest are:
        // a - Expression is Associative.
        // b - Expression is Commutative.
        // c - Left Equal Operator is defined.
        // d - Constant Expression is result type.
        // e - Binary Expression is result type.
        // f - In A + (B-C), (A+B) is the result type.
        //
        // -000-0 100000
        // 0000-- 100000
        // ---00- 100000
        template<typename LhsResultType, typename RhsLhsPolicyType, typename RhsRhsPolicyType,
                 template <typename, typename> class RhsOpType,
                 template <typename, typename> class OpType,
                 typename ResultType>
        class BinaryExpressionEvaluator<ConstantExpressionPolicy<LhsResultType>,
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
                                                >
                                            >
                                         >::type
                                         >
        {
            public:
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
        template<typename LhsResultType, typename RhsLhsPolicyType, typename RhsRhsPolicyType,
                 template <typename, typename> class RhsOpType,
                 template <typename, typename> class OpType,
                 typename ResultType>
        class BinaryExpressionEvaluator<ConstantExpressionPolicy<LhsResultType>,
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
            public:
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
        class BinaryExpressionEvaluator<ConstantExpressionPolicy<LhsResultType>,
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
            public:
                static void Eval(const Expression<ConstantExpressionPolicy<LhsResultType> >& lhs,
                                 const Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >& rhs, 
                                 Accumulator<ResultType>& accum)
                {
                    typedef typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>::ResultType RhsResultType;
                    
                    rhs.Evaluate(accum);
                    OpType<LhsResultType, RhsResultType>::ApplyLeftEqual(accum, lhs);
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
        class BinaryExpressionEvaluator<ConstantExpressionPolicy<LhsResultType>,
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
            public:
                static void Eval(const Expression<ConstantExpressionPolicy<LhsResultType> >& lhs,
                                 const Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >& rhs, 
                                 Accumulator<ResultType>& accum)
                {
                    typedef typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>::ResultType RhsResultType;
                    
                    rhs.Evaluate(accum);
                    OpType<RhsResultType, LhsResultType>::Apply(accum, lhs);
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
        // Case 4 - Evaluate lhs, evaluate rhs with parent op type.
        // 1--1-- 000010       
        template<typename LhsResultType, typename RhsLhsPolicyType, typename RhsRhsPolicyType,
                 template <typename, typename> class RhsOpType,
                 template <typename, typename> class OpType,
                 typename ResultType>
        class BinaryExpressionEvaluator<ConstantExpressionPolicy<LhsResultType>,
                                        BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>,
                                        ResultType, OpType, BinaryNullOp,
                                        typename boost::enable_if
                                        <
                                            boost::mpl::and_
                                            <
                                                typename AssociativeTraits<ConstantExpressionPolicy<LhsResultType>, OpType, BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >::IsAssociative,
                                                boost::is_same<ResultType, typename ConstantExpressionPolicy<LhsResultType>::ResultType>
                                            >
                                         >::type
                                         >
        {
            public:
                static void Eval(const Expression<ConstantExpressionPolicy<LhsResultType> >& lhs,
                                 const Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >& rhs, 
                                 Accumulator<ResultType>& accum)
                {
                    typedef typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>::ResultType RhsResultType;
                    
                    lhs.Evaluate(accum);
                    rhs.template Evaluate<OpType>(accum);
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
        // 
        // Case 5 - Evaluate new (A+B) then op= C
        // 100011 000001
         template<typename LhsResultType, typename RhsLhsPolicyType, typename RhsRhsPolicyType,
                 template <typename, typename> class RhsOpType,
                 template <typename, typename> class OpType,
                 typename ResultType>
        class BinaryExpressionEvaluator<ConstantExpressionPolicy<LhsResultType>,
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
            public:
                static void Eval(const Expression<ConstantExpressionPolicy<LhsResultType> >& lhs,
                                 const Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >& rhs, 
                                 Accumulator<ResultType>& accum)
                {
                    typedef typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>::ResultType RhsResultType;
                    
                    // TODO.
                }
                
                #ifdef NEKTAR_UNIT_TESTS
                static const unsigned int ClassNum = 5;
                #endif //NEKTAR_UNIT_TESTS
        };


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
        class BinaryExpressionEvaluator<BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType>,
                                        ConstantExpressionPolicy<RhsResultType>,
                                        ResultType, OpType, BinaryNullOp,
                                        typename boost::enable_if
                                        <
                                            boost::is_same<typename BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType>::ResultType, ResultType>
                                        >::type >
        {
            public:
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
        class BinaryExpressionEvaluator<BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType>,
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
            public:
                static void Eval(const Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType> >& lhs,
                                 const Expression<ConstantExpressionPolicy<RhsResultType> >& rhs,
                                 Accumulator<ResultType>& accum)
                {
                    typedef typename BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType>::ResulType LhsResultType;
                    
                    LhsResultType lhsTemp = lhs.Evaluate();
                    OpType<LhsResultType, RhsResultType>::Apply(accum, lhsTemp, rhs);
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
        class BinaryExpressionEvaluator<BinaryExpressionPolicy<LhsLhsPolicyType, LhsOpType, LhsRhsPolicyType>,
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
            public:
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
                      
        ////////////////////////////////////////////////////////////////////////////////////////////////
        /// Expression where one side is constant and the other is binary.
        /// If the left side is binary, then we always evaluate it first.  If the return type of the 
        /// lhs is the same as the expression, then we can just call it directly, otherwise we need to 
        /// create a temporary.

        // R = (A + B) - C


        
        // R = (A + B) - C
        //
        // Case 1 - typeof(A+B) = typeof(R)
        // Simply pass the accumulator to the left and then assume -= C exists.
        
        
        // R = (A + B) - C
        //
        // Case 2 - typeof(A+B) != typeof(R)
        //          typeof(A+B) - C is associative (with or without operator change)
        //          typeof(B-C) == typeof(R)
        //
        // Case 2.1
        //          typeof(A) + typeof(B-C) is commutative
        //          Pass accumulator to (B-C), then assume += A is available.
        //
        // Case 2.2 
        //          typeof(A) + typeof(B-C) is not commutative.
        //          typeof(A) =+ typeof(B-C) is defined.
        //          Pass accumulator to (B-C), the =+ A.
        //
        // Case 3 - Temporary to (A+B), then op (A+B) - C
        
        
        
        
        // R = A + (B - C)
        //
        // Case 1 - typeof(A) + typeof(B-C) is commutative
        // Evaluate (B-C) + A
        //
        // Case 2 - typeof(A) + typeof(B-C) is not comutative
        // Case 2.1 - typeof(A) + typeof(B-C) is associative
        //            Evalaute (A+B) - C
        //
        // Case 3 - typeof(A) + typeof(B-C) is not associative or commutative.
        // Case 3.1 - typeof(B-C) = typeof(R)
        //            typeof(A) =+ typeof(R) is defined.
        //            Evaluate (B-C), then =+ A.
        //
        // Case 4 - Otherwise, create a temproary for (B-C) and call Op directly.
        
        
        ////////////////////////////////////////////////////////////////////////////////////
        /// 2 op expression, no parent op types.
        /// 
        /// Expressions such as R = A + (B * C) and R = (A + B) * C
        /// Note that the letters above represent types, not values, and the operator + and *
        /// are just placeholders.
        ///
        /// In order of preference:
        /// 1. Don't create a temporary.
        /// 2. Use op= before leftOp=
        ///
        /// 5 boolean variables for the case of R = A + (B*C)
        /// v = typeof(B*C) = R
        /// w = A + typeof(B*C) is commutative
        /// x = A + (B*C) is associative
        /// y = typeof(B*C) +left= A is defined
        /// z = R = A is defined.
        ///
        /// 32 different scenarios
        /// ----------------------------------------------------------------------------------
        /// Case #1 - 8 scenarios
        /// Conditions: R = A + (B * C)
        ///             Typeof (B * C) is R
        ///             A + R is commutative
        ///
        /// 1. R = (B*C)
        /// 2. R += A
//         template<typename A, typename RhsLhsPolicyType, typename RhsRhsPolicyType, typename R,
//                  template <typename, typename> class RhsOpType, 
//                  template <typename, typename> class OpType>
//         class BinaryExpressionEvaluator<ConstantExpressionPolicy<A>,
//                                         BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType>,
//                                         R, OpType, BinaryNullOp,
//                                         typename boost::enable_if_c
//                                         <
//                                                 boost::is_same
//                                                 <
//                                                     R, 
//                                                     typename Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >::ResultType
//                                                 >::value &&
//                                                 CommutativeTraits<A, OpType, R>::IsCommutative
//                                             >::type >
//         {
//             public:
//                 static void Eval(const Expression<ConstantExpressionPolicy<A> >& lhs, 
//                                 const Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >& rhs,
//                                 Accumulator<R>& accum)
//                 {
//                     rhs.template Apply<BinaryNullOp>(accum);
//                     OpType<R, A>::ApplyEqual(accum, *lhs);
//                 }
//         };
        
        /// ----------------------------------------------------------------------------------
        /// Case #2
        /// Conditions: R = A + (B * C)
        ///             [typeof(B*C) != R] OR (16 scenarios)
        ///             [typeof(B*C) = R and A + R is not commutative, and A + (B*C) is associative, and R = A does not exist, and R +left= A does not exist] (1 scenario)
        ///             [typeof(B*C) = R and A + R is not commutative, and A + (B*C) is not associative] (4 scenarios)
        ///
        /// In this scenario, since the binary expression requires an accumulator other than an 
        /// R, we have no choice but to create a temporary.
        ///
        /// 1. Create T, typeof(B*C)
        /// 2. T = (B*C)
        /// 3. Result = A + T
//         template<typename A, typename RhsLhsPolicyType, typename RhsRhsPolicyType, typename R,
//                 template <typename, typename> class RhsOpType, 
//                 template <typename, typename> class OpType>
//         class BinaryExpressionEvaluator<ConstantExpressionPolicy<A>,
//                                         BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType>,
//                                         R, OpType, BinaryNullOp,
//                                         typename boost::enable_if_c
//                                         <
//                                                 !boost::is_same
//                                                 <
//                                                     R, 
//                                                     typename Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >::ResultType
//                                                 >::value ||
//                                                 (!CommutativeTraits<A, OpType, R>::IsCommutative &&
//                                                 AssociativeTraits<A, OpType, typename RhsLhsPolicyType::ResultType, RhsOpType, typename RhsRhsPolicyType::ResultType>::IsAssociative &&
//                                                 !AssignableTraits<R, A>::IsAssignable &&
//                                                 !OpType<R, A>::HasOpLeftEqual)
//                                             >::type >
//         {
//             public:
//                 typedef typename Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >::ResultType RhsDataType;
// 
//                 static void Eval(const Expression<ConstantExpressionPolicy<A> >& lhs, 
//                                 const Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >& rhs,
//                                 Accumulator<R>& accum)
//                 {
//                     RhsDataType temp = rhs;
//                     OpType<A, RhsDataType>::Apply(accum, *lhs, rhs);
//                 }
//         };

        /// ----------------------------------------------------------------------------------
        /// Case #3 (1 scenario)
        /// Conditions: R = A + (B*C)
        ///             typeof(B*C) is R
        ///             A + R is not commutative
        ///             A + (B*C) is associative
        ///             R +left= A is not defined.
        ///             R = A is defined.
        ///             
        ///
        /// 1. R = A
        /// 2. R += B*C
//         template<typename A, typename RhsLhsPolicyType, typename RhsRhsPolicyType, typename R,
//                 template <typename, typename> class RhsOpType, 
//                 template <typename, typename> class OpType>
//         class BinaryExpressionEvaluator<ConstantExpressionPolicy<A>,
//                                         BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType>,
//                                         R, OpType, BinaryNullOp,
//                                         typename boost::enable_if_c
//                                         <
//                                                 boost::is_same
//                                                 <
//                                                     R, 
//                                                     typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType>::ResultType
//                                                 >::value &&
//                                                 !CommutativeTraits<A, OpType, R>::IsCommutative &&
//                                                 AssociativeTraits<A, OpType, typename RhsLhsPolicyType::ResultType, RhsOpType, typename RhsRhsPolicyType::ResultType>::IsAssociative &&
//                                                 !OpType<R, A>::HasOpLeftEqual &&
//                                                 AssignableTraits<R, A>::IsAssignable
//                                             >::type >
//         {
//             public:
//                 static void Eval(const Expression<ConstantExpressionPolicy<A> >& lhs, 
//                                 const Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >& rhs,
//                                 Accumulator<R>& accum)
//                 {
//                     *accum = *lhs;
//                     rhs.template Apply<OpType>(accum);
//                 }
//         };


        /// ----------------------------------------------------------------------------------
        /// Case #4 (2 scenarios)
        /// Conditions: R = A + (B*C)
        ///             typeof(B*C) = R
        ///             A + R is not commutative
        ///             R +left= A is defined
        ///
        /// 1. R = B*C
        /// 2. R +left= A
//         template<typename A, typename RhsLhsPolicyType, typename RhsRhsPolicyType, typename R,
//                 template <typename, typename> class RhsOpType, 
//                 template <typename, typename> class OpType>
//         class BinaryExpressionEvaluator<ConstantExpressionPolicy<A>,
//                                         BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType>,
//                                         R, OpType, BinaryNullOp,
//                                         typename boost::enable_if_c
//                                         <
//                                                 boost::is_same
//                                                 <
//                                                     R, 
//                                                     typename Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >::ResultType
//                                                 >::value &&
//                                                 !CommutativeTraits<A, OpType, R>::IsCommutative &&
//                                                 OpType<R, A>::HasOpLeftEqual,
//                                                 void
//                                             >::type >
//         {
//             public:
//                 static void Eval(const Expression<ConstantExpressionPolicy<A> >& lhs, 
//                                 const Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >& rhs,
//                                 Accumulator<R>& accum)
//                 {
//                     rhs.template Apply<BinaryNullOp>(accum);
//                     OpType<A, R>::ApplyLeftEqual(accum, *lhs);
//                 }
//         };


        /// ----------------------------------------------------------------------------------
        /// R = (A + B) * C
        /// ----------------------------------------------------------------------------------
               
        /// LhsBinaryRhsConstant Specialization 1
        /// Conditions: R = (A+B) * C
        ///             typeof(A+B) = R
        ///
        /// 1. R = A+B
        /// 2. R += C;
//         template<typename C, typename LhsLhsPolicyType, typename LhsRhsPolicyType, typename R,
//                  template <typename, typename> class LhsOpType, 
//                  template <typename, typename> class OpType>
//         class BinaryExpressionEvaluator<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType>,
//                                         ConstantExpressionPolicy<C>,
//                                         R, OpType, BinaryNullOp,
//                                         typename boost::enable_if
//                                         <
//                                             boost::is_same
//                                             <
//                                                 R, 
//                                                 typename BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType>::ResultType
//                                             >
//                                         >::type >
//         {
//             public:
//                 static void Eval(const Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >& lhs,
//                                  const Expression<ConstantExpressionPolicy<C> >& rhs, 
//                                  Accumulator<R>& accum)
//                 {
//                     lhs.template Apply<BinaryNullOp>(accum);
//                     OpType<R, C>::ApplyEqual(accum, *rhs);
//                 }
//         };
        
        /// LhsBinaryRhsConstant Specialization 2
        /// Conditions: R = (A+B) * C
        ///             typeof(A+B) != R
        ///             (A+B) * C is associative
        ///             typeof(B*C) = R
        ///             A + R is Commutative
        ///
        /// 1. R = B*C
        /// 2. R += A
        ///
        /// TODO - This doesn't work with changed op types.
//         template<typename C, typename LhsLhsPolicyType, typename LhsRhsPolicyType, typename R,
//                  template <typename, typename> class LhsOpType, 
//                  template <typename, typename> class OpType>
//         class BinaryExpressionEvaluator<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType>,
//                                         ConstantExpressionPolicy<C>,
//                                         R, OpType, BinaryNullOp,
//                                         typename boost::enable_if_c
//                                         <
//                                             !boost::is_same
//                                             <
//                                                 R,
//                                                 typename BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType>::ResultType
//                                             >::value &&
//                                             AssociativeTraits<typename LhsLhsPolicyType::ResultType, LhsOpType, typename LhsRhsPolicyType::ResultType, OpType, C>::IsStrictlyAssociative &&
//                                             boost::is_same
//                                             <
//                                                 R,
//                                                 typename BinaryExpressionPolicy<LhsRhsPolicyType, ConstantExpressionPolicy<C>, OpType>::ResultType
//                                             >::value &&
//                                             CommutativeTraits<typename LhsLhsPolicyType::ResultType, LhsOpType, R>::IsCommutative
//                                         >::type >
//         {
//             public:
//                 static void Eval(const Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >& lhs,
//                                  const Expression<ConstantExpressionPolicy<C> >& rhs, 
//                                  Accumulator<R>& accum)
//                 {
//                 }
//         };

        /// Conditions: R = (A+B) * C
        ///             typeof(A+B) != C
        ///
        /// 1. temp = lhs
        /// 2. R = temp + C
//         template<typename C, typename LhsLhsPolicyType, typename LhsRhsPolicyType, typename R,
//                 template <typename, typename> class LhsOpType, 
//                 template <typename, typename> class OpType>
//         class BinaryExpressionEvaluator<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType>,
//                                     ConstantExpressionPolicy<C>,
//                                     R, OpType, BinaryNullOp,
//                                     typename boost::enable_if
//                                     <
//                                             boost::mpl::not_<boost::is_same
//                                             <
//                                                 R, 
//                                                 typename Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >::ResultType
//                                             > >
//                                         >::type >
//         {
//             public:
//                 typedef typename Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >::ResultType LhsDataType;
// 
//                 static void Eval(const Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >& lhs,
//                                 const Expression<ConstantExpressionPolicy<C> >& rhs, 
//                                 Accumulator<R>& accum)
//                 {
//                     LhsDataType t = lhs;
//                     OpType<LhsDataType,  C>::Apply(accum, t, *rhs);
//                 }
//         };
        
        
        ///////////////////////////////
        // Now for 2 op expressions with an incoming parent type.
        // 
        // For the lhs binary expression, the enablers are the same, just the logic is slightly different.
        /// ----------------------------------------------------------------------------------
        /// Now for R = (A + B) * C
        ///
        /// This is considerably easier than the previous cases. 
        ///
        /// ----------------------------------------------------------------------------------
        /// Case #5
        /// Conditions: R = (A+B) * C
        ///             typeof(A+B) = R
        ///
        /// 1. R = A+B
        /// 2. R += C;
/*        template<typename C, typename LhsLhsPolicyType, typename LhsRhsPolicyType, typename R,
                template <typename, typename> class LhsOpType, 
                template <typename, typename> class OpType,
                template <typename, typename> class ParentOpType>
        class BinaryExpressionEvaluator<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType>,
                                    ConstantExpressionPolicy<C>,
                                    R, OpType, ParentOpType,
                                    typename boost::enable_if
                                    <
                                            boost::is_same
                                            <
                                                R, 
                                                typename Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >::ResultType
                                            >
                                        >::type >
        {
            public:
                static void Eval(const Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >& lhs,
                                const Expression<ConstantExpressionPolicy<C> >& rhs, 
                                Accumulator<R>& accum)
                {
                    lhs.template Apply
                    *accum = lhs;
                    OpType<R, C>::ApplyEqual(accum, *rhs);
                }
        };
        
        /// 
        /// ----------------------------------------------------------------------------------
        /// Case #6
        /// Conditions: R = (A+B) * C
        ///             typeof(A+B) != C
        ///
        /// 1. temp = lhs
        /// 2. R = temp + C
        ////////////////////////////////////////////////////////////////////////////////////
        template<typename C, typename LhsLhsPolicyType, typename LhsRhsPolicyType, typename R,
                template <typename, typename> class LhsOpType, 
                template <typename, typename> class OpType>
        class BinaryExpressionEvaluator<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType>,
                                    ConstantExpressionPolicy<C>,
                                    R, OpType, BinaryNullOp,
                                    typename boost::enable_if
                                    <
                                            boost::mpl::not_<boost::is_same
                                            <
                                                R, 
                                                typename Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >::ResultType
                                            > >
                                        >::type >
        {
            public:
                typedef typename Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >::ResultType LhsDataType;

                static void Eval(const Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >& lhs,
                                const Expression<ConstantExpressionPolicy<C> >& rhs, 
                                Accumulator<R>& accum)
                {
                    LhsDataType t = lhs;
                    OpType<LhsDataType,  C>::Apply(accum, t, *rhs);
                }
        };*/
        
        ///
        /// Conditions: R = A + (B*C)
        ///                 R + (B*C) is associative.
//         template<typename A, typename RhsLhsPolicyType, typename RhsRhsPolicyType, typename R,
//                  template <typename, typename> class RhsOpType, 
//                  template <typename, typename> class OpType,
//                  template <typename, typename> class ParentOpType>
//         class BinaryExpressionEvaluator<ConstantExpressionPolicy<A>,
//                                        BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType>,
//                                        R, OpType, ParentOpType,
//                                        typename boost::enable_if_c
//                                        <
//                                             AssociativeTraits<R, OpType, typename RhsLhsPolicyType::ResultType, RhsOpType, typename RhsRhsPolicyType>::IsAssociative,
//                                             void
//                                         >::type >
//         {
//             public:
//                 typedef BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> RhsPolicy;
//                 typedef AssociativeTraits<R, ParentOpType, A, OpType, typename RhsPolicy::ResultType>::OpChage FirstOp;
//                 typedef AssociativeTraits<R, OpType, typename RhsLhsPolicyType::ResultType, RhsOpType, typename RhsRhsPolicyType>::OpChange SecondOp;
//                 
//                 static void Eval(const Expression<ConstantExpressionPolicy<A> >& lhs, 
//                                  const Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >& rhs,
//                                  Accumulator<R>& accum)
//                 {
//                     FirstOp::ApplyEqual(accum, *lhs);
//                     
//                     //This isn't right - Op type should be changed.
//                     rhs.template Apply<OpType>(accum);
//                 }
//         };


        /////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////
        // Binary Expression with two binary leaves
        //
        // This is the general purpose approach.  It would be possible to 
        // use these same algorithms when one branch is constant, but that
        // would be sub-optimal.  For example, In case #2 below, we need to 
        // make a temporary of the rhs, whereas our optimization above can 
        // avoid that temporary in many cases.
        //
        //////////////////////////////////////////////////////////////////////
        
        
        /// \brief Specialization for two binary expressions as children.
        ///
        /// - Both expressions have the same type as the accumulator.
        /// - The expresison is associative
        /// Incoming ParentOp [ (A lop B) op (C rop D) ]
//         template<typename LhsLhsPolicyType, template<typename, typename> class LhsOpType, typename LhsRhsPolicyType,
//                 typename RhsLhsPolicyType, template<typename, typename> class RhsOpType, typename RhsRhsPolicyType,
//                 template <typename, typename> class OpType,
//                 template <typename, typename> class ParentOpType, 
//                 typename ResultType>
//         class BinaryExpressionEvaluator<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType>,
//                                     BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType>,
//                                     ResultType, OpType, ParentOpType,
//                                     typename boost::enable_if_c
//                                     <
//                                             boost::is_same
//                                             <
//                                                 ResultType, 
//                                                 typename Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >::ResultType
//                                             >::value &&
//                                             boost::is_same
//                                             <
//                                                 ResultType, 
//                                                 typename Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >::ResultType
//                                             >::value &&
// 
//                                             AssociativeTraits<typename Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >::ResultType,
//                                                             OpType,
//                                                             typename RhsLhsPolicyType::ResultType,
//                                                             RhsOpType,
//                                                             typename RhsRhsPolicyType::ResultType>::IsAssociative
//                                         >::type >
//         {
//             public:
//                 static void Eval(const Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >& lhs, 
//                                 const Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >& rhs,
//                                 Accumulator<ResultType>& accum)
//                 {
//                     lhs.template ApplyEqual<ParentOpType>(accum);
//                     rhs.template ApplyEqual<OpType>(accum);
//                 }
//         };
//         
//         template<typename LhsLhsPolicyType, template<typename, typename> class LhsOpType, typename LhsRhsPolicyType,
//                 typename RhsLhsPolicyType, template<typename, typename> class RhsOpType, typename RhsRhsPolicyType,
//                 template <typename, typename> class OpType,
//                 template <typename, typename> class ParentOpType, 
//                 typename ResultType>
//         class BinaryExpressionEvaluator<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType>,
//                                     BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType>,
//                                     ResultType, OpType, ParentOpType,
//                                     typename boost::enable_if_c
//                                     <
//                                             boost::is_same
//                                             <
//                                                 ResultType, 
//                                                 typename Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >::ResultType
//                                             >::value &&
//                                             boost::is_same
//                                             <
//                                                 ResultType, 
//                                                 typename Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >::ResultType
//                                             >::value &&
// 
//                                             !AssociativeTraits<typename Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >::ResultType,
//                                                             OpType,
//                                                             typename RhsLhsPolicyType::ResultType,
//                                                             RhsOpType,
//                                                             typename RhsRhsPolicyType::ResultType>::IsAssociative
//                                         >::type >
//         {
//             public:
//                 typedef typename Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >::ResultType RhsType;
//                 static void Eval(const Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >& lhs, 
//                                 const Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >& rhs,
//                                 Accumulator<ResultType>& accum)
//                 {
//                     lhs.template Apply<ParentOpType>(accum);
//                     RhsType temp = rhs;
//                     OpType<ResultType, RhsType>::ApplyEqual(accum, temp);
//                 }
//         };
//         
//         template<typename LhsLhsPolicyType, template<typename, typename> class LhsOpType, typename LhsRhsPolicyType,
//                 typename RhsLhsPolicyType, template<typename, typename> class RhsOpType, typename RhsRhsPolicyType,
//                 template <typename, typename> class OpType,
//                 template <typename, typename> class ParentOpType, 
//                 typename ResultType>
//         class BinaryExpressionEvaluator<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType>,
//                                     BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType>,
//                                     ResultType, OpType, ParentOpType,
//                                     typename boost::enable_if_c
//                                     <
//                                             boost::is_same
//                                             <
//                                                 ResultType, 
//                                                 typename Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >::ResultType
//                                             >::value &&
//                                             !boost::is_same
//                                             <
//                                                 ResultType, 
//                                                 typename Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >::ResultType
//                                             >::value &&
// 
//                                             CommutativeTraits<typename Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >::ResultType,
//                                                             OpType, 
//                                                             typename Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >::ResultType>::IsCommutative
//                                         >::type >
//         {
//             public:
//                 typedef typename Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >::ResultType LhsType;
//                 
//                 static void Eval(const Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >& lhs, 
//                                 const Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >& rhs,
//                                 Accumulator<ResultType>& accum)
//                 {
//                     rhs.template ApplyEqual<ParentOpType>(accum);
//                     LhsType temp = lhs;
//                     OpType<ResultType, LhsType>::ApplyEqual(accum, temp);
//                 }
//         };
//         
//         
//         
//         template<typename LhsLhsPolicyType, template<typename, typename> class LhsOpType, typename LhsRhsPolicyType,
//                 typename RhsLhsPolicyType, template<typename, typename> class RhsOpType, typename RhsRhsPolicyType,
//                 template <typename, typename> class OpType,
//                 template <typename, typename> class ParentOpType, 
//                 typename ResultType>
//         class BinaryExpressionEvaluator<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType>,
//                                     BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType>,
//                                     ResultType, OpType, ParentOpType,
//                                     typename boost::enable_if_c
//                                     <
//                                             (!boost::is_same
//                                             <
//                                                 ResultType, 
//                                                 typename Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >::ResultType
//                                             >::value &&
//                                             !boost::is_same
//                                             <
//                                                 ResultType, 
//                                                 typename Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >::ResultType
//                                             >::value) ||
//                                             
//                                             (boost::is_same
//                                             <
//                                                 ResultType, 
//                                                 typename Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >::ResultType
//                                             >::value &&
//                                             !boost::is_same
//                                             <
//                                                 ResultType, 
//                                                 typename Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >::ResultType
//                                             >::value &&
//                                             !CommutativeTraits<typename Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >::ResultType,
//                                                             OpType, 
//                                                             typename Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >::ResultType>::IsCommutative)
//                                         >::type >
//         {
//             public:
//                 typedef typename Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >::ResultType LhsType;
//                 typedef typename Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >::ResultType RhsType;
//                 
//                 static void Eval(const Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >& lhs, 
//                                 const Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >& rhs,
//                                 Accumulator<ResultType>& accum)
//                 {
//                     LhsType lhs_temp = lhs;
//                     RhsType rhs_temp = rhs;
//                     OpType<LhsType, RhsType>::Apply(accum, lhs_temp, rhs_temp);
//                 }
//         };
}

#endif //NEKTAR_USE_EXPRESSION_TEMPLATES
#endif //NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_BINARY_EXPRESSION_EVALUATOR_HPP

