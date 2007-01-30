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


#include <LibUtilities/ExpressionTemplates/BinaryExpressionEvaluatorFwd.hpp>
#include <LibUtilities/ExpressionTemplates/ConstantExpression.hpp>
#include <LibUtilities/ExpressionTemplates/BinaryExpressionPolicyFwd.hpp>
#include <LibUtilities/ExpressionTemplates/AssociativeTraits.hpp>
#include <LibUtilities/ExpressionTemplates/CommutativeTraits.hpp>
#include <LibUtilities/ExpressionTemplates/AssignableTraits.hpp>
#include <LibUtilities/ExpressionTemplates/BinaryExpressionOperators.hpp>

#include <boost/utility/enable_if.hpp>

namespace Nektar
{
    namespace expt
    {
        // Visual Studio 2005 has a hard time with a construct such as
        // OpType<LhsType, RhsType>::Traits::HasOpEqual.  
        // This class is simply a level of indirection to get this to work.
        // This has the unfortunate side effect of forcing the user to define 
        // HasOpEqual and HasOpLeftEqual in each traits class, where in the past
        // they could just ignore it if they didn't have it.
        template<typename TraitsType>
        class GetTraitsValues
        {
            public:
                static const bool HasOpEqual = TraitsType::HasOpEqual;
                static const bool HasOpLeftEqual = TraitsType::HasOpLeftEqual;
        };

        ////////////////////////////////////////////////
        // 1 op Specializations
        ////////////////////////////////////////////////
        /// Whenever possible, we want to follow the following:
        /// A ParentOpEqual Lhs
        /// A OpEqual Rhs
        /// Where OpEqual is chosen appropriately based on the parent op.
        ///
        /// The idea here is that if we can do this, then the user can do very 
        /// little work as we will use the +=, /=, *=, -= expressions they 
        /// have probably already written.
        /// I don't see any way to avoid forcing them to write the appropriate traits
        /// objects and typedefs to inform me of the existance of these operators, but 
        /// in the ideal case this will be enough.
        ///
        /// This approach can be taken in the following circumstances:
        /// Traits<Accum, Lhs>::op_equal exists 
        /// Traits<Accum, Rhs>::op_eqal exists.
         
        /// This specialization is for expressions in which both sides are constant, of arbitrary
        /// type, and I can opequal both of them in order.  This is the easiest one to program for the 
        /// end user.
        
        /// Case #1 - No incoming parent type.
        /// In this case, by definition, the user must have created a function to perform this operation,
        /// so we can simply call the Apply method on the operator regardless of result, lhs, or rhs data 
        /// type.
        template<typename LhsDataType, typename RhsDataType, typename ResultType, 
                 template <typename, typename> class OpType>
        class EvaluateBinaryExpression<ConstantExpressionPolicy<LhsDataType>,
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
        
        /// Case #2 - Incoming Parent Type, Result op= lhs, Result op= rhs is defined.
        /// In this case, I simply apply the op= operations.
        template<typename LhsDataType, typename RhsDataType, typename ResultType, 
                 template <typename, typename> class OpType,
                 template <typename, typename> class ParentOpType>
        class EvaluateBinaryExpression<ConstantExpressionPolicy<LhsDataType>,
                                       ConstantExpressionPolicy<RhsDataType>,
                                       ResultType, OpType, ParentOpType,
                                       typename boost::enable_if_c
                                       <
                                            GetTraitsValues<typename ParentOpType<ResultType, LhsDataType>::TraitsType>::HasOpEqual &&
                                            GetTraitsValues<typename OpType<ResultType, RhsDataType>::TraitsType>::HasOpEqual &&                                
                                            !boost::is_same<BinaryNullOp<void, void>, typename ParentOpType<LhsDataType, RhsDataType>::template Rebind<void, void>::type >::value,
                                            void
                                       >::type >
        {
            public:
                
                static void Eval(const Expression<ConstantExpressionPolicy<LhsDataType> >& lhs, 
                                 const Expression<ConstantExpressionPolicy<RhsDataType> >& rhs,
                                 Accumulator<ResultType>& accum)
                {
                    ParentOpType<ResultType, LhsDataType>::ApplyEqual(accum, *lhs);
                    AssociativeTraits<ResultType, ParentOpType, LhsDataType, OpType, RhsDataType>::OpChange::ApplyEqual(accum, *rhs);
                }
        };
        
        // The only other case is where the OpEquals are not defined, but this is such a bizarre case that
        // I imagine if they aren't defined they are more likely user errors than a true part of the math.
        // So I will leave it undefined to help catch compiler errors.  
        // TODO - It may be worth defining it just so I can put a boost static assert to aid debugging.



    
        ////////////////////////////////////////////////////////////////////////////////////////////////
        /// Expression where one side is constant and the other is binary.
        /// If the left side is binary, then we always evaluate it first.  If the return type of the 
        /// lhs is the same as the expression, then we can just call it directly, otherwise we need to 
        /// create a temporary.
        
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
        template<typename A, typename RhsLhsPolicyType, typename RhsRhsPolicyType, typename R,
                 template <typename, typename> class RhsOpType, 
                 template <typename, typename> class OpType>
        class EvaluateBinaryExpression<ConstantExpressionPolicy<A>,
                                       BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType>,
                                       R, OpType, BinaryNullOp,
                                       typename boost::enable_if_c
                                       <
                                            boost::is_same
                                            <
                                                R, 
                                                typename Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >::ResultType
                                            >::value &&
                                            CommutativeTraits<A, OpType, R>::IsCommutative
                                        >::type >
        {
            public:
                static void Eval(const Expression<ConstantExpressionPolicy<A> >& lhs, 
                                 const Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >& rhs,
                                 Accumulator<R>& accum)
                {
                    rhs.template Apply<BinaryNullOp>(accum);
                    OpType<R, A>::ApplyEqual(accum, *lhs);
                }
        };
        
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
        template<typename A, typename RhsLhsPolicyType, typename RhsRhsPolicyType, typename R,
                 template <typename, typename> class RhsOpType, 
                 template <typename, typename> class OpType>
        class EvaluateBinaryExpression<ConstantExpressionPolicy<A>,
                                       BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType>,
                                       R, OpType, BinaryNullOp,
                                       typename boost::enable_if_c
                                       <
                                            !boost::is_same
                                            <
                                                R, 
                                                typename Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >::ResultType
                                            >::value ||
                                            (!CommutativeTraits<A, OpType, R>::IsCommutative &&
                                            AssociativeTraits<A, OpType, typename RhsLhsPolicyType::ResultType, RhsOpType, typename RhsRhsPolicyType::ResultType>::IsAssociative &&
                                            !AssignableTraits<R, A>::IsAssignable &&
                                            !GetTraitsValues<typename OpType<R, A>::TraitsType>::HasOpLeftEqual)
                                        >::type >
        {
            public:
                typedef typename Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >::ResultType RhsDataType;

                static void Eval(const Expression<ConstantExpressionPolicy<A> >& lhs, 
                                 const Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >& rhs,
                                 Accumulator<R>& accum)
                {
                    RhsDataType temp = rhs;
                    OpType<A, RhsDataType>::Apply(accum, *lhs, rhs);
                }
        };

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
        template<typename A, typename RhsLhsPolicyType, typename RhsRhsPolicyType, typename R,
                 template <typename, typename> class RhsOpType, 
                 template <typename, typename> class OpType>
        class EvaluateBinaryExpression<ConstantExpressionPolicy<A>,
                                       BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType>,
                                       R, OpType, BinaryNullOp,
                                       typename boost::enable_if_c
                                       <
                                            boost::is_same
                                            <
                                                R, 
                                                typename Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >::ResultType
                                            >::value &&
                                            !CommutativeTraits<A, OpType, R>::IsCommutative &&
                                            AssociativeTraits<A, OpType, typename RhsLhsPolicyType::ResultType, RhsOpType, typename RhsRhsPolicyType::ResultType>::IsAssociative &&
                                            !GetTraitsValues<typename OpType<R, A>::TraitsType>::HasOpLeftEqual &&
                                            AssignableTraits<R, A>::IsAssignable
                                        >::type >
        {
            public:
                static void Eval(const Expression<ConstantExpressionPolicy<A> >& lhs, 
                                 const Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >& rhs,
                                 Accumulator<R>& accum)
                {
                    *accum = *lhs;
                    rhs.template Apply<OpType>(accum);
                }
        };


        /// ----------------------------------------------------------------------------------
        /// Case #4 (2 scenarios)
        /// Conditions: R = A + (B*C)
        ///             typeof(B*C) = R
        ///             A + R is not commutative
        ///             R +left= A is defined
        ///
        /// 1. R = B*C
        /// 2. R +left= A
        template<typename A, typename RhsLhsPolicyType, typename RhsRhsPolicyType, typename R,
                 template <typename, typename> class RhsOpType, 
                 template <typename, typename> class OpType>
        class EvaluateBinaryExpression<ConstantExpressionPolicy<A>,
                                       BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType>,
                                       R, OpType, BinaryNullOp,
                                       typename boost::enable_if_c
                                       <
                                            boost::is_same
                                            <
                                                R, 
                                                typename Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >::ResultType
                                            >::value &&
                                            !CommutativeTraits<A, OpType, R>::IsCommutative &&
                                            GetTraitsValues<typename OpType<R, A>::TraitsType>::HasOpLeftEqual,
                                            void
                                        >::type >
        {
            public:
                static void Eval(const Expression<ConstantExpressionPolicy<A> >& lhs, 
                                 const Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >& rhs,
                                 Accumulator<R>& accum)
                {
                    rhs.template Apply<BinaryNullOp>(accum);
                    OpType<A, R>::ApplyLeftEqual(accum, *lhs);
                }
        };


        ///             
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
        template<typename C, typename LhsLhsPolicyType, typename LhsRhsPolicyType, typename R,
                 template <typename, typename> class LhsOpType, 
                 template <typename, typename> class OpType>
        class EvaluateBinaryExpression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType>,
                                       ConstantExpressionPolicy<C>,
                                       R, OpType, BinaryNullOp,
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
        class EvaluateBinaryExpression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType>,
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
        };
        
        
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
        class EvaluateBinaryExpression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType>,
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
        class EvaluateBinaryExpression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType>,
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
//         class EvaluateBinaryExpression<ConstantExpressionPolicy<A>,
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
        template<typename LhsLhsPolicyType, template<typename, typename> class LhsOpType, typename LhsRhsPolicyType,
                 typename RhsLhsPolicyType, template<typename, typename> class RhsOpType, typename RhsRhsPolicyType,
                 template <typename, typename> class OpType,
                 template <typename, typename> class ParentOpType, 
                 typename ResultType>
        class EvaluateBinaryExpression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType>,
                                       BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType>,
                                       ResultType, OpType, ParentOpType,
                                       typename boost::enable_if_c
                                       <
                                            boost::is_same
                                            <
                                                ResultType, 
                                                typename Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >::ResultType
                                            >::value &&
                                            boost::is_same
                                            <
                                                ResultType, 
                                                typename Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >::ResultType
                                            >::value &&

                                            AssociativeTraits<typename Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >::ResultType,
                                                              OpType,
                                                              typename RhsLhsPolicyType::ResultType,
                                                              RhsOpType,
                                                              typename RhsRhsPolicyType::ResultType>::IsAssociative
                                        >::type >
        {
            public:
                static void Eval(const Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >& lhs, 
                                 const Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >& rhs,
                                 Accumulator<ResultType>& accum)
                {
                    lhs.template ApplyEqual<ParentOpType>(accum);
                    rhs.template ApplyEqual<OpType>(accum);
                }
        };
        
        template<typename LhsLhsPolicyType, template<typename, typename> class LhsOpType, typename LhsRhsPolicyType,
                 typename RhsLhsPolicyType, template<typename, typename> class RhsOpType, typename RhsRhsPolicyType,
                 template <typename, typename> class OpType,
                 template <typename, typename> class ParentOpType, 
                 typename ResultType>
        class EvaluateBinaryExpression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType>,
                                       BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType>,
                                       ResultType, OpType, ParentOpType,
                                       typename boost::enable_if_c
                                       <
                                            boost::is_same
                                            <
                                                ResultType, 
                                                typename Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >::ResultType
                                            >::value &&
                                            boost::is_same
                                            <
                                                ResultType, 
                                                typename Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >::ResultType
                                            >::value &&

                                            !AssociativeTraits<typename Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >::ResultType,
                                                              OpType,
                                                              typename RhsLhsPolicyType::ResultType,
                                                              RhsOpType,
                                                              typename RhsRhsPolicyType::ResultType>::IsAssociative
                                        >::type >
        {
            public:
                typedef typename Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >::ResultType RhsType;
                static void Eval(const Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >& lhs, 
                                 const Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >& rhs,
                                 Accumulator<ResultType>& accum)
                {
                    lhs.template Apply<ParentOpType>(accum);
                    RhsType temp = rhs;
                    OpType<ResultType, RhsType>::ApplyEqual(accum, temp);
                }
        };
        
        template<typename LhsLhsPolicyType, template<typename, typename> class LhsOpType, typename LhsRhsPolicyType,
                 typename RhsLhsPolicyType, template<typename, typename> class RhsOpType, typename RhsRhsPolicyType,
                 template <typename, typename> class OpType,
                 template <typename, typename> class ParentOpType, 
                 typename ResultType>
        class EvaluateBinaryExpression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType>,
                                       BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType>,
                                       ResultType, OpType, ParentOpType,
                                       typename boost::enable_if_c
                                       <
                                            boost::is_same
                                            <
                                                ResultType, 
                                                typename Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >::ResultType
                                            >::value &&
                                            !boost::is_same
                                            <
                                                ResultType, 
                                                typename Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >::ResultType
                                            >::value &&

                                            CommutativeTraits<typename Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >::ResultType,
                                                              OpType, 
                                                              typename Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >::ResultType>::IsCommutative
                                        >::type >
        {
            public:
                typedef typename Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >::ResultType LhsType;
                
                static void Eval(const Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >& lhs, 
                                 const Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >& rhs,
                                 Accumulator<ResultType>& accum)
                {
                    rhs.template ApplyEqual<ParentOpType>(accum);
                    LhsType temp = lhs;
                    OpType<ResultType, LhsType>::ApplyEqual(accum, temp);
                }
        };
        
        
        
        template<typename LhsLhsPolicyType, template<typename, typename> class LhsOpType, typename LhsRhsPolicyType,
                 typename RhsLhsPolicyType, template<typename, typename> class RhsOpType, typename RhsRhsPolicyType,
                 template <typename, typename> class OpType,
                 template <typename, typename> class ParentOpType, 
                 typename ResultType>
        class EvaluateBinaryExpression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType>,
                                       BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType>,
                                       ResultType, OpType, ParentOpType,
                                       typename boost::enable_if_c
                                       <
                                            (!boost::is_same
                                            <
                                                ResultType, 
                                                typename Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >::ResultType
                                            >::value &&
                                            !boost::is_same
                                            <
                                                ResultType, 
                                                typename Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >::ResultType
                                            >::value) ||
                                             
                                            (boost::is_same
                                            <
                                                ResultType, 
                                                typename Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >::ResultType
                                            >::value &&
                                            !boost::is_same
                                            <
                                                ResultType, 
                                                typename Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >::ResultType
                                            >::value &&
                                            !CommutativeTraits<typename Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >::ResultType,
                                                              OpType, 
                                                              typename Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >::ResultType>::IsCommutative)
                                        >::type >
        {
            public:
                typedef typename Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >::ResultType LhsType;
                typedef typename Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >::ResultType RhsType;
                
                static void Eval(const Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >& lhs, 
                                 const Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >& rhs,
                                 Accumulator<ResultType>& accum)
                {
                    LhsType lhs_temp = lhs;
                    RhsType rhs_temp = rhs;
                    OpType<LhsType, RhsType>::Apply(accum, lhs_temp, rhs_temp);
                }
        };
    }
}

#endif //NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_BINARY_EXPRESSION_EVALUATOR_HPP

