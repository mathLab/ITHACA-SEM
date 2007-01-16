////////////////////////////////////////////////////////////////////////////////
//
// BinaryExpression.hpp
// Blake Nelson
//
// Binary expression templates.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_BINARY_EXPRESSION_HPP
#define NEKTAR_LIB_UTILITIES_BINARY_EXPRESSION_HPP

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/logical.hpp>

#include <LibUtilities/ExpressionTemplates/BinaryExpressionOperators.hpp>
#include <LibUtilities/ExpressionTemplates/NullOp.hpp>
#include <LibUtilities/ExpressionTemplates/AssociativeTraits.hpp>
#include <LibUtilities/ExpressionTemplates/CommutativeTraits.hpp>
#include <LibUtilities/ExpressionTemplates/ConstantExpression.hpp>
#include <LibUtilities/ExpressionTemplates/UnaryExpression.hpp>
#include <LibUtilities/ExpressionTemplates/BinaryExpressionPolicy.hpp>
#include <LibUtilities/ExpressionTemplates/BinaryExpressionEvaluator.hpp>
#include <LibUtilities/ExpressionTemplates/AssociativeExpressionTraits.hpp>

#include <iostream>


namespace Nektar
{
    namespace expt
    {
        

        // To aid evaluation, adjust the expression trees as they are being created.
        // Generic constant + binary expression
        
        
        // ReturnTypeGenerator
        // Responsibility is to generate the return type of the given binary operation.
        
        // Responsibilities of the Generator
        // 1. Generate the correct return type for a given expression.
        // 2. Provide an Apply method which creates the appropriate return type.
        // 3. Some future optimization may include re-arranging the tree to optimize evaluation.
        template<typename LhsExpressionPolicyType, 
                 template <typename, typename> class OpType,
                 typename RhsExpressionPolicyType,
                 typename enabled = void>
        class BinaryExpressionGenerator;
        
        template<typename LhsExpressionPolicy, template <typename, typename> class OpType, typename RhsDataType>
        class BinaryExpressionGenerator<LhsExpressionPolicy, OpType, ConstantExpressionPolicy<RhsDataType> >
        {
            public:
                typedef BinaryExpressionPolicy<LhsExpressionPolicy, ConstantExpressionPolicy<RhsDataType>, OpType> ResultPolicy;
                typedef Expression<ResultPolicy> ResultExpression;
                
                static ResultExpression Apply(const Expression<LhsExpressionPolicy>& lhs, typename boost::call_traits<RhsDataType>::const_reference rhs)
                {
                    typename ResultPolicy::DataType d(lhs, Expression<ConstantExpressionPolicy<RhsDataType> >(rhs));
                    return ResultExpression(d);
                }
        };
        
//         // If the expression is associative, then we can re-order 
//         template<typename LhsDataType, template <typename, typename> class OpType, typename RhsLhsDataType,
//                  typename RhsRhsDataType, template <typename, typename> class RhsOpType>
//         class BinaryExpressionGenerator<ConstantExpressionPolicy<LhsDataType>, OpType, 
//                                         BinaryExpressionPolicy<RhsLhsDataType, RhsRhsDataType, RhsOpType>,
//                                         typename boost::enable_if_c
//                                         <
//                                             AssociativeTraits<LhsDataType, OpType, RhsLhsDataType, RhsOpType, RhsRhsDataType>::IsAssociative,
//                                             void
//                                         >::type >
//         {
//             public:
//                 typedef BinaryExpressionPolicy<LhsExpressionPolicy, ConstantExpressionPolicy<RhsDataType>, OpType> ResultPolicy;
//                 typedef Expression<ResultPolicy> ResultExpression;
//                 
//                 static ResultExpression Apply(const Expression<LhsExpressionPolicy>& lhs, typename boost::call_traits<RhsDataType>::const_reference rhs)
//                 {
//                     typename ResultPolicy::DataType d(lhs, Expression<ConstantExpressionPolicy<RhsDataType> >(rhs));
//                     return ResultExpression(d);
//                 }
//         };
        
        
        // Generic binary + constant expression
        
        //////////////////////////////////////////////////////
        /// Addition
        //////////////////////////////////////////////////////

        // Next up - There is a lot of duplication here.  Anything to be done about it?
        // Follow up on the AssociativeExpressionTraits - But for non associative cases and do the 
        // same thing.
        
        template<typename LhsPolicyType, typename RhsPolicyType>
        typename BinaryExpressionTraits<LhsPolicyType, AddOp, RhsPolicyType>::ResultType
        operator+(const Expression<LhsPolicyType>& lhs, const Expression<RhsPolicyType>& rhs)
        {
            typedef BinaryExpressionTraits<LhsPolicyType, AddOp, RhsPolicyType> Traits;
            return Traits::Apply(lhs, rhs);
        }
        
        template<typename LhsDataType, typename RhsPolicyType>
        typename BinaryExpressionTraits<ConstantExpressionPolicy<LhsDataType>, AddOp, RhsPolicyType>::ResultType
        operator+(const LhsDataType& lhs, const Expression<RhsPolicyType>& rhs)
        {
            return Expression<ConstantExpressionPolicy<LhsDataType> >(lhs) + rhs;
        }
        
        template<typename LhsPolicyType, typename RhsDataType>
        typename BinaryExpressionTraits<LhsPolicyType, AddOp, ConstantExpressionPolicy<RhsDataType> >::ResultType
        operator+(const Expression<LhsPolicyType>& lhs, const RhsDataType& rhs)
        {
            return lhs + Expression<ConstantExpressionPolicy<RhsDataType> >(rhs);
        }
        // C + C
/*        template<typename LhsDataType, typename RhsDataType>
        Expression<BinaryExpressionPolicy<ConstantExpresionPolicy<LhsDataType>, ConstantExpressionPolicy<RhsDataType>, AddOp> >
        operator+(const Expression<ConstantExpressionPolicy<LhsDataType> >& lhs, 
                  const Expression<ConstantExpressionPolicy<RhsLhsDataType> >& rhs)
        {
            typedef ConstantExpressionPolicy<LhsDataType> LhsPolicy;
            typedef ConstantExpressionPolicy<RhsDataType> RhsPolicy;
            typedef BinaryExpressionPolicy<LhsPolicy, RhsPolicy> ResultPolicy;
            typedef Expression<ResultPolicy> ResultType;
            
            typename ResultPolicy::DataType d(lhs, rhs);
            return ResultType(d);
        }
        
        // C + U
        template<typename LhsDataType, typename RhsDataType, template <typename> class RhsOpType>
        Expression<BinaryExpressionPolicy<ConstantExpresionPolicy<LhsDataType>, UnaryExpressionPolicy<RhsDataType, RhsOpType>, AddOp> >
        operator+(const Expression<ConstantExpressionPolicy<LhsDataType> >& lhs, 
                  const Expression<UnaryExpressionPolicy<RhsLhsDataType, OpType> >& rhs)
        {
            typedef ConstantExpressionPolicy<LhsDataType> LhsPolicy;
            typedef UnaryExpressionPolicy<RhsDataType, RhsOpType> RhsPolicy;
            typedef BinaryExpressionPolicy<LhsPolicy, RhsPolicy> ResultPolicy;
            typedef Expression<ResultPolicy> ResultType;
            
            typename ResultPolicy::DataType d(lhs, rhs);
            return ResultType(d);
        }
        
        // C + B
        template<typename LhsDataType, typename RhsLhsDataType,
                 template<typename, typename> class RhsOpType,
                 typename RhsRhsDataType>
        typename AssociativeExpressionTraits
                 <
                    ConstantExpressionPolicy<LhsDataType>,
                    AddOp,
                    BinaryExpressionPolicy<RhsLhsDataType, RhsRhsDataType, RhsOpType>
                 >::ResultType
        operator+(const LhsDataType& lhs, const Expression<BinaryExpressionPolicy<RhsLhsDataType, RhsRhsDataType, RhsOpType> >& rhs)
        {
            typedef ConstantExpressionPolicy<LhsDataType> LhsPolicy;
            typedef BinaryExpressionPolicy<RhsLhsDataType, RhsRhsDataType, RhsOpType> RhsPolicy;
            
            return AssociativeExpressionTraits<LhsPolicy, AddOp, RhsPolicy>::Apply(lhs, rhs);
        }
        
        
        // U + C
        template<typename LhsDataType, template<typename> class LhsOpType, typename RhsDataType>
        Expression<BinaryExpressionPolicy<UnaryExpresionPolicy<LhsDataType, LhsOpType>, ConstantExpressionPolicy<RhsDataType>, AddOp> >
        operator+(const Expression<UnaryExpressionPolicy<LhsDataType, LhsOpType> >& lhs, 
                  const Expression<ConstantExpressionPolicy<RhsLhsDataType> >& rhs)
        {
            typedef UnaryExpressionPolicy<LhsDataType, LhsOpType> LhsPolicy;
            typedef ConstantExpressionPolicy<RhsDataType> RhsPolicy;
            typedef BinaryExpressionPolicy<LhsPolicy, RhsPolicy> ResultPolicy;
            typedef Expression<ResultPolicy> ResultType;
            
            typename ResultPolicy::DataType d(lhs, rhs);
            return ResultType(d);
        }
        
        // C + U
        template<typename LhsDataType, typename RhsDataType, template <typename> class OpType>
        Expression<BinaryExpressionPolicy<ConstantExpresionPolicy<LhsDataType>, UnaryExpressionPolicy<RhsDataType, OpType>, AddOp> >
        operator+(const Expression<ConstantExpressionPolicy<LhsDataType> >& lhs, 
                  const Expression<UnaryExpressionPolicy<RhsLhsDataType, OpType> >& rhs)
        {
            typedef ConstantExpressionPolicy<LhsDataType> LhsPolicy;
            typedef UnaryExpressionPolicy<RhsDataType, OpType> RhsPolicy;
            typedef BinaryExpressionPolicy<LhsPolicy, RhsPolicy> ResultPolicy;
            typedef Expression<ResultPolicy> ResultType;
            
            typename ResultPolicy::DataType d(lhs, rhs);
            return ResultType(d);
        }
        
        // C + B
        template<typename LhsDataType, typename RhsLhsDataType,
                 template<typename, typename> class RhsOpType,
                 typename RhsRhsDataType>
        typename AssociativeExpressionTraits
                 <
                    ConstantExpressionPolicy<LhsDataType>,
                    AddOp,
                    BinaryExpressionPolicy<RhsLhsDataType, RhsRhsDataType, RhsOpType>
                 >::ResultType
        operator+(const LhsDataType& lhs, const Expression<BinaryExpressionPolicy<RhsLhsDataType, RhsRhsDataType, RhsOpType> >& rhs)
        {
            typedef ConstantExpressionPolicy<LhsDataType> LhsPolicy;
            typedef BinaryExpressionPolicy<RhsLhsDataType, RhsRhsDataType, RhsOpType> RhsPolicy;
            
            return AssociativeExpressionTraits<LhsPolicy, AddOp, RhsPolicy>::Apply(lhs, rhs);
        }*/
        
/*        
        // Original Addition operators
        template<typename LhsExpressionPolicyType, typename RhsExpressionPolicyType>
        Expression<BinaryExpressionPolicy< LhsExpressionPolicyType, RhsExpressionPolicyType, AddOp > > operator+(
            const Expression<LhsExpressionPolicyType>& lhs, const Expression<RhsExpressionPolicyType>& rhs)
        {
            typedef BinaryExpressionPolicy< LhsExpressionPolicyType, RhsExpressionPolicyType, AddOp> ResultPolicyType;
            typename ResultPolicyType::DataType d(lhs, rhs);
            return Expression<ResultPolicyType>(d);
        }

        template<typename LhsDataType, typename RhsExpressionPolicy>
        Expression<BinaryExpressionPolicy<
            ConstantExpressionPolicy<LhsDataType>,
            RhsExpressionPolicy,
            AddOp > > operator+(
            const LhsDataType& lhs,
            const Expression<RhsExpressionPolicy>& rhs)
        {
            typedef ConstantExpressionPolicy<LhsDataType> LhsExpressionPolicy;
            typedef Expression<LhsExpressionPolicy> LhsExpressionType;
            typedef Expression<RhsExpressionPolicy> RhsExpressionType;

            typedef BinaryExpressionPolicy<LhsExpressionPolicy, RhsExpressionPolicy, AddOp> ResultPolicy;
            typename ResultPolicy::DataType d(Expression<LhsExpressionPolicy>(lhs), rhs);
            return Expression<ResultPolicy>(d);
        }

        template<typename LhsPolicy, typename RhsDataType>
        Expression<BinaryExpressionPolicy<
            LhsPolicy,
            ConstantExpressionPolicy<RhsDataType >,
            AddOp > > operator+(
            const Expression<LhsPolicy>& lhs,
            const RhsDataType& rhs)
        {
            // Even newer
            //return AssociativeExpressionTraits<LhsPolicy, ConstantExpressionPolicy<RhsDataType>, AddOp>::Apply(lhs, rhs);
            // New
            return BinaryExpressionGenerator<LhsPolicy, AddOp, ConstantExpressionPolicy<RhsDataType> >::Apply(lhs, rhs);
        
            // Original
//             typedef Expression<LhsExpressionPolicy> LhsExpressionType;
//             typedef ConstantExpressionPolicy<RhsDataType> RhsExpressionPolicy;
//             typedef Expression<RhsExpressionPolicy> RhsExpressionType;
// 
//             typedef BinaryExpressionPolicy<LhsExpressionPolicy, RhsExpressionPolicy, AddOp> ResultPolicyType;
//             typename ResultPolicyType::DataType d(lhs, Expression<RhsExpressionPolicy>(rhs));
//             return Expression<ResultPolicyType>(d);
        }*/
        
        ////////////////////////////////
        // Multiplication
        ///////////////////////////////
        template<typename LhsPolicyType, typename RhsPolicyType>
        typename BinaryExpressionTraits<LhsPolicyType, MultiplyOp, RhsPolicyType>::ResultType
        operator*(const Expression<LhsPolicyType>& lhs, const Expression<RhsPolicyType>& rhs)
        {
            typedef BinaryExpressionTraits<LhsPolicyType, MultiplyOp, RhsPolicyType> Traits;
            return Traits::Apply(lhs, rhs);
        }
        
        template<typename LhsDataType, typename RhsPolicyType>
        typename BinaryExpressionTraits<ConstantExpressionPolicy<LhsDataType>, MultiplyOp, RhsPolicyType>::ResultType
        operator*(const LhsDataType& lhs, const Expression<RhsPolicyType>& rhs)
        {
            return Expression<ConstantExpressionPolicy<LhsDataType> >(lhs) * rhs;
        }
        
        template<typename LhsPolicyType, typename RhsDataType>
        typename BinaryExpressionTraits<LhsPolicyType, MultiplyOp, ConstantExpressionPolicy<RhsDataType> >::ResultType
        operator*(const Expression<LhsPolicyType>& lhs, const RhsDataType& rhs)
        {
            return lhs * Expression<ConstantExpressionPolicy<RhsDataType> >(rhs);
        }
        
//         template<typename LhsPolicyType, typename RhsPolicyType>
//         typename BinaryExpressionTraits<LhsPolicyType, MultiplyOp, RhsPolicyType>::ResultType
//         operator*(const Expression<LhsPolicyType>& lhs, const Expression<RhsPolicyType>& rhs)
//         {
//             typedef BinaryExpressionTraits<LhsPolicyType, MultiplyOp, RhsPolicyType> Traits;
//             return Traits::Apply(lhs, rhs);
//         }
        
/*        template<typename LhsExpressionPolicyType, typename RhsExpressionPolicyType>
        Expression<BinaryExpressionPolicy< LhsExpressionPolicyType, RhsExpressionPolicyType, MultiplyOp > > operator*(
            const Expression<LhsExpressionPolicyType>& lhs, const Expression<RhsExpressionPolicyType>& rhs)
        {
            typedef BinaryExpressionPolicy< LhsExpressionPolicyType, RhsExpressionPolicyType, MultiplyOp> ResultPolicyType;
            typename ResultPolicyType::DataType d(lhs, rhs);
            return Expression<ResultPolicyType>(d);
            //return Expression<BinaryExpressionPolicy< Expression<LhsExpressionPolicyType>, Expression<RhsExpressionPolicyType>, MultiplyOp > >(lhs, rhs);
        }

        
        template<typename LhsDataType, typename RhsExpressionPolicy>
        Expression<BinaryExpressionPolicy<
            ConstantExpressionPolicy<LhsDataType>,
            RhsExpressionPolicy,
            MultiplyOp > > operator*(const LhsDataType& lhs, const Expression<RhsExpressionPolicy>& rhs)
        {
            typedef ConstantExpressionPolicy<LhsDataType> LhsPolicyType;
            typedef Expression<LhsPolicyType> LhsExpressionType;
            typedef Expression<RhsExpressionPolicy> RhsExpressionType;

            typedef BinaryExpressionPolicy<LhsPolicyType, RhsExpressionPolicy, MultiplyOp> ResultPolicy;
            typename ResultPolicy::DataType d(Expression<LhsPolicyType>(lhs), rhs);
            return Expression<ResultPolicy>(d);
        }

        template<typename LhsExpressionPolicy, typename RhsDataType>
            Expression<BinaryExpressionPolicy<
            LhsExpressionPolicy,
            ConstantExpressionPolicy<RhsDataType>,
            MultiplyOp > > operator*(const Expression<LhsExpressionPolicy>& lhs, const RhsDataType& rhs)
        {
            return BinaryExpressionGenerator<LhsExpressionPolicy, MultiplyOp, ConstantExpressionPolicy<RhsDataType> >::Apply(lhs, rhs);
//             typedef Expression<LhsExpressionPolicy> LhsExpressionType;
//             typedef ConstantExpressionPolicy<RhsDataType> RhsExpressionPolicy;
//             typedef Expression<RhsExpressionPolicy> RhsExpressionType;
// 
//             typedef BinaryExpressionPolicy<LhsExpressionPolicy, RhsExpressionPolicy, MultiplyOp> ResultPolicy;
//             typename ResultPolicy::DataType d(lhs, Expression<RhsExpressionPolicy>(rhs));
//             return Expression<ResultPolicy>(d);
        }*/
        
        //////////////////////////////////////////////////////
        /// Subtraction
        //////////////////////////////////////////////////////
        template<typename LhsPolicyType, typename RhsPolicyType>
        typename BinaryExpressionTraits<LhsPolicyType, SubtractOp, RhsPolicyType>::ResultType
        operator-(const Expression<LhsPolicyType>& lhs, const Expression<RhsPolicyType>& rhs)
        {
            typedef BinaryExpressionTraits<LhsPolicyType, AddOp, RhsPolicyType> Traits;
            return Traits::Apply(lhs, rhs);
        }
        
        template<typename LhsDataType, typename RhsPolicyType>
        typename BinaryExpressionTraits<ConstantExpressionPolicy<LhsDataType>, SubtractOp, RhsPolicyType>::ResultType
        operator-(const LhsDataType& lhs, const Expression<RhsPolicyType>& rhs)
        {
            return Expression<ConstantExpressionPolicy<LhsDataType> >(lhs) - rhs;
        }
        
        template<typename LhsPolicyType, typename RhsDataType>
        typename BinaryExpressionTraits<LhsPolicyType, SubtractOp, ConstantExpressionPolicy<RhsDataType> >::ResultType
        operator-(const Expression<LhsPolicyType>& lhs, const RhsDataType& rhs)
        {
            return lhs - Expression<ConstantExpressionPolicy<RhsDataType> >(rhs);
        }
        
//         template<typename LhsPolicyType, typename RhsPolicyType>
//         typename BinaryExpressionTraits<LhsPolicyType, SubtractOp, RhsPolicyType>::ResultType
//         operator-(const Expression<LhsPolicyType>& lhs, const Expression<RhsPolicyType>& rhs)
//         {
//             typedef BinaryExpressionTraits<LhsPolicyType, SubtractOp, RhsPolicyType> Traits;
//             return Traits::Apply(lhs, rhs);
//         }
        
//         template<typename LhsExpressionPolicyType, typename RhsExpressionPolicyType>
//         Expression<BinaryExpressionPolicy< Expression<LhsExpressionPolicyType>, Expression<RhsExpressionPolicyType>, SubtractOp > > operator-(
//         const Expression<LhsExpressionPolicyType>& lhs, const Expression<RhsExpressionPolicyType>& rhs)
//         {
//             return Expression<BinaryExpressionPolicy< Expression<LhsExpressionPolicyType>, Expression<RhsExpressionPolicyType>, SubtractOp > >(lhs, rhs);
//         }
// 
//         
//         template<typename LhsDataType, typename RhsExpressionPolicy>
//         Expression<BinaryExpressionPolicy<Expression<ConstantExpressionPolicy<LhsDataType> >, Expression<RhsExpressionPolicy>, SubtractOp > > 
//         operator-(const LhsDataType& lhs, const Expression<RhsExpressionPolicy>& rhs)
//         {
//             typedef Expression<ConstantExpressionPolicy<LhsDataType> > LhsExpressionType;
//             typedef Expression<RhsExpressionPolicy> RhsExpressionType;
// 
//             return Expression<BinaryExpressionPolicy<LhsExpressionType, RhsExpressionType, SubtractOp> >(LhsExpressionType(lhs), RhsExpressionType(rhs));
//         }
// 
//         template<typename LhsExpressionPolicy, typename RhsDataType> 
//         Expression<BinaryExpressionPolicy<Expression<LhsExpressionPolicy>, Expression<ConstantExpressionPolicy<RhsDataType> >, SubtractOp > > 
//         operator-(const Expression<LhsExpressionPolicy>& lhs, const RhsDataType& rhs)
//         {
//             typedef Expression<LhsExpressionPolicy> LhsExpressionType;
//             typedef Expression<ConstantExpressionPolicy<RhsDataType> > RhsExpressionType;
// 
//             return Expression<BinaryExpressionPolicy<LhsExpressionType, RhsExpressionType, SubtractOp> >(LhsExpressionType(lhs), RhsExpressionType(rhs));
//         }
        
        //////////////////////////////////////////////////////
        /// Division
        //////////////////////////////////////////////////////
        template<typename LhsPolicyType, typename RhsPolicyType>
        typename BinaryExpressionTraits<LhsPolicyType, DivideOp, RhsPolicyType>::ResultType
        operator/(const Expression<LhsPolicyType>& lhs, const Expression<RhsPolicyType>& rhs)
        {
            typedef BinaryExpressionTraits<LhsPolicyType, DivideOp, RhsPolicyType> Traits;
            return Traits::Apply(lhs, rhs);
        }
        
        template<typename LhsDataType, typename RhsPolicyType>
        typename BinaryExpressionTraits<ConstantExpressionPolicy<LhsDataType>, DivideOp, RhsPolicyType>::ResultType
        operator/(const LhsDataType& lhs, const Expression<RhsPolicyType>& rhs)
        {
            return Expression<ConstantExpressionPolicy<LhsDataType> >(lhs) + rhs;
        }
        
        template<typename LhsPolicyType, typename RhsDataType>
        typename BinaryExpressionTraits<LhsPolicyType, DivideOp, ConstantExpressionPolicy<RhsDataType> >::ResultType
        operator/(const Expression<LhsPolicyType>& lhs, const RhsDataType& rhs)
        {
            return lhs + Expression<ConstantExpressionPolicy<RhsDataType> >(rhs);
        }
        
//         template<typename LhsPolicyType, typename RhsPolicyType>
//         typename BinaryExpressionTraits<LhsPolicyType, DivideOp, RhsPolicyType>::ResultType
//         operator/(const Expression<LhsPolicyType>& lhs, const Expression<RhsPolicyType>& rhs)
//         {
//             typedef BinaryExpressionTraits<LhsPolicyType, DivideOp, RhsPolicyType> Traits;
//             return Traits::Apply(lhs, rhs);
//         }
        
//         template<typename LhsExpressionPolicyType, typename RhsExpressionPolicyType>
//         Expression<BinaryExpressionPolicy< Expression<LhsExpressionPolicyType>, Expression<RhsExpressionPolicyType>, DivideOp > > operator/(
//         const Expression<LhsExpressionPolicyType>& lhs, const Expression<RhsExpressionPolicyType>& rhs)
//         {
//             return Expression<BinaryExpressionPolicy< Expression<LhsExpressionPolicyType>, Expression<RhsExpressionPolicyType>, DivideOp > >(lhs, rhs);
//         }
// 
//         
//         template<typename LhsDataType, typename RhsExpressionPolicy>
//         Expression<BinaryExpressionPolicy<Expression<ConstantExpressionPolicy<LhsDataType> >, Expression<RhsExpressionPolicy>, DivideOp > > 
//         operator/(const LhsDataType& lhs, const Expression<RhsExpressionPolicy>& rhs)
//         {
//             typedef Expression<ConstantExpressionPolicy<LhsDataType> > LhsExpressionType;
//             typedef Expression<RhsExpressionPolicy> RhsExpressionType;
// 
//             return Expression<BinaryExpressionPolicy<LhsExpressionType, RhsExpressionType, DivideOp> >(LhsExpressionType(lhs), RhsExpressionType(rhs));
//         }
// 
//         template<typename LhsExpressionPolicy, typename RhsDataType> 
//         Expression<BinaryExpressionPolicy<Expression<LhsExpressionPolicy>, Expression<ConstantExpressionPolicy<RhsDataType> >, DivideOp > > 
//         operator/(const Expression<LhsExpressionPolicy>& lhs, const RhsDataType& rhs)
//         {
//             typedef Expression<LhsExpressionPolicy> LhsExpressionType;
//             typedef Expression<ConstantExpressionPolicy<RhsDataType> > RhsExpressionType;
// 
//             return Expression<BinaryExpressionPolicy<LhsExpressionType, RhsExpressionType, DivideOp> >(LhsExpressionType(lhs), RhsExpressionType(rhs));
//         }
    }
    

}

#endif //NEKTAR_USE_EXPRESSION_TEMPLATES
#endif // NEKTAR_LIB_UTILITIES_BINARY_EXPRESSION_HPP

/**
    $Log: BinaryExpression.hpp,v $
    Revision 1.9  2007/01/16 05:29:50  bnelson
    Major improvements for expression templates.

    Revision 1.8  2006/11/12 17:58:47  bnelson
    *** empty log message ***

    Revision 1.7  2006/11/08 04:17:32  bnelson
    Made expressions work for complicated, nested expression templates - but suboptimally.

    Revision 1.6  2006/11/06 17:07:18  bnelson
    Continued work on creating temporaries as needed when sub-expression types don't match the type of the accumulator.

    Revision 1.5  2006/10/02 01:16:52  bnelson
    Started working on adding BLAS and LAPACK

    Revision 1.4  2006/09/14 02:08:59  bnelson
    Fixed gcc compile errors

    Revision 1.3  2006/09/11 03:24:24  bnelson
    Updated expression templates so they are all specializations of an Expression object, using policies to differentiate.

    Revision 1.2  2006/08/25 01:33:47  bnelson
    no message

    Revision 1.1  2006/06/01 09:20:55  kirby
    *** empty log message ***

    Revision 1.1  2006/05/04 18:57:40  kirby
    *** empty log message ***

    Revision 1.1  2006/01/10 14:50:30  bnelson
    Initial Revision.

**/
    
/*        template<typename LhsDataType, typename RhsLhsPolicyType, typename RhsRhsPolicyType, typename ResultType,
                 template <typename, typename> class RhsOpType, 
                 template <typename, typename> class OpType>
        class EvaluateBinaryExpression<ConstantExpressionPolicy<LhsDataType>,
                                       BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType>,
                                       ResultType, OpType, BinaryNullOp,
                                       typename boost::enable_if
                                       <
                                            boost::mpl::not_<
                                            boost::is_same
                                            <
                                                ResultType, 
                                                typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType>::ResultType
                                            > >
                                       >::type >
        {
            public:
                typedef typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType>::ResultType RhsDataType;
                
                static void Eval(const Expression<ConstantExpressionPolicy<LhsDataType> >& lhs, 
                                 const Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >& rhs,
                                 Accumulator<ResultType>& accum)
                {
                    ASSERTL0(!accum.IsInitialized(), "Accumulator must not be initialized in the lowest level binary expression evaluator.");
                    lhs.Apply(accum);
                    RhsDataType t = rhs;
                    OpType<LhsDataType, RhsDataType>::ApplyEqual(t, accum);
                }
                
        };
    
        template<typename RhsDataType, typename LhsLhsPolicyType, typename LhsRhsPolicyType, typename ResultType,
                 template <typename, typename> class LhsOpType, 
                 template <typename, typename> class OpType>
        class EvaluateBinaryExpression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType>,
                                       ConstantExpressionPolicy<RhsDataType>,
                                       ResultType, OpType, BinaryNullOp,
                                       typename boost::enable_if
                                       <
                                            boost::mpl::not_<
                                            boost::is_same
                                            <
                                                ResultType, 
                                                typename BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType>::ResultType
                                            > >
                                        >::type >
        {
            public:
            
                typedef typename BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType>::ResultType LhsDataType;
                static void Eval(const Expression<BinaryExpressionPolicy<LhsLhsPolicyType, LhsRhsPolicyType, LhsOpType> >& lhs,
                                 const Expression<ConstantExpressionPolicy<RhsDataType> >& rhs, 
                                 Accumulator<ResultType>& accum)
                {
                    ASSERTL0(!accum.IsInitialized(), "Accumulator must not be initialized in the lowest level binary expression evaluator.");
                    LhsDataType t = lhs;
                    rhs.Apply(accum);
                    OpType<LhsDataType, RhsDataType>::ApplyLeftEqual(t, accum);
                }
        };*/
            // Parent and child priority are different
        // Don't have to worry about this one at the low level - the temporary has already been created.
            
        
        // Conditions for this one.
        // 1. The accumulator and binary expression have the same type.
        // 2. 
//         template<typename LhsDataType, typename RhsLhsPolicyType, typename RhsRhsPolicyType, typename ResultType,
//                  template <typename, typename> class RhsOpType, 
//                  template <typename, typename> class OpType,
//                  template <typename, typename> class ParentOpType>
//         class EvaluateBinaryExpression<ConstantExpressionPolicy<LhsDataType>,
//                                        BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType>,
//                                        ResultType, OpType, ParentOpType,
//                                        typename boost::enable_if
//                                        <
//                                             boost::mpl::and_
//                                             <
//                                                 boost::is_same
//                                                 <
//                                                     ResultType, 
//                                                     typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType>::ResultType
//                                                 >,
//                                                 boost::mpl::not_<boost::is_same<BinaryNullOp<void, void>, ParentOpType<void, void> > >
//                                             >
//                                        >::type >
//         {
//             public:
//                 static void Eval(const Expression<ConstantExpressionPolicy<LhsDataType> >& lhs, 
//                                  const Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsRhsPolicyType, RhsOpType> >& rhs,
//                                  Accumulator<ResultType>& accum)
//                 {
//                     ASSERTL0(accum.IsInitialized(), "Accumulator must be initialized if a parent op type is provided.");
//                     rhs.Apply(accum);
//                     OpType<LhsDataType, ResultType>::ApplyLeftEqual(*lhs, accum);
//                 }
//                 
//         };
