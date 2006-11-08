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

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>

#include <LibUtilities/ExpressionTemplates/BinaryExpressionOperators.hpp>
#include <LibUtilities/ExpressionTemplates/NullOp.hpp>
#include <iostream>

namespace Nektar
{
    namespace expt
    {
        template<typename LhsInputExpressionPolicyType, typename RhsInputExpressionPolicyType, template <typename, typename> class OpType>
        class BinaryExpressionPolicy
        {
        };

        
        // The evaluate expression class is really the Apply method for the BinaryExpression.  I needed to move it outside the class
        // to use the enable_if feature.
        
        /// OpType is the operation that combines the lhs and rhs expression (lhs op rhs)
        
        // One thing to consider - it may not be valid to accumulate the lhs of the expression into the accumulator.
        //
        // Example, A + (B*C).  In this case.  Consider the evaluation of B*C.  The accumulator has the value of A,
        // We need to create a temporary and then add it in.
        //
        // So we need a concept of incoming operation type.
        
        /// LhsExpressionPolicyType - The type of the lhs expression.
        /// RhsExpressionPolicyType - The type of the rhs expression.

                                                
        template<typename LhsExpressionPolicyType, typename RhsExpressionPolicyType, typename ResultType, template <typename, typename> class OpType, 
                   typename ParentOpType = NullOp,
                   typename LhsNeedsTemp = void, typename RhsNeedsTemp = void>
        class EvaluateBinaryExpression
        {
            public:
                static void eval(const Expression<LhsExpressionPolicyType>& lhs, 
                                 const Expression<RhsExpressionPolicyType>& rhs,
                                 typename boost::call_traits<ResultType>::reference result)
                {
//                     typedef typename Expression<LhsExpressionPolicyType>::ResultType LhsType;
//                     typedef typename Expression<RhsExpressionPolicyType>::ResultType RhsType;
//                     lhs.Apply(result);
//                     rhs.template ApplyEqual<OpType<LhsType, RhsType> >(result);
                    typedef typename Expression<LhsExpressionPolicyType>::ResultType LhsType;
                    typedef typename Expression<RhsExpressionPolicyType>::ResultType RhsType;
                    LhsType lhs_temp(lhs);
                    RhsType rhs_temp(rhs);

                    OpType<LhsType, RhsType>::Apply(result, lhs_temp, rhs_temp);
                }
        };
        
        // The requirement for a temporary is either
        // 1. One of the sides returns a data type different than the accumulator type.
        // 2. The priority of the expression's operator is different than the parent operator.
        template<typename LhsExpressionPolicyType, typename RhsExpressionPolicyType, typename ResultType, template <typename, typename> class OpType>
        class EvaluateBinaryExpression<LhsExpressionPolicyType, RhsExpressionPolicyType, ResultType, OpType, 
                                       NullOp,
                                       typename boost::disable_if<boost::is_same<typename Expression<LhsExpressionPolicyType>::ResultType, ResultType> >::type,
                                       typename boost::enable_if<boost::is_same<typename Expression<RhsExpressionPolicyType>::ResultType, ResultType> >::type>
        {
            public:
                static void eval(const Expression<LhsExpressionPolicyType>& lhs, 
                                 const Expression<RhsExpressionPolicyType>& rhs,
                                 typename boost::call_traits<ResultType>::reference result)
                {
//                     typedef typename Expression<LhsExpressionPolicyType>::ResultType LhsType;
//                     typedef typename Expression<RhsExpressionPolicyType>::ResultType RhsType;
//                     LhsType lhs_temp;
//                     lhs.Apply(lhs_temp);
//                     rhs.Apply(result);

                    //OpType<LhsType, RhsType>::ApplyLhsEqual(result, lhs_temp);
                    typedef typename Expression<LhsExpressionPolicyType>::ResultType LhsType;
                    typedef typename Expression<RhsExpressionPolicyType>::ResultType RhsType;
                    LhsType lhs_temp(lhs);
                    RhsType rhs_temp(rhs);

                    OpType<LhsType, RhsType>::Apply(result, lhs_temp, rhs_temp);
                }
        };
        
        template<typename LhsExpressionPolicyType, typename RhsExpressionPolicyType, typename ResultType, template <typename, typename> class OpType, typename ParentOpType>
        class EvaluateBinaryExpression<LhsExpressionPolicyType, RhsExpressionPolicyType, ResultType, OpType, 
                                        ParentOpType,
                                        typename boost::disable_if<boost::is_same<typename Expression<LhsExpressionPolicyType>::ResultType, ResultType> >::type,
                                        typename boost::enable_if<boost::is_same<typename Expression<RhsExpressionPolicyType>::ResultType, ResultType> >::type>
        {
            public:
                static void eval(const Expression<LhsExpressionPolicyType>& lhs, 
                                    const Expression<RhsExpressionPolicyType>& rhs,
                                    typename boost::call_traits<ResultType>::reference result)
                {
                    typedef typename Expression<LhsExpressionPolicyType>::ResultType LhsType;
                    typedef typename Expression<RhsExpressionPolicyType>::ResultType RhsType;
                    LhsType lhs_temp(lhs);
                    RhsType rhs_temp(rhs);

                    OpType<LhsType, RhsType>::Apply(result, lhs_temp, rhs_temp);
                }
        };
        
        template<typename LhsExpressionPolicyType, typename RhsExpressionPolicyType, typename ResultType, template <typename, typename> class OpType>
        class EvaluateBinaryExpression<LhsExpressionPolicyType, RhsExpressionPolicyType, ResultType, OpType, NullOp,
                                       typename boost::enable_if<boost::is_same<typename Expression<LhsExpressionPolicyType>::ResultType, ResultType> >::type,
                                       typename boost::disable_if<boost::is_same<typename Expression<RhsExpressionPolicyType>::ResultType, ResultType> >::type>
        {
            public:
                static void eval(const Expression<LhsExpressionPolicyType>& lhs, 
                                    const Expression<RhsExpressionPolicyType>& rhs,
                                    typename boost::call_traits<ResultType>::reference result)
                {
                    typedef typename Expression<LhsExpressionPolicyType>::ResultType LhsType;
                    typedef typename Expression<RhsExpressionPolicyType>::ResultType RhsType;
                    LhsType lhs_temp(lhs);
                    RhsType rhs_temp(rhs);

                    OpType<LhsType, RhsType>::Apply(result, lhs_temp, rhs_temp);
                }
        };
                

        // OpType - A class with a statis method called Apply that takes a single
        // parameter and returns a result of the same or different type.
        // A template parameter to allow a single OpType templated class to be
        // used for a variety of types.
        //
        // OpType - The operation for the expression.  
        // LhsExpressionType - The type of the lhs side of the expression.  Must be Constant, Unary, or Binary expressions.
        // RhsExpressionType - Similar, but for the rhs.
        template<typename LhsInputExpressionPolicyType, typename RhsInputExpressionPolicyType, template <typename, typename> class OpType>
        class Expression<BinaryExpressionPolicy<Expression<LhsInputExpressionPolicyType>, Expression<RhsInputExpressionPolicyType>, OpType> >
        {
            public:
                // The types of the objects being operated upon.  Not the expression types, but what the 
                // expressions will return.
                typedef Expression<BinaryExpressionPolicy<Expression<LhsInputExpressionPolicyType>, Expression<RhsInputExpressionPolicyType>, OpType> > ThisType;

                typedef Expression<LhsInputExpressionPolicyType> LhsExpressionType;
                typedef Expression<RhsInputExpressionPolicyType> RhsExpressionType;

                typedef typename LhsExpressionType::ResultType LhsResultType;
                typedef typename RhsExpressionType::ResultType RhsResultType;

                typedef typename OpType<LhsResultType, RhsResultType>::ResultType ResultType;
                typedef typename ExpressionMetadataChooser<ResultType>::MetadataType MetadataType;
                typedef typename LhsExpressionType::MetadataType LhsMetadataType;
                typedef typename RhsExpressionType::MetadataType RhsMetadataType;

            public:
                Expression(const LhsExpressionType& lhs, const RhsExpressionType& rhs) :
                    m_lhs(lhs),
                    m_rhs(rhs),
                    m_metadata(OpType<LhsResultType, RhsResultType>::CreateBinaryMetadata(lhs.GetMetadata(), rhs.GetMetadata()))
                {
                }

                Expression(const ThisType& rhs) :
                    m_lhs(rhs.m_lhs),
                    m_rhs(rhs.m_rhs),
                    m_metadata(rhs.m_metadata)
                {
                }


                ~Expression() {}

                //static void EvaluateExpression(typename boost::call_traits<LhsExpressionType>::const_reference lhs, 
                //                               typename boost::call_traits<RhsExpressionType>::const_reference rhs,
                //                               typename boost::call_traits<ResultType>::reference result, 
                //                               typename boost::enable_if<boost::is_same<LhsParameterType, ResultType> >::type* f0 = NULL,
                //                               typename boost::enable_if<boost::is_same<RhsParameterType, ResultType> >::type* f1 = NULL)
                //{
                //    lhs.Apply(result);
                //}

                //static void EvaluateExpression(typename boost::call_traits<LhsExpressionType>::const_reference lhs, 
                //                               typename boost::call_traits<RhsExpressionType>::const_reference rhs,
                //                               typename boost::call_traits<ResultType>::reference result, 
                //                               typename boost::disable_if<boost::is_same<LhsParameterType, ResultType> >::type* f0 = NULL,
                //                               typename boost::enable_if<boost::is_same<RhsParameterType, ResultType> >::type* f1 = NULL)
                //{
                //    LhsParameterType temp(lhs.GetMetadata());
                //    lhs.Apply(temp);
                //}

                // Different cases for apply:
                // All data types (result, lhs, rhs) are the same.
                //   1.  Operations are serial (such as matrix m1 + m2 + m3 + m4).  In this case apply for both sides works.
                //   
                // All same data type, apply can go through as normal (such as matrix m1 + m2 + m3)
                // All same data type, apply must create a temporary (such as matrix m1*m2 + m3*m4*m5)
                //      in this case m3*m4*m5 must have a temporary created.
                // Two cases for the apply method.
                // 1.  Result and Parameter types are the same.
                // 2.  Result and Parameter types are different.
                void Apply(typename boost::call_traits<ResultType>::reference result) const
                {
                    EvaluateBinaryExpression<LhsInputExpressionPolicyType, RhsInputExpressionPolicyType, ResultType, OpType>::eval(m_lhs, m_rhs, result);
                }

                /// ParentOpType is the operation one step higher in the parse tree.
                /// For example, if the expression is "A + BC", then then the OpType for 
                /// this operation is "*" and the ParentOpType is +.
//                 template<typename ParentOpType>
//                 void ApplyEqual(typename boost::call_traits<ResultType>::reference result) const
//                 {
//                     // In this case, we can directly apply the operator to the lhs with the accumulator.
//                     // This will be in cases like "A + B + C".
//                     m_lhs.template  ApplyEqual<ParentOpType>(result);
//                     OpType::ApplyEqual(result, m_value);
//                 }
                
                const MetadataType& GetMetadata() const
                {
                    return m_metadata;
                }

            private:
                ThisType& operator=(const ThisType& rhs);

                LhsExpressionType m_lhs;
                RhsExpressionType m_rhs;

                MetadataType m_metadata;
        };

        
        //////////////////////////////////////////////////////
        /// Addition
        //////////////////////////////////////////////////////
        template<typename LhsExpressionPolicyType, typename RhsExpressionPolicyType>
        Expression<BinaryExpressionPolicy< Expression<LhsExpressionPolicyType>, Expression<RhsExpressionPolicyType>, AddOp > > operator+(
            const Expression<LhsExpressionPolicyType>& lhs, const Expression<RhsExpressionPolicyType>& rhs)
        {
            return Expression<BinaryExpressionPolicy< Expression<LhsExpressionPolicyType>, Expression<RhsExpressionPolicyType>, AddOp > >(lhs, rhs);
        }

        
        template<typename LhsDataType, typename RhsExpressionPolicy>
        Expression<BinaryExpressionPolicy<
            Expression<ConstantExpressionPolicy<LhsDataType> >,
            Expression<RhsExpressionPolicy>,
            AddOp > > operator+(
            const LhsDataType& lhs,
            const Expression<RhsExpressionPolicy>& rhs)
        {
            typedef Expression<ConstantExpressionPolicy<LhsDataType> > LhsExpressionType;
            typedef Expression<RhsExpressionPolicy> RhsExpressionType;

            return Expression<BinaryExpressionPolicy<LhsExpressionType, RhsExpressionType, AddOp> >(LhsExpressionType(lhs), RhsExpressionType(rhs));
        }

        template<typename LhsExpressionPolicy, typename RhsDataType>
        Expression<BinaryExpressionPolicy<
            Expression<LhsExpressionPolicy>,
            Expression<ConstantExpressionPolicy<RhsDataType> >,
            AddOp > > operator+(
            const Expression<LhsExpressionPolicy>& lhs,
            const RhsDataType& rhs)
        {
            typedef Expression<LhsExpressionPolicy> LhsExpressionType;
            typedef Expression<ConstantExpressionPolicy<RhsDataType> > RhsExpressionType;

            return Expression<BinaryExpressionPolicy<LhsExpressionType, RhsExpressionType, AddOp> >(LhsExpressionType(lhs), RhsExpressionType(rhs));
        }
        
        ////////////////////////////////
        // Multiplication
        ///////////////////////////////
        template<typename LhsExpressionPolicyType, typename RhsExpressionPolicyType>
        Expression<BinaryExpressionPolicy< Expression<LhsExpressionPolicyType>, Expression<RhsExpressionPolicyType>, MultiplyOp > > operator*(
            const Expression<LhsExpressionPolicyType>& lhs, const Expression<RhsExpressionPolicyType>& rhs)
        {
            return Expression<BinaryExpressionPolicy< Expression<LhsExpressionPolicyType>, Expression<RhsExpressionPolicyType>, MultiplyOp > >(lhs, rhs);
        }

        
        template<typename LhsDataType, typename RhsExpressionPolicy>
        Expression<BinaryExpressionPolicy<
            Expression<ConstantExpressionPolicy<LhsDataType> >,
            Expression<RhsExpressionPolicy>,
            MultiplyOp > > operator*(const LhsDataType& lhs, const Expression<RhsExpressionPolicy>& rhs)
        {
            typedef Expression<ConstantExpressionPolicy<LhsDataType> > LhsExpressionType;
            typedef Expression<RhsExpressionPolicy> RhsExpressionType;

            return Expression<BinaryExpressionPolicy<LhsExpressionType, RhsExpressionType, MultiplyOp> >(LhsExpressionType(lhs), RhsExpressionType(rhs));
        }

        template<typename LhsExpressionPolicy, typename RhsDataType>
            Expression<BinaryExpressionPolicy<
            Expression<LhsExpressionPolicy>,
            Expression<ConstantExpressionPolicy<RhsDataType> >,
            MultiplyOp > > operator*(const Expression<LhsExpressionPolicy>& lhs, const RhsDataType& rhs)
        {
            typedef Expression<LhsExpressionPolicy> LhsExpressionType;
            typedef Expression<ConstantExpressionPolicy<RhsDataType> > RhsExpressionType;

            return Expression<BinaryExpressionPolicy<LhsExpressionType, RhsExpressionType, MultiplyOp> >(LhsExpressionType(lhs), RhsExpressionType(rhs));
        }
        
        //////////////////////////////////////////////////////
        /// Subtraction
        //////////////////////////////////////////////////////
        template<typename LhsExpressionPolicyType, typename RhsExpressionPolicyType>
        Expression<BinaryExpressionPolicy< Expression<LhsExpressionPolicyType>, Expression<RhsExpressionPolicyType>, SubtractOp > > operator-(
        const Expression<LhsExpressionPolicyType>& lhs, const Expression<RhsExpressionPolicyType>& rhs)
        {
            return Expression<BinaryExpressionPolicy< Expression<LhsExpressionPolicyType>, Expression<RhsExpressionPolicyType>, SubtractOp > >(lhs, rhs);
        }

        
        template<typename LhsDataType, typename RhsExpressionPolicy>
        Expression<BinaryExpressionPolicy<Expression<ConstantExpressionPolicy<LhsDataType> >, Expression<RhsExpressionPolicy>, SubtractOp > > 
        operator-(const LhsDataType& lhs, const Expression<RhsExpressionPolicy>& rhs)
        {
            typedef Expression<ConstantExpressionPolicy<LhsDataType> > LhsExpressionType;
            typedef Expression<RhsExpressionPolicy> RhsExpressionType;

            return Expression<BinaryExpressionPolicy<LhsExpressionType, RhsExpressionType, SubtractOp> >(LhsExpressionType(lhs), RhsExpressionType(rhs));
        }

        template<typename LhsExpressionPolicy, typename RhsDataType> 
        Expression<BinaryExpressionPolicy<Expression<LhsExpressionPolicy>, Expression<ConstantExpressionPolicy<RhsDataType> >, SubtractOp > > 
        operator-(const Expression<LhsExpressionPolicy>& lhs, const RhsDataType& rhs)
        {
            typedef Expression<LhsExpressionPolicy> LhsExpressionType;
            typedef Expression<ConstantExpressionPolicy<RhsDataType> > RhsExpressionType;

            return Expression<BinaryExpressionPolicy<LhsExpressionType, RhsExpressionType, SubtractOp> >(LhsExpressionType(lhs), RhsExpressionType(rhs));
        }
        
        //////////////////////////////////////////////////////
        /// Division
        //////////////////////////////////////////////////////
        template<typename LhsExpressionPolicyType, typename RhsExpressionPolicyType>
        Expression<BinaryExpressionPolicy< Expression<LhsExpressionPolicyType>, Expression<RhsExpressionPolicyType>, DivideOp > > operator/(
        const Expression<LhsExpressionPolicyType>& lhs, const Expression<RhsExpressionPolicyType>& rhs)
        {
            return Expression<BinaryExpressionPolicy< Expression<LhsExpressionPolicyType>, Expression<RhsExpressionPolicyType>, DivideOp > >(lhs, rhs);
        }

        
        template<typename LhsDataType, typename RhsExpressionPolicy>
        Expression<BinaryExpressionPolicy<Expression<ConstantExpressionPolicy<LhsDataType> >, Expression<RhsExpressionPolicy>, DivideOp > > 
        operator/(const LhsDataType& lhs, const Expression<RhsExpressionPolicy>& rhs)
        {
            typedef Expression<ConstantExpressionPolicy<LhsDataType> > LhsExpressionType;
            typedef Expression<RhsExpressionPolicy> RhsExpressionType;

            return Expression<BinaryExpressionPolicy<LhsExpressionType, RhsExpressionType, DivideOp> >(LhsExpressionType(lhs), RhsExpressionType(rhs));
        }

        template<typename LhsExpressionPolicy, typename RhsDataType> 
        Expression<BinaryExpressionPolicy<Expression<LhsExpressionPolicy>, Expression<ConstantExpressionPolicy<RhsDataType> >, DivideOp > > 
        operator/(const Expression<LhsExpressionPolicy>& lhs, const RhsDataType& rhs)
        {
            typedef Expression<LhsExpressionPolicy> LhsExpressionType;
            typedef Expression<ConstantExpressionPolicy<RhsDataType> > RhsExpressionType;

            return Expression<BinaryExpressionPolicy<LhsExpressionType, RhsExpressionType, DivideOp> >(LhsExpressionType(lhs), RhsExpressionType(rhs));
        }
    }
}

#endif // NEKTAR_LIB_UTILITIES_BINARY_EXPRESSION_HPP

/**
    $Log: BinaryExpression.hpp,v $
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
    
