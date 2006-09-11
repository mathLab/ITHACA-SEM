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

namespace Nektar
{
    template<typename LhsInputExpressionPolicyType, typename RhsInputExpressionPolicyType, template <typename, typename> class OpType>
    class BinaryExpressionPolicy
    {
    };

    // Originally these functions were methods inside the binary expression class.
    // Visual Studio 2003 didn't like that, so I had to move them outside the class as static methods.
    template<typename LhsExpressionType, typename RhsExpressionType, typename ResultType, template <typename, typename> class OpType>
    void EvaluateExpression(typename boost::call_traits<LhsExpressionType>::const_reference lhs, 
                                    typename boost::call_traits<RhsExpressionType>::const_reference rhs,
                                    typename boost::call_traits<ResultType>::reference result, 
                                    typename boost::enable_if<boost::is_same<typename LhsExpressionType::ResultType, ResultType> >::type* f0 = NULL,
                                    typename boost::enable_if<boost::is_same<typename RhsExpressionType::ResultType, ResultType> >::type* f1 = NULL)
    {
        typedef typename LhsExpressionType::ResultType LhsType;
        typedef typename RhsExpressionType::ResultType RhsType;
        lhs.Apply(result);
        rhs.template ApplyEqual<OpType<LhsType, RhsType> >(result);
    }

    template<typename LhsExpressionType, typename RhsExpressionType, typename ResultType, template <typename, typename> class OpType>
    void EvaluateExpression(typename boost::call_traits<LhsExpressionType>::const_reference lhs, 
                                    typename boost::call_traits<RhsExpressionType>::const_reference rhs,
                                    typename boost::call_traits<ResultType>::reference result, 
                                    typename boost::disable_if<boost::is_same<typename LhsExpressionType::ResultType, ResultType> >::type* f0 = NULL,
                                    typename boost::enable_if<boost::is_same<typename RhsExpressionType::ResultType, ResultType> >::type* f1 = NULL)
    {
        typedef typename LhsExpressionType::ResultType LhsType;
        typedef typename RhsExpressionType::ResultType RhsType;
        LhsType lhs_value = lhs;
        RhsType rhs_value = rhs;
        

        OpType<LhsType, RhsType>::Apply(result, lhs, rhs);
    }

    template<typename LhsExpressionType, typename RhsExpressionType, typename ResultType, template <typename, typename> class OpType>
    void EvaluateExpression(typename boost::call_traits<LhsExpressionType>::const_reference lhs, 
                                    typename boost::call_traits<RhsExpressionType>::const_reference rhs,
                                    typename boost::call_traits<ResultType>::reference result, 
                                    typename boost::enable_if<boost::is_same<typename LhsExpressionType::ResultType, ResultType> >::type* f0 = NULL,
                                    typename boost::disable_if<boost::is_same<typename RhsExpressionType::ResultType, ResultType> >::type* f1 = NULL)
    {
        typedef typename LhsExpressionType::ResultType LhsType;
        typedef typename RhsExpressionType::ResultType RhsType;
        lhs.Apply(result);
        RhsType rhs_value = rhs;
        OpType<LhsType, RhsType>::ApplyEqual(result, rhs);
    }

    // OpType - A class with a statis method called Apply that takes a single
    // parameter and returns a result of the same or different type.
    // A template parameter to allow a single OpType templated class to be
    // used for a variety of types.
    //
    // OpType - The operation for the expression.  
    // LhsExpressionType - The type of the lhs side of the expression.  Must be Constant, Unary, or Binary expressions.
    // RhsExpressionType - Similar, but for the rhs.
    template<typename LhsInputExpressionPolicyType, typename RhsInputExpressionPolicyType, template <typename, typename> class OpType>
    class Expression<BinaryExpressionPolicy<typename Expression<LhsInputExpressionPolicyType>, typename Expression<RhsInputExpressionPolicyType>, OpType> >
    {
        public:
            // The types of the objects being operated upon.  Not the expression types, but what the 
            // expressions will return.
            typedef Expression<BinaryExpressionPolicy<typename Expression<LhsInputExpressionPolicyType>, typename Expression<RhsInputExpressionPolicyType>, OpType> > ThisType;

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
                EvaluateExpression<LhsExpressionType, RhsExpressionType, ResultType, OpType>(m_lhs, m_rhs, result);
            }

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

}

#endif // NEKTAR_LIB_UTILITIES_BINARY_EXPRESSION_HPP

/**
    $Log: BinaryExpression.hpp,v $
    Revision 1.2  2006/08/25 01:33:47  bnelson
    no message

    Revision 1.1  2006/06/01 09:20:55  kirby
    *** empty log message ***

    Revision 1.1  2006/05/04 18:57:40  kirby
    *** empty log message ***

    Revision 1.1  2006/01/10 14:50:30  bnelson
    Initial Revision.

**/
