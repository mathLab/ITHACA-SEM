////////////////////////////////////////////////////////////////////////////////
// 
// ExpressionTemplateOperators.hpp
// Blake Nelson
//
// Generic expression template operators.
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BinaryExpressionTraits.hpp>
#include <LibUtilities/ExpressionTemplateConcepts.hpp>

#include <boost/concept_check.hpp>

namespace Nektar
{
    template<typename LhsType, typename RhsType>
	class OpAdd
	{
		public:
            BOOST_CLASS_REQUIRE2(LhsType, RhsType, Nektar, OpAddConcept);

            typedef typename GetReturnType<LhsType>::result_type lhs_result_type;
            typedef typename GetReturnType<RhsType>::result_type rhs_result_type;

            typedef typename BinaryExpressionTraits<lhs_result_type, rhs_result_type>::op_add_result_type result_type;

		public:
			static typename boost::call_traits<result_type>::reference 
            apply(typename boost::call_traits<LhsType>::param_type lhs,
				  typename boost::call_traits<RhsType>::param_type rhs, 
                  typename boost::call_traits<result_type>::reference result)
			{
				op_add(lhs, rhs, result);
				return result;
			}
	};



    template<typename LhsType, typename RhsType>
	class OpMultiply
	{
		public:
            BOOST_CLASS_REQUIRE2(LhsType, RhsType, Nektar, OpMultiplyConcept);

            typedef typename GetReturnType<LhsType>::result_type lhs_result_type;
            typedef typename GetReturnType<RhsType>::result_type rhs_result_type;

            typedef typename BinaryExpressionTraits<lhs_result_type, rhs_result_type>::op_multiply_result_type result_type;

		public:
			static typename boost::call_traits<result_type>::reference 
            apply(typename boost::call_traits<LhsType>::param_type lhs,
				  typename boost::call_traits<RhsType>::param_type rhs, 
                  typename boost::call_traits<result_type>::reference result)
			{
				op_multiply(lhs, rhs, result);
				return result;
			}
	};



    template<typename LhsType, typename RhsType>
	class OpDivide
	{
		public:
            BOOST_CLASS_REQUIRE2(LhsType, RhsType, Nektar, OpDivideConcept);

            typedef typename GetReturnType<LhsType>::result_type lhs_result_type;
            typedef typename GetReturnType<RhsType>::result_type rhs_result_type;

            typedef typename BinaryExpressionTraits<lhs_result_type, rhs_result_type>::op_divide_result_type result_type;

		public:
			static typename boost::call_traits<result_type>::reference 
            apply(typename boost::call_traits<LhsType>::param_type lhs,
				  typename boost::call_traits<RhsType>::param_type rhs, 
                  typename boost::call_traits<result_type>::reference result)
			{
				op_divide(lhs, rhs, result);
				return result;
			}
	};


    template<typename LhsType, typename RhsType>
	class OpSubtract
	{
		public:
            BOOST_CLASS_REQUIRE2(LhsType, RhsType, Nektar, OpSubtractConcept);
            typedef typename GetReturnType<LhsType>::result_type lhs_result_type;
            typedef typename GetReturnType<RhsType>::result_type rhs_result_type;

            typedef typename BinaryExpressionTraits<lhs_result_type, rhs_result_type>::op_subtract_result_type result_type;

		public:
			static typename boost::call_traits<result_type>::reference 
            apply(typename boost::call_traits<LhsType>::param_type lhs,
				  typename boost::call_traits<RhsType>::param_type rhs, 
                  typename boost::call_traits<result_type>::reference result)
			{
				op_subtract(lhs, rhs, result);
				return result;
			}
	};

    // Operators - These do the real work in constructing expression templates.
	template<typename LhsType, typename RhsType>
	BinaryExpression<OpAdd, ConstantExpression<LhsType>, ConstantExpression<RhsType> >
	operator+(const ConstantExpression<LhsType>& lhs, const ConstantExpression<RhsType>& rhs)
	{
        return BinaryExpression<OpAdd, ConstantExpression<LhsType>, ConstantExpression<RhsType> >(lhs, rhs);
	}

    template<template <typename,typename> class BinaryOpType, typename BinaryLhsType, typename BinaryRhsType, typename RhsType>
	BinaryExpression<OpAdd, BinaryExpression<BinaryOpType, BinaryLhsType, BinaryRhsType>, ConstantExpression<RhsType> >
	operator+(const BinaryExpression<BinaryOpType, BinaryLhsType, BinaryRhsType>& lhs, const ConstantExpression<RhsType>& rhs)
	{
        return BinaryExpression<OpAdd, BinaryExpression<BinaryOpType, BinaryLhsType, BinaryRhsType>, ConstantExpression<RhsType> >(lhs, rhs);
	}


    template<typename LhsType, typename RhsType>
	BinaryExpression<OpMultiply, ConstantExpression<LhsType>, ConstantExpression<RhsType> >
	operator*(const ConstantExpression<LhsType>& lhs, const ConstantExpression<RhsType>& rhs)
	{
        return BinaryExpression<OpMultiply, ConstantExpression<LhsType>, ConstantExpression<RhsType> >(lhs, rhs);
	}

    template<template <typename,typename> class BinaryOpType, typename BinaryLhsType, typename BinaryRhsType, typename RhsType>
	BinaryExpression<OpMultiply, BinaryExpression<BinaryOpType, BinaryLhsType, BinaryRhsType>, ConstantExpression<RhsType> >
	operator*(const BinaryExpression<BinaryOpType, BinaryLhsType, BinaryRhsType>& lhs, const ConstantExpression<RhsType>& rhs)
	{
        return BinaryExpression<OpMultiply, BinaryExpression<BinaryOpType, BinaryLhsType, BinaryRhsType>, ConstantExpression<RhsType> >(lhs, rhs);
	}





    template<typename LhsType, typename RhsType>
	BinaryExpression<OpSubtract, ConstantExpression<LhsType>, ConstantExpression<RhsType> >
	operator-(const ConstantExpression<LhsType>& lhs, const ConstantExpression<RhsType>& rhs)
	{
        return BinaryExpression<OpSubtract, ConstantExpression<LhsType>, ConstantExpression<RhsType> >(lhs, rhs);
	}

 //   template<template <typename, typename> class BinaryOpType, typename BinaryLhsType, typename BinaryRhsType, typename RhsType>
	//BinaryExpression<OpSubtractEqual, BinaryExpression<BinaryOpType, BinaryLhsType, BinaryRhsType>, ConstantExpression<RhsType> >
	//operator-(const BinaryExpression<BinaryOpType, BinaryLhsType, BinaryRhsType>& lhs, const ConstantExpression<RhsType>& rhs)
	//{
 //       return BinaryExpression<OpSubtractEqual, BinaryExpression<BinaryOpType, BinaryLhsType, BinaryRhsType>, ConstantExpression<RhsType> >(lhs, rhs);
	//}





    template<typename LhsType, typename RhsType>
	BinaryExpression<OpDivide, ConstantExpression<LhsType>, ConstantExpression<RhsType> >
	operator/(const ConstantExpression<LhsType>& lhs, const ConstantExpression<RhsType>& rhs)
	{
        return BinaryExpression<OpDivide, ConstantExpression<LhsType>, ConstantExpression<RhsType> >(lhs, rhs);
	}

 //   template<template <typename,typename> class BinaryOpType, typename BinaryLhsType, typename BinaryRhsType, typename RhsType>
	//BinaryExpression<OpDivideEqual, BinaryExpression<BinaryOpType, BinaryLhsType, BinaryRhsType>, ConstantExpression<RhsType> >
	//operator/(const BinaryExpression<BinaryOpType, BinaryLhsType, BinaryRhsType>& lhs, const ConstantExpression<RhsType>& rhs)
	//{
 //       return BinaryExpression<OpDivideEqual, BinaryExpression<BinaryOpType, BinaryLhsType, BinaryRhsType>, ConstantExpression<RhsType> >(lhs, rhs);
	//}
}

/**
    $Log: ExpressionTemplateOperators.hpp,v $
    Revision 1.1  2006/05/04 18:57:42  kirby
    *** empty log message ***

    Revision 1.2  2006/01/31 13:51:12  bnelson
    Updated for new configure.

    Revision 1.1  2006/01/10 14:50:31  bnelson
    Initial Revision.

**/

