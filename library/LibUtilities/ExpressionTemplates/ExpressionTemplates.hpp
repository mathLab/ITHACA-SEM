////////////////////////////////////////////////////////////////////////////////
// 
// exprt.hpp
// Blake Nelson
//
// Generic expression templates.
//
// Three steps to integrate your classes with the expression templates.
//
// 1.  Provide a copy constructor and an assignment operator for your class which 
//     takes a binary expression as the argument.
//
// 2.  Define a BinaryExpressionTraits object for each possible combination of 
//     operators using your type.
//
// 3.  Define op_add, op_subtract, op_divide, op_multiply.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_HPP
#define NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_HPP

#include <LibUtilities/ExpressionTemplates/ExpressionReturnType.hpp>
#include <LibUtilities/ExpressionTemplates/ExpressionTemplateConcepts.hpp>
#include <LibUtilities/ExpressionTemplates/BinaryExpressionTraits.hpp>
#include <LibUtilities/ExpressionTemplates/ConstantExpression.hpp>
#include <LibUtilities/ExpressionTemplates/BinaryExpression.hpp>
#include <LibUtilities/ExpressionTemplates/ExpressionTemplateOperators.hpp>

#include <boost/call_traits.hpp>
#include <boost/concept_check.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>

#include <iostream>

namespace Nektar
{
	//// Utilities and functors.
	//template<typename LhsType, typename RhsType>
	//class OpAddEqual
	//{
 //       public:
 //           BOOST_CLASS_REQUIRE2(LhsType, RhsType, Nektar, OpAddEqualConcept);
 //           typedef typename GetReturnType<LhsType>::result_type result_type;

	//	public:
	//		typedef typename boost::call_traits<LhsType>::reference reference_type;
	//		typedef typename boost::call_traits<RhsType>::param_type param_type;

	//	public:
	//		static reference_type apply(reference_type lhs, param_type rhs) 
	//		{
	//			return lhs += rhs;
	//		}
	//};

	//template<typename LhsType, typename RhsType>
	//class OpSubtractEqual
	//{
	//	public:
	//		typedef typename boost::call_traits<LhsType>::reference reference_type;
	//		typedef typename boost::call_traits<RhsType>::param_type param_type;
 //           typedef typename GetReturnType<LhsType>::result_type result_type;

	//	public:
	//		static reference_type apply(reference_type lhs, param_type rhs) 
	//		{
	//			return lhs -= rhs;
	//		}
	//};

	//template<typename LhsType, typename RhsType>
	//class OpMultiplyEqual
	//{
	//	public:
	//		typedef typename boost::call_traits<LhsType>::reference reference_type;
	//		typedef typename boost::call_traits<RhsType>::param_type param_type;
 //           typedef typename GetReturnType<LhsType>::result_type result_type;

	//	public:
	//		static reference_type apply(reference_type lhs, param_type rhs) 
	//		{
	//			return lhs *= rhs;
	//		}
	//};

	//template<typename LhsType, typename RhsType>
	//class OpDivideEqual
	//{
	//	public:
	//		typedef typename boost::call_traits<LhsType>::reference reference_type;
	//		typedef typename boost::call_traits<RhsType>::param_type param_type;
 //           typedef typename GetReturnType<LhsType>::result_type result_type;

	//	public:
	//		static reference_type apply(reference_type lhs, param_type rhs) 
	//		{
	//			return lhs /= rhs;
	//		}
	//};





 //   template<typename ExpressionType>
 //   class ExpressionResultType
 //   {
 //       public:
 //           typedef ExpressionType result_type;
 //   };

 //   template<typename DataType>
 //   class ExpressionResultType<ConstantExpression<DataType> >
 //   {
 //       public:
 //           typedef typename ConstantExpression<DataType>::result_type result_type;
 //   };

 //   template<typename ResultType, typename LhsType, typename RhsType, typename OpType>
 //   class ExpressionResultType<BinaryExpression<ResultType, LhsType, RhsType, OpType> >
 //   {
 //       public:
 //           typedef typename BinaryExpression<ResultType, LhsType, RhsType, OpType>::result_type result_type;
 //   };


	



    //// Generic binary expression.
    //template<typename OpType, typename LhsType, typename RhsType>
    //class Expression<2, typename expr_traits<LhsType, RhsType>::result_type, OpType, LhsType, RhsType>
    //{
    //    //public:
    //    //    typedef typename expr_traits<LhsType, RhsType>::result_type result_type;

    //    //public:
    //    //    static void apply(typename boost::call_traits<LhsType>::param_type lhs,
    //    //        typename boost::call_traits<RhsType>::param_type rhs, 
    //    //        typename boost::call_traits<result_type>::reference result)
    //    //    {
    //    //        OpType::apply(lhs, rhs, result);7
    //    //    }
    //};

    //template<template <typename, typename> class OpType, typename LhsType, typename RhsType>
    //class BinaryExpression<OpType, ConstantExpression<LhsType>, RhsType>
    //{
    //    public:
    //        typedef typename GetReturnType<LhsType>::result_type lhs_result_type;
    //        typedef typename GetReturnType<RhsType>::result_type rhs_result_type;

    //        typedef typename OpType<lhs_result_type, rhs_result_type>::result_type result_type;

    //    public:
    //        BinaryExpression(const ConstantExpression<LhsType>& lhs, typename boost::call_traits<RhsType>::param_type rhs) :
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

    //        // Case 1 - lhs result type is the same as the final result type.
    //        void apply(typename boost::call_traits<result_type>::reference result, 
    //            typename boost::enable_if<boost::is_same<lhs_result_type, result_type> >::type* param = 0) const
    //        {
    //            // This assumes the the result type of the lhs is the same as the result type of 
    //            // the rhs.
    //            m_lhs.apply(result);
    //            OpType<typename LhsType::result_type, RhsType>::apply(result, m_rhs.getValue(), result);
    //        }

    //        // Case 2 - rhs result type is the same as the final result type and the lhs is not.
    //        void apply(typename boost::call_traits<result_type>::reference result, 
    //            typename boost::disable_if<boost::is_same<lhs_result_type, result_type> >::type* param1 = 0,
    //            typename boost::enable_if<boost::is_same<rhs_result_type, result_type> >::type* param2 = 0) const
    //        {
    //            // This assumes the the result type of the lhs is the same as the result type of 
    //            // the rhs.
    //            //m_lhs.apply(result);
    //            //m_rhs.apply(result);
    //            //OpType<typename LhsType::result_type, RhsType>::apply(result, m_rhs.getValue(), result);
    //        }

    //    private:
    //        ConstantExpression<LhsType> m_lhs;
    //        RhsType m_rhs;
    //};




	// Utilities.
	template<typename DataType>
	inline ConstantExpression<DataType> make_expr(const DataType& p)
	{
		return ConstantExpression<DataType>(p);
	}


}


/// BACKUP
	//// Utilities and functors.
	//template<typename DataType>
	//class OpPlusEqual
	//{
	//	public:
	//		typedef typename boost::call_traits<DataType>::reference reference_type;
	//		typedef typename boost::call_traits<DataType>::param_type param_type;

	//	public:
	//		static reference_type apply(reference_type lhs, param_type rhs) 
	//		{
	//			return lhs += rhs;
	//		}
	//};

	//template<typename DataType>
	//class OpMinusEqual
	//{
	//	public:
	//		typedef typename boost::call_traits<DataType>::reference reference_type;
	//		typedef typename boost::call_traits<DataType>::param_type param_type;

	//	public:
	//		static reference_type apply(reference_type lhs, param_type rhs) 
	//		{
	//			return lhs -= rhs;
	//		}
	//};

	//template<typename DataType>
	//class OpTimesEqual
	//{
	//	public:
	//		typedef typename boost::call_traits<DataType>::reference reference_type;
	//		typedef typename boost::call_traits<DataType>::param_type param_type;

	//	public:
	//		static reference_type apply(reference_type lhs, param_type rhs) 
	//		{
	//			return lhs *= rhs;
	//		}
	//};

	//template<typename DataType>
	//class OpDivideEqual
	//{
	//	public:
	//		typedef typename boost::call_traits<DataType>::reference reference_type;
	//		typedef typename boost::call_traits<DataType>::param_type param_type;

	//	public:
	//		static reference_type apply(reference_type lhs, param_type rhs) 
	//		{
	//			return lhs /= rhs;
	//		}
	//};

	//template<typename ResultType, typename LhsType, typename RhsType>
	//class OpPlus
	//{
	//	public:

	//	public:
	//		static typename boost::call_traits<ResultType>::reference apply(typename boost::call_traits<LhsType>::param_type lhs,
	//			typename boost::call_traits<RhsType>::param_type rhs, typename boost::call_traits<ResultType>::reference result)
	//		{
	//			// op_plus is a user defined function along the lines of operator+
	//			op_plus(lhs, rhs, result);
	//			return result;
	//		}
	//};


	//template<typename ResultType>
	//class ConstantExpression 
	//{
	//	public:
	//		typedef ResultType result_type;

	//	public:
	//		explicit ConstantExpression(typename boost::call_traits<ResultType>::param_type value) :
	//			m_value(value)
	//		{
	//		}

	//		ConstantExpression(const ConstantExpression<ResultType>& rhs) :
	//			m_value(rhs.m_value)
	//		{
	//		}

	//		~ConstantExpression() {}

	//		ConstantExpression<ResultType>& operator=(const ConstantExpression<ResultType>& rhs)
	//		{
	//			m_value = rhs.m_value;
	//			return *this;
	//		}

	//		void operator()(typename boost::call_traits<ResultType>::reference result) const
	//		{
	//			result = m_value;
	//		}

	//		template<typename OpType>
	//		void evaluate(typename boost::call_traits<ResultType>::reference result) const
	//		{
	//			OpType::apply(result, m_value);
	//		}

	//		const ResultType& getValue() const { return m_value; }

	//	private:
	//		const ResultType& m_value;
	//};


 //   template<typename ResultType, typename LhsExpressionType, typename RhsExpressionType, typename OpType>
 //   class BinaryExpression
 //   {
 //       public:
 //           typedef ResultType result_type;

 //       public:
 //           BinaryExpression(typename boost::call_traits<LhsExpressionType>::param_type lhs, typename boost::call_traits<RhsExpressionType>::param_type rhs) :
 //               m_lhs(lhs),
 //               m_rhs(rhs)
 //           {
 //           }

 //           BinaryExpression(const BinaryExpression& rhs) :
 //               m_lhs(rhs.m_lhs),
 //               m_rhs(rhs.m_rhs)
 //           {
 //           }

 //           ~BinaryExpression() {}

 //           BinaryExpression& operator=(const BinaryExpression& rhs)
 //           {
 //               m_lhs = rhs.m_lhs;
 //               m_rhs = rhs.m_rhs;
 //               return *this;
 //           }

 //           void operator()(typename boost::call_traits<ResultType>::reference result) const
 //           {
 //               // TODO - This is only one of the cases, but good for testing.
 //               m_lhs(result);
 //               OpType::apply(result, m_rhs.getValue(), result);
 //           }

 //       private:
 //           LhsExpressionType m_lhs;
 //           RhsExpressionType m_rhs;
 //   };

 //   template<typename ResultType, typename LhsType, typename RhsType, typename OpType>
 //   class BinaryExpression<ResultType, ConstantExpression<LhsType>, ConstantExpression<RhsType>, OpType >
 //   {
 //       public:
 //           typedef ResultType result_type;

 //       public:
 //           BinaryExpression(const ConstantExpression<LhsType>& lhs, const ConstantExpression<RhsType>& rhs) :
 //               m_lhs(lhs),
 //               m_rhs(rhs)
 //           {
 //           }

 //           BinaryExpression(const BinaryExpression& rhs) :
 //               m_lhs(rhs.m_lhs),
 //               m_rhs(rhs.m_rhs)
 //           {
 //           }

 //           ~BinaryExpression() {}

 //           BinaryExpression& operator=(const BinaryExpression& rhs)
 //           {
 //               m_lhs = rhs.m_lhs;
 //               m_rhs = rhs.m_rhs;
 //               return *this;
 //           }

 //           void operator()(typename boost::call_traits<ResultType>::reference result) const
 //           {
 //               OpType::apply(m_lhs.getValue(), m_rhs.getValue(), result);
 //           }

 //       private:
 //           ConstantExpression<LhsType> m_lhs;
 //           ConstantExpression<RhsType> m_rhs;
 //   };

 //   template<typename ExpressionType>
 //   class ExpressionResultType
 //   {
 //       public:
 //           typedef ExpressionType result_type;
 //   };

 //   template<typename DataType>
 //   class ExpressionResultType<ConstantExpression<DataType> >
 //   {
 //       public:
 //           typedef typename ConstantExpression<DataType>::result_type result_type;
 //   };

 //   template<typename ResultType, typename LhsType, typename RhsType, typename OpType>
 //   class ExpressionResultType<BinaryExpression<ResultType, LhsType, RhsType, OpType> >
 //   {
 //       public:
 //           typedef typename BinaryExpression<ResultType, LhsType, RhsType, OpType>::result_type result_type;
 //   };

 //   // The only requirement on ExpressionType is that it support operator().
 //   template<typename ExpressionType>
	//class Expression
 //   {
 //       public:
 //           typedef typename ExpressionResultType<ExpressionType>::result_type result_type;

 //       public:
 //           Expression(typename boost::call_traits<ExpressionType>::param_type expr) :
 //               m_expression(expr)
 //           {
 //           }

 //           Expression(const Expression<ExpressionType>& rhs) :
 //               m_expression(rhs.m_expression)
 //           {
 //           }

 //           ~Expression()
 //           {
 //           }

 //           Expression<ExpressionType>& operator=(const Expression<ExpressionType>& rhs)
 //           {
 //               m_expression = rhs.m_expression;
 //           }

 //           void operator()(typename boost::call_traits<result_type>::reference result) const
 //           {
 //               m_expression(result);
 //           }

 //       private:
 //           ExpressionType m_expression;
 //   };

	//
	//// Expression template traits.
	//template<typename LhsType, typename RhsType>
	//class expr_traits_base
	//{
	//};

	//template<typename LhsType, typename RhsType>
	//class expr_traits
	//{
	//};

	//// Operators - These do the real work in constructing expression templates.
	//template<typename LhsType, typename RhsType>
	//Expression<BinaryExpression<typename expr_traits<LhsType, RhsType>::op_plus_result_type, ConstantExpression<LhsType>, ConstantExpression<RhsType>, OpPlus<typename expr_traits<LhsType, RhsType>::op_plus_result_type , LhsType, RhsType> > >
	//operator+(const ConstantExpression<LhsType>& lhs, const ConstantExpression<RhsType>& rhs)
	//{
	//	return BinaryExpression<typename expr_traits<LhsType, RhsType>::op_plus_result_type, ConstantExpression<LhsType>, ConstantExpression<RhsType>, OpPlus<typename expr_traits<LhsType, RhsType>::op_plus_result_type , LhsType, RhsType> >(lhs, rhs);
	//}

 //   template<typename LhsExpressionType, typename RhsType>
 //   Expression<BinaryExpression<typename expr_traits<typename LhsExpressionType::result_type, RhsType>::op_plus_result_type, LhsExpressionType, ConstantExpression<RhsType>, OpPlus<typename expr_traits<typename LhsExpressionType::result_type, RhsType>::op_plus_result_type , typename LhsExpressionType::result_type, RhsType> > >
	//operator+(const LhsExpressionType& lhs, const ConstantExpression<RhsType>& rhs)
	//{
	//	return BinaryExpression<typename expr_traits<typename LhsExpressionType::result_type, RhsType>::op_plus_result_type, LhsExpressionType, ConstantExpression<RhsType>, OpPlus<typename expr_traits<typename LhsExpressionType::result_type, RhsType>::op_plus_result_type , typename LhsExpressionType::result_type, RhsType> >(lhs, rhs);
	//}

	//// Utilities.
	//template<typename DataType>
	//ConstantExpression<DataType> make_expr(const DataType& p)
	//{
	//	return ConstantExpression<DataType>(p);
	//}
// }

#endif // NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_HPP

/**
    $Log: ExpressionTemplates.hpp,v $
    Revision 1.1  2006/06/01 09:20:56  kirby
    *** empty log message ***

    Revision 1.1  2006/05/04 18:57:42  kirby
    *** empty log message ***

    Revision 1.2  2006/01/31 13:51:12  bnelson
    Updated for new configure.

    Revision 1.1  2006/01/10 14:50:31  bnelson
    Initial Revision.

**/
