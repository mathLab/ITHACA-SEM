//////////////////////////////////////////////////////////////////////////////////
////
//// BinaryExpression.hpp
//// Blake Nelson
////
//// Binary expression templates.
////
//////////////////////////////////////////////////////////////////////////////////
//
//#ifndef NEKTAR_LIB_UTILITIES_BINARY_EXPRESSION_HPP
//#define NEKTAR_LIB_UTILITIES_BINARY_EXPRESSION_HPP
//
//#include <boost/utility/enable_if.hpp>
//#include <boost/type_traits.hpp>
//
//namespace Nektar
//{
//    // To evaluate an expression, I need to know if the lhs or rhs of the binary tree
//    // matches the overall result type.  This can be done with enable_if and disable_if,
//    // but for some reason the work for global functions but not member functions.  
//    // Therefore, these global functions exist solely so the binary functions can forward 
//    // the evaluation.
//    
//    template<template <typename, typename> class OpType, typename LhsType, typename RhsType>
//    class BinaryExpression
//    {
//        
//    };
//
//    template<template <typename, typename> class OpType, typename LhsType, typename RhsType>
//    class BinaryExpression<OpType, ConstantExpression<LhsType>, ConstantExpression<RhsType> >
//    {
//        public:
//            typedef typename GetReturnType<LhsType>::result_type lhs_result_type;
//            typedef typename GetReturnType<RhsType>::result_type rhs_result_type;
//
//            typedef typename OpType<lhs_result_type, rhs_result_type>::result_type result_type;
//
//        public:
//            BinaryExpression(const ConstantExpression<LhsType>& lhs, const ConstantExpression<RhsType>& rhs) :
//                m_lhs(lhs),
//                m_rhs(rhs)
//            {
//            }
//
//            BinaryExpression(const BinaryExpression& rhs) :
//                m_lhs(rhs.m_lhs),
//                m_rhs(rhs.m_rhs)
//            {
//            }
//
//            ~BinaryExpression() {}
//
//            BinaryExpression& operator=(const BinaryExpression& rhs)
//            {
//                m_lhs = rhs.m_lhs;
//                m_rhs = rhs.m_rhs;
//            }
//
//            void apply(typename boost::call_traits<result_type>::reference result) const
//            {
//                OpType<lhs_result_type, rhs_result_type>::apply(m_lhs.getValue(), m_rhs.getValue(), result);
//            }
//
//        private:
//            ConstantExpression<LhsType> m_lhs;
//            ConstantExpression<RhsType> m_rhs;
//    };
//
//    // Case 1 - The lhs does match the result type.
//    template<template <typename, typename> class OpType, typename LhsType, typename RhsType, typename ResultType>
//    void evaluate_expression(typename boost::call_traits<LhsType>::param_type lhs,
//                             typename boost::call_traits<RhsType>::param_type rhs,
//                             typename boost::call_traits<ResultType>::reference result,
//                             typename boost::enable_if<boost::is_same<typename GetReturnType<LhsType>::result_type, ResultType> >::type* v = 0)
//    {
//        typedef typename GetReturnType<LhsType>::result_type lhs_result_type;
//        typedef typename GetReturnType<RhsType>::result_type rhs_result_type;
//
//        lhs.apply(result);
//        OpType<lhs_result_type, rhs_result_type>::apply(result, rhs.getValue(), result);
//    }
//
//    template<template <typename, typename> class OpType, typename LhsType, typename RhsType, typename ResultType>
//    void evaluate_expression(typename boost::call_traits<LhsType>::param_type lhs,
//                             typename boost::call_traits<RhsType>::param_type rhs,
//                             typename boost::call_traits<ResultType>::reference result,
//                             typename boost::disable_if<boost::is_same<typename GetReturnType<LhsType>::result_type, ResultType> >::type* v = 0)
//    {
//        typedef typename GetReturnType<LhsType>::result_type LhsResultType;
//        typedef typename GetReturnType<RhsType>::result_type RhsResultType;
//
//        LhsResultType temp;
//        lhs.apply(temp);
//        OpType<LhsResultType, RhsResultType>::apply(temp, rhs.getValue(), result);
//    }
//
//    template<template <typename, typename> class OpType, typename LhsType, typename RhsDataType>
//    class BinaryExpression<OpType, LhsType, ConstantExpression<RhsDataType> >
//    {
//        public:
//            typedef ConstantExpression<RhsDataType> RhsType;
//            typedef typename GetReturnType<LhsType>::result_type LhsResultType;
//            typedef typename GetReturnType<RhsType>::result_type RhsResultType;
//
//            typedef typename OpType<LhsResultType, RhsResultType>::result_type ResultType;
//
//        public:
//            BinaryExpression(typename boost::call_traits<LhsType>::param_type lhs, 
//                             typename boost::call_traits<RhsType>::param_type rhs) :
//                m_lhs(lhs),
//                m_rhs(rhs)
//            {
//            }
//
//            BinaryExpression(const BinaryExpression& rhs) :
//                m_lhs(rhs.m_lhs),
//                m_rhs(rhs.m_rhs)
//            {
//            }
//
//            ~BinaryExpression() {}
//
//            BinaryExpression& operator=(const BinaryExpression& rhs)
//            {
//                m_lhs = rhs.m_lhs;
//                m_rhs = rhs.m_rhs;
//            }
//
//            // Case 1 - The return type of the lhs is the same as the result, so we can accumulate.
//            void apply(typename boost::call_traits<ResultType>::reference result) const
//            {
//                evaluate_expression<OpType, LhsType, RhsType, ResultType>(m_lhs, m_rhs, result);
//            }
//
//            //// Case 2 - rhs result type is the same as the final result type and the lhs is not.
//            //void apply(typename boost::call_traits<result_type>::reference result, 
//            //    typename boost::disable_if<boost::is_same<int, int> >::type* param1 = 0,
//            //    double t = 0) const
//            //{
//            //    // This assumes the the result type of the lhs is the same as the result type of 
//            //    // the rhs.
//            //    //m_lhs.apply(result);
//            //    //m_rhs.apply(result);
//            //    //OpType<typename LhsType::result_type, RhsType>::apply(result, m_rhs.getValue(), result);
//            //}
//
//        private:
//            LhsType m_lhs;
//            RhsType m_rhs;
//    };
//
//    template<template <typename, typename> class OpType, typename LhsType, typename RhsType>
//    class GetReturnType<BinaryExpression<OpType, LhsType, RhsType> >
//    {
//        public:
//            typedef typename BinaryExpression<OpType, LhsType, RhsType>::result_type result_type;
//    };
//}
//
//#endif // NEKTAR_LIB_UTILITIES_BINARY_EXPRESSION_HPP
//
///**
//    $Log: BinaryExpression.hpp,v $
//    Revision 1.1  2006/06/01 09:20:55  kirby
//    *** empty log message ***
//
//    Revision 1.1  2006/05/04 18:57:40  kirby
//    *** empty log message ***
//
//    Revision 1.1  2006/01/10 14:50:30  bnelson
//    Initial Revision.
//
//**/
