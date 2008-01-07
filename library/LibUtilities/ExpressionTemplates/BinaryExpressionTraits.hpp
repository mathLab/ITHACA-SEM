////////////////////////////////////////////////////////////////////////////////
// 
// exprt.hpp
// Blake Nelson
//
// Generic expression templates.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_BINARY_EXPRESSION_TEMPLATES_EXPRESSION_TRAITS_H
#define NEKTAR_LIB_UTILITIES_BINARY_EXPRESSION_TEMPLATES_EXPRESSION_TRAITS_H
#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

#include <LibUtilities/ExpressionTemplates/AssociativeTraits.hpp>
#include <LibUtilities/ExpressionTemplates/CommutativeTraits.hpp>
#include <LibUtilities/ExpressionTemplates/BinaryExpressionPolicyFwd.hpp>
#include <LibUtilities/ExpressionTemplates/AssociativeExpressionTraits.hpp>
#include <LibUtilities/ExpressionTemplates/CommutativeExpressionTraits.hpp>

namespace Nektar
{
    // Aggregate information about types.
    
    template<typename LhsPolicy, template<typename, typename> class OpType, typename RhsPolicy, typename enabled = void>
    class BinaryExpressionTraits : 
        public IsCommutativeExpression<LhsPolicy, OpType, RhsPolicy>//, 
        //public AssociativeExpressionTraits<LhsPolicy, OpType, RhsPolicy>
    {
        public:
            typedef BinaryExpressionPolicy<LhsPolicy, OpType, RhsPolicy> ResultPolicy;
            typedef Expression<ResultPolicy> ResultType;
            typedef IsCommutative<LhsPolicy, OpType, RhsPolicy> IsCommutative;
            
            static ResultType Apply(const Expression<LhsPolicy>& lhs, const Expression<RhsPolicy>& rhs)
            {
                typename ResultPolicy::DataType d(lhs, rhs);
                return ResultType(d);
            }
            
    };

    
 //   // Just need two more specializations.        
//    // Binary rhs, non-constant rhs rhs
//    // Binary rhs, constant rhs rhs
//    template<typename LhsInputPolicy, typename RhsLhsPolicy, typename RhsRhsDataType> 
//    class BinaryExpressionTraits<LhsInputPolicy, SubtractOp, 
//                                    BinaryExpressionPolicy<RhsLhsPolicy, ConstantExpressionPolicy<RhsRhsDataType>, AddOp>,
//                                    typename boost::enable_if_c
//                                    <
//                                    AssociativeExpressionTraits<LhsInputPolicy, SubtractOp, BinaryExpressionPolicy<RhsLhsPolicy, ConstantExpressionPolicy<RhsRhsDataType>, AddOp> >::IsAssociativeWithOpChange,
//                                    void
//                                    >::type > :
//        public IsCommutativeExpression<LhsInputPolicy, SubtractOp, BinaryExpressionPolicy<RhsLhsPolicy, ConstantExpressionPolicy<RhsRhsDataType>, AddOp> >,
//        public AssociativeExpressionTraits<LhsInputPolicy, SubtractOp, BinaryExpressionPolicy<RhsLhsPolicy, ConstantExpressionPolicy<RhsRhsDataType>, AddOp> >
//    {
//        public:
//            typedef BinaryExpressionPolicy<RhsLhsPolicy, ConstantExpressionPolicy<RhsRhsDataType>, AddOp> RhsInputPolicy;
//            
//            typedef BinaryExpressionPolicy<LhsInputPolicy, RhsLhsPolicy, SubtractOp> LhsOutputPolicy;
//            typedef ConstantExpressionPolicy<RhsRhsDataType> RhsOutputPolicy;
//            typedef BinaryExpressionPolicy<LhsOutputPolicy, RhsOutputPolicy, SubtractOp> ResultPolicy;
//            typedef Expression<ResultPolicy> ResultType;
//            
//            static ResultType Apply(const Expression<LhsInputPolicy>& lhs, const Expression<RhsInputPolicy>& rhs)
//            {
//                typename LhsOutputPolicy::DataType lhs_data(lhs, (*rhs).first);
//                typename RhsOutputPolicy::DataType rhs_data((*rhs).second);
//                
//                typename ResultPolicy::DataType d(lhs_data, rhs_data);
//                return ResultType(d);
//            }
//    };
//    
//    
//    template<typename LhsInputPolicy, typename RhsLhsPolicy, typename RhsRhsDataType> 
//    class BinaryExpressionTraits<LhsInputPolicy, SubtractOp, 
//                                    BinaryExpressionPolicy<RhsLhsPolicy, ConstantExpressionPolicy<RhsRhsDataType>, SubtractOp>,
//                                    typename boost::enable_if_c
//                                    <
//                                    AssociativeExpressionTraits<LhsInputPolicy, SubtractOp, BinaryExpressionPolicy<RhsLhsPolicy, ConstantExpressionPolicy<RhsRhsDataType>, SubtractOp> >::IsAssociativeWithOpChange,
//                                    void
//                                    >::type > :
//        public IsCommutativeExpression<LhsInputPolicy, SubtractOp, BinaryExpressionPolicy<RhsLhsPolicy, ConstantExpressionPolicy<RhsRhsDataType>, SubtractOp> >,
//        public AssociativeExpressionTraits<LhsInputPolicy, SubtractOp, BinaryExpressionPolicy<RhsLhsPolicy, ConstantExpressionPolicy<RhsRhsDataType>, SubtractOp> >
//    {
//        public:
//            typedef BinaryExpressionPolicy<RhsLhsPolicy, ConstantExpressionPolicy<RhsRhsDataType>, SubtractOp> RhsInputPolicy;
//            
//            typedef BinaryExpressionPolicy<LhsInputPolicy, RhsLhsPolicy, AddOp> LhsOutputPolicy;
//            typedef ConstantExpressionPolicy<RhsRhsDataType> RhsOutputPolicy;
//            typedef BinaryExpressionPolicy<LhsOutputPolicy, RhsOutputPolicy, SubtractOp> ResultPolicy;
//            typedef Expression<ResultPolicy> ResultType;
//            
//            static ResultType Apply(const Expression<LhsInputPolicy>& lhs, const Expression<RhsInputPolicy>& rhs)
//            {
//                typename LhsOutputPolicy::DataType lhs_data(lhs, (*rhs).first);
//                typename RhsOutputPolicy::DataType rhs_data((*rhs).second);
//                
//                typename ResultPolicy::DataType d(lhs_data, rhs_data);
//                return ResultType(d);
//            }
//    };
//    
//    template<typename LhsInputPolicy, typename RhsLhsPolicy, typename RhsRhsDataType> 
//    class BinaryExpressionTraits<LhsInputPolicy, DivideOp, 
//                                    BinaryExpressionPolicy<RhsLhsPolicy, ConstantExpressionPolicy<RhsRhsDataType>, MultiplyOp>,
//                                    typename boost::enable_if_c
//                                    <
//                                    AssociativeExpressionTraits<LhsInputPolicy, DivideOp, BinaryExpressionPolicy<RhsLhsPolicy, ConstantExpressionPolicy<RhsRhsDataType>, MultiplyOp> >::IsAssociativeWithOpChange,
//                                    void
//                                    >::type > :
//        public IsCommutativeExpression<LhsInputPolicy, DivideOp, BinaryExpressionPolicy<RhsLhsPolicy, ConstantExpressionPolicy<RhsRhsDataType>, MultiplyOp> >,
//        public AssociativeExpressionTraits<LhsInputPolicy, DivideOp, BinaryExpressionPolicy<RhsLhsPolicy, ConstantExpressionPolicy<RhsRhsDataType>, MultiplyOp> >
//    {
//        public:
//            typedef BinaryExpressionPolicy<RhsLhsPolicy, ConstantExpressionPolicy<RhsRhsDataType>, MultiplyOp> RhsInputPolicy;
//            
//            typedef BinaryExpressionPolicy<LhsInputPolicy, RhsLhsPolicy, DivideOp> LhsOutputPolicy;
//            typedef ConstantExpressionPolicy<RhsRhsDataType> RhsOutputPolicy;
//            typedef BinaryExpressionPolicy<LhsOutputPolicy, RhsOutputPolicy, DivideOp> ResultPolicy;
//            typedef Expression<ResultPolicy> ResultType;
//            
//            static ResultType Apply(const Expression<LhsInputPolicy>& lhs, const Expression<RhsInputPolicy>& rhs)
//            {
//                typename LhsOutputPolicy::DataType lhs_data(lhs, (*rhs).first);
//                typename RhsOutputPolicy::DataType rhs_data((*rhs).second);
//                
//                typename ResultPolicy::DataType d(lhs_data, rhs_data);
//                return ResultType(d);
//            }
//    };
//    
//    template<typename LhsInputPolicy, typename RhsLhsPolicy, typename RhsRhsDataType> 
//    class BinaryExpressionTraits<LhsInputPolicy, DivideOp, 
//                                    BinaryExpressionPolicy<RhsLhsPolicy, ConstantExpressionPolicy<RhsRhsDataType>, DivideOp>,
//                                    typename boost::enable_if_c
//                                    <
//                                    AssociativeExpressionTraits<LhsInputPolicy, DivideOp, BinaryExpressionPolicy<RhsLhsPolicy, ConstantExpressionPolicy<RhsRhsDataType>, DivideOp> >::IsAssociativeWithOpChange,
//                                    void
//                                    >::type > :
//        public IsCommutativeExpression<LhsInputPolicy, DivideOp, BinaryExpressionPolicy<RhsLhsPolicy, ConstantExpressionPolicy<RhsRhsDataType>, DivideOp> >,
//        public AssociativeExpressionTraits<LhsInputPolicy, DivideOp, BinaryExpressionPolicy<RhsLhsPolicy, ConstantExpressionPolicy<RhsRhsDataType>, DivideOp> >
//    {
//        public:
//            typedef BinaryExpressionPolicy<RhsLhsPolicy, ConstantExpressionPolicy<RhsRhsDataType>, DivideOp> RhsInputPolicy;
//            
//            typedef BinaryExpressionPolicy<LhsInputPolicy, RhsLhsPolicy, MultiplyOp> LhsOutputPolicy;
//            typedef ConstantExpressionPolicy<RhsRhsDataType> RhsOutputPolicy;
//            typedef BinaryExpressionPolicy<LhsOutputPolicy, RhsOutputPolicy, DivideOp> ResultPolicy;
//            typedef Expression<ResultPolicy> ResultType;
//            
//            static ResultType Apply(const Expression<LhsInputPolicy>& lhs, const Expression<RhsInputPolicy>& rhs)
//            {
//                typename LhsOutputPolicy::DataType lhs_data(lhs, (*rhs).first);
//                typename RhsOutputPolicy::DataType rhs_data((*rhs).second);
//                
//                typename ResultPolicy::DataType d(lhs_data, rhs_data);
//                return ResultType(d);
//            }
//    };
//    
//    
    template<typename LhsType, template <typename, typename> class OpType, typename RhsType>
    class BinaryExpressionType
    {
        public:
            typedef Expression<BinaryExpressionPolicy<ConstantExpressionPolicy<LhsType>, OpType, ConstantExpressionPolicy<RhsType> > > Type;
    };
    
    template<template <typename, typename> class OpType, typename LhsType, typename RhsType>
    typename BinaryExpressionType<LhsType, OpType, RhsType>::Type 
    CreateBinaryExpression(const LhsType& lhs, const RhsType& rhs)
    {
        typedef ConstantExpressionPolicy<LhsType> LhsPolicy;
        typedef ConstantExpressionPolicy<RhsType> RhsPolicy;
        typedef BinaryExpressionPolicy<LhsPolicy, OpType, RhsPolicy> ResultPolicy; 
        typedef Expression<ResultPolicy> ResultExpressionType; 
        
        // Note the () around the first parameter.  Failure to do this will result in a compiler error.
        // For more information, lookup Scott Myers discussion on C++'s most vexing parse.
        Expression<LhsPolicy> lhsTemp(lhs);
        Expression<RhsPolicy> rhsTemp(rhs);
        typename ResultPolicy::DataType d(lhsTemp, rhsTemp); 
        //typename ResultPolicy::DataType d((Expression<LhsPolicy>(lhs)), Expression<RhsPolicy>(rhs)); 
        return ResultExpressionType(d); 
    }
}

#endif //NEKTAR_USE_EXPRESSION_TEMPLATES
#endif // NEKTAR_LIB_UTILITIES_BINARY_EXPRESSION_TEMPLATES_EXPRESSION_TRAITS_H

/**
    $Log: BinaryExpressionTraits.hpp,v $
    Revision 1.13  2007/10/04 03:48:54  bnelson
    *** empty log message ***

    Revision 1.12  2007/10/03 02:58:03  bnelson
    *** empty log message ***

    Revision 1.11  2007/08/16 02:14:21  bnelson
    Moved expression templates to the Nektar namespace.

    Revision 1.10  2007/07/20 00:36:05  bnelson
    *** empty log message ***

    Revision 1.9  2007/03/14 21:24:07  sherwin
    Update for working version of MultiRegions up to ExpList1D

    Revision 1.8  2007/01/29 01:34:05  bnelson
    Updates to compile in windows.

    Revision 1.7  2007/01/16 17:37:55  bnelson
    Wrapped everything with #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

    Revision 1.6  2007/01/16 05:29:50  bnelson
    Major improvements for expression templates.

    Revision 1.5  2006/11/06 17:07:19  bnelson
    Continued work on creating temporaries as needed when sub-expression types don't match the type of the accumulator.

    Revision 1.4  2006/09/14 02:08:59  bnelson
    Fixed gcc compile errors

    Revision 1.3  2006/09/11 03:24:24  bnelson
    Updated expression templates so they are all specializations of an Expression object, using policies to differentiate.

    Revision 1.2  2006/08/25 01:33:47  bnelson
    no message

    Revision 1.1  2006/06/01 09:20:55  kirby
    *** empty log message ***

    Revision 1.1  2006/05/04 18:57:41  kirby
    *** empty log message ***

    Revision 1.1  2006/01/10 14:50:30  bnelson
    Initial Revision.

**/

