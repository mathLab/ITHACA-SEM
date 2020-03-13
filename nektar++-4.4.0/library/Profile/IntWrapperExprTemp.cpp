//
//#define NEKTAR_USE_EXPRESSION_TEMPLATES
//
//#include "IntWrapperExprTemp.h"
//#include <LibUtilities/ExpressionTemplates/Expression.hpp>
//#include <LibUtilities/ExpressionTemplates/ExpressionTemplates.hpp>
//#include <LibUtilities/ExpressionTemplates/ExpressionMetadata.hpp>
//#include <LibUtilities/ExpressionTemplates/ConstantExpressionTraits.hpp>
//#include <LibUtilities/ExpressionTemplates/BinaryOperators.hpp>
//#include <LibUtilities/ExpressionTemplates/BinaryExpressionTraits.hpp>
//#include <LibUtilities/BasicUtils/OperatorGenerators.hpp>
//#include <iostream>
//
//using namespace std;
//using namespace Nektar;
//
//void NekAdd(expr::IntWrapper& result, const expr::IntWrapper& lhs, const expr::IntWrapper& rhs)
//{
//    result.SetValue(lhs.GetValue() + rhs.GetValue());
//}
//    
//void NekAddEqual(expr::IntWrapper& result, const expr::IntWrapper& rhs)
//{
//    result += rhs;
//}
//
//expr::IntWrapper NekAdd(const expr::IntWrapper& lhs, const expr::IntWrapper& rhs)
//{
//    return expr::IntWrapper(lhs.GetValue() + rhs.GetValue());
//}
//        
//namespace expr
//{
//        
//    
//    GENERATE_ADDITION_OPERATOR(IntWrapper, 0, IntWrapper, 0);
//
//
//    void AddIntWrapperExprTemp(IntWrapper& result, const IntWrapper& r1, const IntWrapper& r2)
//    {
//        Expression<BinaryExpressionPolicy<ConstantExpressionPolicy<IntWrapper>, AddOp, ConstantExpressionPolicy<IntWrapper> > > exp =
//                r1 + r2;
//        result = exp;
//    }
//
//    void AddIntWrapperExprTemp(IntWrapper& result, const IntWrapper& r1, const IntWrapper& r2,
//                      const IntWrapper& r3)
//    {
//        result = r1 + r2 + r3;
//    }
//
//    void AddIntWrapperExprTemp(IntWrapper& result, const IntWrapper& r1, const IntWrapper& r2,
//                      const IntWrapper& r3,
//                      const IntWrapper& r4)
//    {
//        result = r1 + r2 + r3 + r4;
//    }
//
//    void AddIntWrapperExprTemp(IntWrapper& result, const IntWrapper& r1, const IntWrapper& r2,
//                      const IntWrapper& r3,
//                      const IntWrapper& r4,
//                      const IntWrapper& r5)
//    {
//        result = r1 + r2 + r3 + r4 + r5;
//    }
//
//    void AddIntWrapperExprTemp(IntWrapper& result, const IntWrapper& r1, const IntWrapper& r2,
//                      const IntWrapper& r3,
//                      const IntWrapper& r4,
//                      const IntWrapper& r5,
//                      const IntWrapper& r6)
//    {
//        result = r1 + r2 + r3 + r4 + r5 + r6;
//    }
//
//    void AddIntWrapperExprTemp(IntWrapper& result, const IntWrapper& r1, const IntWrapper& r2,
//                      const IntWrapper& r3,
//                      const IntWrapper& r4,
//                      const IntWrapper& r5,
//                      const IntWrapper& r6,
//                      const IntWrapper& r7)
//    {
//        result = r1 + r2 + r3 + r4 + r5 + r6 + r7;
//    }
//}
