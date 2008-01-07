
#define NEKTAR_USE_EXPRESSION_TEMPLATES

#include "IntWrapper.h"
#include <LibUtilities/ExpressionTemplates/Expression.hpp>
#include <LibUtilities/ExpressionTemplates/ExpressionTemplates.hpp>
#include <LibUtilities/ExpressionTemplates/ExpressionMetadata.hpp>
#include <LibUtilities/ExpressionTemplates/ConstantExpressionTraits.hpp>
#include <LibUtilities/ExpressionTemplates/BinaryOperators.hpp>
#include <LibUtilities/ExpressionTemplates/BinaryExpressionTraits.hpp>
#include <iostream>

using namespace std;
using namespace Nektar;

void NekAdd(IntWrapper& result, const IntWrapper& lhs, const IntWrapper& rhs)
{
    result.SetValue(lhs.GetValue() + rhs.GetValue());
}
    
void NekAddEqual(IntWrapper& result, const IntWrapper& rhs)
{
    result += rhs;
}

IntWrapper NekAdd(const IntWrapper& lhs, const IntWrapper& rhs)
{
    return IntWrapper(lhs.GetValue() + rhs.GetValue());
}

// Expression<BinaryExpressionPolicy<ConstantExpressionPolicy<IntWrapper>, ConstantExpressionPolicy<IntWrapper>, AddOp > >
// operator+(const IntWrapper& lhs, const IntWrapper& rhs)
// {
//     return CreateBinaryExpression<AddOp>(lhs, rhs);
// }
       

// template<>
// class ConstantExpressionTraits<IntWrapper>
// {
//     public:
//         typedef IntWrapper result_type;
//         typedef StringMetadata MetadataType;
// };

// template<>
// class BinaryExpressionMetadataTraits<IntWrapper, IntWrapper, Nektar::AddOp>
// {
//     public:
//         typedef StringMetadata MetadataType;
// };



    
// IntWrapper operator+(const IntWrapper& lhs, const IntWrapper& rhs)
// {
//     return IntWrapper(lhs.GetValue() + rhs.GetValue());
// }

void AddIntWrapperExprTemp(IntWrapper& result, const IntWrapper& r1, const IntWrapper& r2)
{
    Expression<BinaryExpressionPolicy<ConstantExpressionPolicy<IntWrapper>, AddOp, ConstantExpressionPolicy<IntWrapper> > > exp =
            MakeExpr(r1) + MakeExpr(r2);
    result = exp;
}

void AddIntWrapperExprTemp(IntWrapper& result, const IntWrapper& r1, const IntWrapper& r2,
                  const IntWrapper& r3)
{
    result = MakeExpr(r1) + MakeExpr(r2) + MakeExpr(r3);
}

void AddIntWrapperExprTemp(IntWrapper& result, const IntWrapper& r1, const IntWrapper& r2,
                  const IntWrapper& r3,
                  const IntWrapper& r4)
{
    result = MakeExpr(r1) + MakeExpr(r2) + MakeExpr(r3) + MakeExpr(r4);
}

void AddIntWrapperExprTemp(IntWrapper& result, const IntWrapper& r1, const IntWrapper& r2,
                  const IntWrapper& r3,
                  const IntWrapper& r4,
                  const IntWrapper& r5)
{
    result = MakeExpr(r1) + MakeExpr(r2) + MakeExpr(r3) + MakeExpr(r4) + MakeExpr(r5);
}

void AddIntWrapperExprTemp(IntWrapper& result, const IntWrapper& r1, const IntWrapper& r2,
                  const IntWrapper& r3,
                  const IntWrapper& r4,
                  const IntWrapper& r5,
                  const IntWrapper& r6)
{
    result = MakeExpr(r1) + MakeExpr(r2) + MakeExpr(r3) + MakeExpr(r4) + MakeExpr(r5) + MakeExpr(r6);
}

void AddIntWrapperExprTemp(IntWrapper& result, const IntWrapper& r1, const IntWrapper& r2,
                  const IntWrapper& r3,
                  const IntWrapper& r4,
                  const IntWrapper& r5,
                  const IntWrapper& r6,
                  const IntWrapper& r7)
{
    result = MakeExpr(r1) + MakeExpr(r2) + MakeExpr(r3) + MakeExpr(r4) + MakeExpr(r5) + MakeExpr(r6) + MakeExpr(r7);
}

