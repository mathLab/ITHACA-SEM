
#define NEKTAR_USE_EXPRESSION_TEMPLATES

#include <LibUtilities/ExpressionTemplates/Expression.hpp>
#include <LibUtilities/ExpressionTemplates/ExpressionTemplates.hpp>
#include <LibUtilities/ExpressionTemplates/ExpressionMetadata.hpp>
#include <LibUtilities/ExpressionTemplates/ConstantExpressionTraits.hpp>
#include <LibUtilities/ExpressionTemplates/BinaryOperators.hpp>
#include <LibUtilities/ExpressionTemplates/BinaryExpressionTraits.hpp>
#include <iostream>

#include "MatrixOpsExprTemp.h"
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>

namespace Nektar
{
//    template<typename MatrixType>
//    class BinaryExpressionEvaluator<BinaryExpressionPolicy<ConstantExpressionPolicy<MatrixType>, AddOp, ConstantExpressionPolicy<MatrixType> >,
//                                    BinaryExpressionPolicy<ConstantExpressionPolicy<MatrixType>, AddOp, ConstantExpressionPolicy<MatrixType> >,
//                                    MatrixType, AddOp, BinaryNullOp>
//    {
//        public:
//            typedef BinaryExpressionPolicy<ConstantExpressionPolicy<MatrixType>, AddOp, ConstantExpressionPolicy<MatrixType> > PolicyType;
//            static void Eval(const Expression<PolicyType>& lhs, 
//                             const Expression<PolicyType>& rhs,
//                             Accumulator<MatrixType>& result)
//            {
//                const MatrixType& m1 = (*lhs).first;
//                const MatrixType& m2 = (*lhs).second;
//                const MatrixType& m3 = (*rhs).first;
//                const MatrixType& m4 = (*rhs).second;
//                
//                const typename MatrixType::NumberType* b1 = m1.GetRawPtr();
//                const typename MatrixType::NumberType* b2 = m2.GetRawPtr();
//                const typename MatrixType::NumberType* b3 = m3.GetRawPtr();
//                const typename MatrixType::NumberType* b4 = m4.GetRawPtr();
//                
//                unsigned int size = m1.GetStorageSize();
//
//                typename MatrixType::NumberType* r = (*result).GetRawPtr();
//                typename MatrixType::NumberType* end = r+size;
//                while( r < end )
//                {
//                    *r = *b1 + *b2 + *b3 + *b4;
//                    ++b1;
//                    ++b2;
//                    ++b3;
//                    ++b4;
//                    ++r;
//                } 
//            }
//    };
}

void AddMatricesExprTemp(Nektar::NekMatrix<double>& result,
                const Nektar::NekMatrix<double>& v1,
                const Nektar::NekMatrix<double>& v2)
{
    result = MakeExpr(v1) + MakeExpr(v2);
}

void AddMatricesExprTemp(Nektar::NekMatrix<double>& result,
                const Nektar::NekMatrix<double>& v1,
                const Nektar::NekMatrix<double>& v2,
                const Nektar::NekMatrix<double>& v3)
{
    result = MakeExpr(v1) + MakeExpr(v2) + MakeExpr(v3);
}

void AddMatricesExprTemp(Nektar::NekMatrix<double>& result,
                const Nektar::NekMatrix<double>& v1,
                const Nektar::NekMatrix<double>& v2,
                const Nektar::NekMatrix<double>& v3,
                const Nektar::NekMatrix<double>& v4)
{
    result = MakeExpr(v1) + MakeExpr(v2) + MakeExpr(v3) + MakeExpr(v4);
}

void AddMatricesExprTempResultAlloc(
                const Nektar::NekMatrix<double>& v1,
                const Nektar::NekMatrix<double>& v2)
{
    Nektar::NekMatrix<double> result = MakeExpr(v1) + MakeExpr(v2);
}

void AddMatricesExprTempResultAlloc(
                const Nektar::NekMatrix<double>& v1,
                const Nektar::NekMatrix<double>& v2,
                const Nektar::NekMatrix<double>& v3)
{
    Nektar::NekMatrix<double> result = MakeExpr(v1) + MakeExpr(v2) + MakeExpr(v3);
}

void AddMatricesExprTempResultAlloc(
                const Nektar::NekMatrix<double>& v1,
                const Nektar::NekMatrix<double>& v2,
                const Nektar::NekMatrix<double>& v3,
                const Nektar::NekMatrix<double>& v4)
{
    Nektar::NekMatrix<double> result = (MakeExpr(v1) + MakeExpr(v2)) + (MakeExpr(v3) + MakeExpr(v4));
}
