
#ifndef NEKTAR_PROFILE_MATRIX_OPS_EXPR_TEMP_H
#define NEKTAR_PROFILE_MATRIX_OPS_EXPR_TEMP_H

#include <LibUtilities/LinearAlgebra/NekMatrixFwd.hpp>


void AddMatricesExprTemp(Nektar::NekMatrix<double>& result,
                const Nektar::NekMatrix<double>& v1,
                const Nektar::NekMatrix<double>& v2);

void AddMatricesExprTemp(Nektar::NekMatrix<double>& result,
                const Nektar::NekMatrix<double>& v1,
                const Nektar::NekMatrix<double>& v2,
                const Nektar::NekMatrix<double>& v3);

void AddMatricesExprTemp(Nektar::NekMatrix<double>& result,
                const Nektar::NekMatrix<double>& v1,
                const Nektar::NekMatrix<double>& v2,
                const Nektar::NekMatrix<double>& v3,
                const Nektar::NekMatrix<double>& v4);

void AddMatricesExprTempResultAlloc(
                const Nektar::NekMatrix<double>& v1,
                const Nektar::NekMatrix<double>& v2);

void AddMatricesExprTempResultAlloc(
                const Nektar::NekMatrix<double>& v1,
                const Nektar::NekMatrix<double>& v2,
                const Nektar::NekMatrix<double>& v3);

void AddMatricesExprTempResultAlloc(
                const Nektar::NekMatrix<double>& v1,
                const Nektar::NekMatrix<double>& v2,
                const Nektar::NekMatrix<double>& v3,
                const Nektar::NekMatrix<double>& v4);
                
#endif //NEKTAR_PROFILE_MATRIX_OPS_EXPR_TEMP_H
