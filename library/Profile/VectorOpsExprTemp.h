
#ifndef NEKTAR_PROFILE_VECTOR_OPS_EXPR_TEMP_H
#define NEKTAR_PROFILE_VECTOR_OPS_EXPR_TEMP_H

#include <LibUtilities/LinearAlgebra/NekVectorFwd.hpp>


void AddVectorsExprTemp(Nektar::NekVector<double>& result,
                const Nektar::NekVector<double>& v1,
                const Nektar::NekVector<double>& v2);

                
#endif //NEKTAR_PROFILE_VECTOR_OPS_EXPR_TEMP_H
