
#ifndef NEKTAR_PROFILE_VECTOR_OPS_H
#define NEKTAR_PROFILE_VECTOR_OPS_H

#include <LibUtilities/LinearAlgebra/NekVector.hpp>

void AddVectors(Nektar::NekVector<double>& result,
                const Nektar::NekVector<double>& v1,
                const Nektar::NekVector<double>& v2)
{
    result = v1 + v2;
}
                
void AddVectorsAccum(Nektar::NekVector<double>& result,
                const Nektar::NekVector<double>& v1,
                const Nektar::NekVector<double>& v2)
{
    result = v1;
    result += v2;
}
                
#endif //NEKTAR_PROFILE_VECTOR_OPS_H
