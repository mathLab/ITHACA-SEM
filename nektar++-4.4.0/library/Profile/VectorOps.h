
#ifndef NEKTAR_PROFILE_VECTOR_OPS_H
#define NEKTAR_PROFILE_VECTOR_OPS_H

#include <LibUtilities/LinearAlgebra/NekVector.hpp>

void AddVectors(Nektar::NekVector<double>& result,
                const Nektar::NekVector<double>& v1,
                const Nektar::NekVector<double>& v2)
{
    result = v1 + v2;
}

void AddVectors(Nektar::NekVector<double>& result,
                const Nektar::NekVector<double>& v1,
                const Nektar::NekVector<double>& v2,
                const Nektar::NekVector<double>& v3)
{
    result = v1 + v2 + v3;
}

void AddVectors(Nektar::NekVector<double>& result,
                const Nektar::NekVector<double>& v1,
                const Nektar::NekVector<double>& v2,
                const Nektar::NekVector<double>& v3,
                const Nektar::NekVector<double>& v4)
{
    result = v1 + v2 + v3 + v4;
}

void AddVectorsAccum(Nektar::NekVector<double>& result,
                const Nektar::NekVector<double>& v1,
                const Nektar::NekVector<double>& v2)
{
    result = v1;
    result += v2;
}

void AddVectorsAccum(Nektar::NekVector<double>& result,
                const Nektar::NekVector<double>& v1,
                const Nektar::NekVector<double>& v2,
                const Nektar::NekVector<double>& v3)
{
    result = v1;
    result += v2;
    result += v3;
}

void AddVectorsAccum(Nektar::NekVector<double>& result,
                const Nektar::NekVector<double>& v1,
                const Nektar::NekVector<double>& v2,
                const Nektar::NekVector<double>& v3,
                const Nektar::NekVector<double>& v4)
{
    result = v1;
    result += v2;
    result += v3;
    result += v4;
}

#endif //NEKTAR_PROFILE_VECTOR_OPS_H
