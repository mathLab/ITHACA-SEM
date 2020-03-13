
#ifndef NEKTAR_PROFILE_MATRIX_OPS_H
#define NEKTAR_PROFILE_MATRIX_OPS_H

#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>

void AddMatrices(Nektar::NekMatrix<double>& result,
                const Nektar::NekMatrix<double>& v1,
                const Nektar::NekMatrix<double>& v2)
{
    result = v1 + v2;
}

void AddMatrices(Nektar::NekMatrix<double>& result,
                const Nektar::NekMatrix<double>& v1,
                const Nektar::NekMatrix<double>& v2,
                const Nektar::NekMatrix<double>& v3)
{
    result = v1 + v2 + v3;
}

void AddMatrices(Nektar::NekMatrix<double>& result,
                const Nektar::NekMatrix<double>& v1,
                const Nektar::NekMatrix<double>& v2,
                const Nektar::NekMatrix<double>& v3,
                const Nektar::NekMatrix<double>& v4)
{
    result = v1 + v2 + v3 + v4;
}

void AddMatricesResultAlloc(
                const Nektar::NekMatrix<double>& v1,
                const Nektar::NekMatrix<double>& v2)
{
    Nektar::NekMatrix<double> result = v1 + v2;
}

void AddMatricesResultAlloc(
                const Nektar::NekMatrix<double>& v1,
                const Nektar::NekMatrix<double>& v2,
                const Nektar::NekMatrix<double>& v3)
{
    Nektar::NekMatrix<double> result = v1 + v2 + v3;
}

void AddMatricesResultAlloc(
                const Nektar::NekMatrix<double>& v1,
                const Nektar::NekMatrix<double>& v2,
                const Nektar::NekMatrix<double>& v3,
                const Nektar::NekMatrix<double>& v4)
{
    Nektar::NekMatrix<double> result = v1 + v2 + v3 + v4;
}

                
void AddMatricesAccum(Nektar::NekMatrix<double>& result,
                const Nektar::NekMatrix<double>& v1,
                const Nektar::NekMatrix<double>& v2)
{
    result = v1;
    NekAddEqual(result, v2);
}

void AddMatricesAccum(Nektar::NekMatrix<double>& result,
                const Nektar::NekMatrix<double>& v1,
                const Nektar::NekMatrix<double>& v2,
                const Nektar::NekMatrix<double>& v3)
{
    result = v1;
    NekAddEqual(result, v2);
    NekAddEqual(result, v3);
}

void AddMatricesAccum(Nektar::NekMatrix<double>& result,
                const Nektar::NekMatrix<double>& v1,
                const Nektar::NekMatrix<double>& v2,
                const Nektar::NekMatrix<double>& v3,
                const Nektar::NekMatrix<double>& v4)
{
    result = v1;
    NekAddEqual(result, v2);
    NekAddEqual(result, v3);
    NekAddEqual(result, v4);
}

void AddMatricesAccumResultAlloc(
                const Nektar::NekMatrix<double>& v1,
                const Nektar::NekMatrix<double>& v2)
{
    Nektar::NekMatrix<double> result = v1;
    NekAddEqual(result, v2);
}

void AddMatricesAccumResultAlloc(
                const Nektar::NekMatrix<double>& v1,
                const Nektar::NekMatrix<double>& v2,
                const Nektar::NekMatrix<double>& v3)
{
    Nektar::NekMatrix<double> result = v1;
    NekAddEqual(result, v2);
    NekAddEqual(result, v3);
}

void AddMatricesAccumResultAlloc(
                const Nektar::NekMatrix<double>& v1,
                const Nektar::NekMatrix<double>& v2,
                const Nektar::NekMatrix<double>& v3,
                const Nektar::NekMatrix<double>& v4)
{
    Nektar::NekMatrix<double> result = v1;
    NekAddEqual(result, v2);
    NekAddEqual(result, v3);
    NekAddEqual(result, v4);
}

void AddMatricesHandCoded(Nektar::NekMatrix<double>& result,
                const Nektar::NekMatrix<double>& v1,
                const Nektar::NekMatrix<double>& v2)
{
    double* r = result.GetRawPtr();
    const double* b1 = v1.GetRawPtr();
    const double* b2 = v2.GetRawPtr();
    
    unsigned int numElements = result.GetRows()*result.GetColumns();
    for(unsigned int i = 0; i < numElements; ++i)
    {
        r[i] = b1[i] + b2[i];
    }
}

void AddMatricesHandCoded(Nektar::NekMatrix<double>& result,
                const Nektar::NekMatrix<double>& v1,
                const Nektar::NekMatrix<double>& v2,
                const Nektar::NekMatrix<double>& v3)
{
    double* r = result.GetRawPtr();
    const double* b1 = v1.GetRawPtr();
    const double* b2 = v2.GetRawPtr();
    const double* b3 = v3.GetRawPtr();

    unsigned int numElements = result.GetRows()*result.GetColumns();
    for(unsigned int i = 0; i < numElements; ++i)
    {
        r[i] = b1[i] + b2[i] + b3[i];
    }
}

void AddMatricesHandCoded(Nektar::NekMatrix<double>& result,
                const Nektar::NekMatrix<double>& v1,
                const Nektar::NekMatrix<double>& v2,
                const Nektar::NekMatrix<double>& v3,
                const Nektar::NekMatrix<double>& v4)
{
    double* r = result.GetRawPtr();
    const double* b1 = v1.GetRawPtr();
    const double* b2 = v2.GetRawPtr();
    const double* b3 = v3.GetRawPtr();
    const double* b4 = v4.GetRawPtr();

    unsigned int numElements = result.GetRows()*result.GetColumns();
    for(unsigned int i = 0; i < numElements; ++i)
    {
        r[i] = b1[i] + b2[i] + b3[i] + b4[i];
    }
}


void AddMatricesHandCodedResultAlloc(
                const Nektar::NekMatrix<double>& v1,
                const Nektar::NekMatrix<double>& v2)
{
    double* r = new double[v1.GetStorageSize()];
    const double* b1 = v1.GetRawPtr();
    const double* b2 = v2.GetRawPtr();
    
    unsigned int numElements = v1.GetStorageSize();
    for(unsigned int i = 0; i < numElements; ++i)
    {
        r[i] = b1[i] + b2[i];
    }
    delete [] r;
}

void AddMatricesHandCodedResultAlloc(
                const Nektar::NekMatrix<double>& v1,
                const Nektar::NekMatrix<double>& v2,
                const Nektar::NekMatrix<double>& v3)
{
    double* r = new double[v1.GetStorageSize()];
    const double* b1 = v1.GetRawPtr();
    const double* b2 = v2.GetRawPtr();
    const double* b3 = v3.GetRawPtr();

    unsigned int numElements = v1.GetStorageSize();
    for(unsigned int i = 0; i < numElements; ++i)
    {
        r[i] = b1[i] + b2[i] + b3[i];
    }

    delete [] r;
}

void AddMatricesHandCodedResultAlloc(
                const Nektar::NekMatrix<double>& v1,
                const Nektar::NekMatrix<double>& v2,
                const Nektar::NekMatrix<double>& v3,
                const Nektar::NekMatrix<double>& v4)
{
    double* r = new double[v1.GetStorageSize()];
    const double* b1 = v1.GetRawPtr();
    const double* b2 = v2.GetRawPtr();
    const double* b3 = v3.GetRawPtr();
    const double* b4 = v4.GetRawPtr();

    unsigned int numElements = v1.GetStorageSize();
    for(unsigned int i = 0; i < numElements; ++i)
    {
        r[i] = b1[i] + b2[i] + b3[i] + b4[i];
    }

    delete [] r;
}
#endif //NEKTAR_PROFILE_VECTOR_OPS_H

