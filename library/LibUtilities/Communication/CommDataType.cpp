///////////////////////////////////////////////////////////////////////////////
//
// File CommDataType.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Define static members for the different data types
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/Communication/CommDataType.h>

#ifdef NEKTAR_USING_PETSC
#include "petscsys.h"
#endif

namespace Nektar
{
namespace LibUtilities
{

/**
 * @brief Return the size in bytes of a data type @p dt.
 *
 * @param dt  Data type
 *
 * @return Size of @p dt in bytes.
 */
int CommDataTypeGetSize(CommDataType dt)
{
#ifdef NEKTAR_USE_MPI
    int size;
    MPI_Type_size(dt, &size);
    return size;
#elif NEKTAR_USING_PETSC
    if (dt == MPI_CHAR)
       return sizeof(char);
    else if (dt == MPI_INT)
       return sizeof(int);
    else if (dt == MPI_UNSIGNED)
        return sizeof(unsigned);
    else if (dt == MPI_LONG)
        return sizeof(long);
    else if (dt == MPI_UNSIGNED_LONG)
        return sizeof(unsigned long);
    else if (dt == MPI_LONG_LONG)
        return sizeof(long long);
    else if (dt == MPI_UNSIGNED_LONG_LONG)
        return sizeof(unsigned long long);
    else if (dt == MPI_FLOAT)
        return sizeof(float);
    else if (dt == MPI_DOUBLE)
        return sizeof(double);
    else if (dt == MPI_LONG_DOUBLE)
        return sizeof(long double);
    return sizeof(int);
#else
    switch (dt)
    {
        case MPI_CHAR:
            return sizeof(char);
        case MPI_INT:
            return sizeof(int);
        case MPI_UNSIGNED:
            return sizeof(unsigned);
        case MPI_LONG:
            return sizeof(long);
        case MPI_UNSIGNED_LONG:
            return sizeof(unsigned long);
        case MPI_LONG_LONG:
            return sizeof(long long);
        case MPI_UNSIGNED_LONG_LONG:
            return sizeof(unsigned long long);
        case MPI_FLOAT:
            return sizeof(float);
        case MPI_DOUBLE:
            return sizeof(double);
        case MPI_LONG_DOUBLE:
            return sizeof(long double);
        default:
            ASSERTL0(false, "Unrecognised datatype!");
    }
    return sizeof(int);
#endif
}

template<> CommDataType &CommDataTypeTraits<char>::GetDataType()
{
    static CommDataType type = MPI_CHAR;
    return type;
}

template<> CommDataType &CommDataTypeTraits<int>::GetDataType()
{
    static CommDataType type = MPI_INT;
    return type;
}

template<> CommDataType &CommDataTypeTraits<unsigned>::GetDataType()
{
    static CommDataType type = MPI_UNSIGNED;
    return type;
}

template<> CommDataType &CommDataTypeTraits<long>::GetDataType()
{
    static CommDataType type = MPI_LONG;
    return type;
}

template<> CommDataType &CommDataTypeTraits<unsigned long>::GetDataType()
{
    static CommDataType type = MPI_UNSIGNED_LONG;
    return type;
}

template<> CommDataType &CommDataTypeTraits<long long>::GetDataType()
{
    static CommDataType type = MPI_LONG_LONG;
    return type;
}

template<> CommDataType &CommDataTypeTraits<unsigned long long>::GetDataType()
{
    static CommDataType type = MPI_UNSIGNED_LONG_LONG;
    return type;
}

template<> CommDataType &CommDataTypeTraits<float>::GetDataType()
{
    static CommDataType type = MPI_FLOAT;
    return type;
}

template<> CommDataType &CommDataTypeTraits<double>::GetDataType()
{
    static CommDataType type = MPI_DOUBLE;
    return type;
}

template<> CommDataType &CommDataTypeTraits<long double>::GetDataType()
{
    static CommDataType type = MPI_LONG_DOUBLE;
    return type;
}

}
}
