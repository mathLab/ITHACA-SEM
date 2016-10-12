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
// License for the specific language governing rights and limitations under
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
#else
    switch (dt)
    {
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

/// Type trait mapping an int to MPI_INT
template <> CommDataType CommDataTypeTraits<int>::type      = MPI_INT;
/// Type trait mapping an unsigned int to MPI_UNSIGNED
template <> CommDataType CommDataTypeTraits<unsigned>::type = MPI_UNSIGNED;
/// Type trait mapping a long int to MPI_LONG
template <> CommDataType CommDataTypeTraits<long>::type     = MPI_LONG;
/// Type trait mapping an unsigned long int to MPI_UNSIGNED_LONG
template <>
CommDataType CommDataTypeTraits<unsigned long>::type = MPI_UNSIGNED_LONG;
/// Type trait mapping a long long int to MPI_LONG_LONG
template <> CommDataType CommDataTypeTraits<long long>::type = MPI_LONG_LONG;
/// Type trait mapping an unsigned long long int to MPI_UNSIGNED_LONG_LONG
template <>
CommDataType CommDataTypeTraits<unsigned long long>::type =
    MPI_UNSIGNED_LONG_LONG;
/// Type trait mapping a float to MPI_FLOAT
template <> CommDataType CommDataTypeTraits<float>::type  = MPI_FLOAT;
/// Type trait mapping a double to MPI_DOUBLE
template <> CommDataType CommDataTypeTraits<double>::type = MPI_DOUBLE;
/// Type trait mapping a long double to MPI_LONG_DOUBLE
template <>
CommDataType CommDataTypeTraits<long double>::type = MPI_LONG_DOUBLE;
}
}
