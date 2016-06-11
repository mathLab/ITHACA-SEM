///////////////////////////////////////////////////////////////////////////////
//
// File CommDataType.h
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
// Description: Describes data types (using MPI_Datatype if available)
//
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_LIB_UTILITIES_COMMDATATYPE_H
#define NEKTAR_LIB_UTILITIES_COMMDATATYPE_H

#include <vector>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>

#ifdef NEKTAR_USE_MPI
#include <mpi.h>

namespace Nektar
{
    namespace LibUtilities
    {
        typedef MPI_Datatype CommDataType;
    }
}
#else

namespace Nektar
{
    namespace LibUtilities
    {
        enum CommDataType
        {
            MPI_INT, MPI_UNSIGNED,
            MPI_LONG, MPI_UNSIGNED_LONG,
            MPI_LONG_LONG, MPI_UNSIGNED_LONG_LONG,
            MPI_FLOAT, MPI_DOUBLE, MPI_LONG_DOUBLE
        };
    }
}
#endif

namespace Nektar
{
    template<typename Dim, typename DataType> class Array;

    namespace LibUtilities
    {
        int CommDataTypeGetSize(CommDataType);

        template<class T>
        class CommDataTypeTraits
        {
                static CommDataType type;
            public:
                static CommDataType& GetDataType()
                {
                    return type;
                }
                static void* GetPointer(T& val)
                {
                    return &val;
                }
                static const void* GetPointer(const T& val)
                {
                    return &val;
                }
                static int GetCount(const T& val)
                {
                    return 1;
                }

                const static bool IsVector = false;
        };

        /**
         * Partial specialisation for vectors
         */
        template<class elemT>
        class CommDataTypeTraits<std::vector<elemT> >
        {
            public:
                static CommDataType& GetDataType()
                {
                    return CommDataTypeTraits<elemT>::GetDataType();
                }
                static void* GetPointer(std::vector<elemT>& val)
                {
                    return &val[0];
                }
                static const void* GetPointer(const std::vector<elemT>& val)
                {
                    return &val[0];
                }
                static int GetCount(const std::vector<elemT>& val)
                {
                    return val.size();
                }
                const static bool IsVector = true;
        };

        /**
         * Partial specialisation for vectors
         */
        template<class elemT>
        class CommDataTypeTraits<Array<OneD, elemT> >
        {
            public:
                static CommDataType& GetDataType()
                {
                    return CommDataTypeTraits<elemT>::GetDataType();
                }
                static void* GetPointer(Array<OneD, elemT>& val)
                {
                    return val.get();
                }
                static const void* GetPointer(const Array<OneD, elemT>& val)
                {
                    return val.get();
                }
                static int GetCount(const Array<OneD, elemT>& val)
                {
                    return val.num_elements();
                }
                const static bool IsVector = true;
        };
    }
}

#endif
