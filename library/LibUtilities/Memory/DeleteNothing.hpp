///////////////////////////////////////////////////////////////////////////////
//
// File: DeleteNothing.hpp
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
// Description: Utility class to allow points to exist in a boost
//              shared array without being automatically deleted.
//
// 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_MEMORY_DELETE_NOTHING_HPP
#define NEKTAR_LIB_UTILITIES_MEMORY_DELETE_NOTHING_HPP

#include <LibUtilities/LinearAlgebra/MatrixBlockType.h>

namespace Nektar
{
    template<typename DataType>
    class DeleteNothing
    {
        public:
            void operator()(DataType* p) {}
            void operator()(const DataType* p) const {}
    };
            
};
    
#endif //NEKTAR_LIB_UTILITIES_MEMORY_DELETE_NOTHING_HPP

/**
    $Log: DeleteNothing.hpp,v $
    Revision 1.2  2007/02/04 04:27:49  bnelson
    Updated linear systems to work with normal objects instead of only shared pointers.

    Revision 1.1  2006/10/30 05:10:13  bnelson
    Object to use when creating a shared array and you don't want to delete the pointer.

 **/

