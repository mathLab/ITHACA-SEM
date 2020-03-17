///////////////////////////////////////////////////////////////////////////////
//
// File: ConsistentObjectAccess.hpp
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
// Description: A wrapper around objects or pointers to objects which give
// the same interface for both types (i.e., -> will work for both).
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_BASIC_UTILS_CONSISTENT_ACCESS_OBJECT_HPP
#define NEKTAR_LIB_UTILITIES_BASIC_UTILS_CONSISTENT_ACCESS_OBJECT_HPP

#include <boost/core/ignore_unused.hpp>

#include <memory>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>

namespace Nektar
{    
    template<typename DataType>
    struct ConsistentObjectAccess
    {     
        static DataType& reference(DataType& o) { return o; }
        static const DataType& const_reference(const DataType& o) { return o; }
        static DataType* pointer(DataType& o) { return &o; }
        static const DataType* const_pointer(const DataType& o) { return &o; }
        
        static bool ReferencesObject(const DataType& o) {
            boost::ignore_unused(o);
            return true;
        }
    };
    
    template<typename DataType>
    struct ConsistentObjectAccess<DataType*> 
    {
        static const DataType& const_reference(DataType* o) { ASSERTL1(o != 0, "Can't dereference null pointer."); return *o; }
        static const DataType* const_pointer(DataType* o) { return o; }
        static bool ReferencesObject(DataType* o) { return o != 0; }

        static DataType& reference(DataType* o) { ASSERTL1(o != 0, "Can't dereference null pointer."); return *o; }
        static DataType* pointer(DataType* o) { return o; }
    };    
    

    template<typename DataType>
    struct ConsistentObjectAccess<std::shared_ptr<DataType> >
    {
        static const DataType& const_reference(const std::shared_ptr<DataType>& o) { ASSERTL1(o, "Can't dereference null pointer."); return *o; }
        static const DataType* const_pointer(const std::shared_ptr<DataType>& o) { return o.get(); }
        static DataType& reference(const std::shared_ptr<DataType>& o) { ASSERTL1(o, "Can't dereference null pointer."); return *o; }
        static DataType* pointer(const std::shared_ptr<DataType>& o) { return o.get(); }
        static bool ReferencesObject(const std::shared_ptr<DataType>& o) { return o.get(); }
    };
}
    
#endif //NEKTAR_LIB_UTILITIES_BASIC_UTILS_CONSISTENT_ACCESS_OBJECT_HPP

