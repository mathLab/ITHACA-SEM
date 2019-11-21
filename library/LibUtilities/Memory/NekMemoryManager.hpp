////////////////////////////////////////////////////////////////////////////////
//
//  File: NekMemoryManager.hpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description:
//
//  A memory manager that allocates memory from thread specific pools.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_NEK_MEMORY_MANAGER_H
#define NEKTAR_LIB_UTILITIES_NEK_MEMORY_MANAGER_H

#include <memory>
#include <type_traits>

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/Memory/ThreadSpecificPool.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>

#include <vector>

#ifdef max
#undef max
#endif

namespace Nektar
{

/// @brief General purpose memory allocation routines with the ability
///        to allocate from thread specific memory pools.
///
/// If compiled with NEKTAR_MEMORY_POOL_ENABLED, the MemoryManager
/// allocates from thread specific memory pools for small objects.
/// Large objects are managed with the system supplied new/delete.
/// These memory pools provide faster allocation and deallocation
/// of small objects (particularly useful for shared pointers which
/// allocate many 4 byte objects).
///
/// @warning All memory allocated from the memory manager must be returned
/// to the memory manager.  Calling delete on memory allocated from the
/// manager will likely cause undefined behavior.  A particularly subtle
/// violation of this rule occurs when giving memory allocated from the
/// manager to a shared pointer.
/// @code
/// std::shared_ptr<Obj> f(MemoryManager<Obj>::Allocate());
/// @endcode
/// Shared pointers call delete when they go out of scope, so this line of
/// code will cause problems.  Instead, you should call the
/// AllocateSharedPtr method:
/// @code
/// std::shared_ptr<Obj> f = MemoryManager<Obj>::AllocateSharedPtr();
/// @endcode
template<typename DataType>
class MemoryManager
{
public:
    /// @brief Deallocate a pointer allocated by
    /// MemoryManager::Allocate.
    /// @note Results are undefined if called with a pointer to
    /// something that was not allocated with the memory manager.
    ///
    /// Use this method to deallocate a pointer you have allocated from
    /// the MemoryManager using the Allocate method.
    ///
    /// Example:
    /// @code
    /// CustObj* c = MemoryManager::Allocate<CustObj>();
    /// MemoryManager::Deallocate(c);
    /// @endcode
    static void Deallocate(DataType*& data)
    {
#ifdef NEKTAR_MEMORY_POOL_ENABLED
        data->~DataType();
        GetMemoryPool().Deallocate(data, sizeof(DataType));
#else
        delete data;
#endif

        data = NULL;
    }

#ifdef NEKTAR_MEMORY_POOL_ENABLED
    /// @brief Allocates a single object from the memory pool.
    /// @throws unknown If the object throws an exception during
    /// construction, this method will catch it, release the memory
    /// back to the pool, then rethrow it.
    ///
    /// The allocated object must be returned to the memory pool
    /// via Deallocate.
    template<typename... Args>
    static DataType* Allocate(const Args &...args)
    {
        DataType* result = static_cast<DataType*>(
            GetMemoryPool().Allocate(sizeof(DataType)));

        if (result)
        {
            try
            {
                new (result) DataType(args...);
            }
            catch(...)
            {
                GetMemoryPool().Deallocate(result, sizeof(DataType));
                throw;
            }
        }

        return result;
    }

#else //NEKTAR_MEMORY_POOL_ENABLED
    /// @brief Allocates a single object from the memory pool.
    /// @throws unknown Any exception thrown by DataType's default
    /// constructor will propogate through this method.
    ///
    /// The allocated object must be returned to the memory pool
    /// via Deallocate.
    template<typename... Args>
    static DataType* Allocate(const Args &...args)
    {
        return new DataType(args...);
    }
#endif //NEKTAR_MEMORY_POOL_ENABLED

    /// @brief Allocate a shared pointer from the memory pool.
    ///
    /// The shared pointer does not need to be returned to the memory
    /// pool. When the reference count to this object reaches 0, the
    /// shared pointer will automatically return the memory.
    template<typename... Args>
    static std::shared_ptr<DataType> AllocateSharedPtr(const Args &...args)
    {
        return AllocateSharedPtrD( [](DataType *){}, args...);
    }

    template<typename DeallocatorType, typename... Args>
    static std::shared_ptr<DataType> AllocateSharedPtrD(
        const DeallocatorType& d, const Args &...args)
    {
        DataType* data = Allocate(args...);
        return std::shared_ptr<DataType>(
            data, [=](DataType *ptr){
                d(ptr);
                MemoryManager<DataType>::Deallocate(ptr);
            });
    }

    /// \brief Allocates a chunk of raw, uninitialized memory, capable of
    /// holding NumberOfElements objects.
    ///
    /// \param NumberOfElements The number of elements the array should be
    ///                         capable of holding.
    ///
    /// This method is not meant to be called by client code.  Use Array
    /// instead.  Any memory allocated from this method must be returned to the
    /// memory pool via RawDeallocate.  Failure to do so will result in memory
    /// leaks and undefined behavior.
    static DataType* RawAllocate(size_t NumberOfElements)
    {
#ifdef NEKTAR_MEMORY_POOL_ENABLED
        return static_cast<DataType*>(GetMemoryPool().Allocate(sizeof(DataType)*NumberOfElements));
#else //NEKTAR_MEMORY_POOL_ENABLED
        return static_cast<DataType*>(::operator new(NumberOfElements * sizeof(DataType)));
#endif //NEKTAR_MEMORY_POOL_ENABLED
    }


    /// \brief Deallocates memory allocated from RawAllocate.
    /// \param array A pointer to the memory returned from RawAllocate.
    /// \param NumberOfElements The number of object held in the array.
    ///
    /// This method is not meant to be called by client code.  Use Array instead.
    /// Only memory allocated via RawAllocate should be returned to the pool here.
    static void RawDeallocate(DataType* array, size_t NumberOfElements)
    {
#ifdef NEKTAR_MEMORY_POOL_ENABLED
        GetMemoryPool().Deallocate(array, sizeof(DataType)*NumberOfElements);
#else //NEKTAR_MEMORY_POOL_ENABLED
        ::operator delete(array);
#endif //NEKTAR_MEMORY_POOL_ENABLED
    }

    /////////////////////////////////////////////////////////////////
    ///\name Allocator Interface
    /// The allocator interface allows a MemoryManager object to be used
    /// in any object that allows an allocator parameter, such as STL
    /// containers.
    /////////////////////////////////////////////////////////////////
    typedef DataType value_type;
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;
    typedef DataType* pointer;
    typedef const DataType* const_pointer;
    typedef DataType& reference;
    typedef const DataType& const_reference;

    MemoryManager() {}
    template<typename T>
    MemoryManager(const MemoryManager<T>& rhs)
    {
        boost::ignore_unused(rhs);
    }
    ~MemoryManager() {}

    pointer address(reference r) const { return &r; }
    const_pointer address(const_reference r) const { return &r; }

    pointer allocate(size_type n, std::allocator<void>::const_pointer hint = 0)//typename MemoryManager<void>::pointer hint = 0)
    {
        boost::ignore_unused(hint);
        return RawAllocate(n);
    }

    void deallocate(pointer p, size_type n)
    {
        return RawDeallocate(p, n);
    }

    void construct(pointer p, const_reference val)
    {
        new(p) DataType(val);
    }

    void destroy(pointer p)
    {
        p->~DataType();
    }

    size_type max_size()
    {
        return std::numeric_limits<size_type>::max()/sizeof(DataType);
    }

    template<typename U>
    struct rebind
    {
        typedef MemoryManager<U> other;
    };

    /////////////////////////////////////////////////////////////////
    ///@}
    /////////////////////////////////////////////////////////////////

private:

};

template<typename DataType>
bool operator==(const MemoryManager<DataType>& lhs, const MemoryManager<DataType>& rhs)
{
    boost::ignore_unused(lhs,rhs);
    return true;
}

template<typename DataType>
bool operator!=(const MemoryManager<DataType>& lhs, const MemoryManager<DataType>& rhs)
{
    return !(lhs == rhs);
}

}

#endif //NEKTAR_LIB_UTILITIES_NEK_MEMORY_MANAGER_H


