///////////////////////////////////////////////////////////////////////////////
//
// File: MutexTypeDefs.h
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
// Description: Type definitions of mutex objects.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_BASICUTILS_MUTEXTYPEDEFS_H
#define NEKTAR_LIB_UTILITIES_BASICUTILS_MUTEXTYPEDEFS_H

#include <boost/thread/shared_mutex.hpp>
#include <boost/thread/locks.hpp>

#include <LibUtilities/LibUtilitiesDeclspec.h>

namespace Nektar
{
namespace LibUtilities
{
#ifdef NEKTAR_USING_THREADS
    typedef boost::unique_lock<boost::shared_mutex> WriteLock;
    typedef boost::shared_lock<boost::shared_mutex> ReadLock;
    typedef boost::mutex::scoped_lock               ScopedLock;
    typedef boost::shared_mutex                     SharedMutex;
    typedef boost::mutex                            Mutex;
#else
    class Lock
    {
        public:
            LIB_UTILITIES_EXPORT Lock(){}

            LIB_UTILITIES_EXPORT Lock(int i){}

            LIB_UTILITIES_EXPORT ~Lock(){}

            LIB_UTILITIES_EXPORT inline void lock(){}

            LIB_UTILITIES_EXPORT inline void unlock(){}
    };
    typedef Lock WriteLock;
    typedef Lock ReadLock;
    typedef Lock ScopedLock;
    typedef int SharedMutex;
    typedef int Mutex;
#endif
}
}

#endif //NEKTAR_LIB_UTILITIES_BASICUTILS_MUTEXTYPEDEFS_H
