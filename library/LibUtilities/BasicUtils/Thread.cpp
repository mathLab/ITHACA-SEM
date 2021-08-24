///////////////////////////////////////////////////////////////////////////////
//
// File Thread.cpp
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
// Description: Thread manager
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <boost/core/ignore_unused.hpp>

#include "LibUtilities/BasicUtils/Thread.h"

namespace Nektar
{
namespace Thread
{


/**
 *
 */
ThreadManagerFactory& GetThreadManagerFactory()
{
    static ThreadManagerFactory instance;
    return instance;
}

/**
 * @brief ThreadJob implementation
 */
ThreadJob::ThreadJob()
{
    // empty
}


/**
 *
 */
ThreadJob::~ThreadJob()
{
    // empty
}


/**
 * Part of the thread interface.  Do not use unless you're
 * an implementation of ThreadManager.
 * @warning Do not use: needs to be public so thread implementation can call it.
 */
void ThreadJob::SetWorkerNum(unsigned int num)
{
    m_workerNum = num;
}


/**
 *
 */
unsigned int ThreadJob::GetWorkerNum()
{
    return m_workerNum;
}


/**
 * Implementations should not override this function, since they
 * are initialised.
 *
 * @returns True if this ThreadManager has been initialised.
 */
bool ThreadManager::IsInitialised()
{
    return true;
}


/**
 * Should stop worker threads and clean up.
 * This shouldn't be called until the program is exiting
 * anyway, some implementations may need to unlock threads
 * to allow clean exit.
 */
ThreadManager::~ThreadManager()
{
    // empty
}


/**
 * ThreadMaster implementation
 */
ThreadMaster::ThreadMaster() : m_threadManagers(THREADMANAGER_MAX),
    m_mutex(), m_threadingType()
{
    // empty
}


/**
 * Will clear the list of ThreadManagers, destructing them.
 */
ThreadMaster::~ThreadMaster()
{
    // Locking is a bit pointless, since the map is empty after this call.
    m_threadManagers.clear();
}


/**
 *
 */
ThreadMaster& GetThreadMaster()
{
    static ThreadMaster instance;
    return instance;
}


/**
 * @param p_type String to be passed to the ThreadManagerFactory
 *
 * Subsequent CreateInstance calls will pass this string to the
 * ThreadManagerFactory to create a ThreadManager.
 *
 * It is an error to call this more than once (since having different kinds of
 * ThreadManager active is probably a stupendously bad idea).
 */
void ThreadMaster::SetThreadingType(const std::string &p_type)
{
    ASSERTL0(m_threadingType.empty(),
             "Tried to SetThreadingType when it was already set");
    m_threadingType = p_type;
}


/**
 * @returns a shared pointer to /em either a ThreadStartupManager
 * or whatever ThreadManager has been created for the string s
 * with CreateInstance.
 *
 * Calling code may store the result if it is sure the call to
 * GetInstance(s) has occurred after the call to CreateInstance(s).
 * This cannot be before threadedcommReader::StartThreads(), as that's
 * where SetThreadingType is called.
 *
 * @warning Static initialisation may want to access a ThreadManager.
 * Such code must be able to cope with the temporary ThreadStartupManager.
 */
ThreadManagerSharedPtr& ThreadMaster::GetInstance(const ThreadManagerName t)
{
    if ( !m_threadManagers[t] )
    {
        m_threadManagers[t] = ThreadManagerSharedPtr(
                                            new ThreadStartupManager());
        return m_threadManagers[t];
    }
    return m_threadManagers[t];
}


/**
 * @return Shared pointer to the created ThreadManager.
 *
 * An error occurs if this is called before SetThreadingType.
 */
ThreadManagerSharedPtr ThreadMaster::CreateInstance(const ThreadManagerName t,
    unsigned int nThr)
{
    ASSERTL0(!m_threadingType.empty(),
             "Trying to create a ThreadManager before SetThreadingType called");
    return m_threadManagers[t] =
        Thread::GetThreadManagerFactory().CreateInstance(m_threadingType, nThr);
}


/**
 * @brief ThreadDefaultManager
 */
ThreadStartupManager::ThreadStartupManager() : m_type("Threading starting up")
{
    // empty
}


/**
 *
 */
ThreadStartupManager::~ThreadStartupManager()
{
    // empty
}


/**
 *
 */
void ThreadStartupManager::QueueJobs(std::vector<ThreadJob*>& joblist)
{
    boost::ignore_unused(joblist);
    NEKERROR(ErrorUtil::efatal,
             "Attempted to QueueJobs in ThreadDefaultManager");
}


/**
 *
 */
void ThreadStartupManager::QueueJob(ThreadJob* job)
{
    boost::ignore_unused(job);
    NEKERROR(ErrorUtil::efatal,
             "Attempted to QueueJob in ThreadDefaultManager");
}


/**
 *
 */
unsigned int ThreadStartupManager::GetNumWorkers()
{
    return 1;
}


/**
 *
 */
unsigned int ThreadStartupManager::GetWorkerNum()
{
    return 0;
}


/**
 *
 */
void ThreadStartupManager::SetNumWorkers(const unsigned int num)
{
    ASSERTL0(num==1,
             "Attempted to SetNumWorkers to != 1 in ThreadDefaultManager");
}


/**
 *
 */
void ThreadStartupManager::SetNumWorkers()
{
    return;
}


/**
 *
 */
unsigned int ThreadStartupManager::GetMaxNumWorkers()
{
    return 1;
}


/**
 *
 */
void ThreadStartupManager::Wait()
{
    return;
}


/**
 *
 */
void ThreadStartupManager::SetChunkSize(unsigned int chnk)
{
    boost::ignore_unused(chnk);
    NEKERROR(ErrorUtil::efatal,
             "Attempted to SetChunkSize in ThreadDefaultManager");
}


/**
 *
 */
void ThreadStartupManager::SetSchedType(SchedType s)
{
    boost::ignore_unused(s);
    NEKERROR(ErrorUtil::efatal,
             "Attempted to SetSchedType in ThreadDefaultManager");
}


/**
 *
 */
bool ThreadStartupManager::InThread()
{
    return false;
}


/**
 *
 */
void ThreadStartupManager::Hold()
{
    return;
}


/**
 *
 */
bool ThreadStartupManager::IsInitialised()
{
    return false;
}


/**
 *
 */
const std::string& ThreadStartupManager::GetType() const
{
    return m_type;
}


/**
 * @brief ThreadDefaultManager copy constructor
 */
ThreadStartupManager& ThreadStartupManager::operator=(
        const ThreadStartupManager& src)
{
    boost::ignore_unused(src);
    return *this;
}

} // Thread
} /* namespace Nektar */
