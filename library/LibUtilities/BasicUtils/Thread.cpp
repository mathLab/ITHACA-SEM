/*
 * Thread.cpp
 *
 */
#include <iostream>

#include "LibUtilities/BasicUtils/Thread.h"

namespace Nektar
{
    namespace Thread
    {
//        ThreadManagerSharedPtr ThreadManager::m_instance;
        ThreadManagerFactory& GetThreadManagerFactory()
        {
            typedef Loki::SingletonHolder<ThreadManagerFactory,
                Loki::CreateUsingNew,
                Loki::NoDestroy,
                Loki::SingleThreaded> Type;
            return Type::Instance();
        }

        // ThreadJob implementation
        ThreadJob::ThreadJob()
        {
            // empty
        }
        
        ThreadJob::~ThreadJob()
        {
            // empty
        }
        
        void ThreadJob::SetWorkerNum(unsigned int num)
        {
            m_workerNum = num;
        }
        
        unsigned int ThreadJob::GetWorkerNum()
        {
            return m_workerNum;
        }

        // ThreadManager implementation.
        bool ThreadManager::IsInitialised()
        {
            return true;
        }

        ThreadManager::~ThreadManager()
        {
            // empty
        }
        
        // ThreadMaster implementation
        ThreadMaster::ThreadMaster() : m_threadManagers(1), m_mutex(), m_threadingType()
        {
            // empty
        }

        ThreadMaster::~ThreadMaster()
        {
            // Locking is a bit pointless, since the map is empty after this call.
            m_threadManagers.clear();
        }

        ThreadMaster& GetThreadMaster()
        {
            typedef Loki::SingletonHolder<ThreadMaster,
                    Loki::CreateUsingNew,
                    Loki::NoDestroy,
                    Loki::SingleThreaded> Type;
            return Type::Instance();
        }

        void ThreadMaster::SetThreadingType(const std::string &p_type)
        {
            ASSERTL0(m_threadingType.empty(), "Tried to SetThreadingType when it was already set");
            m_threadingType = p_type;
        }

        ThreadManagerSharedPtr& ThreadMaster::GetInstance(const ThreadManagerName t)
        {
            if ( !m_threadManagers[t] )
            {
                m_threadManagers[t] = ThreadManagerSharedPtr(new ThreadStartupManager());
                return m_threadManagers[t];
            }
            return m_threadManagers[t];
        }

        ThreadManagerSharedPtr ThreadMaster::CreateInstance(const ThreadManagerName t,
            unsigned int nThr)
        {
            ASSERTL0(!m_threadingType.empty(), "Trying to create a ThreadManager before SetThreadingType called");
            return m_threadManagers[t] =
                Thread::GetThreadManagerFactory().CreateInstance(m_threadingType, nThr);
        }

        // ThreadDefaultManager
        ThreadStartupManager::ThreadStartupManager() : m_type("Threading starting up")
        {
            // empty
        }
        ThreadStartupManager::~ThreadStartupManager()
        {
            // empty
        }
        void ThreadStartupManager::QueueJobs(std::vector<ThreadJob*>& joblist)
        {
            NEKERROR(ErrorUtil::efatal, "Attempted to QueueJobs in ThreadDefaultManager");
        }
        void ThreadStartupManager::QueueJob(ThreadJob* job)
        {
            NEKERROR(ErrorUtil::efatal, "Attempted to QueueJob in ThreadDefaultManager");
        }
        unsigned int ThreadStartupManager::GetNumWorkers()
        {
            return 1;
        }
        unsigned int ThreadStartupManager::GetWorkerNum()
        {
            return 0;
        }
        void ThreadStartupManager::SetNumWorkers(const unsigned int num)
        {
            ASSERTL0(num==1, "Attempted to SetNumWorkers to != 1 in ThreadDefaultManager");
        }
        void ThreadStartupManager::SetNumWorkers()
        {
            return;
        }
        unsigned int ThreadStartupManager::GetMaxNumWorkers()
        {
            return 1;
        }
        void ThreadStartupManager::Wait()
        {
            return;
        }
        void ThreadStartupManager::SetChunkSize(unsigned int chnk)
        {
            NEKERROR(ErrorUtil::efatal, "Attempted to SetChunkSize in ThreadDefaultManager");
        }
        void ThreadStartupManager::SetSchedType(SchedType s)
        {
            NEKERROR(ErrorUtil::efatal, "Attempted to SetSchedType in ThreadDefaultManager");
        }
        bool ThreadStartupManager::InThread()
        {
            return false;
        }
        void ThreadStartupManager::Hold()
        {
            return;
        }
        bool ThreadStartupManager::IsInitialised()
        {
            return false;
        }
        const std::string& ThreadStartupManager::GetType() const
        {
            return m_type;
        }

    } // Thread
} /* namespace Nektar */
