/*
 * Thread.cpp
 *
 */
#include <iostream>

#include "Thread.h"
#include <loki/Singleton.h>

namespace Nektar
{
    namespace Thread
    {
        ThreadManagerSharedPtr ThreadManager::m_instance;
        ThreadManagerFactory& GetThreadManager()
        {
            typedef Loki::SingletonHolder<ThreadManagerFactory,
                Loki::CreateUsingNew,
                Loki::NoDestroy > Type;
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
        ThreadManager::~ThreadManager()
        {
        	// empty
        }
        
        // ThreadHandle implementation.
        ThreadHandle::ThreadHandle(SchedType sched)
        {
        	Setup(sched, 1);
        }

        ThreadHandle::ThreadHandle(SchedType sched, unsigned int chnk)
        {
        	Setup(sched, chnk);
        }

        void ThreadHandle::Setup(SchedType sched, unsigned int chnk)
        {
        	m_tm = ThreadManager::GetInstance();
        	if (!m_tm)
        	{
        		std::cerr << "Attempted to construct a ThreadHandle before a ThreadManager has been created." << std::endl;
        		std::abort();
        	}
        	m_tm->SetSchedType(sched);
        	m_tm->SetChunkSize(chnk);
        	m_tm->SetNumWorkers();
        }

        ThreadHandle::~ThreadHandle()
        {
        	m_tm->Wait();
        }

    } // Thread
} /* namespace Nektar */
