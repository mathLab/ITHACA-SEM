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
        ThreadManagerSharedPtr ThreadManager::instance;
        ThreadManagerFactory& GetThreadManager()
        {
            typedef Loki::SingletonHolder<ThreadManagerFactory,
                Loki::CreateUsingNew,
                Loki::NoDestroy > Type;
            return Type::Instance();
        }

        
        ThreadJob::ThreadJob()
        {
            // empty
        }
        
        ThreadJob::~ThreadJob()
        {
            // empty
        }
        
        unsigned int ThreadJob::getWorkerNum()
        {
            return m_workerNum;
        }
        
        void ThreadJob::setWorkerNum(unsigned int num)
        {
            m_workerNum = num;
        }
        
        ThreadManager::~ThreadManager()
        {
        	// empty
        }
        
        ThreadHandle::ThreadHandle(SchedType sched)
        {
        	setup(sched, 1);
        }

        ThreadHandle::ThreadHandle(SchedType sched, unsigned int chnk)
        {
        	setup(sched, chnk);
        }

        void ThreadHandle::setup(SchedType sched, unsigned int chnk)
        {
        	//FIXME: check if instance is initialised.  If not, throw.
        	tm = ThreadManager::getInstance();
        	tm->setSchedType(sched);
        	tm->setChunkSize(chnk);
        	tm->setNumWorkers();
        }

        ThreadHandle::~ThreadHandle()
        {
        	tm->wait();
        }

    } // Thread
} /* namespace Nektar */
