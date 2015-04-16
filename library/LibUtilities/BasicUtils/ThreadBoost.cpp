/*
 * Thread.cpp
 *
 */

#include "LibUtilities/BasicUtils/ThreadBoost.h"
#include <iostream>

namespace Nektar
{
    namespace Thread
    {

        std::string ThreadManagerBoost::className =
            GetThreadManagerFactory().RegisterCreatorFunction("ThreadManagerBoost",
                    ThreadManagerBoost::Create, "Threading using Boost.");

        ThreadManagerBoost::~ThreadManagerBoost()
        {
            // This is an immediate teardown.  We attempt to kill everything.
            // we daren't lock anything as we may cause a deadlock
            for (unsigned int i=0; i<m_numThreads; i++)
            {
                m_threadList[i]->Stop();
            }

            m_masterQueueCondVar.notify_all();
            m_masterActiveCondVar.notify_all();
            for (unsigned int i=0; i<m_numThreads; i++)
            {
                m_threadThreadList[i]->join();
                delete m_threadThreadList[i];
                delete m_threadList[i];
            }
        
            delete []m_threadList;
            delete []m_threadThreadList;
            delete []m_threadActiveList;
            delete []m_threadBusyList;
            delete m_barrier;
        }
        
        ThreadManagerBoost::ThreadManagerBoost(unsigned int numT) :
                m_numThreads(numT), m_numWorkers(numT-1), m_masterQueue(), m_masterQueueMutex(),
                m_masterActiveMutex(), m_masterQueueCondVar(), m_masterActiveCondVar(),
                m_chunkSize(1), m_schedType(e_dynamic),
                m_threadMap()
        {
            using namespace std;
            try {
                m_threadList = new ThreadWorkerBoost*[m_numThreads];
                m_threadThreadList = new boost::thread*[m_numThreads];
                m_threadBusyList = new bool[m_numThreads];
                m_threadActiveList = new bool[m_numThreads];
            } catch (exception &e) {
                cerr << "Exception while allocating thread storage: "
                        << e.what() << endl;
                abort();
            }
            unsigned int i = 0;
            while (i < m_numThreads)
            {
                ThreadWorkerBoost *tw;
                try {
                    tw = new ThreadWorkerBoost(this, i);
                } catch (exception &e) {
                    cerr << "Exception while allocating worker threads: "
                            << e.what() << endl;
                    abort();
                }
        
                m_threadList[i] = tw;
                m_threadBusyList[i] = false;
                m_threadActiveList[i] = true;
        
                try {
                    m_threadThreadList[i] = new boost::thread(boost::ref(*tw));
                    boost::thread::id id = m_threadThreadList[i]->get_id();
                    m_threadMap[id] = i;
                } catch (...) {
                    std::cerr << "Exception while creating worker threads" << std::endl;
                    abort();
                }
                i++;
            }
            m_threadActiveList[m_numThreads-1] = false;
            m_masterThreadId = boost::this_thread::get_id();
            m_barrier = new boost::barrier(m_numWorkers > 0 ? m_numWorkers : 1);
            m_type = "Threading with Boost";
        }
        
        void ThreadManagerBoost::QueueJobs(std::vector<ThreadJob*> &joblist)
        {
            std::vector<ThreadJob *>::iterator it;
            for (it=joblist.begin(); it<joblist.end(); ++it)
            {
                QueueJob(*it);
            }
        }
        
        void ThreadManagerBoost::QueueJob(ThreadJob *job)
        {
            Lock masterQueueLock(m_masterQueueMutex); // locks the queue
            m_masterQueue.push(job);
            m_masterQueueCondVar.notify_all(); // alert a waiting thread.
        }   // queue unlocked
        
        bool ThreadManagerBoost::IsWorking()
        {
            bool working = false;
            Lock masterActiveLock(m_masterActiveMutex);
            for (unsigned int i = 0; i < m_numWorkers; i++)
            {
                working = working || m_threadBusyList[i];
            }
            return working;
        }
        
        void ThreadManagerBoost::SetChunkSize(unsigned int chnk)
        {
            Lock masterQueueLock(m_masterQueueMutex); // locks the queue
            m_chunkSize = std::max(chnk, 1U);
        }

        void ThreadManagerBoost::SetSchedType(SchedType s)
        {
            Lock masterQueueLock(m_masterQueueMutex); // locks the queue
            m_schedType = s;
        }

        bool ThreadManagerBoost::InThread()
        {
            boost::thread::id id = boost::this_thread::get_id();
            return (id != m_masterThreadId);
        }

        void ThreadManagerBoost::Wait()
        {
            bool working;
            Lock masterQueueLock(m_masterQueueMutex); // locks the queue
            unsigned int nw = m_numWorkers;
            SetNumWorkersImpl(nw+1);
            working = IsWorking();
            while (!m_masterQueue.empty() || working)
            {
                // while waiting, master queue is unlocked
                m_masterQueueCondVar.wait(masterQueueLock);
                // on exiting wait master queue is locked again
                working = IsWorking();
            }
            SetNumWorkersImpl(nw-1);
        }
        
        unsigned int ThreadManagerBoost::GetNumWorkers()
        {
            return m_numWorkers;
        }
        
        unsigned int ThreadManagerBoost::GetWorkerNum()
        {
            boost::thread::id id = boost::this_thread::get_id();
            return m_threadMap[id];
        }

        void ThreadManagerBoost::SetNumWorkersImpl(const unsigned int num)
        {
            Lock masterActiveLock(m_masterActiveMutex); // locks the active

            if (m_numWorkers == num)
            {
                return;
            }

            delete m_barrier;
            m_barrier = new boost::barrier(num > 0 ? num : 1);

            m_numWorkers = num;
            for (unsigned int i = 0; i < m_numThreads; i++)
            {
                m_threadActiveList[i] = i < m_numWorkers ? true : false;
            }
            m_masterActiveCondVar.notify_all();
        } // Lock on active released here

        void ThreadManagerBoost::SetNumWorkers(unsigned int num)
        {
            num = std::min(num, m_numThreads);
            num = std::max(num, static_cast<unsigned int>(0));
            --num;
            SetNumWorkersImpl(num);
        }
        
        void ThreadManagerBoost::SetNumWorkers()
        {
            SetNumWorkersImpl(m_numThreads-1);
        }
        
        unsigned int ThreadManagerBoost::GetMaxNumWorkers()
        {
            return m_numThreads;
        }
        
        void ThreadManagerBoost::Hold()
        {
                m_barrier->wait();
        }

        const std::string& ThreadManagerBoost::GetType() const
        {
            return m_type;
        }


        ThreadWorkerBoost::ThreadWorkerBoost(ThreadManagerBoost *tm, unsigned int workerNum) :
                m_threadManager(tm), m_workerQueue(),
                m_keepgoing(true), m_threadNum(workerNum)
        {
            // Nothing to see here
        
        }
        
        ThreadWorkerBoost::~ThreadWorkerBoost()
        {
            if (m_keepgoing)
            {
                std::cerr << "Warning: ThreadWorker: " << m_threadNum
                        << "destroyed while running!" << std::endl;
            }
            // on destuction the m_workerQueue will be destructed and that
            // will destruct any ThreadJobs still in there.
        }
        
        void ThreadWorkerBoost::LoadJobs()
        {
            // Lock the master queue
            Lock masterQueueLock(m_threadManager->m_masterQueueMutex);
            m_threadManager->m_threadBusyList[m_threadNum] = false;
            m_threadManager->m_masterQueueCondVar.notify_all();
            while (m_threadManager->m_masterQueue.empty()
                && m_keepgoing)
            {
                // while waiting, master queue is unlocked
                m_threadManager->m_masterQueueCondVar.wait(masterQueueLock);
                // on exiting wait master queue is locked again
            }
            bool active;
            {
                Lock masterActiveLock(m_threadManager->m_masterActiveMutex);
                active = m_threadManager->m_threadActiveList[m_threadNum];
            }
            if (active && m_keepgoing)
            {
                unsigned int numToLoad = GetNumToLoad();
                while (m_workerQueue.size() < numToLoad
                        && !m_threadManager->m_masterQueue.empty())
                {
                    ThreadJob *tj = m_threadManager->m_masterQueue.front();
                    m_workerQueue.push(tj);
                    m_threadManager->m_masterQueue.pop();
                }
            }
            //std::cerr << "Loaded thr " << m_threadNum << " with " << m_workerQueue.size() << " jobs" << std::endl;
            m_threadManager->m_threadBusyList[m_threadNum] = !m_workerQueue.empty();
        } // lock on master queue released here


        unsigned int ThreadWorkerBoost::GetNumToLoad()
        {
            unsigned int numToLoad;
            switch (m_threadManager->m_schedType)
            {
                case e_guided:
                    numToLoad = std::max(static_cast<unsigned long>(m_threadManager->m_chunkSize),
                            m_threadManager->m_masterQueue.size() / (2*m_threadManager->m_numWorkers +1));
                    break;

                case e_dynamic:
                    numToLoad = m_threadManager->m_chunkSize;
                    break;
            }
            return numToLoad;
        }
        
        void ThreadWorkerBoost::WaitForActive()
        {
            Lock masterActiveLock(m_threadManager->m_masterActiveMutex);
        
            while (!m_threadManager->m_threadActiveList[m_threadNum]
                && m_keepgoing)
            {
                // while waiting, master active is unlocked
                m_threadManager->m_masterActiveCondVar.wait(masterActiveLock);
                // on exiting wait master active is locked again
            }
        }
        
        void ThreadWorkerBoost::MainLoop()
        {
            while (m_keepgoing)
            {
                WaitForActive();
                LoadJobs();
                RunJobs();
            }
        } // exiting here should terminate the thread
        
        void ThreadWorkerBoost::RunJobs()
        {
            while (!m_workerQueue.empty() && m_keepgoing)
            {
                ThreadJob * tj;
                try
                {
                    tj = m_workerQueue.front();
                    tj->SetWorkerNum(m_threadNum);
                    tj->Run();
                    m_workerQueue.pop();
                    delete tj;
                } catch(...)
                {
                    // something bad happened, probably time to die
                    // maybe signal ThreadManager
                    throw;
                }
            }
        }
        
    } // Thread
} /* namespace Nektar */
