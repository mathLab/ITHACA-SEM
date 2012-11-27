/*
 * Thread.cpp
 *
 */

#include "Thread.h"
#include <iostream>

namespace Nektar
{
    namespace Thread
    {
        ThreadManager *ThreadManager::instance = 0;
        
        ThreadJob::ThreadJob() :
                m_workerNum(-1)
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
            // This is an immediate teardown.  We attempt to kill everything.
            // we daren't lock anything as we may cause a deadlock
            std::cerr << "exiting threadmanager" << std::endl;
            for (unsigned int i=0; i<m_numThreads; i++)
            {
                m_threadList[i]->stop();
            }
        
            for (unsigned int i=0; i<m_numThreads; i++)
            {
                m_threadThreadList[i]->join();
                std::cerr << "joined " << i << std::endl;
                delete m_threadThreadList[i];
                delete m_threadList[i];
            }
        
            delete []m_threadList;
            delete []m_threadThreadList;
            delete []m_threadActiveList;
            delete []m_threadBusyList;
        }
        
        ThreadManager::ThreadManager(unsigned int numT, unsigned int chnk) :
                m_numThreads(numT), m_numWorkers(numT), m_masterQueue(), m_masterQueueMutex(),
                m_masterActiveMutex(), m_masterQueueCondVar(), m_masterActiveCondVar(),
                m_schedType(e_static)
        {
            using namespace std;
            try {
                m_threadList = new ThreadWorker*[m_numThreads];
                m_threadThreadList = new boost::thread*[m_numThreads];
                m_threadBusyList = new bool[m_numThreads];
                m_threadActiveList = new bool[m_numThreads];
            } catch (exception &e) {
                cerr << "Exception while allocating thread storage: "
                        << e.what() << endl;
                abort();
            }
            m_chunkSize = std::max(chnk, static_cast<unsigned int>(1)); // really, C++?
            unsigned int i = 0;
            while (i < m_numThreads)
            {
                ThreadWorker *tw;
                try {
                    tw = new ThreadWorker(this, i);
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
                } catch (...) {
                    std::cerr << "Exception while creating worker threads" << std::endl;
                    abort();
                }
        
                i++;
            }
        }
        
        void ThreadManager::queueJobs(std::vector<ThreadJob*> &joblist)
        {
            std::vector<ThreadJob *>::iterator it;
            for (it=joblist.begin(); it<joblist.end(); ++it)
            {
                queueJob(*it);
            }
        }
        
        void ThreadManager::queueJob(ThreadJob *job)
        {
            Lock masterQueueLock(m_masterQueueMutex); // locks the queue
            m_masterQueue.push(job);
            m_masterQueueCondVar.notify_all(); // alert a waiting thread.
        }   // queue unlocked
        
        bool ThreadManager::isWorking()
        {
            bool working = false;
            Lock masterActiveLock(m_masterActiveMutex);
            for (unsigned int i = 0; i < m_numWorkers; i++)
            {
                working = working || m_threadBusyList[i];
            }
            return working;
        }
        
        void ThreadManager::setChunkSize(unsigned int chnk)
        {
            Lock masterQueueLock(m_masterQueueMutex); // locks the queue
            m_chunkSize = std::max(chnk, static_cast<unsigned int>(1)); // really, C++?
        }

        void ThreadManager::setSchedType(SchedType s)
        {
            Lock masterQueueLock(m_masterQueueMutex); // locks the queue
            m_schedType = s;
        }

        void ThreadManager::wait()
        {
            bool working;
            Lock masterQueueLock(m_masterQueueMutex); // locks the queue
            working = isWorking();
            while (!m_masterQueue.empty() || working)
            {
                // while waiting, master queue is unlocked
                m_masterQueueCondVar.wait(masterQueueLock);
                // on exiting wait master queue is locked again
                working = isWorking();
            }
        }
        
        unsigned int ThreadManager::getNumWorkers()
        {
            return m_numWorkers;
        }
        
        void ThreadManager::setNumWorkers(const unsigned int num)
        {
            Lock masterActiveLock(m_masterActiveMutex); // locks the active
            m_numWorkers = std::min(num, m_numThreads);
            m_numWorkers = std::max(m_numWorkers, static_cast<unsigned int>(0));
            for (unsigned int i=0; i<m_numThreads; i++)
            {
                m_threadActiveList[i] = i<m_numWorkers ? true: false;
            }
            m_masterActiveCondVar.notify_all();
        } // lock on active released
        
        void ThreadManager::setNumWorkers()
        {
            setNumWorkers(m_numThreads);
        }
        
        unsigned int ThreadManager::getMaxNumWorkers()
        {
            return m_numThreads;
        }
        
        ThreadWorker::ThreadWorker(ThreadManager *tm, int workerNum) :
                m_threadManager(tm), m_workerQueue(),
                m_keepgoing(true), m_threadNum(workerNum)
        {
            // Nothing to see here
        
        }
        
        ThreadWorker::~ThreadWorker()
        {
            if (m_keepgoing)
            {
                std::cerr << "Warning: ThreadWorker: " << m_threadNum
                        << "destroyed while running!" << std::endl;
            }
        }
        
        void ThreadWorker::loadJobs()
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
                unsigned int numToLoad = getNumToLoad();
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


        unsigned int ThreadWorker::getNumToLoad()
        {
            unsigned int numToLoad;
            switch (m_threadManager->m_schedType)
            {
                case e_guided:
                    numToLoad = std::max(static_cast<unsigned long>(m_threadManager->m_chunkSize),
                            m_threadManager->m_masterQueue.size() / (2*m_threadManager->m_numWorkers +1));
                    break;

                case e_static:
                    numToLoad = m_threadManager->m_chunkSize;
                    break;
            }
            return numToLoad;
        }
        
        void ThreadWorker::waitForActive()
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
        
        void ThreadWorker::mainLoop()
        {
            while (m_keepgoing)
            {
                waitForActive();
                loadJobs();
                runJobs();
            }
        } // exiting here should terminate the thread
        
        void ThreadWorker::runJobs()
        {
            while (!m_workerQueue.empty())
            {
                ThreadJob * tj;
                try
                {
                    tj = m_workerQueue.front();
                    tj->setWorkerNum(m_threadNum);
                    tj->run();
                    m_workerQueue.pop();
                    delete tj;
                } catch(...)
                {
                    // something bad happened, probably time to die
                    // maybe signal ThreadManager
                }
            }
        }
        
    } // Thread
} /* namespace Nektar */
