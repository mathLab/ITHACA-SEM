/*
 * Thread.h
 *      Author: simon
 */

#ifndef THREADBOOST_H_
#define THREADBOOST_H_

#include "Thread.h"
#include <queue>
#include <vector>
#include <boost/thread/mutex.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/condition_variable.hpp>
#include <boost/thread/thread.hpp>

#include "LibUtilities/Memory/NekMemoryManager.hpp"

namespace Nektar
{
    namespace Thread
    {

    	typedef boost::unique_lock<boost::mutex> Lock;
        
        class ThreadWorkerBoost;
        
        class ThreadManagerBoost: public ThreadManager
        {
        
        friend class ThreadWorkerBoost;
        
        public:
    		ThreadManagerBoost(unsigned int numWorkers);
    		~ThreadManagerBoost();
        	void queueJobs(std::vector<ThreadJob*>& joblist);
        	void queueJob(ThreadJob* job);
        	unsigned int getNumWorkers();
        	void setNumWorkers(const unsigned int num);
        	void setNumWorkers();
        	unsigned int getMaxNumWorkers();
        	void wait();
        	void setChunkSize(unsigned int chnk);
        	void setSchedType(SchedType s);

        	static ThreadManagerSharedPtr create(unsigned int numT)
        	{
        		if (!instance)
        		{
        			instance = MemoryManager<ThreadManagerBoost>::AllocateSharedPtr(numT);
        			//instance = new ThreadManagerBoost(numT);
        		}
        		return instance;
        	}

        private:
        	ThreadManagerBoost(const ThreadManagerBoost&);
        	ThreadManagerBoost();
        	bool isWorking();
        	void setNumWorkersImpl(const unsigned int num);
        	const unsigned int m_numThreads;
        	unsigned int m_numWorkers;
        	std::queue<ThreadJob*> m_masterQueue;
        	boost::mutex m_masterQueueMutex;
        	boost::mutex m_masterActiveMutex;
        	boost::condition_variable m_masterQueueCondVar;
        	boost::condition_variable m_masterActiveCondVar;
        	ThreadWorkerBoost** m_threadList;
        	boost::thread** m_threadThreadList;
        	bool* m_threadBusyList;
        	bool* m_threadActiveList;
        	unsigned int m_chunkSize;
        	SchedType m_schedType;
        	static std::string className;
        };


        class ThreadWorkerBoost
        {

        public:
        	ThreadWorkerBoost(ThreadManagerBoost *threadManager, int workerNum);
        	~ThreadWorkerBoost();
        	void operator()() { mainLoop(); };
        	unsigned int getWorkerNum() { return m_threadNum; };
        	void stop() { m_keepgoing = false;} ;

        private:
        	ThreadWorkerBoost();
        	ThreadWorkerBoost(const ThreadWorkerBoost &);
        	void mainLoop();
        	void loadJobs();
        	inline unsigned int getNumToLoad();
        	void waitForActive();
        	void runJobs();

        	// Member variables
        	ThreadManagerBoost *m_threadManager;
        	std::queue<ThreadJob *> m_workerQueue;
        	bool m_keepgoing;
        	unsigned int m_threadNum;

        };


    }
}
#endif /* THREADBOOST_H_ */
