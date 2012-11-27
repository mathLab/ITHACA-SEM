/*
 * Thread.h
 *      Author: simon
 */

#ifndef THREAD_H_
#define THREAD_H_

#include <queue>
#include <vector>
#include <boost/thread/mutex.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/condition_variable.hpp>
#include <boost/thread/thread.hpp>

namespace Nektar
{
    namespace Thread
    {
        enum SchedType
        {
            e_guided,
            e_static
        };

        typedef boost::unique_lock<boost::mutex> Lock;
        
        class ThreadWorker;
        
        /**
         * Base class for tasks to be sent to the ThreadManager
         * to run.  Override the run() method with the desired task.
         */
        class ThreadJob {
        
        friend class ThreadWorker;
        
        public:
            /**
             * Base constructor
             */
        	ThreadJob();
            /**
             * Base destructor.
             */
        	virtual ~ThreadJob();
            /**
             * This method will be called when the task is loaded
             * onto a worker thread and is ready to run.
             */
        	virtual void run() = 0;
        
        protected:
            /**
             * Returns an integer identifying the worker thread the
             * job is running on.  Value will be 0...N, where N is
             * the number of active worker threads.
             */
            unsigned int getWorkerNum();
        
        private:
        	void setWorkerNum(unsigned int num);
        
        	// Member variables
        	unsigned int m_workerNum;
        };
        
        /**
         * The controller for the worker threads and jobs.
         *
         * There is only a single instance of this class.  When it instantiates
         * (through the static method createThreadManager) it spawns the
         * required number of worker threads.
         *
         * Jobs are sent to the worker threads by placing them in the master
         * queue with queueJob() and queueJobs().  Jobs will be run in
         * indeterminate order, however when threads are permitted to take more
         * than one job onto their local queues they will take and run them in
         * the order they were queued.  Jobs are permitted to run *immediately*
         * they are queued.  If this is not desired set the number of active
         * workers to 0 until all jobs are queued.
         *
         * Jobs are loaded from the master queue by active worker threads.
         * Currently there are two strategies to determine the number of jobs.
         * e_static always loads m_chunkSize jobs onto each worker,
         * e_guided uses a simple algorithm that loads a larger proportion of jobs
         * onto the worker as the size of the queue list is large, but scales down
         * to no fewer than m_chunkSize as the queue empties.
         *
         */ 
        class ThreadManager
        {
        
        friend class ThreadWorker;
        
        public:
        	~ThreadManager();
            /**
             * Copy jobs into the master queue.
             * The ThreadJob pointers are copied into the master queue.
             * After jobs are run they are destructed.
             */
        	void queueJobs(std::vector<ThreadJob *> &joblist);
            /**
             * Copy a job into the master queue.
             * The ThreadJob pointer is copied into the master queue.
             * After jobs are run they are destructed.
             */
        	void queueJob(ThreadJob *job);
            /**
             * Returns the number of currently active workers.
             */
        	unsigned int getNumWorkers();
            /**
             * Sets the number of currently active workers.
             */
        	void setNumWorkers(const unsigned int num);
            /**
             * Sets the number of currently active workers to the maximum.
             */
        	void setNumWorkers();
            /**
             * Gets the maximum number of workers.
             */
        	unsigned int getMaxNumWorkers();
            /**
             * The calling thread (which should be the master thread) blocks
             * until all jobs submitted have finished.  This means on return
             * the queue will be empty and all workers idle.
             *
             * Currently this will deadlock if called when there are no active
             * workers.
             */
        	void wait();
            /**
             * Sets the m_chunkSize.
             * This is used by the scheduling type to determine how many jobs
             * are loaded onto each worker at once.
             */
            void setChunkSize(unsigned int chnk);
            /**
             * Sets the scheduling type
             */
            void setSchedType(SchedType s);
            /**
             * Creates the ThreadManager instance.
             * If called more than once, the parameters are ignored and
             * a pointer to the original instance is returned.
             *
             * @param numT The number of worker threads to create
             * @param chnksz The initial m_chunkSize.
             */
        	static ThreadManager *createThreadManager(unsigned int numT,
        			unsigned int chnksz=1)
        	{
                if (instance == 0)
                {
                    instance = new ThreadManager(numT, chnksz);
                }
        		return instance;
        	}
            /**
             * Returns a pointer to the ThreadManager
             */
            static ThreadManager *getInstance()
            {
                return instance;
            }

        
        
        private:
        	// Member functions
        	ThreadManager(unsigned int numWorkers, unsigned int chnk);
        	ThreadManager(const ThreadManager &); //not defined
        	ThreadManager(); // not defined
        	bool isWorking();
        
        	// Member variables (NB remember constructor)
        	const unsigned int m_numThreads;
        	unsigned int m_numWorkers;
        	std::queue<ThreadJob *> m_masterQueue;
        	boost::mutex m_masterQueueMutex;
        	boost::mutex m_masterActiveMutex;
        	boost::condition_variable m_masterQueueCondVar;
        	boost::condition_variable m_masterActiveCondVar;
        	ThreadWorker **m_threadList;
        	boost::thread **m_threadThreadList;
        	bool *m_threadBusyList;
        	bool *m_threadActiveList;
        	unsigned int m_chunkSize;
        	SchedType m_schedType;
        
        	// Static Member variables
        	static ThreadManager *instance;
        
        };
        
        
        class ThreadWorker {
        
        public:
        	ThreadWorker(ThreadManager *threadManager, int workerNum);
        	virtual ~ThreadWorker();
        	void operator()() { mainLoop(); };
        	unsigned int getWorkerNum() { return m_threadNum; };
        	void stop() { m_keepgoing = false;} ;
        
        private:
        	ThreadWorker();
        	ThreadWorker(const ThreadWorker &);
        	void mainLoop();
        	void loadJobs();
            inline unsigned int getNumToLoad();
        	void waitForActive();
        	void runJobs();
        
        	// Member variables
        	ThreadManager *m_threadManager;
        	std::queue<ThreadJob *> m_workerQueue;
        	bool m_keepgoing;
        	unsigned int m_threadNum;
        
        };
        
        
    }
}
#endif /* THREAD_H_ */
