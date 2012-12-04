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

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

namespace Nektar
{
    namespace Thread
    {

    	enum SchedType
    	{
    		e_guided,
    		e_dynamic
    	};

    	class ThreadManager;
    	typedef boost::shared_ptr<ThreadManager> ThreadManagerSharedPtr;
    	typedef LibUtilities::NekFactory< std::string, ThreadManager, unsigned int> ThreadManagerFactory;
    	LIB_UTILITIES_EXPORT ThreadManagerFactory& GetThreadManager();

        //typedef boost::unique_lock<boost::mutex> Lock;

        /**
         * Base class for tasks to be sent to the ThreadManager
         * to run.  Override the run() method with the desired task.
         */
        class ThreadJob {

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

        	void setWorkerNum(unsigned int num);
        
        protected:
            /**
             * Returns an integer identifying the worker thread the
             * job is running on.  Value will be 0...N, where N is
             * the number of active worker threads.
             */
            unsigned int getWorkerNum();
        
        private:
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
        class ThreadManager : public boost::enable_shared_from_this<ThreadManager>
        {

        public:
        	virtual ~ThreadManager();
        	virtual void queueJobs(std::vector<ThreadJob*>& joblist) = 0;
        	virtual void queueJob(ThreadJob* job) = 0;
        	virtual unsigned int getNumWorkers() = 0;
        	virtual void setNumWorkers(const unsigned int num) = 0;
        	virtual void setNumWorkers() = 0;
        	virtual unsigned int getMaxNumWorkers() = 0;
        	virtual void wait() = 0;
        	virtual void setChunkSize(unsigned int chnk) = 0;
        	virtual void setSchedType(SchedType s) = 0;

/*
        	static ThreadManager* createThreadManager(unsigned int numT,
        			unsigned int chnksz = 1)
        	{
        		if (instance == 0)
        		{
        			instance = new ThreadManager(numT, chnksz);
        		}
        		return instance;
        	}
*/

        	static ThreadManagerSharedPtr getInstance()
        	{
        		return instance;
        	}

        protected:
        	static ThreadManagerSharedPtr instance;

        };


        class ThreadHandle : public ThreadManager
        	{
        	public:
        		ThreadHandle(SchedType sched=e_dynamic);
        		ThreadHandle(SchedType sched, unsigned int chnk=1);
        		~ThreadHandle();

        		/**
        		 * Copy jobs into the master queue.
        		 * The ThreadJob pointers are copied into the master queue.
        		 * After jobs are run they are destructed.
        		 */
        		void queueJobs(std::vector<ThreadJob *> &joblist)
        		{
        			tm->queueJobs(joblist);
        		}
        		/**
        		 * Copy a job into the master queue.
        		 * The ThreadJob pointer is copied into the master queue.
        		 * After jobs are run they are destructed.
        		 */
        		void queueJob(ThreadJob *job)
        		{
        			tm->queueJob(job);
        		}
        		/**
        		 * Returns the number of currently active workers.
        		 */
        		unsigned int getNumWorkers()
        		{
        			return tm->getNumWorkers();
        		}
        		/**
        		 * Sets the number of currently active workers.
        		 */
        		void setNumWorkers(const unsigned int num)
        		{
        			tm->setNumWorkers(num);
        		}
        		/**
        		 * Sets the number of currently active workers to the maximum.
        		 */
        		void setNumWorkers()
        		{
        			tm->setNumWorkers();
        		}
        		/**
        		 * Gets the maximum number of workers.
        		 */
        		unsigned int getMaxNumWorkers()
        		{
        			return tm->getMaxNumWorkers();
        		}
        		/**
        		 * The calling thread (which should be the master thread) blocks
        		 * until all jobs submitted have finished.  This means on return
        		 * the queue will be empty and all workers idle.
        		 *
        		 * Currently this will deadlock if called when there are no active
        		 * workers.
        		 */
        		void wait()
        		{
        			tm->wait();
        		}
        		/**
        		 * Sets the m_chunkSize.
        		 * This is used by the scheduling type to determine how many jobs
        		 * are loaded onto each worker at once.
        		 */
        		void setChunkSize(unsigned int chnk)
        		{
        			tm->setChunkSize(chnk);
        		}
        		/**
        		 * Sets the scheduling type
        		 */
        		void setSchedType(SchedType sched)
        		{
        			tm->setSchedType(sched);
        		}

        	private:
        		ThreadManagerSharedPtr tm;

        		void setup(SchedType sched, unsigned int chnk);

        	};
    }
}
#endif /* THREAD_H_ */
