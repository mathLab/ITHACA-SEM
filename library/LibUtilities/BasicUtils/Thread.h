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

    	/**
    	 * @brief Identifies the algorithm for scheduling.
    	 *
    	 * Currently there are two strategies to determine the number of jobs.
         * e_static always loads m_chunkSize jobs onto each worker,
         * e_guided uses a simple algorithm that loads a larger proportion of jobs
         * onto the worker as the size of the queue list is large, but scales down
         * to no fewer than m_chunkSize as the queue empties.
         *
     	 * @see ThreadManager::SetChunkSize()
         *
    	 */
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
         * @brief Base class for tasks to be sent to the ThreadManager to run.
         *
         * For a parallel region the overall task should be divided into tasklets, each
         * capable of running independently of other tasklets, in any order.  Each of
         * these tasklets should be represented by an instance of a subclass of ThreadJob.
         * Exactly how to partition the overall task will be problem dependent, but the
         * ideal case is to have as many tasklets as there are worker threads (available
         * through ThreadManager::GetNumWorkers() ) and for each tasklet to take the same
         * amount of computational effort.  Clearly, the tasklets should not overwrite
         * one another's data, and data locations they are writing to should be
         * sufficiently far apart to avoid cache ping pong.
         *
         * Subclasses may define as many variables and functions as are necessary.
         * Instance of the subclass will be destructed once the ThreadManager has
         * finished the Run() method.
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
             * onto a worker thread and is ready to run.  When Run
             * has finished this instance will be destructed.
             */
        	virtual void Run() = 0;

        	/**
        	 * Part of the thread interface.  Do not use unless you're
        	 * an implementation of ThreadManager.
        	 * @warning Do not use: needs to be public so thread implementation can call it.
        	 */
        	void SetWorkerNum(unsigned int num);
        
        protected:
            /**
             * Returns an integer identifying the worker thread the
             * job is running on.  Value will be 0...N, where N is
             * the number of active worker threads.
             */
            unsigned int GetWorkerNum();
        
        private:
            unsigned int m_workerNum;

        };
        
        /**
         * @brief The controller for the worker threads and jobs.
         *
         * There is only a single instance of this class.  When it instantiates
         * (through the static method createThreadManager) it spawns the
         * required number of worker threads.  Use GetInstance() to get a pointer
         * to the instance.
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
         *
         * We have the concept of the *master thread*.  This is the thread
         * that instantiates the ThreadManager, calls QueueJob(), Wait()
         * and so on.  This thread must be included in the number of threads
         * started when the ThreadManager is created, so that no more than
         * this number of threads are ever executing work at any time.  The
         * master thread is assumed to be working unless it is in Wait().
         *
         * Only the master thread should call ThreadManager methods.
         *
         * The simplest and preferred way to interact with the ThreadManager is to
         * create a ThreadHandle instance, and use that.
         *
         *
         */ 
        class ThreadManager : public boost::enable_shared_from_this<ThreadManager>
        {

        public:
        	/**
        	 * @brief Destructor.
        	 *
        	 * Should stop worker threads and clean up.
        	 * This shouldn't be called until the program is exiting
        	 * anyway, some implementations may need to unlock threads
        	 * to allow clean exit.
        	 */
        	virtual ~ThreadManager();
        	/**
        	 * @brief Pass a list of tasklets to the master queue.
        	 * @param joblist Vector of ThreadJob pointers.
        	 *
        	 * The list of jobs is copied into the master queue.
        	 * Jobs may be available for running *immediately*,
        	 * even before the list has been fully copied.  This can
        	 * have consequences for the scheduling.
        	 * If this is an issue then suspend the workers with SetNumWorkers(0)
        	 * until the jobs are queued.
        	 *
        	 * @see SchedType
        	 */
        	virtual void QueueJobs(std::vector<ThreadJob*>& joblist) = 0;
        	/**
        	 * @brief Pass a single job to the master queue.
        	 * @param job A pointer to a ThreadJob subclass.
        	 *
        	 * The job may become available for running immediately.
        	 */
        	virtual void QueueJob(ThreadJob* job) = 0;
        	/**
        	 * @brief Return the number of active workers.
        	 *
        	 * Active workers are threads that are either running jobs
        	 * or are waiting for jobs to be queued.
        	 */
        	virtual unsigned int GetNumWorkers() = 0;
        	/**
        	 * @brief Sets the number of active workers.
        	 * @param num The number of active workers.
        	 *
        	 * Active workers are threads that are either running jobs
        	 * or are waiting for jobs to be queued.
        	 *
        	 * If num is greater than the maximum allowed number
        	 * of active workers, then the maximum value will be
        	 * used instead.
        	 */
        	virtual void SetNumWorkers(const unsigned int num) = 0;
        	/**
        	 * @brief Sets the number of active workers to the maximum.
        	 *
        	 * Sets the number of active workers to the maximum available.
        	 */
        	virtual void SetNumWorkers() = 0;
        	/**
        	 * @brief Gets the maximum available number of threads.
        	 * @return The maximum number of workers.
        	 */
        	virtual unsigned int GetMaxNumWorkers() = 0;
        	/**
        	 * @brief Waits until all queued jobs are finished.
        	 *
        	 * If there are no jobs running or queued this method returns
        	 * immediately.  Otherwise it blocks until the queue is empty
        	 * and the worker threads are idle.
        	 *
        	 * Implementations *must*
        	 * ensure that trivial deadlocks are not possible from this method,
        	 * that is, that this code:
        	 * @code
        	 *  // assume ThreadManager* tm
        	 *  // assume SomeJob is subclass of ThreadJob
        	 *  // assume SomeJob job
        	 *
        	 *  tm->SetNumWorkers(0);
        	 *  tm->QueueJob(job);
        	 *  tm->Wait();
        	 * @endcode
        	 * does not wait forever.  Since the master thread is counted in
        	 * the number of worker threads, implementations should increase
        	 * the number of active workers by 1 on entering Wait().
        	 */
        	virtual void Wait() = 0;
        	/**
        	 * @brief Controls how many jobs are sent to each worker at a time.
        	 *
        	 * The exact meaning of this parameter depends on the current scheduling algorithm.
        	 *
        	 * @see SchedType
        	 * @see SetSchedType()
        	 */
        	virtual void SetChunkSize(unsigned int chnk) = 0;
        	/**
        	 * @brief Sets the current scheduling algorithm.
        	 * @see SetChunkSize()
        	 */
        	virtual void SetSchedType(SchedType s) = 0;
        	/**
        	 * @brief Indicates whether the code is in a worker thread or not.
        	 * @return True if the caller is in a worker thread.
        	 */
        	virtual bool InThread() = 0;
        	/**
        	 * @brief Returns pointer to the single instance.
        	 * @return Pointer to the single instance.
        	 */
        	static ThreadManagerSharedPtr GetInstance()
        	{
        		return m_instance;
        	}

        protected:
        	/**
        	 * @brief Pointer to an instance of a ThreadManager.
        	 *
        	 * Implementations must set this variable when they are instantiated.
        	 */
        	static ThreadManagerSharedPtr m_instance;

        };

        /**
         * @brief Proxy class for ThreadManager.
         *
         * This is the preferred way to use threading.  It provides a slightly
         * simpler interface than the ThreadManager class.  One can create
         * an instance of ThreadHandle on the stack, use it, and then allow it
         * to go out of scope without affecting the underlying ThreadManager.
         * On construction the number of active workers is set to the maximum.
         * The constructors also provide defaults for the schedule type and
         * chunksize.
         *
         * When the destructor is called, such as when a stack instance of
         * this class goes out of scope, Wait() is called.  This means one can
         * use the scope of the instance to act as a barrier for parallel activity.
         * e.g.
         * @code
         *   //assume SomeJob is subclass of ThreadJob
         *   //assume job1, job2, ... are instances of SomeJob
         *   //assume some implementation of ThreadManager has been constructed.
         *
         *   void foo(job1, job2, job3, ...)
         *   {
         *     ThreadHandle THandle; // sets everything to default
         *     THandle.QueueJob(job1);
         *     THandle.QueueJob(job2);
         *     THandle.QueueJob(job3);
         *     //...
         *   } // return from foo will only happen when all jobs are finished.
         * @endcode
         */
        class ThreadHandle : public ThreadManager
        	{
        	public:
        		/**
        		 * @brief Constructor
        		 * @param sched The desired scheduling algorithm.
        		 *
        		 * Construct a handle for the ThreadManager, set the number
        		 * of active workers to the maximum, set the scheduling algorithm
        		 * to sched, set the chunksize to 1.
        		 */
        		ThreadHandle(SchedType sched=e_dynamic);
        		/**
        		 * @brief Constructor
        		 * @param sched The desired scheduling algorithm.
        		 * @param chnk The desired chunksize.
        		 *
        		 * Construct a handle for the ThreadManager, set the number
        		 * of active workers to the maximum, set the scheduling algorithm
        		 * to sched, set the chunksize to chnk.
        		 */
            	ThreadHandle(SchedType sched, unsigned int chnk=1);
            	/**
            	 * @brief Remove instance from scope when all jobs have finished.
            	 *
            	 * The ThreadHandle will be destructed when the job queue is empty
            	 * and all running jobs are finished.  The underlying ThreadManager
            	 * is not destructed.
            	 */
        		virtual ~ThreadHandle();
        		void QueueJobs(std::vector<ThreadJob *> &joblist)
        		{
        			m_tm->QueueJobs(joblist);
        		}
         		void QueueJob(ThreadJob *job)
        		{
        			m_tm->QueueJob(job);
        		}
         		unsigned int GetNumWorkers()
        		{
        			return m_tm->GetNumWorkers();
        		}
         		void SetNumWorkers(const unsigned int num)
        		{
        			m_tm->SetNumWorkers(num);
        		}
         		void SetNumWorkers()
        		{
        			m_tm->SetNumWorkers();
        		}
        		unsigned int GetMaxNumWorkers()
        		{
        			return m_tm->GetMaxNumWorkers();
        		}
         		void Wait()
        		{
        			m_tm->Wait();
        		}
        		void SetChunkSize(unsigned int chnk)
        		{
        			m_tm->SetChunkSize(chnk);
        		}
         		void SetSchedType(SchedType sched)
        		{
        			m_tm->SetSchedType(sched);
        		}
        		bool InThread()
        		{
        			return m_tm->InThread();
        		}

        	private:
        		ThreadManagerSharedPtr m_tm;

        		void Setup(SchedType sched, unsigned int chnk);

        	};
    }
}
#endif /* THREAD_H_ */
