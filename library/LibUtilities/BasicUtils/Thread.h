///////////////////////////////////////////////////////////////////////////////
//
// File Thread.h
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
// Description:
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_THREAD_H_
#define NEKTAR_LIB_UTILITIES_THREAD_H_

#include <queue>
#include <vector>
#include <memory>
#include <boost/thread/mutex.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/condition_variable.hpp>
#include <boost/thread/thread.hpp>

#include <LibUtilities/BasicUtils/NekFactory.hpp>

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
typedef std::shared_ptr<ThreadManager> ThreadManagerSharedPtr;
typedef LibUtilities::NekFactory<std::string, ThreadManager, unsigned int>
    ThreadManagerFactory;
LIB_UTILITIES_EXPORT ThreadManagerFactory& GetThreadManagerFactory();

/**
 * @brief Base class for tasks to be sent to the ThreadManager to run.
 *
 * For a parallel region the overall task should be divided into tasklets, each
 * capable of running independently of other tasklets, in any order.  Each of
 * these tasklets should be represented by an instance of a subclass of
 * ThreadJob.
 *
 * Exactly how to partition the overall task will be problem dependent, but the
 * ideal case is to have as many tasklets as there are worker threads (available
 * through ThreadManager::GetNumWorkers() ) and for each tasklet to take the
 * same amount of computational effort.  Clearly, the tasklets should not
 * overwrite one another's data, and data locations they are writing to should
 * be sufficiently far apart to avoid cache ping pong.
 *
 * Subclasses may define as many variables and functions as are necessary.
 * Instance of the subclass will be destructed once the ThreadManager has
 * finished the Run() method.
 */
class ThreadJob
{
    public:
        /// Base constructor
        LIB_UTILITIES_EXPORT ThreadJob();
        /// Base destructor.
        LIB_UTILITIES_EXPORT virtual ~ThreadJob();

        /**
         * This method will be called when the task is loaded
         * onto a worker thread and is ready to run.  When Run
         * has finished this instance will be destructed.
         */
        LIB_UTILITIES_EXPORT virtual void Run() = 0;

        /// Set number of worker threads.
        LIB_UTILITIES_EXPORT void SetWorkerNum(unsigned int num);

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
 * @brief The interface class for the controller for worker threads and jobs.
 *
 * There may be multiple instances of implementations of this class.
 * They are instantiated by ThreadMaster::CreateInstance(string, int).
 * Use ThreadMaster::GetInstance(string)
 * to get a pointer to the instance.  Each ThreadManager is associated
 * in ThreadMaster with a string.
 *
 * Jobs are sent to the worker threads by placing them in the
 * queue with QueueJob() and QueueJobs().  Jobs may be run in
 * indeterminate order.  Jobs are permitted to run *immediately*
 * they are queued.  If this is not desired set the number of active
 * workers to 0 until all jobs are queued.
 *
 * Jobs are taken from the master queue by active worker threads.
 *
 * We have the concept of the *master thread*.  This is the thread
 * that instantiates the ThreadManager, calls QueueJob(), Wait()
 * and so on.  This thread must be included in the number of threads
 * started when the ThreadManager is created, so that no more than
 * this number of threads are ever executing work at any time.  The
 * master thread is assumed to be working unless it is in Wait().
 *
 * Thus if the ThreadManager is called with CreateInstance(string s, int N)
 * there should be N threads /em including the master thread.  If the master
 * thread calls Wait() (which causes the master thread to sleep until
 * the job queue is empty) then the master thread should become available
 * to run jobs from the queue.
 *
 * Only the master thread should call ThreadManager methods.  Although
 * it should always be possible for code to identify whether it's the
 * master thread (through design) there is a InThread() method that returns
 * true if the thread is a worker.
 *
 */
class ThreadManager : public std::enable_shared_from_this<ThreadManager>
{
    public:
        /// Destructor.
        LIB_UTILITIES_EXPORT virtual ~ThreadManager();
        /**
         * @brief Pass a list of tasklets to the master queue.
         * @param joblist Vector of ThreadJob pointers.
         *
         * The list of jobs is copied into the master queue.  Jobs may be
         * available for running *immediately*, even before the list has been
         * fully copied.  This can have consequences for the scheduling.  If
         * this is an issue then suspend the workers with SetNumWorkers(0) until
         * the jobs are queued.
         *
         * @see SchedType
         */
        LIB_UTILITIES_EXPORT virtual void QueueJobs(std::vector<ThreadJob*>& joblist) = 0;
        /**
         * @brief Pass a single job to the master queue.
         * @param job A pointer to a ThreadJob subclass.
         *
         * The job may become available for running immediately.  If this is an
         * issue then suspend the workers with SetNumWorkers(0) until the jobs
         * are queued.
         */
        LIB_UTILITIES_EXPORT virtual void QueueJob(ThreadJob* job) = 0;
        /**
         * @brief Return the number of active workers.
         *
         * Active workers are threads that are either running jobs
         * or are waiting for jobs to be queued.
         */
        LIB_UTILITIES_EXPORT virtual unsigned int GetNumWorkers() = 0;
        /**
         * @brief Returns the worker number of the executing thread.
         *
         * Returns an unsigned int between 0 and N-1 where N is the number of
         * active worker threads.  Repeated calls from within this thread will
         * always return the same value and the value will be the same as
         * returned from ThreadJob.GetWorkerNum().  The same thread will run a
         * job until it finishes.
         *
         * Although if there are active threads then thread 0 is always one of
         * them, it is possible that thread 0 does not run for a given set of
         * jobs.  For example, if there are 4 active threads and 3 jobs are
         * submitted with a e_static scheduling strategy and a chunksize of 1,
         * then it is possible that threads 1,2, and 3 pick up the jobs and
         * thread 0 remains idle.
         *
         * Returns 0 if called by non-thread.
         */
        LIB_UTILITIES_EXPORT virtual unsigned int GetWorkerNum() = 0;
        /**
         * @brief Sets the number of active workers.
         * @param num The number of active workers.
         *
         * Active workers are threads that are either running jobs or are
         * waiting for jobs to be queued.
         *
         * If num is greater than the maximum allowed number of active workers,
         * then the maximum value will be used instead.
         */
        LIB_UTILITIES_EXPORT virtual void SetNumWorkers(const unsigned int num) = 0;
        /**
         * @brief Sets the number of active workers to the maximum.
         *
         * Sets the number of active workers to the maximum available.
         */
        LIB_UTILITIES_EXPORT virtual void SetNumWorkers() = 0;
        /**
         * @brief Gets the maximum available number of threads.
         * @return The maximum number of workers.
         */
        LIB_UTILITIES_EXPORT virtual unsigned int GetMaxNumWorkers() = 0;
        /**
         * @brief Waits until all queued jobs are finished.
         *
         * If there are no jobs running or queued this method returns
         * immediately.  Otherwise it blocks until the queue is empty and the
         * worker threads are idle.
         *
         * Implementations *must* ensure that trivial deadlocks are not possible
         * from this method, that is, that this code:
         * @code
         *  // assume ThreadManager* tm
         *  // assume SomeJob is subclass of ThreadJob
         *  // assume SomeJob job
         *
         *  tm->SetNumWorkers(0);
         *  tm->QueueJob(job);
         *  tm->Wait();
         * @endcode
         * does not wait forever.  Since the master thread is counted in the
         * number of worker threads, implementations should increase the number
         * of active workers by 1 on entering Wait().
         */
        LIB_UTILITIES_EXPORT virtual void Wait() = 0;
        /**
         * @brief Controls how many jobs are sent to each worker at a time.
         *
         * The exact meaning of this parameter depends on the current scheduling
         * algorithm.
         *
         * @see SchedType
         * @see SetSchedType()
         */
        LIB_UTILITIES_EXPORT virtual void SetChunkSize(unsigned int chnk) = 0;
        /**
         * @brief Sets the current scheduling algorithm.
         * @see SetChunkSize()
         */
        LIB_UTILITIES_EXPORT virtual void SetSchedType(SchedType s) = 0;
        /**
         * @brief Indicates whether the code is in a worker thread or not.
         * @return True if the caller is in a worker thread.
         */
        LIB_UTILITIES_EXPORT virtual bool InThread() = 0;
        /**
         * @brief A calling threads holds until all active threads call this
         * method.
         *
         * When called, the calling thread will sleep until all active workers
         * have called this method.  Once all have done so all threads awake and
         * continue execution.
         *
         * @note Behaviour is likely undefined if the number of active workers
         * is altered after a thread has called this method.  It is only safe to
         * call SetNumWorkers() when no threads are holding.
         */
        LIB_UTILITIES_EXPORT virtual void Hold() = 0;
        /**
         * @brief Returns a description of the type of threading.
         *
         * E.g. "Threading with Boost"
         */
        LIB_UTILITIES_EXPORT virtual const std::string& GetType() const = 0;

        /// ThreadManager implementation.
        LIB_UTILITIES_EXPORT virtual bool IsInitialised();

        inline int GetThrFromPartition(int pPartition)
        {
            return pPartition % GetMaxNumWorkers();
        }

        inline int GetRankFromPartition(int pPartition)
        {
            return pPartition / GetMaxNumWorkers();
        }

        inline int GetPartitionFromRankThr(int pRank, unsigned int pThr)
        {
            return pRank * GetMaxNumWorkers() + pThr;
        }
};

typedef boost::unique_lock<boost::shared_mutex> WriteLock;
typedef boost::shared_lock<boost::shared_mutex> ReadLock;

/**
 * A class to manage multiple ThreadManagers.  It also acts as a cut-out during
 * static initialisation, where code attempts to grab a ThreadManager before any
 * can have been initialised.  When this happens the ThreadMaster creates a
 * simple implementation of ThreadManager, ThreadStartupManager, which behaves
 * as a single-threaded ThreadManager that fails if anything non-trivial is
 * attempted. This means code should be cautious about caching the value of
 * GetInstance(string), as it might change once that particular threading system
 * is initialised.
 *
 * Multiple ThreadManagers should now be possible.  Each is identified with a
 * string and can be recovered by calling
 * Thread::GetThreadMaster().GetInstance(string).
 *
 * ThreadMaster is thread-safe apart from CreateInstance.  Call CreateInstance
 * in single threaded code only.
 */
class ThreadMaster
{
    private:
        std::vector<ThreadManagerSharedPtr> m_threadManagers;
        boost::shared_mutex m_mutex;
        std::string m_threadingType;

    public:
        enum ThreadManagerName {
                                SessionJob,
                                THREADMANAGER_MAX
                                };
        /// Constructor
        LIB_UTILITIES_EXPORT ThreadMaster();
        /// Destructor
        LIB_UTILITIES_EXPORT ~ThreadMaster();
        /// Sets what ThreadManagers will be created in CreateInstance.
        LIB_UTILITIES_EXPORT void SetThreadingType(const std::string &p_type);
        /// Gets the ThreadManager associated with string s.
        LIB_UTILITIES_EXPORT ThreadManagerSharedPtr& GetInstance(const ThreadManagerName t);
        /// Creates an instance of a ThreadManager (which one is determined by
        /// a previous call to SetThreadingType) and associates it with
        /// the string s.
        LIB_UTILITIES_EXPORT ThreadManagerSharedPtr CreateInstance(const ThreadManagerName t,
            unsigned int nThr);

};
LIB_UTILITIES_EXPORT ThreadMaster& GetThreadMaster();


/**
 * @brief A default ThreadManager.
 *
 * This will be returned by ThreadMaster if a ThreadManager has not been
 * initialised, such as if the code is still in static initialisation.
 *
 * This manager pretends to be a ThreadManager with 1 thread.  It will
 * cause an error if anything more than trivial functions are called.
 */
class ThreadStartupManager: public ThreadManager
{
    public:
        ThreadStartupManager();
        ThreadStartupManager(const ThreadStartupManager& src) = default;
        virtual ~ThreadStartupManager();
        virtual void QueueJobs(std::vector<ThreadJob*>& joblist);
        virtual void QueueJob(ThreadJob* job);
        virtual unsigned int GetNumWorkers();
        virtual unsigned int GetWorkerNum();
        virtual void SetNumWorkers(const unsigned int num);
        virtual void SetNumWorkers();
        virtual unsigned int GetMaxNumWorkers();
        virtual void Wait();
        virtual void SetChunkSize(unsigned int chnk);
        virtual void SetSchedType(SchedType s);
        virtual bool InThread();
        virtual void Hold();
        virtual bool IsInitialised();
        virtual const std::string& GetType() const;
    private:
        // Do not allow assignment as m_type is const
        ThreadStartupManager& operator=(const ThreadStartupManager& src);

        const std::string m_type;
};

}

}
#endif /* THREAD_H_ */
