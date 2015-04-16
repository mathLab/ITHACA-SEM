/*
 * Thread.h
 *      Author: simon
 */

#ifndef THREADBOOST_H_
#define THREADBOOST_H_

#include "LibUtilities/BasicUtils/Thread.h"
#include <queue>
#include <vector>
#include <map>
#include <boost/thread/barrier.hpp>

#include "LibUtilities/Memory/NekMemoryManager.hpp"

namespace Nektar
{
    namespace Thread
    {

        typedef boost::unique_lock<boost::mutex> Lock;
        
        class ThreadWorkerBoost;

        /**
         * @brief Implementation of ThreadManager using Boost threads.
         */
        class ThreadManagerBoost: public ThreadManager
        {
        
        /**
         * So the workers can access the master queue and locks.
         */
        friend class ThreadWorkerBoost;

        public:
            /**
             * @brief Constructs a ThreadManagerBoost.
             * @param numWorkers The number of threads to start (including master thread).
             * @note Do not use, use factory instead.
             */
            ThreadManagerBoost(unsigned int numWorkers);
            /**
             * @brief Shuts down threading.
             *
             * Terminates all running threads (they will finish their current job),
             * releases resources and destructs.
             */
            virtual ~ThreadManagerBoost();
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
            virtual const std::string& GetType() const;


            /**
             * @brief Called by the factory method.
             */
            static ThreadManagerSharedPtr Create(unsigned int numT)
            {
//                if (!m_instance)
//                {
//                    m_instance = boost::shared_ptr<ThreadManager>(new ThreadManagerBoost(numT));
////                    m_instance = MemoryManager<ThreadManagerBoost>::AllocateSharedPtr(numT);
//                }
//                return m_instance;
                return boost::shared_ptr<ThreadManager>(new ThreadManagerBoost(numT));
            }

        private:
            ThreadManagerBoost();
            ThreadManagerBoost(const ThreadManagerBoost&);
            bool IsWorking();
            void SetNumWorkersImpl(const unsigned int num);

            // Member variables
            const unsigned int m_numThreads;
            unsigned int m_numWorkers;
            std::queue<ThreadJob*> m_masterQueue;
            boost::mutex m_masterQueueMutex;
            boost::mutex m_masterActiveMutex;
            boost::condition_variable m_masterQueueCondVar;
            boost::condition_variable m_masterActiveCondVar;
            ThreadWorkerBoost** m_threadList;
            boost::thread** m_threadThreadList;
            boost::thread::id m_masterThreadId;
            bool* m_threadBusyList;
            bool* m_threadActiveList;
            unsigned int m_chunkSize;
            SchedType m_schedType;
            boost::barrier *m_barrier;
            std::map<boost::thread::id, unsigned int> m_threadMap;
            static std::string className;
            std::string m_type;
        };

        /**
         * @brief Implementation class for ThreadManagerBoost.
         *
         * Each instance of this class corresponds to a worker thread.
         * Instances manage their own queue of jobs to run, grabbing new
         * jobs from the master queue when it is exhausted.
         *
         */
        class ThreadWorkerBoost
        {

        public:
            /**
             * @brief Constructor
             * @param threadManager Pointer to the ThreadManagerBoost that is controlling this worker.
             * @param workerNum Unique number from 0..(number_of_threads - 1)
             *
             * Called by the ThreadManagerBoost instance.
             */
            ThreadWorkerBoost(ThreadManagerBoost *threadManager, unsigned int workerNum);
            /**
             * @brief Destructor.
             *
             * Winds up this thread's execution.  Jobs in its queue are lost.
             */
            ~ThreadWorkerBoost();
            /**
             * @brief This provides the interface that boost::thread uses to start the worker.
             */
            void operator()() { MainLoop(); };
            /**
             * @brief Return the index of the worker thread.
             * @returns Index of worker thread, an integer between 0 and (number_of_threads - 1)
             */
            unsigned int GetWorkerNum() { return m_threadNum; };
            /**
             * @brief A signal to shut down.
             *
             * If this method is called the worker will shut down.  Used by
             * the ThreadManagerBoost to stop threading.
             */
            void Stop() { m_keepgoing = false;} ;

        private:
            ThreadWorkerBoost();
            ThreadWorkerBoost(const ThreadWorkerBoost &);
            void MainLoop();
            void LoadJobs();
            unsigned int GetNumToLoad();
            void WaitForActive();
            void RunJobs();

            // Member variables
            ThreadManagerBoost *m_threadManager;
            std::queue<ThreadJob *> m_workerQueue;
            bool m_keepgoing;
            unsigned int m_threadNum;

        };


    }
}
#endif /* THREADBOOST_H_ */
