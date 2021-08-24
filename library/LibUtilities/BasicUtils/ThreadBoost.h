///////////////////////////////////////////////////////////////////////////////
//
// File ThreadBoost.h
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

#ifndef NEKTAR_LIBUTILITIES_THREADBOOST_H_
#define NEKTAR_LIBUTILITIES_THREADBOOST_H_

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
        /// Constructs a ThreadManagerBoost.
        ThreadManagerBoost(unsigned int numWorkers);
        /// Shuts down threading.
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


        /// Called by the factory method.
        static ThreadManagerSharedPtr Create(unsigned int numT)
        {
            return std::shared_ptr<ThreadManager>(
                new ThreadManagerBoost(numT));
        }

    private:
        ThreadManagerBoost();
        ThreadManagerBoost(const ThreadManagerBoost&);
        bool IsWorking();
        void SetNumWorkersImpl(const unsigned int num);

        // Member variables
        const unsigned int                         m_numThreads;
        unsigned int                               m_numWorkers;
        std::queue<ThreadJob*>                     m_masterQueue;
        boost::mutex                               m_masterQueueMutex;
        boost::mutex                               m_masterActiveMutex;
        boost::condition_variable                  m_masterQueueCondVar;
        boost::condition_variable                  m_masterActiveCondVar;
        ThreadWorkerBoost**                        m_threadList;
        boost::thread**                            m_threadThreadList;
        boost::thread::id                          m_masterThreadId;
        bool*                                      m_threadBusyList;
        bool*                                      m_threadActiveList;
        unsigned int                               m_chunkSize;
        SchedType                                  m_schedType;
        boost::barrier*                            m_barrier;
        std::map<boost::thread::id, unsigned int>  m_threadMap;
        static std::string                         className;
        std::string                                m_type;
};

/**
 * @brief Implementation class for ThreadManagerBoost.
 *
 * Each instance of this class corresponds to a worker thread.
 * Instances manage their own queue of jobs to run, grabbing new
 * jobs from the master queue when it is exhausted.
 */
class ThreadWorkerBoost
{
    public:
        /// Constructor
        ThreadWorkerBoost(ThreadManagerBoost *threadManager,
                          unsigned int workerNum);
        /// Destructor.
        ~ThreadWorkerBoost();
        /// This provides the interface that boost::thread uses to start the
        /// worker.
        void operator()() { MainLoop(); };
        /**
         * @brief Return the index of the worker thread.
         * @return Index of worker thread, an integer between 0 and
         *        (number_of_threads - 1)
         */
        unsigned int GetWorkerNum() { return m_threadNum; };
        /**
         * @brief A signal to shut down.
         *
         * If this method is called the worker will shut down.  Used by the
         * ThreadManagerBoost to stop threading.
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
        ThreadManagerBoost*      m_threadManager;
        std::queue<ThreadJob *>  m_workerQueue;
        bool                     m_keepgoing;
        unsigned int             m_threadNum;
};


}
}
#endif /* THREADBOOST_H_ */
