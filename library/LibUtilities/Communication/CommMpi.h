///////////////////////////////////////////////////////////////////////////////
//
// File CommMpi.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
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
// Description: CommMpi header
//
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_LIB_UTILITIES_COMMMPI_H
#define NEKTAR_LIB_UTILITIES_COMMMPI_H

#include <string>
#include <mpi.h>

#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

#ifndef MPI_SYNC
#define MPISYNC 0
#else
#define MPISYNC 1
#endif

namespace Nektar
{
    namespace LibUtilities
    {
        // Forward declarations
        class CommMpi;

        /// Pointer to a Communicator object.
        typedef boost::shared_ptr<CommMpi> CommMpiSharedPtr;

        /// A global linear system.
        class CommMpi : public Comm
        {
        public:
            /// Creates an instance of this class
            static CommSharedPtr create(int narg, char* arg[])
            {
                return MemoryManager<CommMpi>::AllocateSharedPtr(narg,arg);
            }

            /// Name of class
            static std::string className;

            CommMpi(int narg, char* arg[]);
            virtual ~CommMpi();

            MPI_Comm GetComm();

        protected:
            virtual void v_Finalise();
            virtual int  v_GetRank();
            virtual void v_Block();
            virtual bool v_TreatAsRankZero(void);
            virtual void v_Send(int pProc, Array<OneD, NekDouble>& pData);
            virtual void v_Send(int pProc, Array<OneD, int>& pData);
            virtual void v_Send(int pProc, std::vector<unsigned int>& pData);
            virtual void v_Recv(int pProc, Array<OneD, NekDouble>& pData);
            virtual void v_Recv(int pProc, Array<OneD, int>& pData);
            virtual void v_Recv(int pProc, std::vector<unsigned int>& pData);
            virtual void v_SendRecv(int pSendProc,
                                    Array<OneD, NekDouble>& pSendData,
                                    int pRecvProc,
                                    Array<OneD, NekDouble>& pRecvData);
            virtual void v_SendRecv(int pSendProc,
                                    Array<OneD, int>& pSendData,
                                    int pRecvProc,
                                    Array<OneD, int>& pRecvData);
			virtual void v_SendRecvReplace(int pSendProc,
										   int pRecvProc,
										   Array<OneD, NekDouble>& pSendData);
			virtual void v_SendRecvReplace(int pSendProc,
										   int pRecvProc,
										   Array<OneD, int>& pSendData);
            virtual void v_AllReduce(NekDouble& pData,
                                     enum ReduceOperator pOp);
            virtual void v_AllReduce(int& pData,
                                     enum ReduceOperator pOp);
            virtual void v_AllReduce(Array<OneD, NekDouble>& pData,
                                     enum ReduceOperator pOp);
            virtual void v_AllReduce(Array<OneD, int      >& pData,
                                     enum ReduceOperator pOp);
            virtual void v_AllReduce(std::vector<unsigned int>& pData,
                                     enum ReduceOperator pOp);
			virtual void v_AlltoAll(Array<OneD, NekDouble>& pSendData,
									Array<OneD, NekDouble>& pRecvData);
            virtual void v_AlltoAll(Array<OneD, int>& pSendData,
									Array<OneD, int>& pRecvData);
			virtual void v_AlltoAllv(Array<OneD, NekDouble>& pSendData,
									Array<OneD, int>& pSendDataSizeMap,
									Array<OneD, int>& pSendDataOffsetMap,
									Array<OneD, NekDouble>& pRecvData,
									Array<OneD, int>& pRecvDataSizeMap,
									Array<OneD, int>& pRecvDataOffsetMap);
			virtual void v_AlltoAllv(Array<OneD, int>& pSendData,
									Array<OneD, int>& pSendDataSizeMap,
									Array<OneD, int>& pSendDataOffsetMap,
									Array<OneD, int>& pRecvData,
									Array<OneD, int>& pRecvDataSizeMap,
									Array<OneD, int>& pRecvDataOffsetMap);
            virtual void v_SplitComm(int pRows, int pColumns);

        private:
            MPI_Comm m_comm;
            int m_rank;

            CommMpi(MPI_Comm pComm);
        };
    }
}

#endif
