///////////////////////////////////////////////////////////////////////////////
//
// File CommSerial.h
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
// Description: CommSerial header
//
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_LIB_UTILITIES_COMMSERIAL_H
#define NEKTAR_LIB_UTILITIES_COMMSERIAL_H

#include <string>

#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>

namespace Nektar
{
    namespace LibUtilities
    {
        // Forward declarations
        class CommSerial;

        /// Pointer to a Communicator object.
        typedef boost::shared_ptr<CommSerial> CommSerialSharedPtr;

        /// A global linear system.
        class CommSerial : public Comm
        {
        public:
            /// Creates an instance of this class
            LIB_UTILITIES_EXPORT static CommSharedPtr create(int narg, char* arg[])
            {
                return MemoryManager<CommSerial>::AllocateSharedPtr(narg,arg);
            }

            /// Name of class
            LIB_UTILITIES_EXPORT static std::string className;

            LIB_UTILITIES_EXPORT CommSerial(int argc, char* argv[]);
            LIB_UTILITIES_EXPORT virtual ~CommSerial();

        protected:
            LIB_UTILITIES_EXPORT virtual void v_Finalise();
            LIB_UTILITIES_EXPORT virtual int  v_GetRank();
            LIB_UTILITIES_EXPORT virtual bool v_TreatAsRankZero(void);

            LIB_UTILITIES_EXPORT virtual void v_Block();
            LIB_UTILITIES_EXPORT virtual void v_Send(int pProc, Array<OneD, NekDouble>& pData);
            LIB_UTILITIES_EXPORT virtual void v_Send(int pProc, Array<OneD, int>& pData);
            LIB_UTILITIES_EXPORT virtual void v_Send(int pProc, std::vector<unsigned int>& pData);
            LIB_UTILITIES_EXPORT virtual void v_Recv(int pProc, Array<OneD, NekDouble>& pData);
            LIB_UTILITIES_EXPORT virtual void v_Recv(int pProc, Array<OneD, int>& pData);
            LIB_UTILITIES_EXPORT virtual void v_Recv(int pProc, std::vector<unsigned int>& pData);
            LIB_UTILITIES_EXPORT virtual void v_SendRecv(int pSendProc,
                                    Array<OneD, NekDouble>& pSendData,
                                    int pRecvProc,
                                    Array<OneD, NekDouble>& pRecvData);
            LIB_UTILITIES_EXPORT virtual void v_SendRecv(int pSendProc,
                                    Array<OneD, int>& pSendData,
                                    int pRecvProc,
                                    Array<OneD, int>& pRecvData);
            LIB_UTILITIES_EXPORT virtual void v_SendRecvReplace(int pSendProc,
                                                                int pRecvProc,
                                                                Array<OneD, NekDouble>& pSendData);
            LIB_UTILITIES_EXPORT virtual void v_SendRecvReplace(int pSendProc,
                                                                int pRecvProc,
                                                                Array<OneD, int>& pSendData);
            LIB_UTILITIES_EXPORT virtual void v_AllReduce(NekDouble& pData,
                                                          enum ReduceOperator pOp);
            LIB_UTILITIES_EXPORT virtual void v_AllReduce(int& pData,
                                                          enum ReduceOperator pOp);
            LIB_UTILITIES_EXPORT virtual void v_AllReduce(Array<OneD, NekDouble>& pData,
                                                          enum ReduceOperator pOp);
            LIB_UTILITIES_EXPORT virtual void v_AllReduce(Array<OneD, int      >& pData,
                                                          enum ReduceOperator pOp);
            LIB_UTILITIES_EXPORT virtual void v_AllReduce(std::vector<unsigned int>& pData,
                                                          enum ReduceOperator pOp);
            LIB_UTILITIES_EXPORT virtual void v_AlltoAll(Array<OneD, NekDouble>& pSendData,
                                                         Array<OneD, NekDouble>& pRecvData);
            LIB_UTILITIES_EXPORT virtual void v_AlltoAll(Array<OneD, int>& pSendData,
                                                         Array<OneD, int>& pRecvData);
            LIB_UTILITIES_EXPORT virtual void v_AlltoAllv(Array<OneD, NekDouble>& pSendData,
                                                          Array<OneD, int>& pSendDataSizeMap,
                                                          Array<OneD, int>& pSendDataOffsetMap,
                                                          Array<OneD, NekDouble>& pRecvData,
                                                          Array<OneD, int>& pRecvDataSizeMap,
                                                          Array<OneD, int>& pRecvDataOffsetMap);
            LIB_UTILITIES_EXPORT virtual void v_AlltoAllv(Array<OneD, int>& pSendData,
                                                          Array<OneD, int>& pSendDataSizeMap,
                                                          Array<OneD, int>& pSendDataOffsetMap,
                                                          Array<OneD, int>& pRecvData,
                                                          Array<OneD, int>& pRecvDataSizeMap,
                                                          Array<OneD, int>& pRecvDataOffsetMap);
            LIB_UTILITIES_EXPORT virtual void v_SplitComm(int pRows, int pColumns);
            
        };

    }
}

#endif
