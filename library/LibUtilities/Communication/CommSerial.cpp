///////////////////////////////////////////////////////////////////////////////
//
// File CommSerial.cpp
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
// Description: Serial (= no) communication implementation
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Communication/CommSerial.h>

namespace Nektar
{
    namespace LibUtilities
    {
        std::string CommSerial::className
            = GetCommFactory().RegisterCreatorFunction(
                    "Serial",
                    CommSerial::create,
                    "Single-process serial communication.");

        CommSerial::CommSerial(int argc, char* argv[]) :
                Comm(argc, argv)
        {
            m_size = 1;
            m_type = "Serial";
        }

        CommSerial::~CommSerial()
        {

        }

        /**
         *
         */
        void CommSerial::v_Finalise()
        {

        }

        /**
         *
         */
        int CommSerial::v_GetRank()
        {
            return 0;
        }


        /**
         *
         */
        bool CommSerial::v_TreatAsRankZero(void)
        {
            return true;
        }


        /**
         *
         */
        void CommSerial::v_Block()
        {
        }

        /**
         *
         */
        double CommSerial::v_Wtime()
        {
	    return 0;
        }

        /**
         *
         */
        void CommSerial::v_Send(const void* buf, int count, CommDataType dt, int dest)
        {
        }


        /**
         *
         */
        void CommSerial::v_Recv(void* buf, int count, CommDataType dt, int source)
        {
        }

        /**
         *
         */
        void CommSerial::v_Sendrecv(const void *sendbuf, int sendcount, CommDataType sendtype, int dest,
                void *recvbuf, int recvcount, CommDataType recvtype, int source)
        {
        }

		/**
         *
         */
        void CommSerial::v_SendRecvReplace(void* buf, int count, CommDataType dt,
                int pSendProc, int pRecvProc)
		{
		}
		
		
        /**
         *
         */
        void CommSerial::v_AllReduce(void* buf, int count, CommDataType dt, enum ReduceOperator pOp)
        {

        }

        /**
         *
         */
        void CommSerial::v_AlltoAll(const void* sendbuf, int sendcount, CommDataType sendtype,
                                    void* recvbuf, int recvcount, CommDataType recvtype)
        {

        }
		
		/**
         *
         */
        void CommSerial::v_AlltoAllv(const void *sendbuf, const int sendcounts[], const int sensdispls[], CommDataType sendtype,
                void *recvbuf, const int recvcounts[], const int rdispls[], CommDataType recvtype)
        {

        }

		void CommSerial::v_Bcast(void* buffer, int count, CommDataType dt, int root)
		{

		}

        void CommSerial::v_Exscan(const Array<OneD, unsigned long long>& pData, const enum ReduceOperator pOp, Array<OneD, unsigned long long>& ans)
        {

        }

        void CommSerial::v_Gather(const void* sendbuf, int sendcount, CommDataType sendtype,
                void *recvbuf, int recvcount, CommDataType recvtype, int root)
        {
            std::memcpy(recvbuf, sendbuf, sendcount*CommDataTypeGetSize(sendtype));
        }

        void CommSerial::v_Scatter(const void *sendbuf, int sendcount, CommDataType sendtype,
                void *recvbuf, int recvcount, CommDataType recvtype, int root)
        {
            std::memcpy(recvbuf, sendbuf, sendcount*CommDataTypeGetSize(sendtype));
        }
        /**
         *
         */
        void CommSerial::v_SplitComm(int pRows, int pColumns)
        {
            ASSERTL0(false, "Cannot split a serial process.");
        }
        /**
         *
         */
        CommSharedPtr CommSerial::v_CommCreateIf(int flag)
        {
            ASSERTL0(false, "Cannot split a serial process.");
        }

    }
}
