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
        void CommSerial::v_Send(int pProc, Array<OneD, NekDouble>& pData)
        {
        }


        /**
         *
         */
        void CommSerial::v_Recv(int pProc, Array<OneD, NekDouble>& pData)
        {
        }


        /**
         *
         */
        void CommSerial::v_Send(int pProc, Array<OneD, int>& pData)
        {
        }


        /**
         *
         */
        void CommSerial::v_Recv(int pProc, Array<OneD, int>& pData)
        {
        }


        /**
         *
         */
        void CommSerial::v_Send(int pProc, std::vector<unsigned int>& pData)
        {
        }


        /**
         *
         */
        void CommSerial::v_Recv(int pProc, std::vector<unsigned int>& pData)
        {
        }


        /**
         *
         */
        void CommSerial::v_SendRecv(int pSendProc,
                                Array<OneD, NekDouble>& pSendData,
                                int pRecvProc,
                                Array<OneD, NekDouble>& pRecvData)
        {
        }


        /**
         *
         */
        void CommSerial::v_SendRecv(int pSendProc,
                                Array<OneD, int>& pSendData,
                                int pRecvProc,
                                Array<OneD, int>& pRecvData)
        {
        }
		
		
		/**
         *
         */
        void CommSerial::v_SendRecvReplace(int pSendProc,
										   int pRecvProc,
										   Array<OneD, NekDouble>& pSendData)
		{
		}
		
		
		/**
         *
         */
        void CommSerial::v_SendRecvReplace(int pSendProc,
										   int pRecvProc,
										   Array<OneD, int>& pSendData)
		{
		}


        /**
         *
         */
        void CommSerial::v_AllReduce(NekDouble& pData, enum ReduceOperator pOp)
        {

        }


        /**
         *
         */
        void CommSerial::v_AllReduce(int& pData, enum ReduceOperator pOp)
        {

        }


        /**
         *
         */
        void CommSerial::v_AllReduce(Array<OneD, NekDouble>& pData, enum ReduceOperator pOp)
        {

        }


        /**
         *
         */
        void CommSerial::v_AllReduce(Array<OneD, int      >& pData, enum ReduceOperator pOp)
        {

        }
		
		
        /**
         *
         */
        void CommSerial::v_AllReduce(std::vector<unsigned int>& pData, enum ReduceOperator pOp)
        {

        }


		/**
         *
         */
		void CommSerial::v_AlltoAll(Array<OneD, NekDouble>& pSendData,Array<OneD, NekDouble>& pRecvData)
		{
			
        }
		
		
		/**
         *
         */
		void CommSerial::v_AlltoAll(Array<OneD, int>& pSendData,Array<OneD, int>& pRecvData)
		{
			
        }
		
		
		/**
         *
         */
		void CommSerial::v_AlltoAllv(Array<OneD, NekDouble>& pSendData,
									 Array<OneD, int>& pSendDataSizeMap,
									 Array<OneD, int>& pSendDataOffsetMap,
									 Array<OneD, NekDouble>& pRecvData,
									 Array<OneD, int>& pRecvDataSizeMap,
									 Array<OneD, int>& pRecvDataOffsetMap)
		{
			
        }
		
		
		/**
         *
         */
		void CommSerial::v_AlltoAllv(Array<OneD, int>& pSendData,
									 Array<OneD, int>& pSendDataSizeMap,
									 Array<OneD, int>& pSendDataOffsetMap,
									 Array<OneD, int>& pRecvData,
									 Array<OneD, int>& pRecvDataSizeMap,
									 Array<OneD, int>& pRecvDataOffsetMap)
		{
			
        }


        /**
         *
         */
        void CommSerial::v_SplitComm(int pRows, int pColumns)
        {
            ASSERTL0(false, "Cannot split a serial process.");
        }
    }
}
