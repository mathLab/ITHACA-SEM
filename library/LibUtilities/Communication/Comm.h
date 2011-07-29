///////////////////////////////////////////////////////////////////////////////
//
// File GlobalLinSys.h
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
// Description: GlobalLinSys header
//
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_LIB_UTILITIES_COMM_H
#define NEKTAR_LIB_UTILITIES_COMM_H

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>

namespace Nektar
{
    namespace LibUtilities
    {
        // Forward declarations
        class Comm;

        /// Pointer to a Communicator object.
        typedef boost::shared_ptr<Comm> CommSharedPtr;

        /// Datatype of the NekFactory used to instantiate classes derived from
        /// the EquationSystem class.
        typedef LibUtilities::NekFactory< std::string, Comm, int, char** > CommFactory;

        LIB_UTILITIES_EXPORT CommFactory& GetCommFactory();

		/// Type of operation to perform in AllReduce.
        enum ReduceOperator
        {
            ReduceSum,
            ReduceMax,
            ReduceMin
        };

        /// Base communications class
        class Comm
        {
            public:
                LIB_UTILITIES_EXPORT Comm(int narg, char* arg[]);
                LIB_UTILITIES_EXPORT virtual ~Comm();

                LIB_UTILITIES_EXPORT inline void Finalise();

                /// Returns number of processes
                LIB_UTILITIES_EXPORT inline int GetSize();
                LIB_UTILITIES_EXPORT inline int GetRank();

                /// Block execution until all processes reach this point
                LIB_UTILITIES_EXPORT inline void Block();
                LIB_UTILITIES_EXPORT inline void Send(int pProc, Array<OneD, NekDouble>& pData);
                LIB_UTILITIES_EXPORT inline void Send(int pProc, Array<OneD, int>& pData);
                LIB_UTILITIES_EXPORT inline void Recv(int pProc, Array<OneD, NekDouble>& pData);
                LIB_UTILITIES_EXPORT inline void Recv(int pProc, Array<OneD, int>& pData);
                LIB_UTILITIES_EXPORT inline void SendRecv(int pSendProc,
                                     Array<OneD, NekDouble>& pSendData,
                                     int pRecvProc,
                                     Array<OneD, NekDouble>& pRecvData);
                LIB_UTILITIES_EXPORT inline void SendRecv(int pSendProc,
                                     Array<OneD, int>& pSendData,
                                     int pRecvProc,
                                     Array<OneD, int>& pRecvData);
                LIB_UTILITIES_EXPORT inline void AllReduce(NekDouble& pData, enum ReduceOperator pOp);
                LIB_UTILITIES_EXPORT inline void AllReduce(int& pData, enum ReduceOperator pOp);
                LIB_UTILITIES_EXPORT inline void AllReduce(Array<OneD, NekDouble>& pData,
                                         enum ReduceOperator pOp);
                LIB_UTILITIES_EXPORT inline void AllReduce(Array<OneD, int      >& pData,
                                         enum ReduceOperator pOp);

            protected:
                int m_size;     ///< Number of processes

                virtual void v_Finalise() = 0;
                virtual int  v_GetRank() = 0;
                virtual void v_Block() = 0;
                virtual void v_Send(int pProc, Array<OneD, NekDouble>& pData) = 0;
                virtual void v_Send(int pProc, Array<OneD, int>& pData) = 0;
                virtual void v_Recv(int pProc, Array<OneD, NekDouble>& pData) = 0;
                virtual void v_Recv(int pProc, Array<OneD, int>& pData) = 0;
                virtual void v_SendRecv(int pSendProc,
                                        Array<OneD, NekDouble>& pSendData,
                                        int pRecvProc,
                                        Array<OneD, NekDouble>& pRecvData) = 0;
                virtual void v_SendRecv(int pSendProc,
                                        Array<OneD, int>& pSendData,
                                        int pRecvProc,
                                        Array<OneD, int>& pRecvData) = 0;
                virtual void v_AllReduce(NekDouble& pData,
                                         enum ReduceOperator pOp) = 0;
                virtual void v_AllReduce(int& pData,
                                         enum ReduceOperator pOp) = 0;
                virtual void v_AllReduce(Array<OneD, NekDouble>& pData,
                                         enum ReduceOperator pOp) = 0;
                virtual void v_AllReduce(Array<OneD, int      >& pData,
                                         enum ReduceOperator pOp) = 0;
        };


        /**
         *
         */
        inline void Comm::Finalise()
        {
            v_Finalise();
        }

        /**
         *
         */
        inline int Comm::GetSize()
        {
            return m_size;
        }

        /**
         *
         */
        inline int Comm::GetRank()
        {
            return v_GetRank();
        }

        /**
         *
         */
        inline void Comm::Block()
        {
            v_Block();
        }

        /**
         *
         */
        inline void Comm::Send(int pProc, Array<OneD, NekDouble>& pData)
        {
            v_Send(pProc, pData);
        }

        /**
         *
         */
        inline void Comm::Recv(int pProc, Array<OneD, NekDouble>& pData)
        {
            v_Recv(pProc, pData);
        }


        /**
         *
         */
        inline void Comm::Send(int pProc, Array<OneD, int>& pData)
        {
            v_Send(pProc, pData);
        }

        /**
         *
         */
        inline void Comm::Recv(int pProc, Array<OneD, int>& pData)
        {
            v_Recv(pProc, pData);
        }

        /**
         *
         */
        inline void Comm::SendRecv(int pSendProc,
                             Array<OneD, NekDouble>& pSendData,
                             int pRecvProc,
                             Array<OneD, NekDouble>& pRecvData)
        {
            v_SendRecv(pSendProc, pSendData, pRecvProc, pRecvData);
        }


        /**
         *
         */
        inline void Comm::SendRecv(int pSendProc,
                             Array<OneD, int>& pSendData,
                             int pRecvProc,
                             Array<OneD, int>& pRecvData)
        {
            v_SendRecv(pSendProc, pSendData, pRecvProc, pRecvData);
        }


        /**
         *
         */
        inline void Comm::AllReduce(NekDouble& pData, enum ReduceOperator pOp)
        {
            v_AllReduce(pData, pOp);
        }


        /**
         *
         */
        inline void Comm::AllReduce(int& pData, enum ReduceOperator pOp)
        {
            v_AllReduce(pData, pOp);
        }


        /**
         *
         */
        inline void Comm::AllReduce(Array<OneD, NekDouble>& pData, enum ReduceOperator pOp)
        {
            v_AllReduce(pData, pOp);
        }


        /**
         *
         */
        inline void Comm::AllReduce(Array<OneD, int>& pData, enum ReduceOperator pOp)
        {
            v_AllReduce(pData, pOp);
        }

    }
}

#endif
