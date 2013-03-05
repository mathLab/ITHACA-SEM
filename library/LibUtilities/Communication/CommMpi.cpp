///////////////////////////////////////////////////////////////////////////////
//
// File CommMpi.cpp
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
// Description: MPI communication implementation
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Communication/CommMpi.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

namespace Nektar
{
    namespace LibUtilities
    {
        std::string CommMpi::className
            = GetCommFactory().RegisterCreatorFunction(
                "ParallelMPI",
                CommMpi::create,
                "Parallel communication using MPI.");

        /**
         *
         */
        CommMpi::CommMpi(int narg, char* arg[])
                : Comm(narg,arg)
        {
            int retval = MPI_Init(&narg, &arg);
            if (retval != MPI_SUCCESS)
            {
                ASSERTL0(false, "Failed to initialise MPI");
            }
            m_comm = MPI_COMM_WORLD;
            MPI_Comm_size( m_comm, &m_size );
            MPI_Comm_rank( m_comm, &m_rank );

            m_type = "Parallel MPI";
        }


        /**
         *
         */
        CommMpi::CommMpi(MPI_Comm pComm)
                : Comm()
        {
            m_comm = pComm;
            MPI_Comm_size( m_comm, &m_size );
            MPI_Comm_rank( m_comm, &m_rank );

            m_type = "Parallel MPI";
        }


        /**
         *
         */
        CommMpi::~CommMpi()
        {

        }


        /**
         *
         */
        MPI_Comm CommMpi::GetComm()
        {
            return m_comm;
        }


        /**
         *
         */
        void CommMpi::v_Finalise()
        {
            MPI_Finalize();
        }


        /**
         *
         */
        int CommMpi::v_GetRank()
        {
            return m_rank;
        }


        /**
         *
         */
        void CommMpi::v_Block()
        {
            MPI_Barrier(m_comm);
        }


        /**
         *
         */
        void CommMpi::v_Send(int pProc, Array<OneD, NekDouble>& pData)
        {
            if (MPISYNC)
            {
                MPI_Ssend( pData.get(),
                          (int) pData.num_elements(),
                          MPI_DOUBLE,
                          pProc,
                          0,
                          m_comm);
            }
            else
            {
                MPI_Send( pData.get(),
                          (int) pData.num_elements(),
                          MPI_DOUBLE,
                          pProc,
                          0,
                          m_comm);
            }
        }


        /**
         *
         */
        void CommMpi::v_Recv(int pProc, Array<OneD, NekDouble>& pData)
        {
            MPI_Status status;
            MPI_Recv( pData.get(),
                      (int) pData.num_elements(),
                      MPI_DOUBLE,
                      pProc,
                      0,
                      m_comm,
                      &status);

            //ASSERTL0(status.MPI_ERROR == MPI_SUCCESS,
            //         "MPI error receiving data.");
        }


        /**
         *
         */
        void CommMpi::v_Send(int pProc, Array<OneD, int>& pData)
        {
            if (MPISYNC)
            {
                MPI_Ssend( pData.get(),
                          (int) pData.num_elements(),
                          MPI_INT,
                          pProc,
                          0,
                          m_comm);
            }
            else
            {
                MPI_Send( pData.get(),
                          (int) pData.num_elements(),
                          MPI_INT,
                          pProc,
                          0,
                          m_comm);
            }
        }


        /**
         *
         */
        void CommMpi::v_Recv(int pProc, Array<OneD, int>& pData)
        {
            MPI_Status status;
            MPI_Recv( pData.get(),
                      (int) pData.num_elements(),
                      MPI_INT,
                      pProc,
                      0,
                      m_comm,
                      &status);

            //ASSERTL0(status.MPI_ERROR == MPI_SUCCESS,
            //         "MPI error receiving data.");
        }


        /**
         *
         */
        void CommMpi::v_SendRecv(int pSendProc,
                                Array<OneD, NekDouble>& pSendData,
                                int pRecvProc,
                                Array<OneD, NekDouble>& pRecvData)
        {
            MPI_Status status;
            int retval = MPI_Sendrecv(pSendData.get(),
                         (int) pSendData.num_elements(),
                         MPI_DOUBLE,
                         pRecvProc,
                         0,
                         pRecvData.get(),
                         (int) pRecvData.num_elements(),
                         MPI_DOUBLE,
                         pSendProc,
                         0,
                         m_comm,
                         &status);

            ASSERTL0(retval == MPI_SUCCESS,
                     "MPI error performing send-receive of data.");
        }


        /**
         *
         */
        void CommMpi::v_SendRecv(int pSendProc,
                                Array<OneD, int>& pSendData,
                                int pRecvProc,
                                Array<OneD, int>& pRecvData)
        {
            MPI_Status status;
            int retval = MPI_Sendrecv(pSendData.get(),
                         (int) pSendData.num_elements(),
                         MPI_INT,
                         pRecvProc,
                         0,
                         pRecvData.get(),
                         (int) pRecvData.num_elements(),
                         MPI_INT,
                         pSendProc,
                         0,
                         m_comm,
                         &status);

            ASSERTL0(retval == MPI_SUCCESS,
                     "MPI error performing send-receive of data.");
        }
		
		/**
         *
         */
        void CommMpi::v_SendRecvReplace(int pSendProc,
								        int pRecvProc,
										Array<OneD, NekDouble>& pSendData)
        {
            MPI_Status status;
            int retval = MPI_Sendrecv_replace(pSendData.get(),
									         (int) pSendData.num_elements(),
									          MPI_DOUBLE,
									          pRecvProc,
									          0,
									          pSendProc,
									          0,
									          m_comm,
									          &status);
			
            ASSERTL0(retval == MPI_SUCCESS,
                     "MPI error performing Send-Receive-Replace of data.");
        }
		
		
        /**
         *
         */
        void CommMpi::v_SendRecvReplace(int pSendProc,
										int pRecvProc,
								         Array<OneD, int>& pSendData)
							
        {
            MPI_Status status;
            int retval = MPI_Sendrecv_replace(pSendData.get(),
											  (int) pSendData.num_elements(),
									          MPI_INT,
									          pRecvProc,
									          0,
									          pSendProc,
									          0,
									          m_comm,
									          &status);
			
            ASSERTL0(retval == MPI_SUCCESS,
                     "MPI error performing Send-Receive-Replace of data.");
        }


        /**
         *
         */
        void CommMpi::v_AllReduce(NekDouble& pData, enum ReduceOperator pOp)
        {
            if (GetSize() == 1)
            {
                return;
            }

            MPI_Op vOp;
            switch (pOp)
            {
            case ReduceMax: vOp = MPI_MAX; break;
            case ReduceMin: vOp = MPI_MIN; break;
            case ReduceSum:
            default:        vOp = MPI_SUM; break;
            }
            int retval = MPI_Allreduce( MPI_IN_PLACE,
                                        &pData,
                                        1,
                                        MPI_DOUBLE,
                                        vOp,
                                        m_comm);

            ASSERTL0(retval == MPI_SUCCESS,
                     "MPI error performing All-reduce.");
        }


        /**
         *
         */
        void CommMpi::v_AllReduce(int& pData, enum ReduceOperator pOp)
        {
            if (GetSize() == 1)
            {
                return;
            }

            MPI_Op vOp;
            switch (pOp)
            {
            case ReduceMax: vOp = MPI_MAX; break;
            case ReduceMin: vOp = MPI_MIN; break;
            case ReduceSum:
            default:        vOp = MPI_SUM; break;
            }
            int retval = MPI_Allreduce( MPI_IN_PLACE,
                                        &pData,
                                        1,
                                        MPI_INT,
                                        vOp,
                                        m_comm);

            ASSERTL0(retval == MPI_SUCCESS,
                     "MPI error performing All-reduce.");
        }


        /**
         *
         */
        void CommMpi::v_AllReduce(Array<OneD, NekDouble>& pData, enum ReduceOperator pOp)
        {
            if (GetSize() == 1)
            {
                return;
            }

            MPI_Op vOp;
            switch (pOp)
            {
            case ReduceMax: vOp = MPI_MAX; break;
            case ReduceMin: vOp = MPI_MIN; break;
            case ReduceSum:
            default:        vOp = MPI_SUM; break;
            }
            int retval = MPI_Allreduce( MPI_IN_PLACE,
                                        pData.get(),
                                        (int) pData.num_elements(),
                                        MPI_DOUBLE,
                                        vOp,
                                        m_comm);

            ASSERTL0(retval == MPI_SUCCESS,
                     "MPI error performing All-reduce.");
        }


        /**
         *
         */
        void CommMpi::v_AllReduce(Array<OneD, int>& pData, enum ReduceOperator pOp)
        {
            if (GetSize() == 1)
            {
                return;
            }

            MPI_Op vOp;
            switch (pOp)
            {
            case ReduceMax: vOp = MPI_MAX; break;
            case ReduceMin: vOp = MPI_MIN; break;
            case ReduceSum:
            default:        vOp = MPI_SUM; break;
            }
            int retval = MPI_Allreduce( MPI_IN_PLACE,
                                        pData.get(),
                                        (int) pData.num_elements(),
                                        MPI_INT,
                                        vOp,
                                        m_comm);

            ASSERTL0(retval == MPI_SUCCESS,
                     "MPI error performing All-reduce.");
        }
		
		
		/**
         *
         */
		void CommMpi::v_AlltoAll(Array<OneD, NekDouble>& pSendData,Array<OneD, NekDouble>& pRecvData)
		{
			int retval = MPI_Alltoall(pSendData.get(),
									  (int) pSendData.num_elements()/GetSize(),
									  MPI_DOUBLE,
									  pRecvData.get(),
									  (int) pRecvData.num_elements()/GetSize(),
									  MPI_DOUBLE,
									  m_comm);
									  
			ASSERTL0(retval == MPI_SUCCESS,
                     "MPI error performing All-to-All.");
		}
		
		
		/**
         *
         */
		void CommMpi::v_AlltoAll(Array<OneD, int>& pSendData,Array<OneD, int>& pRecvData)
		{
			int retval = MPI_Alltoall(pSendData.get(),
									  (int) pSendData.num_elements()/GetSize(),
									  MPI_INT,
									  pRecvData.get(),
									  (int) pRecvData.num_elements()/GetSize(),
									  MPI_INT,
									  m_comm);
			
			ASSERTL0(retval == MPI_SUCCESS,
                     "MPI error performing All-to-All.");
		}
		
		
		/**
         *
         */
		void CommMpi::v_AlltoAllv(Array<OneD, NekDouble>& pSendData,
								 Array<OneD, int>& pSendDataSizeMap,
								 Array<OneD, int>& pSendDataOffsetMap,
								 Array<OneD, NekDouble>& pRecvData,
								 Array<OneD, int>& pRecvDataSizeMap,
								 Array<OneD, int>& pRecvDataOffsetMap)
		{
			int retval = MPI_Alltoallv(pSendData.get(),
									   pSendDataSizeMap.get(),
									   pSendDataOffsetMap.get(),
									   MPI_DOUBLE,
									   pRecvData.get(),
									   pRecvDataSizeMap.get(),
									   pRecvDataOffsetMap.get(),
									   MPI_DOUBLE,
									   m_comm);
									   
		    ASSERTL0(retval == MPI_SUCCESS,
					 "MPI error performing All-to-All-v.");
		}
		
		/**
         *
         */
		void CommMpi::v_AlltoAllv(Array<OneD, int>& pSendData,
								  Array<OneD, int>& pSendDataSizeMap,
								  Array<OneD, int>& pSendDataOffsetMap,
								  Array<OneD, int>& pRecvData,
								  Array<OneD, int>& pRecvDataSizeMap,
								  Array<OneD, int>& pRecvDataOffsetMap)
		{
			int retval = MPI_Alltoallv(pSendData.get(),
									   pSendDataSizeMap.get(),
									   pSendDataOffsetMap.get(),
									   MPI_INT,
									   pRecvData.get(),
									   pRecvDataSizeMap.get(),
									   pRecvDataOffsetMap.get(),
									   MPI_INT,
									   m_comm);
									   
			ASSERTL0(retval == MPI_SUCCESS,
					 "MPI error performing All-to-All-v.");
		}


        /**
         * Processes are considered as a grid of size pRows*pColumns. Comm
         * objects are created corresponding to the rows and columns of this
         * grid. The row and column to which this process belongs is stored in
         * #m_commRow and #m_commColumn.
         */
        void CommMpi::v_SplitComm(int pRows, int pColumns)
        {
            ASSERTL0(pRows*pColumns == m_size,
                    "Rows/Columns do not match comm size.");

            MPI_Comm newComm;

            // Compute row and column in grid.
            int myCol = m_rank % pColumns;
            int myRow = (m_rank - myCol) / pColumns;

            // Split Comm into rows - all processes with same myRow are put in
            // the same communicator. The rank within this communicator is the
            // column index.
            MPI_Comm_split(m_comm, myRow, myCol, &newComm);
            m_commRow = boost::shared_ptr<Comm>(new CommMpi(newComm));

            // Split Comm into columns - all processes with same myCol are put
            // in the same communicator. The rank within this communicator is
            // the row index.
            MPI_Comm_split(m_comm, myCol, myRow, &newComm);
            m_commColumn = boost::shared_ptr<Comm>(new CommMpi(newComm));
        }
    }
}
