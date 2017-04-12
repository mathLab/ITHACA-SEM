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

#ifdef NEKTAR_USING_PETSC
#include "petscsys.h"
#endif

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Communication/CommMpi.h>

namespace Nektar
{
namespace LibUtilities
{
std::string CommMpi::className = GetCommFactory().RegisterCreatorFunction(
    "ParallelMPI", CommMpi::create, "Parallel communication using MPI.");

/**
 *
 */
CommMpi::CommMpi(int narg, char *arg[]) : Comm(narg, arg)
{
    int init = 0;
    MPI_Initialized(&init);
    ASSERTL0(!init, "MPI has already been initialised.");

    int retval = MPI_Init(&narg, &arg);
    if (retval != MPI_SUCCESS)
    {
        ASSERTL0(false, "Failed to initialise MPI");
    }

    m_comm = MPI_COMM_WORLD;
    MPI_Comm_size(m_comm, &m_size);
    MPI_Comm_rank(m_comm, &m_rank);

#ifdef NEKTAR_USING_PETSC
    PetscInitializeNoArguments();
#endif

    m_type = "Parallel MPI";
}

/**
 *
 */
CommMpi::CommMpi(MPI_Comm pComm) : Comm()
{
    m_comm = pComm;
    MPI_Comm_size(m_comm, &m_size);
    MPI_Comm_rank(m_comm, &m_rank);

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
#ifdef NEKTAR_USING_PETSC
    PetscFinalize();
#endif
    int flag;
    MPI_Finalized(&flag);
    if (!flag)
    {
        MPI_Finalize();
    }
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
bool CommMpi::v_TreatAsRankZero(void)
{
    if (m_rank == 0)
    {
        return true;
    }
    else
    {
        return false;
    }
    return true;
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
double CommMpi::v_Wtime()
{
    return MPI_Wtime();
}

/**
 *
 */
void CommMpi::v_Send(void *buf, int count, CommDataType dt, int dest)
{
    if (MPISYNC)
    {
        MPI_Ssend(buf, count, dt, dest, 0, m_comm);
    }
    else
    {
        MPI_Send(buf, count, dt, dest, 0, m_comm);
    }
}

/**
 *
 */
void CommMpi::v_Recv(void *buf, int count, CommDataType dt, int source)
{
    MPI_Recv(buf, count, dt, source, 0, m_comm, MPI_STATUS_IGNORE);
    // ASSERTL0(status.MPI_ERROR == MPI_SUCCESS,
    //         "MPI error receiving data.");
}

/**
 *
 */
void CommMpi::v_SendRecv(void *sendbuf, int sendcount, CommDataType sendtype,
                         int dest, void *recvbuf, int recvcount,
                         CommDataType recvtype, int source)
{
    MPI_Status status;
    int retval = MPI_Sendrecv(sendbuf, sendcount, sendtype, dest, 0, recvbuf,
                              recvcount, recvtype, source, 0, m_comm, &status);

    ASSERTL0(retval == MPI_SUCCESS,
             "MPI error performing send-receive of data.");
}

/**
*
*/
void CommMpi::v_SendRecvReplace(void *buf, int count, CommDataType dt,
                                int pSendProc, int pRecvProc)
{
    MPI_Status status;
    int retval = MPI_Sendrecv_replace(buf, count, dt, pRecvProc, 0, pSendProc,
                                      0, m_comm, &status);

    ASSERTL0(retval == MPI_SUCCESS,
             "MPI error performing Send-Receive-Replace of data.");
}

/**
 *
 */
void CommMpi::v_AllReduce(void *buf, int count, CommDataType dt,
                          enum ReduceOperator pOp)
{
    if (GetSize() == 1)
    {
        return;
    }

    MPI_Op vOp;
    switch (pOp)
    {
        case ReduceMax:
            vOp = MPI_MAX;
            break;
        case ReduceMin:
            vOp = MPI_MIN;
            break;
        case ReduceSum:
        default:
            vOp = MPI_SUM;
            break;
    }
    int retval = MPI_Allreduce(MPI_IN_PLACE, buf, count, dt, vOp, m_comm);

    ASSERTL0(retval == MPI_SUCCESS, "MPI error performing All-reduce.");
}

/**
 *
 */
void CommMpi::v_AlltoAll(void *sendbuf, int sendcount, CommDataType sendtype,
                         void *recvbuf, int recvcount, CommDataType recvtype)
{
    int retval = MPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                              recvtype, m_comm);

    ASSERTL0(retval == MPI_SUCCESS, "MPI error performing All-to-All.");
}

/**
 *
 */
void CommMpi::v_AlltoAllv(void *sendbuf, int sendcounts[], int sdispls[],
                          CommDataType sendtype, void *recvbuf,
                          int recvcounts[], int rdispls[],
                          CommDataType recvtype)
{
    int retval = MPI_Alltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf,
                               recvcounts, rdispls, recvtype, m_comm);

    ASSERTL0(retval == MPI_SUCCESS, "MPI error performing All-to-All-v.");
}

void CommMpi::v_Bcast(void *buffer, int count, CommDataType dt, int root)
{
    int retval = MPI_Bcast(buffer, count, dt, root, m_comm);
    ASSERTL0(retval == MPI_SUCCESS, "MPI error performing Bcast-v.");
}

void CommMpi::v_Exscan(Array<OneD, unsigned long long> &pData,
                       const enum ReduceOperator pOp,
                       Array<OneD, unsigned long long> &ans)
{
    int n = pData.num_elements();
    ASSERTL0(n == ans.num_elements(), "Array sizes differ in Exscan");

    MPI_Op vOp;
    switch (pOp)
    {
        case ReduceMax:
            vOp = MPI_MAX;
            break;
        case ReduceMin:
            vOp = MPI_MIN;
            break;
        case ReduceSum:
        default:
            vOp = MPI_SUM;
            break;
    }

    int retval = MPI_Exscan(pData.get(), ans.get(), n, MPI_UNSIGNED_LONG_LONG,
                            vOp, m_comm);
    ASSERTL0(retval == MPI_SUCCESS, "MPI error performing Exscan-v.");
}

void CommMpi::v_Gather(void *sendbuf, int sendcount, CommDataType sendtype,
                       void *recvbuf, int recvcount, CommDataType recvtype,
                       int root)
{
    int retval = MPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                            recvtype, root, m_comm);

    ASSERTL0(retval == MPI_SUCCESS, "MPI error performing Gather.");
}

void CommMpi::v_Scatter(void *sendbuf, int sendcount, CommDataType sendtype,
                        void *recvbuf, int recvcount, CommDataType recvtype,
                        int root)
{
    int retval = MPI_Scatter(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                             recvtype, root, m_comm);
    ASSERTL0(retval == MPI_SUCCESS, "MPI error performing Scatter.");
}

/**
 * Processes are considered as a grid of size pRows*pColumns. Comm
 * objects are created corresponding to the rows and columns of this
 * grid. The row and column to which this process belongs is stored in
 * #m_commRow and #m_commColumn.
 */
void CommMpi::v_SplitComm(int pRows, int pColumns)
{
    ASSERTL0(pRows * pColumns == m_size,
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

/**
 * Create a new communicator if the flag is non-zero.
 */
CommSharedPtr CommMpi::v_CommCreateIf(int flag)
{
    MPI_Comm newComm;
    // color == MPI_UNDEF => not in the new communicator
    // key == 0 on all => use rank to order them. OpenMPI, at least,
    // implies this is faster than ordering them ourselves.
    MPI_Comm_split(m_comm, flag ? 0 : MPI_UNDEFINED, 0, &newComm);

    if (flag == 0)
    {
        // flag == 0 => get back MPI_COMM_NULL, return a null ptr instead.
        return boost::shared_ptr<Comm>();
    }
    else
    {
        // Return a real communicator
        return boost::shared_ptr<Comm>(new CommMpi(newComm));
    }
}
}
}
