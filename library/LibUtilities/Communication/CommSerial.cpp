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

#include <boost/core/ignore_unused.hpp>

#ifdef NEKTAR_USING_PETSC
#include "petscsys.h"
#endif

#include <LibUtilities/Communication/CommSerial.h>

namespace Nektar
{
namespace LibUtilities
{
std::string CommSerial::className = GetCommFactory().RegisterCreatorFunction(
    "Serial", CommSerial::create, "Single-process serial communication.");

CommSerial::CommSerial(int argc, char *argv[]) : Comm(argc, argv)
{
#ifdef NEKTAR_USING_PETSC
    PetscInitializeNoArguments();
#endif
    m_size = 1;
    m_type = "Serial";
}

CommSerial::~CommSerial() = default;

/**
 *
 */
void CommSerial::v_Finalise()
{
#ifdef NEKTAR_USING_PETSC
    PetscFinalize();
#endif
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
bool CommSerial::v_TreatAsRankZero()
{
    return true;
}

/**
 *
 */
bool CommSerial::v_IsSerial()
{
    return true;
}

std::tuple<int, int, int> CommSerial::v_GetVersion()
{
    return std::make_tuple(0, 0, 0);
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
NekDouble CommSerial::v_Wtime()
{
    return 0;
}

/**
 *
 */
void CommSerial::v_Send(void *buf, int count, CommDataType dt, int dest)
{
    boost::ignore_unused(buf, count, dt, dest);
}

/**
 *
 */
void CommSerial::v_Recv(void *buf, int count, CommDataType dt, int source)
{
    boost::ignore_unused(buf, count, dt, source);
}

/**
 *
 */
void CommSerial::v_SendRecv(void *sendbuf, int sendcount, CommDataType sendtype,
                            int dest, void *recvbuf, int recvcount,
                            CommDataType recvtype, int source)
{
    boost::ignore_unused(sendbuf, sendcount, sendtype, dest, recvbuf, recvcount,
                         recvtype, source);
}

/**
 *
 */
void CommSerial::v_SendRecvReplace(void *buf, int count, CommDataType dt,
                                   int pSendProc, int pRecvProc)
{
    boost::ignore_unused(buf, count, dt, pSendProc, pRecvProc);
}

/**
 *
 */
void CommSerial::v_AllReduce(void *buf, int count, CommDataType dt,
                             enum ReduceOperator pOp)
{
    boost::ignore_unused(buf, count, dt, pOp);
}

/**
 *
 */
void CommSerial::v_AlltoAll(void *sendbuf, int sendcount, CommDataType sendtype,
                            void *recvbuf, int recvcount, CommDataType recvtype)
{
    boost::ignore_unused(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                         recvtype);
}

/**
 *
 */
void CommSerial::v_AlltoAllv(void *sendbuf, int sendcounts[], int sdispls[],
                             CommDataType sendtype, void *recvbuf,
                             int recvcounts[], int rdispls[],
                             CommDataType recvtype)
{
    boost::ignore_unused(sendbuf, sendcounts, sdispls, sendtype, recvbuf,
                         recvcounts, rdispls, recvtype);
}

/**
 *
 */
void CommSerial::v_AllGather(void *sendbuf, int sendcount,
                             CommDataType sendtype, void *recvbuf,
                             int recvcount, CommDataType recvtype)
{
    boost::ignore_unused(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                         recvtype);
}

void CommSerial::v_AllGatherv(void *sendbuf, int sendcount,
                              CommDataType sendtype, void *recvbuf,
                              int recvcounts[], int rdispls[],
                              CommDataType recvtype)
{
    boost::ignore_unused(sendbuf, sendcount, sendtype, recvbuf, recvcounts,
                         rdispls, recvtype);
}

void CommSerial::v_AllGatherv(void *recvbuf, int recvcounts[], int rdispls[],
                              CommDataType recvtype)
{
    boost::ignore_unused(recvbuf, recvcounts, rdispls, recvtype);
}

void CommSerial::v_Bcast(void *buffer, int count, CommDataType dt, int root)
{
    boost::ignore_unused(buffer, count, dt, root);
}

void CommSerial::v_Exscan(Array<OneD, unsigned long long> &pData,
                          const enum ReduceOperator pOp,
                          Array<OneD, unsigned long long> &ans)
{
    boost::ignore_unused(pData, pOp, ans);
}

void CommSerial::v_Gather(void *sendbuf, int sendcount, CommDataType sendtype,
                          void *recvbuf, int recvcount, CommDataType recvtype,
                          int root)
{
    boost::ignore_unused(recvcount, recvtype, root);
    std::memcpy(recvbuf, sendbuf, sendcount * CommDataTypeGetSize(sendtype));
}

void CommSerial::v_Scatter(void *sendbuf, int sendcount, CommDataType sendtype,
                           void *recvbuf, int recvcount, CommDataType recvtype,
                           int root)
{
    boost::ignore_unused(recvcount, recvtype, root);
    std::memcpy(recvbuf, sendbuf, sendcount * CommDataTypeGetSize(sendtype));
}

void CommSerial::v_DistGraphCreateAdjacent(int indegree, const int sources[],
                                           const int sourceweights[],
                                           int reorder)
{
    boost::ignore_unused(indegree, sources, sourceweights, reorder);
}

void CommSerial::v_NeighborAlltoAllv(void *sendbuf, int sendcounts[],
                                     int senddispls[], CommDataType sendtype,
                                     void *recvbuf, int recvcounts[],
                                     int rdispls[], CommDataType recvtype)
{
    boost::ignore_unused(sendbuf, sendcounts, senddispls, sendtype, recvbuf,
                         recvcounts, rdispls, recvtype);
}

void CommSerial::v_Irsend(void *buf, int count, CommDataType dt, int dest,
                          CommRequestSharedPtr request, int loc)
{
    boost::ignore_unused(buf, count, dt, dest, request, loc);
}

void CommSerial::v_SendInit(void *buf, int count, CommDataType dt, int dest,
                            CommRequestSharedPtr request, int loc)
{
    boost::ignore_unused(buf, count, dt, dest, request, loc);
}

void CommSerial::v_Irecv(void *buf, int count, CommDataType dt, int source,
                         CommRequestSharedPtr request, int loc)
{
    boost::ignore_unused(buf, count, dt, source, request, loc);
}

void CommSerial::v_RecvInit(void *buf, int count, CommDataType dt, int source,
                            CommRequestSharedPtr request, int loc)
{
    boost::ignore_unused(buf, count, dt, source, request, loc);
}

void CommSerial::v_WaitAll(CommRequestSharedPtr request)
{
    boost::ignore_unused(request);
}

void CommSerial::v_StartAll(CommRequestSharedPtr request)
{
    boost::ignore_unused(request);
}

CommRequestSharedPtr CommSerial::v_CreateRequest(int num)
{
    boost::ignore_unused(num);
    return std::shared_ptr<CommRequest>(new CommRequest);
}

/**
 *
 */
void CommSerial::v_SplitComm(int pRows, int pColumns)
{
    boost::ignore_unused(pRows, pColumns);
    ASSERTL0(false, "Cannot split a serial process.");
}

/**
 *
 */
CommSharedPtr CommSerial::v_CommCreateIf(int flag)
{
    if (flag == 0)
    {
        // flag == 0 => get back MPI_COMM_NULL, return a null ptr instead.
        return std::shared_ptr<Comm>();
    }
    else
    {
        // Return a real communicator
        return shared_from_this();
    }
}
} // namespace LibUtilities
} // namespace Nektar
