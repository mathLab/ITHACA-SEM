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

#include <mpi.h>
#include <string>

#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/LibUtilitiesDeclspec.h>
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
typedef std::shared_ptr<CommMpi> CommMpiSharedPtr;

/// Class for communicator request type
class CommRequestMpi final : public CommRequest
{
public:
    /// Creates an instance of this class
    inline explicit CommRequestMpi(int num) : m_num(num)
    {
        m_request = Array<OneD, MPI_Request>(num, MPI_REQUEST_NULL);
    }

    /// Default destructor
    inline ~CommRequestMpi() final = default;

    inline MPI_Request *GetRequest(int i)
    {
        return &m_request[i];
    }

    inline int &GetNumRequest()
    {
        return m_num;
    }

private:
    int m_num;
    Array<OneD, MPI_Request> m_request;
};

typedef std::shared_ptr<CommRequestMpi> CommRequestMpiSharedPtr;

/// A global linear system.
class CommMpi : public Comm
{
public:
    /// Creates an instance of this class
    static CommSharedPtr create(int narg, char *arg[])
    {
        return MemoryManager<CommMpi>::AllocateSharedPtr(narg, arg);
    }

    /// Name of class
    LIB_UTILITIES_EXPORT static std::string className;

    LIB_UTILITIES_EXPORT CommMpi(int narg, char *arg[]);
    LIB_UTILITIES_EXPORT virtual ~CommMpi() override;

    LIB_UTILITIES_EXPORT MPI_Comm GetComm();

protected:
    MPI_Comm m_comm;
    int m_rank{};

    explicit CommMpi(MPI_Comm pComm);

    virtual void v_Finalise() final;
    virtual int v_GetRank() final;
    virtual void v_Block() final;
    virtual double v_Wtime() final;
    virtual bool v_TreatAsRankZero() final;
    virtual bool v_IsSerial() final;
    virtual std::tuple<int, int, int> v_GetVersion() final;
    virtual void v_Send(void *buf, int count, CommDataType dt, int dest) final;
    virtual void v_Recv(void *buf, int count, CommDataType dt,
                        int source) final;
    virtual void v_SendRecv(void *sendbuf, int sendcount, CommDataType sendtype,
                            int dest, void *recvbuf, int recvcount,
                            CommDataType recvtype, int source) final;
    virtual void v_SendRecvReplace(void *buf, int count, CommDataType dt,
                                   int pSendProc, int pRecvProc) final;
    virtual void v_AllReduce(void *buf, int count, CommDataType dt,
                             enum ReduceOperator pOp) final;
    virtual void v_AlltoAll(void *sendbuf, int sendcount, CommDataType sendtype,
                            void *recvbuf, int recvcount,
                            CommDataType recvtype) final;
    virtual void v_AlltoAllv(void *sendbuf, int sendcounts[], int sensdispls[],
                             CommDataType sendtype, void *recvbuf,
                             int recvcounts[], int rdispls[],
                             CommDataType recvtype) final;
    virtual void v_AllGather(void *sendbuf, int sendcount,
                             CommDataType sendtype, void *recvbuf,
                             int recvcount, CommDataType recvtype) final;
    virtual void v_AllGatherv(void *sendbuf, int sendcount,
                              CommDataType sendtype, void *recvbuf,
                              int recvcounts[], int rdispls[],
                              CommDataType recvtype) final;
    virtual void v_AllGatherv(void *recvbuf, int recvcounts[], int rdispls[],
                              CommDataType recvtype) final;
    virtual void v_Bcast(void *buffer, int count, CommDataType dt,
                         int root) final;
    virtual void v_Exscan(Array<OneD, unsigned long long> &pData,
                          enum ReduceOperator pOp,
                          Array<OneD, unsigned long long> &ans) final;

    virtual void v_Gather(void *sendbuf, int sendcount, CommDataType sendtype,
                          void *recvbuf, int recvcount, CommDataType recvtype,
                          int root) final;
    virtual void v_Scatter(void *sendbuf, int sendcount, CommDataType sendtype,
                           void *recvbuf, int recvcount, CommDataType recvtype,
                           int root) final;

    virtual void v_DistGraphCreateAdjacent(int indegree, const int sources[],
                                           const int sourceweights[],
                                           int reorder) final;
    virtual void v_NeighborAlltoAllv(void *sendbuf, int sendcounts[],
                                     int sensdispls[], CommDataType sendtype,
                                     void *recvbuf, int recvcounts[],
                                     int rdispls[],
                                     CommDataType recvtype) final;
    virtual void v_Irsend(void *buf, int count, CommDataType dt, int dest,
                          CommRequestSharedPtr request, int loc) final;
    virtual void v_SendInit(void *buf, int count, CommDataType dt, int dest,
                            CommRequestSharedPtr request, int loc) final;
    virtual void v_Irecv(void *buf, int count, CommDataType dt, int source,
                         CommRequestSharedPtr request, int loc) final;
    virtual void v_RecvInit(void *buf, int count, CommDataType dt, int source,
                            CommRequestSharedPtr request, int loc) final;
    virtual void v_StartAll(CommRequestSharedPtr request) final;
    virtual void v_WaitAll(CommRequestSharedPtr request) final;
    virtual CommRequestSharedPtr v_CreateRequest(int num) final;

    virtual void v_SplitComm(int pRows, int pColumns) override;
    virtual CommSharedPtr v_CommCreateIf(int flag) final;
};
} // namespace LibUtilities
} // namespace Nektar

#endif
