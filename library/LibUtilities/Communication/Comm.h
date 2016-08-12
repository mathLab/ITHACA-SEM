///////////////////////////////////////////////////////////////////////////////
//
// File Comm.h
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
// Description: Base communication class
//
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_LIB_UTILITIES_COMM_H
#define NEKTAR_LIB_UTILITIES_COMM_H

#include <vector>

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <boost/enable_shared_from_this.hpp>
#include <boost/static_assert.hpp>

#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
// namespace Nektar { template <typename Dim, typename DataType> class Array; }
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Communication/CommDataType.h>

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
typedef LibUtilities::NekFactory<std::string, Comm, int, char **> CommFactory;

LIB_UTILITIES_EXPORT CommFactory &GetCommFactory();

/// Type of operation to perform in AllReduce.
enum ReduceOperator
{
    ReduceSum,
    ReduceMax,
    ReduceMin
};

/// Base communications class
class Comm : public boost::enable_shared_from_this<Comm>
{
public:
    LIB_UTILITIES_EXPORT Comm(int narg, char *arg[]);
    LIB_UTILITIES_EXPORT virtual ~Comm();

    LIB_UTILITIES_EXPORT inline void Finalise();

    /// Returns number of processes
    LIB_UTILITIES_EXPORT inline int GetSize();
    LIB_UTILITIES_EXPORT inline int GetRank();
    LIB_UTILITIES_EXPORT inline const std::string &GetType() const;

    /// Block execution until all processes reach this point
    LIB_UTILITIES_EXPORT inline void Block();

    /// Return the time in seconds
    LIB_UTILITIES_EXPORT inline NekDouble Wtime();

    template <class T> void Send(int pProc, T &pData);
    template <class T> void Recv(int pProc, T &pData);
    template <class T>
    void SendRecv(int pSendProc, T &pSendData, int pRecvProc, T &pRecvData);
    template <class T>
    void SendRecvReplace(int pSendProc, int pRecvProc, T &pData);

    template <class T> void AllReduce(T &pData, enum ReduceOperator pOp);

    template <class T> void AlltoAll(T &pSendData, T &pRecvData);
    template <class T>
    void AlltoAllv(Array<OneD, T> &pSendData,
                   Array<OneD, int> &pSendDataSizeMap,
                   Array<OneD, int> &pSendDataOffsetMap,
                   Array<OneD, T> &pRecvData,
                   Array<OneD, int> &pRecvDataSizeMap,
                   Array<OneD, int> &pRecvDataOffsetMap);

    template <class T> void Bcast(T &data, int rootProc);

    template <class T>
    void Exscan(T &pData, const enum ReduceOperator pOp, T &ans);

    template <class T> T Gather(const int rootProc, T &val);
    template <class T> T Scatter(const int rootProc, T &pData);

    LIB_UTILITIES_EXPORT inline CommSharedPtr CommCreateIf(int flag);

    LIB_UTILITIES_EXPORT inline void SplitComm(int pRows, int pColumns);
    LIB_UTILITIES_EXPORT inline CommSharedPtr GetRowComm();
    LIB_UTILITIES_EXPORT inline CommSharedPtr GetColumnComm();

    LIB_UTILITIES_EXPORT inline bool TreatAsRankZero(void);
    LIB_UTILITIES_EXPORT inline bool RemoveExistingFiles(void);

protected:
    int m_size;                 ///< Number of processes
    std::string m_type;         ///< Type of communication
    CommSharedPtr m_commRow;    ///< Row communicator
    CommSharedPtr m_commColumn; ///< Column communicator

    Comm();

    virtual void v_Finalise()   = 0;
    virtual int v_GetRank()     = 0;
    virtual void v_Block()      = 0;
    virtual NekDouble v_Wtime() = 0;
    virtual void v_Send(void *buf, int count, CommDataType dt, int dest)   = 0;
    virtual void v_Recv(void *buf, int count, CommDataType dt, int source) = 0;
    virtual void v_SendRecv(void *sendbuf, int sendcount, CommDataType sendtype,
                            int dest, void *recvbuf, int recvcount,
                            CommDataType recvtype, int source) = 0;
    virtual void v_SendRecvReplace(void *buf, int count, CommDataType dt,
                                   int pSendProc, int pRecvProc) = 0;
    virtual void v_AllReduce(void *buf, int count, CommDataType dt,
                             enum ReduceOperator pOp) = 0;
    virtual void v_AlltoAll(void *sendbuf, int sendcount, CommDataType sendtype,
                            void *recvbuf, int recvcount,
                            CommDataType recvtype) = 0;
    virtual void v_AlltoAllv(void *sendbuf, int sendcounts[], int sensdispls[],
                             CommDataType sendtype, void *recvbuf,
                             int recvcounts[], int rdispls[],
                             CommDataType recvtype) = 0;
    virtual void v_Bcast(void *buffer, int count, CommDataType dt,
                         int root) = 0;

    virtual void v_Exscan(Array<OneD, unsigned long long> &pData,
                          const enum ReduceOperator pOp,
                          Array<OneD, unsigned long long> &ans) = 0;

    virtual void v_Gather(void *sendbuf, int sendcount, CommDataType sendtype,
                          void *recvbuf, int recvcount, CommDataType recvtype,
                          int root) = 0;
    virtual void v_Scatter(void *sendbuf, int sendcount, CommDataType sendtype,
                           void *recvbuf, int recvcount, CommDataType recvtype,
                           int root) = 0;

    virtual CommSharedPtr v_CommCreateIf(int flag) = 0;
    virtual void v_SplitComm(int pRows, int pColumns) = 0;
    virtual bool v_TreatAsRankZero(void) = 0;
    LIB_UTILITIES_EXPORT virtual bool v_RemoveExistingFiles(void);
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
inline const std::string &Comm::GetType() const
{
    return m_type;
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
inline double Comm::Wtime()
{
    return v_Wtime();
}

template <class T> void Comm::Send(int pProc, T &pData)
{
    v_Send(CommDataTypeTraits<T>::GetPointer(pData),
           CommDataTypeTraits<T>::GetCount(pData),
           CommDataTypeTraits<T>::GetDataType(), pProc);
}

template <class T> void Comm::Recv(int pProc, T &pData)
{
    v_Recv(CommDataTypeTraits<T>::GetPointer(pData),
           CommDataTypeTraits<T>::GetCount(pData),
           CommDataTypeTraits<T>::GetDataType(), pProc);
}

/**
 *
 */
template <class T>
void Comm::SendRecv(int pSendProc, T &pSendData, int pRecvProc, T &pRecvData)
{
    v_SendRecv(CommDataTypeTraits<T>::GetPointer(pSendData),
               CommDataTypeTraits<T>::GetCount(pSendData),
               CommDataTypeTraits<T>::GetDataType(), pSendProc,
               CommDataTypeTraits<T>::GetPointer(pRecvData),
               CommDataTypeTraits<T>::GetCount(pRecvData),
               CommDataTypeTraits<T>::GetDataType(), pRecvProc);
}

/**
 *
 */
template <class T>
void Comm::SendRecvReplace(int pSendProc, int pRecvProc, T &pData)
{
    v_SendRecvReplace(CommDataTypeTraits<T>::GetPointer(pData),
                      CommDataTypeTraits<T>::GetCount(pData),
                      CommDataTypeTraits<T>::GetDataType(), pSendProc,
                      pRecvProc);
}

/**
 *
 */
template <class T> void Comm::AllReduce(T &pData, enum ReduceOperator pOp)
{
    v_AllReduce(CommDataTypeTraits<T>::GetPointer(pData),
                CommDataTypeTraits<T>::GetCount(pData),
                CommDataTypeTraits<T>::GetDataType(), pOp);
}

template <class T> void Comm::AlltoAll(T &pSendData, T &pRecvData)
{
    BOOST_STATIC_ASSERT_MSG(
        CommDataTypeTraits<T>::IsVector,
        "AlltoAll only valid with Array or vector arguments.");
    int sendSize = CommDataTypeTraits<T>::GetCount(pSendData);
    int recvSize = CommDataTypeTraits<T>::GetCount(pRecvData);
    ASSERTL0(sendSize == recvSize,
             "Send and Recv arrays have incompatible sizes in AlltoAll");

    int count = sendSize / GetSize();
    ASSERTL0(count * GetSize() == sendSize,
             "Array size incompatible with size of communicator");

    v_AlltoAll(CommDataTypeTraits<T>::GetPointer(pSendData), count,
               CommDataTypeTraits<T>::GetDataType(),
               CommDataTypeTraits<T>::GetPointer(pRecvData), count,
               CommDataTypeTraits<T>::GetDataType());
}

/**
 *
 */
template <class T>
void Comm::AlltoAllv(Array<OneD, T> &pSendData,
                     Array<OneD, int> &pSendDataSizeMap,
                     Array<OneD, int> &pSendDataOffsetMap,
                     Array<OneD, T> &pRecvData,
                     Array<OneD, int> &pRecvDataSizeMap,
                     Array<OneD, int> &pRecvDataOffsetMap)
{
    v_AlltoAllv(pSendData.get(), pSendDataSizeMap.get(),
                pSendDataOffsetMap.get(), CommDataTypeTraits<T>::GetDataType(),
                pRecvData.get(), pRecvDataSizeMap.get(),
                pRecvDataOffsetMap.get(), CommDataTypeTraits<T>::GetDataType());
}

/**
 *
 */
template <class T> void Comm::Bcast(T &pData, int pRoot)
{
    v_Bcast(CommDataTypeTraits<T>::GetPointer(pData),
            CommDataTypeTraits<T>::GetCount(pData),
            CommDataTypeTraits<T>::GetDataType(), pRoot);
}

template <class T>
void Comm::Exscan(T &pData, const enum ReduceOperator pOp, T &ans)
{
    ASSERTL0(CommDataTypeTraits<T>::GetCount(pData) ==
                 CommDataTypeTraits<T>::GetCount(ans),
             "Input and output array sizes don't match");
    v_Exscan(CommDataTypeTraits<T>::GetPointer(pData),
             CommDataTypeTraits<T>::GetPointer(ans),
             CommDataTypeTraits<T>::GetCount(pData),
             CommDataTypeTraits<T>::GetDataType(), pOp);
}

/**
 * Concatenate all the input arrays, in rank order, onto the process with rank
 * == rootProc
 */
template <class T> T Comm::Gather(const int rootProc, T &val)
{
    BOOST_STATIC_ASSERT_MSG(
        CommDataTypeTraits<T>::IsVector,
        "Gather only valid with Array or vector arguments.");
    bool amRoot  = (GetRank() == rootProc);
    unsigned nEl = CommDataTypeTraits<T>::GetCount(val);

    unsigned nOut = amRoot ? GetSize() * nEl : 0;
    T ans(nOut);
    void *recvbuf = amRoot ? CommDataTypeTraits<T>::GetPointer(ans) : NULL;

    v_Gather(CommDataTypeTraits<T>::GetPointer(val), nEl,
             CommDataTypeTraits<T>::GetDataType(), recvbuf, nEl,
             CommDataTypeTraits<T>::GetDataType(), rootProc);
    return ans;
}
/**
 * Scatter pData across ranks in chunks of len(pData)/num_ranks
 */
template <class T> T Comm::Scatter(const int rootProc, T &pData)
{
    BOOST_STATIC_ASSERT_MSG(
        CommDataTypeTraits<T>::IsVector,
        "Scatter only valid with Array or vector arguments.");

    bool amRoot  = (GetRank() == rootProc);
    unsigned nEl = CommDataTypeTraits<T>::GetCount(pData) / GetSize();

    void *sendbuf = amRoot ? CommDataTypeTraits<T>::GetPointer(pData) : NULL;
    T ans(nEl);

    v_Scatter(sendbuf, nEl, CommDataTypeTraits<T>::GetDataType(),
              CommDataTypeTraits<T>::GetPointer(ans), nEl,
              CommDataTypeTraits<T>::GetDataType(), rootProc);
    return ans;
}

/**
 * @brief If the flag is non-zero create a new communicator.
 */
inline CommSharedPtr Comm::CommCreateIf(int flag)
{
    return v_CommCreateIf(flag);
}

/**
 * @brief Splits this communicator into a grid of size pRows*pColumns
 * and creates row and column communicators. By default the communicator
 * is a single row.
 */
inline void Comm::SplitComm(int pRows, int pColumns)
{
    v_SplitComm(pRows, pColumns);
}

/**
 * @brief Retrieve the row communicator to which this process belongs.
 */
inline CommSharedPtr Comm::GetRowComm()
{
    if (!m_commRow.get())
    {
        return shared_from_this();
    }
    else
    {
        return m_commRow;
    }
}

/**
 * @brief Retrieve the column communicator to which this process
 * belongs.
 */
inline CommSharedPtr Comm::GetColumnComm()
{
    if (!m_commColumn.get())
    {
        return shared_from_this();
    }
    else
    {
        return m_commColumn;
    }
}

inline bool Comm::TreatAsRankZero(void)
{
    return v_TreatAsRankZero();
}

inline bool Comm::RemoveExistingFiles(void)
{
    return v_RemoveExistingFiles();
}
}
}

#endif
