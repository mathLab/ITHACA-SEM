///////////////////////////////////////////////////////////////////////////////
//
// File AssemblyCommDG.h
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
// Description: Parallel communication methods for DG with MPI, header file
//
///////////////////////////////////////////////////////////////////////////////

#ifndef MULTIREGIONS_ASSEMBLY_COMM_DG_H
#define MULTIREGIONS_ASSEMBLY_COMM_DG_H

#include <LibUtilities/BasicUtils/Timer.h>
#include <MultiRegions/ExpList.h>

namespace Nektar
{
namespace MultiRegions
{

/**
 * The ExchangeMethod classes contain the required structure to distribute the
 * Fwd trace of partition edges to the matching locations in the Bwd trace in
 * the corresponding adjacent partitions. This allows for communication between
 * neighbouring partitions in the physical mesh by exchanging quadrature point
 * values.
 */
class ExchangeMethod
{
public:
    /// Default constructor
    MULTI_REGIONS_EXPORT ExchangeMethod() = default;

    /// Default destructor
    MULTI_REGIONS_EXPORT virtual ~ExchangeMethod() = default;

    /**
     * Perform MPI comm exchange taking the Fwd trace and
     * sending partition edge trace values to the matching locations in the
     * Bwd trace of corresponding adjacent partitions.
     *
     * @param[in] testFwd The values to send to adjacent partitions
     * @param[out] testBwd The values received from adjacent partitions
     */
    MULTI_REGIONS_EXPORT virtual void PerformExchange(
        const Array<OneD, NekDouble> &testFwd,
        Array<OneD, NekDouble> &testBwd) = 0;
};

typedef std::shared_ptr<ExchangeMethod> ExchangeMethodSharedPtr;

/**
 * If parallel operation is not indicated then use the Serial subclass which
 * does not perform any exchange.
 */
class Serial final : public ExchangeMethod
{
public:
    /// Default constructor
    MULTI_REGIONS_EXPORT Serial() = default;

    MULTI_REGIONS_EXPORT inline void PerformExchange(
        const Array<OneD, NekDouble> &testFwd,
        Array<OneD, NekDouble> &testBwd) final
    {
        boost::ignore_unused(testFwd, testBwd);
    }
};

/**
 * Uses the MPI_AllToAll collective operation to perform the exchange of
 * quadrature values. This does not allow for varying exchange array sizes so
 * padding is used to ensure all partitions send/receive the same length array.
 * All ranks communicate full array sizes to all other ranks. One collective
 * operation is posted on each rank which requires communication.
 */
class AllToAll final : public ExchangeMethod
{
public:
    /// Default constructor.
    MULTI_REGIONS_EXPORT AllToAll(
        const LibUtilities::CommSharedPtr &comm, const int &maxQuad,
        const int &nRanks,
        const std::map<int, std::vector<int>> &rankSharedEdges,
        const std::map<int, std::vector<int>> &edgeToTrace);

    MULTI_REGIONS_EXPORT void PerformExchange(
        const Array<OneD, NekDouble> &testFwd,
        Array<OneD, NekDouble> &testBwd) final;

private:
    /// Communicator
    LibUtilities::CommSharedPtr m_comm;
    /// Max number of quadrature points in an element
    int m_maxQuad = 0;
    /// Number of ranks/processes/partitions
    int m_nRanks = 0;
    /// List of trace map indices of the quad points to exchange
    std::vector<int> m_allEdgeIndex;
    /// Largest shared partition edge
    int m_maxCount = 0;
};

/**
 * Uses the MPI_AllToAllV collective operation to perform the exchange of
 * quadrature values. This allows for varying exchange array sizes to minimise
 * communication data size. All ranks communicate to all other ranks, however
 * the array size can be 0 to avoid unnecessary data transfer. One collective
 * peration is posted on each rank which requires communication.
 */
class AllToAllV final : public ExchangeMethod
{
public:
    /// Default constructor.
    MULTI_REGIONS_EXPORT AllToAllV(
        const LibUtilities::CommSharedPtr &comm,
        const std::map<int, std::vector<int>> &rankSharedEdges,
        const std::map<int, std::vector<int>> &edgeToTrace, const int &nRanks);

    MULTI_REGIONS_EXPORT void PerformExchange(
        const Array<OneD, NekDouble> &testFwd,
        Array<OneD, NekDouble> &testBwd) final;

private:
    /// Communicator
    LibUtilities::CommSharedPtr m_comm;
    /// List of trace map indices of the quad points to exchange
    std::vector<int> m_allVEdgeIndex;
    /// List of counts for MPI_alltoallv
    Array<OneD, int> m_allVSendCount;
    /// List of displacements for MPI_alltoallv
    Array<OneD, int> m_allVSendDisp;
};

/**
 * Uses the MPI_NeighborAllToAllV collective operation to perform the exchange
 * of quadrature values. This allows for varying exchange array sizes to
 * minimise communication data size. Ranks only communicate with ranks with
 * which they need to exchange data, i.e. are adjacent in the mesh or share a
 * periodic boundary condition, this further minimises unnecessary data transfer
 * over just reducing array sizes to 0 such as in MPI_AllToAllV. One collective
 * operation is posted on each rank which requires communication.
 */
class NeighborAllToAllV final : public ExchangeMethod
{
public:
    /// Default constructor.
    MULTI_REGIONS_EXPORT NeighborAllToAllV(
        const LibUtilities::CommSharedPtr &comm,
        const std::map<int, std::vector<int>> &rankSharedEdges,
        const std::map<int, std::vector<int>> &edgeToTrace);

    MULTI_REGIONS_EXPORT void PerformExchange(
        const Array<OneD, NekDouble> &testFwd,
        Array<OneD, NekDouble> &testBwd) final;

private:
    /// Communicator
    LibUtilities::CommSharedPtr m_comm;
    /// List of displacements
    Array<OneD, int> m_sendDisp;
    /// List of trace map indices of the quad points to exchange
    std::vector<int> m_edgeTraceIndex;
    /// List of counts
    Array<OneD, int> m_sendCount;
};

/**
 * Uses persistent MPI_Irecv and MPI_Isend operations to perform the exchange of
 * quadrature values. This allows for varying exchange array sizes to minimise
 * communication data size. Ranks only communicate with ranks with which they
 * need to exchange data, i.e. are adjacent in the mesh or share a periodic
 * boundary condition. On each rank there are 'n' receives and 'n' sends posted
 * where 'n' is the number of other ranks with which communication is needed.
 * We use persistent communication methods to reduce overhead.
 */
class Pairwise final : public ExchangeMethod
{
public:
    MULTI_REGIONS_EXPORT Pairwise(
        const LibUtilities::CommSharedPtr &comm,
        const std::map<int, std::vector<int>> &rankSharedEdges,
        const std::map<int, std::vector<int>> &edgeToTrace);

    MULTI_REGIONS_EXPORT void PerformExchange(
        const Array<OneD, NekDouble> &testFwd,
        Array<OneD, NekDouble> &testBwd) final;

private:
    /// Communicator
    LibUtilities::CommSharedPtr m_comm;
    /// List of partition to trace map indices of the quad points to exchange
    std::vector<std::pair<int, std::vector<int>>> m_vecPairPartitionTrace;
    /// List of displacements
    Array<OneD, int> m_sendDisp;
    /// Receive buffer for exchange
    Array<OneD, NekDouble> m_recvBuff;
    /// Send buffer for exchange
    Array<OneD, NekDouble> m_sendBuff;
    /// List of receive requests
    LibUtilities::CommRequestSharedPtr m_recvRequest;
    /// List of send requests
    LibUtilities::CommRequestSharedPtr m_sendRequest;
};

/**
 * @brief Implements communication for populating forward and backwards spaces
 * across processors in the discontinuous Galerkin routines.
 *
 * The AssemblyCommDG class constructs various exchange methods for performing
 * the action of communicating trace data from the forwards space of one
 * processor to the backwards space of the corresponding neighbour element, and
 * vice versa.
 *
 * This class initialises the structure for all exchange methods and then times
 * to determine the fastest method for the particular system configuration, if
 * running in serial configuration it assigns the #Serial exchange method. It
 * then acts as a pass through to the chosen exchange method for the
 * #PerformExchange function.
 */
class AssemblyCommDG
{
public:
    /// Default destructor
    MULTI_REGIONS_EXPORT ~AssemblyCommDG() = default;

    // Constructor for MPI communication methods
    MULTI_REGIONS_EXPORT AssemblyCommDG(
        const ExpList &locExp, const ExpListSharedPtr &trace,
        const Array<OneD, Array<OneD, LocalRegions::ExpansionSharedPtr>>
            &elmtToTrace,
        const Array<OneD, const ExpListSharedPtr> &bndCondExp,
        const Array<OneD, const SpatialDomains::BoundaryConditionShPtr>
            &bndCond,
        const PeriodicMap &perMap);

    /**
     * @brief Perform the trace exchange between processors, given the forwards
     * and backwards spaces.
     *
     * @param testFwd  Local forwards space of the trace (which will be sent)
     * @param testBwd  Local bacwards space of the trace (which will receive
     *                 contributions)
     */
    MULTI_REGIONS_EXPORT inline void PerformExchange(
        const Array<OneD, NekDouble> &testFwd, Array<OneD, NekDouble> &testBwd)
    {
        m_exchange->PerformExchange(testFwd, testBwd);
    }

private:
    /// Chosen exchange method (either fastest parallel or serial)
    ExchangeMethodSharedPtr m_exchange;
    /// Max number of quadrature points in an element
    int m_maxQuad = 0;
    /// Number of ranks/processes/partitions
    int m_nRanks = 0;
    /// Map of process to shared edge IDs
    std::map<int, std::vector<int>> m_rankSharedEdges;
    /// Map of edge ID to quad point trace indices
    std::map<int, std::vector<int>> m_edgeToTrace;

    /// Initalises the structure for the MPI communication
    void InitialiseStructure(
        const ExpList &locExp, const ExpListSharedPtr &trace,
        const Array<OneD, Array<OneD, LocalRegions::ExpansionSharedPtr>>
            &elmtToTrace,
        const Array<OneD, const ExpListSharedPtr> &bndCondExp,
        const Array<OneD, const SpatialDomains::BoundaryConditionShPtr>
            &bndCond,
        const PeriodicMap &perMap, const LibUtilities::CommSharedPtr &comm);

    /// Timing of the MPI exchange method.
    static std::tuple<NekDouble, NekDouble, NekDouble> Timing(
        const LibUtilities::CommSharedPtr &comm, const int &count,
        const int &num, const ExchangeMethodSharedPtr& f);
};

typedef std::shared_ptr<AssemblyCommDG> AssemblyCommDGSharedPtr;

} // namespace MultiRegions
} // namespace Nektar

#endif
