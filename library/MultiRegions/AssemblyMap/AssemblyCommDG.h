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
// Description: Local to Global DG mapping routines, header file
//
///////////////////////////////////////////////////////////////////////////////

#ifndef MULTIREGIONS_ASSEMBLY_COMM_DG_H
#define MULTIREGIONS_ASSEMBLY_COMM_DG_H

#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <LibUtilities/BasicUtils/Timer.h>

namespace Nektar
{
namespace MultiRegions
{
class AssemblyCommDG : public AssemblyMapDG
{
    public:
        /// Default constructor.
        MULTI_REGIONS_EXPORT AssemblyCommDG();

        /// Destructor.
        MULTI_REGIONS_EXPORT ~AssemblyCommDG();

        // Constructor for MPI communication methods
        MULTI_REGIONS_EXPORT AssemblyCommDG(
            const ExpList &locExp, const ExpListSharedPtr &trace,
            const Array<OneD, const MultiRegions::ExpListSharedPtr> &bndCondExp,
            const Array<OneD, const SpatialDomains::BoundaryConditionShPtr>
                &bndCond,
            const PeriodicMap &perMap);

        MULTI_REGIONS_EXPORT void PerformExchange(const Array<OneD, double> &testFwd, Array<OneD, double> &testBwd);

    protected:
        /// Max number of quadrature points in an element
        int m_maxQuad = 0;
        /// Number of ranks/processes/partitions
        int m_nRanks = 0;
        /// Map of process to shared edge IDs
        std::map<int, std::vector<int>> m_rankSharedEdges;
        /// Map of edge ID to quad point trace indices
        std::map<int, std::vector<int>> m_edgeToTrace;

    private:
        constexpr static const char* const MPITypeMap[] =
            {
                "AllToAll",
                "AllToAllV",
                "NeighborAllToAllV",
                "PairwiseSendRecv"
            };

        MULTI_REGIONS_EXPORT std::function<void(
            const Array<OneD, NekDouble> &, Array<OneD, NekDouble> &)> m_MPITraceAssemble;

    private:
        void InitialiseStructure(
            const ExpList &locExp, const ExpListSharedPtr &trace,
            const Array<OneD, const MultiRegions::ExpListSharedPtr> &bndCondExp,
            const Array<OneD, const SpatialDomains::BoundaryConditionShPtr>
            &bndCond,
            const PeriodicMap &perMap);
};

class AllToAll: public AssemblyCommDG
{
    public:
        /// Default constructor.
        AllToAll();

        void PerformExchange(const Array<OneD, double> &testFwd, Array<OneD, double> &testBwd);

    private:
        /// List of trace map indices of the quad points to exchange
        std::vector<int> m_allEdgeIndex;
        /// Largest shared partition edge
        int m_maxCount = 0;
}; //  class

class AllToAllV: public AssemblyCommDG
{
    public:
        /// Default constructor.
        AllToAllV();

        void PerformExchange(const Array<OneD, double> &testFwd, Array<OneD, double> &testBwd);

    private:
        /// List of trace map indices of the quad points to exchange
        std::vector<int> m_allVEdgeIndex;
        /// List of counts for MPI_alltoallv
        Array<OneD, int> m_allVSendCount;
        /// List of displacements for MPI_alltoallv
        Array<OneD, int> m_allVSendDisp;
};

class NeighborAllToAllV: public AssemblyCommDG
{
    public:
        /// Default constructor.
        NeighborAllToAllV();

        void PerformExchange(const Array<OneD, double> &testFwd, Array<OneD, double> &testBwd);

    private:
        /// List of displacements
        Array<OneD, int> m_sendDisp;
        ///Distributed graph communicator
        MPI_Comm m_commGraph;
        /// List of trace map indices of the quad points to exchange
        std::vector<int> m_edgeTraceIndex;
        /// List of counts
        Array<OneD,int> m_sendCount;
};

class Pairwise: public AssemblyCommDG
{
    public:
        Pairwise();
        void PerformExchange(const Array<OneD, double> &testFwd, Array<OneD, double> &testBwd);

    private:
        /// List of partition to trace map indices of the quad points to exchange
        std::vector<std::pair<int, std::vector<int>>> m_vecPairPartitionTrace;
        /// Total quadrature points to send/recv
        int m_totSends = 0;
        /// List of displacements
        Array<OneD,int> m_sendDisp;
};

typedef std::function<void(const Array<OneD, NekDouble> &, Array<OneD, NekDouble> &)> func_t;

static inline std::tuple<NekDouble, NekDouble, NekDouble>  MPITiming(
    const LibUtilities::CommSharedPtr &comm,
    const int &count,
    const int &num,
    const func_t &f)
{

    Array<OneD, double> testFwd(num, 1);
    Array<OneD, double> testBwd(num, -2);

    LibUtilities::Timer t;
    t.Start();

    for (size_t i = 0; i < count; ++i)
    {
        f(testFwd, testBwd);
    }

    t.Stop();

    // These can just be 'reduce' but need to setup the wrapper in comm.h
    Array<OneD, NekDouble> minTime(1, t.TimePerTest(count));
    comm->AllReduce(minTime, LibUtilities::ReduceMin);

    Array<OneD, NekDouble> maxTime(1, t.TimePerTest(count));
    comm->AllReduce(maxTime, LibUtilities::ReduceMax);

    Array<OneD, NekDouble> sumTime(1, t.TimePerTest(count));
    comm->AllReduce(sumTime, LibUtilities::ReduceSum);

    return {sumTime[0]/comm->GetSize(), minTime[0], maxTime[0]};
}
} // end of namespace
} // end of namespace
#endif // NEKTAR_ASSEMBLYCOMMDG_H
