///////////////////////////////////////////////////////////////////////////////
//
// File AssemblyMapDG.h
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

#ifndef MULTIREGIONS_ASSEMBLY_MAP_DG_H
#define MULTIREGIONS_ASSEMBLY_MAP_DG_H

#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/AssemblyMap/AssemblyMap.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList0D.h>

#include <LibUtilities/Communication/CommMpi.h>  // TODO: Implement the graph constructor using virtual functions / wrappers so can remove this include
#include <LibUtilities/BasicUtils/Timer.h>
namespace Nektar
{
    namespace MultiRegions
    {
        class AssemblyMapDG;
        typedef std::shared_ptr<AssemblyMapDG>  AssemblyMapDGSharedPtr;

        ///
        class AssemblyMapDG: public AssemblyMap
        {
        public:
            NekDouble m_mpiTime;
            /// Default constructor.
            MULTI_REGIONS_EXPORT AssemblyMapDG();

            /// Constructor for trace map for one-dimensional expansion.
            MULTI_REGIONS_EXPORT AssemblyMapDG( 
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const SpatialDomains::MeshGraphSharedPtr &graph1D,
                const ExpListSharedPtr &trace,
                const ExpList &locExp,
                const Array<OneD, const MultiRegions::ExpListSharedPtr>
                                                                &bndConstraint,
                const Array<OneD, const SpatialDomains::BoundaryConditionShPtr>
                                                                &bndCond,
                const PeriodicMap &periodicTrace,
                const std::string variable = "DefaultVar");

            /// Destructor.
            MULTI_REGIONS_EXPORT virtual ~AssemblyMapDG();

            /// Return the number of boundary segments on which Dirichlet
            /// boundary conditions are imposed.
            MULTI_REGIONS_EXPORT int GetNumDirichletBndPhys();

            MULTI_REGIONS_EXPORT Array<OneD, LocalRegions::ExpansionSharedPtr>
                &GetElmtToTrace(const int i);

            MULTI_REGIONS_EXPORT 
                Array<OneD,Array<OneD,LocalRegions::ExpansionSharedPtr> >
                &GetElmtToTrace();

            MULTI_REGIONS_EXPORT int GetTraceToUniversalMap(int i);

            MULTI_REGIONS_EXPORT int GetTraceToUniversalMapUnique(int i);

            MULTI_REGIONS_EXPORT void UniversalTraceAssemble(
                    Array<OneD, NekDouble>& Fwd, Array<OneD, NekDouble>& Bwd);

            MULTI_REGIONS_EXPORT void UniversalTraceAssembleGS(
                    Array<OneD, NekDouble> &pGlobal) const;

            MULTI_REGIONS_EXPORT std::function<void(const Array<OneD, NekDouble> &, Array<OneD, NekDouble> &)> MPITraceAssemble;


        protected:
            Gs::gs_data * m_traceGsh;
            
            /// Number of physical dirichlet boundary values in trace
            int m_numDirichletBndPhys;

            /// list of edge expansions for a given element
            Array<OneD, Array<OneD, LocalRegions::ExpansionSharedPtr> > m_elmtToTrace;
            /// Integer map of process trace space quadrature points to
            /// universal space.
            Array<OneD,int> m_traceToUniversalMap;
            /// Integer map of unique process trace space quadrature points to
            /// universal space (signed).
            Array<OneD,int> m_traceToUniversalMapUnique;

            /// Map of ID to quad point trace indices
            std::map<int, std::vector<int>> m_edgeToTrace;

            void MPIInitialiseStructure(
                    const LocalRegions::ExpansionVector &locExpVector,
                    const Array<OneD, const MultiRegions::ExpListSharedPtr> &bndCondExp,
                    const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndCond,
                    const PeriodicMap &perMap);
            /// Map of process to shared edge IDs
            std::map<int, std::vector<int>> m_rankSharedEdges;

            void MPISetupAllToAll();
            /// List of trace map indices of the quad points to exchange
            std::vector<int> m_allEdgeIndex;
            /// Max number of quadrature points in an element
            int m_maxQuad = 0;
            /// Largest shared partition edge
            int m_maxCount = 0;
            /// Number of ranks/processes/partitions
            int m_nRanks = 0;
            void MPIPerformAllToAll(const Array<OneD, double> &testFwd, Array<OneD, double> &testBwd);

            void MPISetupAllToAllV();
            /// List of trace map indices of the quad points to exchange
            std::vector<int> m_allVEdgeIndex;
            /// List of counts for MPI_alltoallv
            Array<OneD, int> m_allVSendCount;
            /// List of displacements for MPI_alltoallv
            Array<OneD, int> m_allVSendDisp;
            void MPIPerformAllToAllV(const Array<OneD, double> &testFwd, Array<OneD, double> &testBwd);

            void MPISetupNeighborAllToAllV();
            ///Distributed graph communicator
            MPI_Comm m_commGraph;
            /// List of trace map indices of the quad points to exchange
            std::vector<int> m_edgeTraceIndex;
            /// List of counts
            Array<OneD,int> m_sendCount;
            /// List of displacements
            Array<OneD, int> m_sendDisp;
            void MPIPerformNeighborAllToAllV(const Array<OneD, double> &testFwd, Array<OneD, double> &testBwd);

            void MPISetupPairwise();
            /// List of partition to trace map indices of the quad points to exchange
            std::vector<std::pair<int, std::vector<int>>> m_vecPairPartitionTrace;
            /// Total quadrature points to send/recv
            int m_totSends = 0;
            void MPIPerformPairwise(const Array<OneD, double> &testFwd, Array<OneD, double> &testBwd);

            void SetUpUniversalDGMap(const ExpList &locExp);

            void SetUpUniversalTraceMap(
                const ExpList         &locExp,
                const ExpListSharedPtr trace,
                const Array<OneD, const MultiRegions::ExpListSharedPtr>         &bndCondExp,
                const Array<OneD, const SpatialDomains::BoundaryConditionShPtr>  &bndCond,
                const PeriodicMap     &perMap = NullPeriodicMap);

            virtual int v_GetLocalToGlobalMap(const int i) const;

            virtual int v_GetGlobalToUniversalMap(const int i) const;

            virtual int v_GetGlobalToUniversalMapUnique(const int i) const;

            virtual const Array<OneD,const int>&  v_GetLocalToGlobalMap();

            virtual const Array<OneD, const int>& v_GetGlobalToUniversalMap();

            virtual const Array<OneD, const int>& v_GetGlobalToUniversalMapUnique();

            virtual NekDouble v_GetLocalToGlobalSign(const int i) const;

            virtual void v_LocalToGlobal(
                    const Array<OneD, const NekDouble>& loc,
                    Array<OneD,       NekDouble>& global,
                    bool useComm) const;

            virtual void v_LocalToGlobal(
                    const NekVector<NekDouble>& loc,
                    NekVector<      NekDouble>& global,
                    bool useComm) const;

            virtual void v_GlobalToLocal(
                    const Array<OneD, const NekDouble>& global,
                          Array<OneD,       NekDouble>& loc) const;

            virtual void v_GlobalToLocal(
                    const NekVector<NekDouble>& global,
                          NekVector<      NekDouble>& loc) const;

            virtual void v_Assemble(
                    const Array<OneD, const NekDouble> &loc,
                          Array<OneD,       NekDouble> &global) const;

            virtual void v_Assemble(
                    const NekVector<NekDouble>& loc,
                          NekVector<      NekDouble>& global) const;

            virtual void v_UniversalAssemble(
                          Array<OneD,     NekDouble>& pGlobal) const;

            virtual void v_UniversalAssemble(
                          NekVector<      NekDouble>& pGlobal) const;

            virtual int v_GetFullSystemBandWidth() const;

            void RealignTraceElement(
                Array<OneD, int>        &toAlign,
                StdRegions::Orientation  orient,
                int                      nquad1,
                int                      nquad2 = 0);
        }; // class

        using func_t = std::function<void
                (const Array<OneD, NekDouble> &, Array<OneD, NekDouble> &)>;

        static inline std::tuple<NekDouble, NekDouble, NekDouble>  MPITiming(
                const LibUtilities::CommSharedPtr &comm,
                const int &count,
                const Array<OneD, double> &testFwd,
                Array<OneD, double> &testBwd,
                const func_t &f)
        {
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

        enum MPIType
        {
            eAllToAll ,
            eAllToAllV,
            eNeighborAllToAll,
            ePairwise
        };

        const char* const MPITypeMap[] =
        {
                "AllToAll",
                "AllToAllV",
                "NeighborAllToAllV",
                "PairwiseSendRecv"
        };

    } // end of namespace
} // end of namespace

#endif //MULTIREGIONS_ASSEMBLY_MAP_DG_H
