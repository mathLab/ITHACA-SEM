///////////////////////////////////////////////////////////////////////////////
//
// File SubStructuredGraph.h
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
// Description: a collection of classes that facilitates the implementation
//              of the multi-level static condensation routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef MULTIREGIONS_SUBSTRUCTUREDGRAPH_H
#define MULTIREGIONS_SUBSTRUCTUREDGRAPH_H

#include <MultiRegions/MultiRegionsDeclspec.h>

#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include <vector>
#include <set>
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>

namespace Nektar
{
    namespace MultiRegions
    {
        class BottomUpSubStructuredGraph;
        class SubGraph;
        class MultiLevelBisectedGraph;
        class PatchMap;

        typedef std::shared_ptr<BottomUpSubStructuredGraph> 
            BottomUpSubStructuredGraphSharedPtr;
        typedef std::shared_ptr<SubGraph>
            SubGraphSharedPtr;
        typedef std::shared_ptr<MultiLevelBisectedGraph>
            MultiLevelBisectedGraphSharedPtr;
        typedef std::shared_ptr<PatchMap>
            PatchMapSharedPtr;

        class PatchMap
        {
        public:
            MULTI_REGIONS_EXPORT      PatchMap(void);
            MULTI_REGIONS_EXPORT      PatchMap(const int vals);
            MULTI_REGIONS_EXPORT     ~PatchMap(void);
            MULTI_REGIONS_EXPORT void SetPatchMap(
                const int          n, 
                const int          patchId, 
                const int          dofId,
                const unsigned int bndPatch,
                const NekDouble    sign);

            void SetNewLevelMap(Array<OneD, const unsigned int> numLocalBndCondPerPatch,
                                Array<OneD, const unsigned int> numLocalIntCondPerPatch);


            inline Array<OneD, const int> GetPatchId() const 
            {
                return m_patchId;
            }

            inline Array<OneD, const int>  GetDofId() const
            {
                return m_dofId;
            }
            
            inline Array<OneD, const int>  GetNewLevelMap() const
            {
                return m_newLevelMap;
            }
            
            inline Array<OneD, const unsigned int> IsBndDof() const
            {
                return m_bndPatch;
            }
            
            inline Array<OneD, const NekDouble> GetSign() const
            {
                return m_sign;
            }
            
        protected:
            Array<OneD, int>          m_patchId;
            Array<OneD, int>          m_dofId;
            Array<OneD, int>          m_newLevelMap;
            Array<OneD, unsigned int> m_bndPatch; 
            Array<OneD, NekDouble>    m_sign; 
        };

        class SubGraph
        {
        public:
            
            MULTI_REGIONS_EXPORT SubGraph(
                const int nVerts, const int idOffset = 0) :
                m_nVerts(nVerts),
                m_idOffset(idOffset)
            {
            }
            
            MULTI_REGIONS_EXPORT ~SubGraph(void)
            {
            }
            
            inline int GetNverts(void) const
            {
                return m_nVerts;
            }

            inline void SetNverts(const int i) 
            {
                m_nVerts = i;
            }

            inline int GetIdOffset(void) const
            {
                return m_idOffset;
            }

            inline void SetIdOffset(const int i) 
            {
                m_idOffset = i;
            }

        protected:
            int m_nVerts;
            int m_idOffset;
        };

        bool SubGraphWithoutVerts(const SubGraphSharedPtr g);

        class MultiLevelBisectedGraph
        {
        public:
            MULTI_REGIONS_EXPORT MultiLevelBisectedGraph(
                MultiLevelBisectedGraphSharedPtr oldLevel,
                const int                        nPartition);
            MULTI_REGIONS_EXPORT MultiLevelBisectedGraph(
                const int nBndDofs);
            MULTI_REGIONS_EXPORT ~MultiLevelBisectedGraph(void);

            MULTI_REGIONS_EXPORT int  GetTotDofs() const;
            MULTI_REGIONS_EXPORT void SetGlobalNumberingOffset();
            MULTI_REGIONS_EXPORT void DumpNBndDofs(void) const;
            MULTI_REGIONS_EXPORT void CollectLeaves(
                std::vector<SubGraphSharedPtr>& leaves) const;
            MULTI_REGIONS_EXPORT int  CutLeaves();
            MULTI_REGIONS_EXPORT int  CutEmptyLeaves();
            MULTI_REGIONS_EXPORT inline int GetNdaughterGraphs() const
            {
                return m_daughterGraphs.size();
            }

            inline const SubGraphSharedPtr GetBndDofsGraph() const
            {
                return m_BndDofs;
            }

            inline
            std::vector<MultiLevelBisectedGraphSharedPtr> &GetDaughterGraphs()
            {
                return m_daughterGraphs;
            }

        protected:
            SubGraphSharedPtr                             m_BndDofs;
            std::vector<MultiLevelBisectedGraphSharedPtr> m_daughterGraphs;
        };


        class BottomUpSubStructuredGraph
        {
        public:
            MULTI_REGIONS_EXPORT BottomUpSubStructuredGraph(
                MultiLevelBisectedGraphSharedPtr graph,
                int nPartition = 0,
                bool globaloffset = false);
            MULTI_REGIONS_EXPORT BottomUpSubStructuredGraph(
                const int nVerts);
            MULTI_REGIONS_EXPORT ~BottomUpSubStructuredGraph(void);

            MULTI_REGIONS_EXPORT int GetTotDofs() const;

            MULTI_REGIONS_EXPORT void UpdateBottomUpReordering(
                Array<OneD, int>& perm,  
                Array<OneD, int>& iperm) const;

            MULTI_REGIONS_EXPORT void ExpandGraphWithVertexWeights(
                const Array<OneD, const int>& wgts);

            MULTI_REGIONS_EXPORT void MaskPatches(
                const int               leveltomask, 
                Array<OneD, NekDouble>& maskarray) const;
            
            MULTI_REGIONS_EXPORT int GetNpatchesWithInterior(
                const int whichlevel) const;

            MULTI_REGIONS_EXPORT void GetNintDofsPerPatch(
                const int                  whichlevel, 
                Array<OneD, unsigned int> &outarray) const;
            
            MULTI_REGIONS_EXPORT int GetInteriorOffset(
                const int whichlevel, 
                const int patch = 0) const;

            MULTI_REGIONS_EXPORT std::vector<SubGraphSharedPtr> 
                GetInteriorBlocks(const int whichlevel) const;

            MULTI_REGIONS_EXPORT int GetNumGlobalDofs(
                const int whichlevel) const;

            MULTI_REGIONS_EXPORT int GetNlevels() const;

            MULTI_REGIONS_EXPORT void Dump() const;

        protected:
            std::vector<SubGraphSharedPtr> m_IntBlocks;
            BottomUpSubStructuredGraphSharedPtr m_daughterGraph;

        private:
            void SetBottomUpReordering(Array<OneD, int>& iperm) const;

            inline BottomUpSubStructuredGraphSharedPtr GetDaughterGraph() const
            {
                return m_daughterGraph;
            }
            inline std::vector<SubGraphSharedPtr> GetInteriorBlocks() const
            {
                return m_IntBlocks;
            }
        };

        namespace 
        {
            typedef boost::adjacency_list<
                boost::setS, boost::vecS, boost::undirectedS> BoostGraph;
        }

        MULTI_REGIONS_EXPORT void CuthillMckeeReordering(
            const BoostGraph& graph,
            Array<OneD, int>& perm,
            Array<OneD, int>& iperm);

        MULTI_REGIONS_EXPORT void MultiLevelBisectionReordering(
            const BoostGraph                    &graph,
            Array<OneD, int>                    &perm,
            Array<OneD, int>                    &iperm,
            BottomUpSubStructuredGraphSharedPtr &substructgraph,
            std::set<int>                        partVerts = std::set<int>(),
            int                                  mdswitch  = 1);
        
        // The parameter MDSWITCH.
        //
        // This parameters defines the maximal size of the smallest patches.  If
        // at a certain level, a patch bundles less than MDSWITCH
        // graph-vertices, metis is not going to partition this subgraph any
        // further.  Some quick and basis test have shown that 30 seems to be
        // good value.  However, this optimal value will probably depend on the
        // polynomial order of the expansion and there is still room for
        // optimisation here.

        MULTI_REGIONS_EXPORT void NoReordering(const BoostGraph& graph,
                          Array<OneD, int>& perm,
                          Array<OneD, int>& iperm);
    } // end of namespace
} // end of namespace

#endif // MULTIREGIONS_SUBSTRUCTUREDGRAPH_H
