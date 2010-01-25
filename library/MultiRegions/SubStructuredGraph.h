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
// Description: a collection of classes that facilitates the implementation
//              of the multi-level static condensation routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef MULTIREGIONS_SUBSTRUCTUREDGRAPH_H
#define MULTIREGIONS_SUBSTRUCTUREDGRAPH_H

#include <MultiRegions/MultiRegions.hpp>

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

        typedef boost::shared_ptr<BottomUpSubStructuredGraph> BottomUpSubStructuredGraphSharedPtr;
        typedef boost::shared_ptr<SubGraph>                   SubGraphSharedPtr;
        typedef boost::shared_ptr<MultiLevelBisectedGraph>    MultiLevelBisectedGraphSharedPtr;
        typedef boost::shared_ptr<PatchMap>                   PatchMapSharedPtr;


        class PatchMap
        {
        public:

        PatchMap(const int patchId,
                 const int dofId,
                 const bool bndPatch,
                 const NekDouble sign):
            m_patchId(patchId),
                m_dofId(dofId),
                m_bndPatch(bndPatch),
                m_sign(sign)
                {
                }
            
            inline int GetPatchId() const 
            {
                return m_patchId;
            }

            inline int GetDofId() const
            {
                return m_dofId;
            }

            inline bool IsBndDof() const
            {
                return m_bndPatch;
            }

            inline NekDouble GetSign() const
            {
                return m_sign;
            }

        protected:
            int       m_patchId;
            int       m_dofId;
            bool      m_bndPatch; 
            NekDouble m_sign; 
        };



        class SubGraph
        {
        public:

            SubGraph(const int nVerts, const int idOffset = 0):
                m_nVerts(nVerts),
                m_idOffset(idOffset)
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
            MultiLevelBisectedGraph(const Array<OneD, const int> sepTree);
            MultiLevelBisectedGraph(const int nBndDofs);

            int  GetTotDofs() const;

            void SetGlobalNumberingOffset();

            void DumpNBndDofs(void) const;

            void CollectLeaves(vector<SubGraphSharedPtr>& leaves) const;

            inline int  GetNdaughterGraphs() const;

            int CutLeaves();

            int CutEmptyLeaves();

            inline const SubGraphSharedPtr GetBndDofsGraph() const
            {
                return m_BndDofs;
            }

        protected:
            SubGraphSharedPtr m_BndDofs;
            MultiLevelBisectedGraphSharedPtr m_leftDaughterGraph;
            MultiLevelBisectedGraphSharedPtr m_rightDaughterGraph;        
        };


        class BottomUpSubStructuredGraph
        {
        public:
            BottomUpSubStructuredGraph(const Array<OneD, const int> septree);
            BottomUpSubStructuredGraph(const MultiLevelBisectedGraphSharedPtr& graph);
            BottomUpSubStructuredGraph(const int nVerts);

            int GetTotDofs() const;

            void UpdateBottomUpReordering(Array<OneD,       int>& perm, 
                                          Array<OneD,       int>& iperm) const;

            void ExpandGraphWithVertexWeights(const Array<OneD, const int>& wgts);

            void MaskPatches(const int leveltomask, Array<OneD, NekDouble>& maskarray) const;
            
            int GetNpatchesWithInterior(const int whichlevel) const;

            void GetNintDofsPerPatch(const int whichlevel, Array<OneD, unsigned int>& outarray) const;
            
            int GetInteriorOffset(const int whichlevel, const int patch = 0) const;

            int GetNumGlobalDofs(const int whichlevel) const;

            int GetNlevels() const;

            void Dump() const;

        protected:
            vector<SubGraphSharedPtr> m_IntBlocks;
            BottomUpSubStructuredGraphSharedPtr m_daughterGraph;

        private:
            void SetBottomUpReordering(Array<OneD, int>& iperm) const;

            inline BottomUpSubStructuredGraphSharedPtr GetDaughterGraph() const
            {
                return m_daughterGraph;
            }
            inline vector<SubGraphSharedPtr> GetInteriorBlocks() const
            {
                return m_IntBlocks;
            }
        };


        namespace 
        {
            typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> BoostGraph;
        }

        void CuthillMckeeReordering(const BoostGraph& graph,
                                    Array<OneD, int>& perm,
                                    Array<OneD, int>& iperm);

        void MultiLevelBisectionReordering(const BoostGraph& graph,
                                           const Array<OneD, const int>& vwgts,
                                           Array<OneD, int>& perm,
                                           Array<OneD, int>& iperm,
                                           BottomUpSubStructuredGraphSharedPtr& substructgraph,
                                           const int mdswitch = 1);
        // The parameter MDSWITCH.
        // This parameters defines the maximal size of the smallest patches.
        // If at a certain level, a patch bundles less than MDSWITCH graph-vertices,
        // metis is not going to partition this subgraph any further.
        // Some quick and basis test have shown that 30 seems to be good value.
        // However, this optimal value will probably depend on the polynomial order
        // of the expansion and there is still room for optimisation here.

        void NoReordering(const BoostGraph& graph,
                          Array<OneD, int>& perm,
                          Array<OneD, int>& iperm);



        
    } // end of namespace
} // end of namespace

#endif // MULTIREGIONS_SUBSTRUCTUREDGRAPH_H

/**
* $Log: SubStructuredGraph.h,v $
* Revision 1.4  2009/11/19 11:41:07  pvos
* Fixed various bugs
*
* Revision 1.3  2009/11/09 15:57:11  pvos
* multi-level recursion bug fixes
*
* Revision 1.2  2009/11/02 11:19:44  pvos
* Fixed a bug for reordering a graph without edges
*
* Revision 1.1  2009/10/30 14:02:55  pvos
* Multi-level static condensation updates
*
**/
