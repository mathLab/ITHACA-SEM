///////////////////////////////////////////////////////////////////////////////
//
// File AssemblyMapCG.cpp
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
// Description: C0-continuous Local to Global mapping routines, base class
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/AssemblyMap/AssemblyMapCG.h>
#include <MultiRegions/ExpList.h>

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class AssemblyMapCG
         * Mappings are created for three possible global solution types:
         *  - Direct full matrix
         *  - Direct static condensation
         *  - Direct multi-level static condensation
         * In the latter case, mappings are created recursively for the
         * different levels of static condensation.
         *
         * These mappings are used by GlobalLinSys to generate the global
         * system.
         */

        /**
         *
         */
        AssemblyMapCG::AssemblyMapCG(
                const LibUtilities::SessionReaderSharedPtr &pSession):
            AssemblyMap(pSession)
        {
            pSession->LoadParameter("MaxStaticCondLevel",m_maxStaticCondLevel,100);
        }


        /**
         *
         */
        AssemblyMapCG::AssemblyMapCG(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const int numLocalCoeffs,
                const ExpList &locExp):
            AssemblyMap(pSession)
        {
            ASSERTL0(false,"AssemblyMapCG: you need to instantiate dimension-specific derived class.");
        }



        /**
         *
         */
        AssemblyMapCG::~AssemblyMapCG()
        {
        }



        /**
         * Sets up the global to universal mapping of degrees of freedom across
         * processors.
         */
        void AssemblyMapCG::SetUpUniversalC0ContMap(
                const ExpList &locExp)
        {
            StdRegions::StdExpansionSharedPtr locExpansion;
            int nDim = 0;
            int nVert = 0;
            int nEdge = 0;
            int nFace = 0;
            int maxEdgeDof = 0;
            int maxFaceDof = 0;
            int maxIntDof = 0;
            int dof = 0;
            int cnt;
            int i,j,k;
            int meshVertId;
            int meshEdgeId;
            int meshFaceId;
            int elementId;
            int vGlobalId;
            int maxBndGlobalId = 0;
            StdRegions::Orientation         edgeOrient;
            StdRegions::Orientation         faceOrient;
            Array<OneD, unsigned int>           edgeInteriorMap;
            Array<OneD, int>                    edgeInteriorSign;
            Array<OneD, unsigned int>           faceInteriorMap;
            Array<OneD, int>                    faceInteriorSign;
            Array<OneD, unsigned int>           interiorMap;

            const StdRegions::StdExpansionVector &locExpVector = *(locExp.GetExp());
            LibUtilities::CommSharedPtr vCommRow = m_comm->GetRowComm();

            m_globalToUniversalMap = Nektar::Array<OneD, int>(m_numGlobalCoeffs, -1);
            m_globalToUniversalMapUnique = Nektar::Array<OneD, int>(m_numGlobalCoeffs, -1);
            m_globalToUniversalBndMap = Nektar::Array<OneD, int>(m_numGlobalBndCoeffs, -1);
            m_globalToUniversalBndMapUnique = Nektar::Array<OneD, int>(m_numGlobalBndCoeffs, -1);

            // Loop over all the elements in the domain to gather mesh data
            for(i = 0; i < locExpVector.size(); ++i)
            {
                locExpansion = boost::dynamic_pointer_cast<StdRegions::StdExpansion>(locExpVector[i]);
                nVert += locExpansion->GetNverts();
                nEdge += locExpansion->GetNedges();
                nFace += locExpansion->GetNfaces();
                // Loop over all edges (and vertices) of element i
                for(j = 0; j < locExpansion->GetNedges(); ++j)
                {
                    dof = locExpansion->GetEdgeNcoeffs(j)-2;
                    maxEdgeDof = (dof > maxEdgeDof ? dof : maxEdgeDof);
                }
                for(j = 0; j < locExpansion->GetNfaces(); ++j)
                {
                    dof = locExpansion->GetFaceIntNcoeffs(j);
                    maxFaceDof = (dof > maxFaceDof ? dof : maxFaceDof);
                }
                locExpansion->GetInteriorMap(interiorMap);
                dof = interiorMap.num_elements();
                maxIntDof = (dof > maxIntDof ? dof : maxIntDof);
            }

            // Tell other processes about how many dof we have
            vCommRow->AllReduce(nVert, LibUtilities::ReduceSum);
            vCommRow->AllReduce(nEdge, LibUtilities::ReduceSum);
            vCommRow->AllReduce(nFace, LibUtilities::ReduceSum);
            vCommRow->AllReduce(maxEdgeDof, LibUtilities::ReduceMax);
            vCommRow->AllReduce(maxFaceDof, LibUtilities::ReduceMax);
            vCommRow->AllReduce(maxIntDof, LibUtilities::ReduceMax);

            // Assemble global to universal mapping for this process
            for(i = 0; i < locExpVector.size(); ++i)
            {
                locExpansion = boost::dynamic_pointer_cast<StdRegions::StdExpansion>(locExpVector[i]);
                nDim = locExpansion->GetShapeDimension();
                cnt = locExp.GetCoeff_Offset(i);

                // Loop over all vertices of element i
                for(j = 0; j < locExpansion->GetNverts(); ++j)
                {
                    meshVertId   = (locExpansion->GetGeom())->GetVid(j);
                    vGlobalId    = m_localToGlobalMap[cnt+locExpansion->GetVertexMap(j)];
                    m_globalToUniversalMap[vGlobalId] = meshVertId + 1;
                    m_globalToUniversalBndMap[vGlobalId]=m_globalToUniversalMap[vGlobalId];
                    maxBndGlobalId = (vGlobalId > maxBndGlobalId ? vGlobalId : maxBndGlobalId);
                }

                // Loop over all edges of element i
                for(j = 0; j < locExpansion->GetNedges(); ++j)
                {
                    edgeOrient = (locExpansion->GetGeom())->GetEorient(j);
                    locExpansion->GetEdgeInteriorMap(j,edgeOrient,edgeInteriorMap,edgeInteriorSign);
                    dof = locExpansion->GetEdgeNcoeffs(j)-2;
                    if (nDim == 2)
                    {
                        meshEdgeId   = (locExpansion->GetGeom2D())->GetEid(j);
                    }
                    else
                    {
                        meshEdgeId   = (locExpansion->GetGeom3D())->GetEid(j);
                    }

                    // Set the global DOF's for the interior modes of edge j
                    for(k = 0; k < dof; ++k)
                    {
                        vGlobalId = m_localToGlobalMap[cnt+edgeInteriorMap[k]];
                        m_globalToUniversalMap[vGlobalId]
                           = nVert + meshEdgeId * maxEdgeDof + k + 1;
                        m_globalToUniversalBndMap[vGlobalId]=m_globalToUniversalMap[vGlobalId];
                        maxBndGlobalId = (vGlobalId > maxBndGlobalId ? vGlobalId : maxBndGlobalId);
                    }
                }

                // Loop over all faces of element i
                for(j = 0; j < locExpansion->GetNfaces(); ++j)
                {
                    faceOrient          = (locExpansion->GetGeom3D())->GetFaceOrient(j);

                    locExpansion->GetFaceInteriorMap(j,faceOrient,faceInteriorMap,faceInteriorSign);
                    dof = locExpansion->GetFaceIntNcoeffs(j);
                    meshFaceId = (locExpansion->GetGeom3D())->GetFid(j);

                    for(k = 0; k < dof; ++k)
                    {
                        vGlobalId = m_localToGlobalMap[cnt+faceInteriorMap[k]];
                        m_globalToUniversalMap[vGlobalId]
                           = nVert + nEdge*maxEdgeDof + meshFaceId * maxFaceDof
                                   + k + 1;
                        m_globalToUniversalBndMap[vGlobalId]=m_globalToUniversalMap[vGlobalId];

                        maxBndGlobalId = (vGlobalId > maxBndGlobalId ? vGlobalId : maxBndGlobalId);
                    }
                }

                // Add interior DOFs to complete universal numbering
                locExpansion->GetInteriorMap(interiorMap);
                dof = interiorMap.num_elements();
                elementId = (locExpansion->GetGeom())->GetGlobalID();
                for (k = 0; k < dof; ++k)
                {
                    vGlobalId = m_localToGlobalMap[cnt+interiorMap[k]];
                    m_globalToUniversalMap[vGlobalId]
                           = nVert + nEdge*maxEdgeDof + nFace*maxFaceDof + elementId*maxIntDof + k + 1;
                }
            }

            // Set up the GSLib universal assemble mapping
            // Internal DOF do not participate in any data
            // exchange, so we keep these set to the special GSLib id=0 so
            // they are ignored.
            Nektar::Array<OneD, long> tmp(m_numGlobalCoeffs);
            Vmath::Zero(m_numGlobalCoeffs, tmp, 1);
            Nektar::Array<OneD, long> tmp2(m_numGlobalBndCoeffs, tmp);
            for (unsigned int i = 0; i < m_numGlobalBndCoeffs; ++i)
            {
                tmp[i] = m_globalToUniversalMap[i];
            }

            m_gsh = Gs::Init(tmp, vCommRow);
            m_bndGsh = Gs::Init(tmp2, vCommRow);
            Gs::Unique(tmp, vCommRow);
            for (unsigned int i = 0; i < m_numGlobalCoeffs; ++i)
            {
                m_globalToUniversalMapUnique[i] = (tmp[i] >= 0 ? 1 : 0);
            }
            for (unsigned int i = 0; i < m_numGlobalBndCoeffs; ++i)
            {
                m_globalToUniversalBndMapUnique[i] = (tmp2[i] >= 0 ? 1 : 0);
            }
        }

        /**
         * The bandwidth calculated here corresponds to what is referred to as
         * half-bandwidth.  If the elements of the matrix are designated as
         * a_ij, it corresponds to the maximum value of |i-j| for non-zero
         * a_ij.  As a result, the value also corresponds to the number of
         * sub- or super-diagonals.
         *
         * The bandwith can be calculated elementally as it corresponds to the
         * maximal elemental bandwith (i.e. the maximal difference in global
         * DOF index for every element).
         *
         * We caluclate here the bandwith of the full global system.
         */
        void AssemblyMapCG::CalculateFullSystemBandWidth()
        {
            int i,j;
            int cnt = 0;
            int locSize;
            int maxId;
            int minId;
            int bwidth = -1;
            for(i = 0; i < m_numPatches; ++i)
            {
                locSize = m_numLocalBndCoeffsPerPatch[i]+m_numLocalIntCoeffsPerPatch[i];
                maxId = -1;
                minId = m_numLocalCoeffs+1;
                for(j = 0; j < locSize; j++)
                {
                    if(m_localToGlobalMap[cnt+j] >= m_numGlobalDirBndCoeffs)
                    {
                        if(m_localToGlobalMap[cnt+j] > maxId)
                        {
                            maxId = m_localToGlobalMap[cnt+j];
                        }

                        if(m_localToGlobalMap[cnt+j] < minId)
                        {
                            minId = m_localToGlobalMap[cnt+j];
                        }
                    }
                }
                bwidth = (bwidth>(maxId-minId))?bwidth:(maxId-minId);

                cnt+=locSize;
            }

            m_fullSystemBandWidth = bwidth;
        }


        int AssemblyMapCG::v_GetLocalToGlobalMap(const int i) const
        {
            return m_localToGlobalMap[i];
        }

        int AssemblyMapCG::v_GetGlobalToUniversalMap(const int i) const
        {
            return m_globalToUniversalMap[i];
        }

        int AssemblyMapCG::v_GetGlobalToUniversalMapUnique(const int i) const
        {
            return m_globalToUniversalMapUnique[i];
        }

        const Array<OneD,const int>&
                    AssemblyMapCG::v_GetLocalToGlobalMap(void)
        {
            return m_localToGlobalMap;
        }

        const Array<OneD,const int>&
                    AssemblyMapCG::v_GetGlobalToUniversalMap(void)
        {
            return m_globalToUniversalMap;
        }

        const Array<OneD,const int>&
                    AssemblyMapCG::v_GetGlobalToUniversalMapUnique(void)
        {
            return m_globalToUniversalMapUnique;
        }

        NekDouble AssemblyMapCG::v_GetLocalToGlobalSign(
                    const int i) const
        {
            if(m_signChange)
            {
                return m_localToGlobalSign[i];
            }
            else
            {
                return 1.0;
            }
        }

        const Array<OneD, NekDouble>& AssemblyMapCG::v_GetLocalToGlobalSign() const
        {
            return m_localToGlobalSign;
        }

        void AssemblyMapCG::v_LocalToGlobal(
                    const Array<OneD, const NekDouble>& loc,
                          Array<OneD,       NekDouble>& global) const
        {
            Array<OneD, const NekDouble> local;
            if(global.data() == loc.data())
            {
                local = Array<OneD, NekDouble>(loc.num_elements(),loc.data());
            }
            else
            {
                local = loc; // create reference
            }


            if(m_signChange)
            {
                Vmath::Scatr(m_numLocalCoeffs, m_localToGlobalSign.get(), local.get(), m_localToGlobalMap.get(), global.get());
            }
            else
            {
                Vmath::Scatr(m_numLocalCoeffs, local.get(), m_localToGlobalMap.get(), global.get());
            }
        }

        void AssemblyMapCG::v_LocalToGlobal(
                    const NekVector<NekDouble>& loc,
                          NekVector<      NekDouble>& global) const
        {
            LocalToGlobal(loc.GetPtr(),global.GetPtr());
        }

        void AssemblyMapCG::v_GlobalToLocal(
                    const Array<OneD, const NekDouble>& global,
                          Array<OneD,       NekDouble>& loc) const
        {
            Array<OneD, const NekDouble> glo;
            if(global.data() == loc.data())
            {
                glo = Array<OneD, NekDouble>(global.num_elements(),global.data());
            }
            else
            {
                glo = global; // create reference
            }
            
                
            if(m_signChange)
            {
                Vmath::Gathr(m_numLocalCoeffs, m_localToGlobalSign.get(), glo.get(), m_localToGlobalMap.get(), loc.get());
            }
            else
            {
                Vmath::Gathr(m_numLocalCoeffs, glo.get(), m_localToGlobalMap.get(), loc.get());
            }
        }

        void AssemblyMapCG::v_GlobalToLocal(
                    const NekVector<NekDouble>& global,
                          NekVector<      NekDouble>& loc) const
        {
            GlobalToLocal(global.GetPtr(),loc.GetPtr());
        }

        void AssemblyMapCG::v_Assemble(
                    const Array<OneD, const NekDouble> &loc,
                          Array<OneD,       NekDouble> &global) const
        {
            Array<OneD, const NekDouble> local;
            if(global.data() == loc.data())
            {
                local = Array<OneD, NekDouble>(local.num_elements(),local.data());
            }
            else
            {
                local = loc; // create reference
            }
            //ASSERTL1(loc.get() != global.get(),"Local and Global Arrays cannot be the same");

            Vmath::Zero(m_numGlobalCoeffs, global.get(), 1);

            if(m_signChange)
            {
                Vmath::Assmb(m_numLocalCoeffs, m_localToGlobalSign.get(), local.get(), m_localToGlobalMap.get(), global.get());
            }
            else
            {
                Vmath::Assmb(m_numLocalCoeffs, local.get(), m_localToGlobalMap.get(), global.get());
            }
            UniversalAssemble(global);
        }

        void AssemblyMapCG::v_Assemble(
                    const NekVector<NekDouble>& loc,
                          NekVector<      NekDouble>& global) const
        {
            Assemble(loc.GetPtr(),global.GetPtr());
        }

        void AssemblyMapCG::v_UniversalAssemble(
                      Array<OneD,     NekDouble>& pGlobal) const
        {
            Gs::Gather(pGlobal, Gs::gs_add, m_gsh);
        }

        void AssemblyMapCG::v_UniversalAssemble(
                      NekVector<      NekDouble>& pGlobal) const
        {
            UniversalAssemble(pGlobal.GetPtr());
        }

        void AssemblyMapCG::v_UniversalAssemble(
                      Array<OneD,     NekDouble>& pGlobal,
                      int                         offset) const
        {
            Array<OneD, NekDouble> tmp(offset);
            Vmath::Vcopy(offset, pGlobal, 1, tmp, 1);
            UniversalAssemble(pGlobal);
            Vmath::Vcopy(offset, tmp, 1, pGlobal, 1);
        }

        int AssemblyMapCG::v_GetFullSystemBandWidth() const
        {
            return m_fullSystemBandWidth;
        }

        int AssemblyMapCG::v_GetNumNonDirVertexModes() const
        {
            return m_numNonDirVertexModes;
        }

        int AssemblyMapCG::v_GetNumNonDirEdgeModes() const
        {
            return m_numNonDirEdgeModes;
        }

        int AssemblyMapCG::v_GetNumNonDirFaceModes() const
        {
            return m_numNonDirFaceModes;
        }

    } // namespace
} // namespace
