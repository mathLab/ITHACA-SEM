///////////////////////////////////////////////////////////////////////////////
//
// File CoupledLcoalToGlobalC0ContMap.cpp
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
// Description: Wrapper class around the library
// LocalToGlobalC0ContMap class for use in the Couplied Linearised NS
// solver.
///////////////////////////////////////////////////////////////////////////////

#include <LinearElasticSolver/EquationSystems/CoupledAssemblyMap.h>
#include <SpatialDomains/MeshGraph.h>
#include <LocalRegions/SegExp.h>
#include <LocalRegions/Expansion1D.h>
#include <LocalRegions/Expansion2D.h>
#include <MultiRegions/GlobalLinSysDirectStaticCond.h>

namespace Nektar
{    
    /** 
     * This is an vector extension of
     * MultiRegion::LocalToGlobalC0BaseMap::SetUp2DExpansionC0ContMap
     * related to the Linearised Navier Stokes problem
     */
    CoupledAssemblyMap::CoupledAssemblyMap(
        const LibUtilities::SessionReaderSharedPtr        &pSession,
        const SpatialDomains::MeshGraphSharedPtr          &graph,
        const SpatialDomains::BoundaryConditionsSharedPtr &boundaryConditions,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const bool CheckforSingularSys):
        AssemblyMapCG2D(pSession)
    {

        int i,j,k,n;
        int cnt = 0,offset=0;
        int meshVertId;
        int meshEdgeId;
        int bndEdgeCnt;
        int globalId;
        int nEdgeCoeffs;
        int nEdgeInteriorCoeffs;
        int firstNonDirGraphVertId;
        int nLocBndCondDofs = 0;
        int nLocDirBndCondDofs = 0;
        int nExtraDirichlet = 0;
        LocalRegions::ExpansionSharedPtr  locExpansion;
        LocalRegions::SegExpSharedPtr     bndSegExp;
        LibUtilities::BasisType           bType;
        StdRegions::Orientation           edgeOrient;
        Array<OneD, unsigned int>         edgeInteriorMap;
        Array<OneD, int>                  edgeInteriorSign;
        int nvel = fields.num_elements();
        
        const LocalRegions::ExpansionVector &locExpVector = *(fields[0]->GetExp());
        int eid, id, diff;
        int nel = fields[0]->GetNumElmts();

        MultiRegions::PeriodicMap periodicVerts;
        MultiRegions::PeriodicMap periodicEdges;
        Array<OneD, map<int,int> > ReorderedGraphVertId(2);
        MultiRegions::BottomUpSubStructuredGraphSharedPtr bottomUpGraph;
        int staticCondLevel = 0;

        if(CheckforSingularSys)
        {
            m_systemSingular = true;
        }
        else
        {
            m_systemSingular = false;
        }
        
        /**
         * STEP 1: Wrap boundary conditions vector in an array
         * (since routine is set up for multiple fields) and call
         * the graph re-odering subroutine to obtain the reordered
         * values
         */

        // Obtain any periodic information and allocate default mapping array
        fields[0]->GetPeriodicEntities(periodicVerts,periodicEdges);

        const Array<OneD, const MultiRegions::ExpListSharedPtr> &bndCondExp
            = fields[0]->GetBndCondExpansions();

        Array<OneD, Array<OneD, const SpatialDomains::BoundaryConditionShPtr> > bndConditionsVec(nvel);
        for(i = 0; i < nvel; ++i)
        {
            bndConditionsVec[i] = fields[i]->GetBndConditions(); 
        }

        for(j = 0; j < bndCondExp.num_elements(); ++j)
        {
            // Make sure that only Dirichlet boundary conditions have been set.
            for (i = 0; i < nvel; ++i)
            {
                ASSERTL0(bndConditionsVec[i][j]->GetBoundaryConditionType() ==
                             SpatialDomains::eDirichlet,
                         "Boundary conditions must be Dirichlet; others not "
                         "supported yet.");
            }
        }

        Array<OneD, map<int, int> > Dofs(2);
        int edgeId, vertId;

        ASSERTL0(!m_systemSingular, "System should not be singular...");

        for(i = 0; i < nel; ++i)
        {
            locExpansion = locExpVector[eid];
            for(j = 0; j < locExpVector[i]->GetNverts(); ++j)
            {
                // Vert Dofs
                Dofs[0][locExpansion->GetGeom()->GetVid(j)] = nvel;
                // Edge Dofs
                Dofs[1][locExpansion->GetGeom()->GetEid(j)] =
                    nvel*(exp2d->GetEdgeNcoeffs(j)-2);
            }
        }

        set<int> extraDir;
        SetUp2DGraphC0ContMap(*fields[0],
                              bndCondExp,
                              bndConditionsVec,
                              periodicVerts,          periodicEdges,
                              Dofs,                   ReorderedGraphVertId,
                              firstNonDirGraphVertId, nExtraDirichlet,
                              bottomUpGraph, extraDir,  false,  4);

        /**
         * STEP 2: Count out the number of Dirichlet vertices and edges first
         */
        for(i = 0; i < bndCondExp.num_elements(); i++)
        {
            for(j = 0; j < bndCondExp[i]->GetNumElmts(); j++)
            {
                bndSegExp = boost::dynamic_pointer_cast<
                    LocalRegions::SegExp>(bndCondExp[i]->GetExp(j));
                nLocDirBndCondDofs += bndSegExp->GetNcoeffs() * nvel;
                nLocBndCondDofs    += bndSegExp->GetNcoeffs() * nvel;
            }
        }

        m_numLocalDirBndCoeffs = nLocDirBndCondDofs + nExtraDirichlet;

        /**
         * STEP 3: Set up an array which contains the offset information of the
         * different graph vertices.
         *
         * This basically means to identify how many global degrees of freedom
         * the individual graph vertices correspond. Obviously, the graph
         * vertices corresponding to the mesh-vertices account for a single
         * global DOF. However, the graph vertices corresponding to the element
         * edges correspond to 2*(N-2) global DOF where N is equal to the number
         * of boundary modes on this edge.
         */
        Array<OneD, int> graphVertOffset(
            nvel * (ReorderedGraphVertId[0].size() +
                    ReorderedGraphVertId[1].size()), 0);
        graphVertOffset[0] = 0;
        
        m_signChange = false;

        // Initially set up graphVertOffset exactly in the single-field case.
        for(i = 0; i < nel; ++i)
        {
            eid = fields[0]->GetOffset_Elmt_Id(i);
            locExpansion = locExpVector[eid];

            for(j = 0; j < locExpansion->GetNedges(); ++j)
            {
                nEdgeCoeffs = locExpansion->GetEdgeNcoeffs(j);
                meshEdgeId  = locExpansion->GetGeom()->GetEid(j);
                meshVertId  = locExpansion->GetGeom()->GetVid(j);
                bType       = locExpansion->GetEdgeBasisType(j);

                graphVertOffset[ReorderedGraphVertId[0][meshVertId]+1] = 1;
                graphVertOffset[ReorderedGraphVertId[1][meshEdgeId]+1] = nEdgeCoeffs-2;

                // need a sign vector for modal expansions if nEdgeCoeffs
                // >=4
                if (nEdgeCoeffs >= 4                     &&
                    (bType == LibUtilities::eModified_A  ||
                     bType == LibUtilities::eModified_B))
                {
                    m_signChange = true;
                }
            }
        }

        for(i = 1; i < graphVertOffset.num_elements(); i++)
        {
            graphVertOffset[i] += graphVertOffset[i-1];
        }

        // Determine number of global Dirichlet degrees of freedom for single
        // field.
        m_numGlobalDirBndCoeffs = graphVertOffset[firstNonDirGraphVertId];

        graphVertOffset[0] = 0;

        // Armed with this information, go back and reset graphVertOffset
        // correctly.
        for (nv = 0; nv < nvel; ++nv)
        {
            for(i = 0; i < nel; ++i)
            {
                eid = fields[0]->GetOffset_Elmt_Id(i);
                locExpansion = locExpVector[eid];

                for(j = 0; j < locExpansion->GetNedges(); ++j)
                {
                    nEdgeCoeffs = locExpansion->GetEdgeNcoeffs(j);
                    meshEdgeId  = locExpansion->GetGeom()->GetEid(j);
                    meshVertId  = locExpansion->GetGeom()->GetVid(j);
                    bType       = locExpansion->GetEdgeBasisType(j);

                    graphVertOffset[ReorderedGraphVertId[0][meshVertId]*nv+1] = 1;
                    graphVertOffset[ReorderedGraphVertId[1][meshEdgeId]*nv+1] = nEdgeCoeffs-2;

                    // need a sign vector for modal expansions if nEdgeCoeffs
                    // >=4
                    if (nEdgeCoeffs >= 4                     &&
                        (bType == LibUtilities::eModified_A  ||
                         bType == LibUtilities::eModified_B))
                    {
                        m_signChange = true;
                    }
                }
            }
        }

        m_numGlobalDirBndCoeffs *= 2;

        // Allocate the proper amount of space for the class-data and fill
        // information that is already known
        cnt = 0; 
        m_numLocalBndCoeffs = 0;
        m_numLocalCoeffs = 0;

        for(i = 0; i < nel; ++i)
        {
            m_numLocalBndCoeffs += (nvel*locExpVector[i]->NumBndryCoeffs() + 1);
        }

        m_numLocalCoeffs += m_numLocalBndCoeffs; 


        m_localToGlobalMap    = Array<OneD, int>(m_numLocalCoeffs,-1);
        m_localToGlobalBndMap = Array<OneD, int>(m_numLocalBndCoeffs,-1);
        m_bndCondCoeffsToGlobalCoeffsMap = Array<OneD, int>(nLocBndCondDofs,-1);
        

        // Set default sign array. 
        m_localToGlobalSign    = Array<OneD, NekDouble>(m_numLocalCoeffs,1.0);
        m_localToGlobalBndSign = Array<OneD, NekDouble>(m_numLocalBndCoeffs,1.0);

        m_staticCondLevel = staticCondLevel;
        m_numPatches = nel;
        
        m_numLocalBndCoeffsPerPatch = Array<OneD, unsigned int>(nel);
        m_numLocalIntCoeffsPerPatch = Array<OneD, unsigned int>(nel);

        for(i = 0; i < nel; ++i)
        {
            m_numLocalBndCoeffsPerPatch[i] = (unsigned int) (nvel*locExpVector[fields[0]->GetOffset_Elmt_Id(i)]->NumBndryCoeffs() + 1);
            //m_numLocalIntCoeffsPerPatch[i] = (unsigned int) nz_loc*(pressure->GetExp(eid)->GetNcoeffs()-1);
        }
        
        /**
         * STEP 4: Now, all ingredients are ready to set up the actual
         * local to global mapping.
         *
         * The remainder of the map consists of the element-interior
         * degrees of freedom. This leads to the block-diagonal submatrix
         * as each element-interior mode is globally orthogonal to modes
         * in all other elements.
         */
        cnt = 0;
        int nv,velnbndry;
        Array<OneD, unsigned int> bmap;
        
       
        // Loop over all the elements in the domain in shuffled
        // ordering (element type consistency)
        for(i = 0; i < nel; ++i)
        {
            locExpansion = boost::dynamic_pointer_cast<StdRegions::StdExpansion2D>(locExpVector[eid]);

            velnbndry = locExpansion->NumBndryCoeffs();

            // require an inverse ordering of the bmap system to store
            // local numbering system which takes matrix these
            // matrices. Therefore get hold of elemental bmap and set
            // up an inverse map
            map<int,int> inv_bmap;
            locExpansion->GetBoundaryMap(bmap);
            for(j = 0; j < bmap.num_elements(); ++j)
            {
                inv_bmap[bmap[j]] = j;
            }
            
            // Loop over all edges (and vertices) of element i
            for(j = 0; j < locExpansion->GetNedges(); ++j)
            {
                nEdgeInteriorCoeffs = locExpansion->GetEdgeNcoeffs(j)-2;
                edgeOrient          = locExpansion->GetGeom()->GetEorient(j);
                meshEdgeId          = locExpansion->GetGeom()->GetEid(j);
                meshVertId          = locExpansion->GetGeom()->GetVid(j);

                locExpansion->GetEdgeInteriorMap(j,edgeOrient,edgeInteriorMap,edgeInteriorSign);

                // Set the global DOF for vertex j of element i
                for(nv = 0; nv < nvel; ++nv)
                {
                    m_localToGlobalMap[cnt+nv*velnbndry+inv_bmap[locExpansion->GetVertexMap(j)]] = graphVertOffset[ReorderedGraphVertId[0][meshVertId]*nvel+ nv];
                    // Set the global DOF's for the interior modes of edge j
                    for(k = 0; k < nEdgeInteriorCoeffs; ++k)
                    {
                        m_localToGlobalMap[cnt+nv*velnbndry+inv_bmap[edgeInteriorMap[k]]] =  graphVertOffset[ReorderedGraphVertId[1][meshEdgeId]*nvel+nv]+k;
                    }
                }

                // Fill the sign vector if required
                if(m_signChange)
                {
                    for(nv = 0; nv < nvel; ++nv)
                    {
                        for(k = 0; k < nEdgeInteriorCoeffs; ++k)
                        {
                            m_localToGlobalSign[cnt+nv*velnbndry + inv_bmap[edgeInteriorMap[k]]] = (NekDouble) edgeInteriorSign[k];
                        }
                    }
                }
            }

            // use difference between two edges of the AddMeanPressureEdgeId to det nEdgeInteriorCoeffs. 
            //nEdgeInteriorCoeffs = graphVertOffset[(ReorderedGraphVertId[1][eid])*nvel+1] - graphVertOffset[(ReorderedGraphVertId[1][eid])*nvel];
            cnt += (velnbndry*nvel+ psize);
        }
       
        // Set up the mapping for the boundary conditions
        offset = cnt = 0;
        for(nv = 0; nv < nvel; ++nv)
        {
            for(i = 0; i < bndCondExp.num_elements(); i++)
            {
                int ncoeffcnt = 0;
                for(j = 0; j < bndCondExp[i]->GetNumElmts(); j++)
                {
                    bndSegExp  = boost::dynamic_pointer_cast<LocalRegions::SegExp>(bndCondExp[i]->GetExp(j));

                    cnt = offset + bndCondExp[i]->GetCoeff_Offset(j);
                    for(k = 0; k < 2; k++)
                    {
                        meshVertId = (bndSegExp->GetGeom1D())->GetVid(k);
                        m_bndCondCoeffsToGlobalCoeffsMap[cnt+bndSegExp->GetVertexMap(k)] = graphVertOffset[ReorderedGraphVertId[0][meshVertId]*nvel+nv];
                    }
                    
                    meshEdgeId = (bndSegExp->GetGeom1D())->GetEid();
                    bndEdgeCnt = 0;
                    nEdgeCoeffs = bndSegExp->GetNcoeffs();
                    for(k = 0; k < nEdgeCoeffs; k++)
                    {
                        if(m_bndCondCoeffsToGlobalCoeffsMap[cnt+k] == -1)
                        {
                            m_bndCondCoeffsToGlobalCoeffsMap[cnt+k] =
                                graphVertOffset[ReorderedGraphVertId[1][meshEdgeId]*nvel+nv]+bndEdgeCnt;
                            bndEdgeCnt++;
                        }
                    }
                    ncoeffcnt += nEdgeCoeffs; 
                }
                // Note: Can not use bndCondExp[i]->GetNcoeffs()
                // due to homogeneous extension not returning just
                // the value per plane
                offset += ncoeffcnt; 
            }
        }
        
        globalId = Vmath::Vmax(m_numLocalCoeffs,&m_localToGlobalMap[0],1)+1;
        m_numGlobalBndCoeffs = globalId; 
        
        /**
         * STEP 5: The boundary condition mapping is generated from the
         * same vertex renumbering and fill in a unique interior map.
         */
        cnt=0;
        for(i = 0; i < m_numLocalCoeffs; ++i)
        {
            if(m_localToGlobalMap[i] == -1)
            {
                m_localToGlobalMap[i] = globalId++;
            }
            else
            {
                if(m_signChange)
                {
                    m_localToGlobalBndSign[cnt]=m_localToGlobalSign[i];
                }
                m_localToGlobalBndMap[cnt++]=m_localToGlobalMap[i];
            }
        }
        m_numGlobalCoeffs = globalId;

        // Set up the local to global map for the next level when using
        // multi-level static condensation
        if( m_session->MatchSolverInfoAsEnum("GlobalSysSoln", MultiRegions::eDirectMultiLevelStaticCond) )
        {
            if (m_staticCondLevel < (bottomUpGraph->GetNlevels()-1))
            {
                Array<OneD, int> vwgts_perm(
                    Dofs[0].size()+Dofs[1].size()-firstNonDirGraphVertId);
                for(i = 0; i < locExpVector.size(); ++i)
                {
                    int eid = fields[0]->GetOffset_Elmt_Id(i);
                    locExpansion = boost::dynamic_pointer_cast<
                        StdRegions::StdExpansion2D>(locExpVector[eid]);
                    for(j = 0; j < locExpansion->GetNverts(); ++j)
                    {
                        meshEdgeId = (LocalRegions::Expansion2D::FromStdExp(locExpansion)->GetGeom2D())->GetEid(j);
                        meshVertId = (LocalRegions::Expansion2D::FromStdExp(locExpansion)->GetGeom2D())->GetVid(j);
                            
                        if(ReorderedGraphVertId[0][meshVertId] >= 
                           firstNonDirGraphVertId)
                        {
                            vwgts_perm[ReorderedGraphVertId[0][meshVertId]-
                                       firstNonDirGraphVertId] = 
                                Dofs[0][meshVertId];
                        }
                            
                        if(ReorderedGraphVertId[1][meshEdgeId] >= 
                           firstNonDirGraphVertId)
                        {
                            vwgts_perm[ReorderedGraphVertId[1][meshEdgeId]-
                                       firstNonDirGraphVertId] = 
                                Dofs[1][meshEdgeId];
                        }
                    }
                }

                bottomUpGraph->ExpandGraphWithVertexWeights(vwgts_perm);
                    
                m_nextLevelLocalToGlobalMap = MemoryManager<AssemblyMap>::
                    AllocateSharedPtr(this,bottomUpGraph);
            }
        }
    }
}
