///////////////////////////////////////////////////////////////////////////////
//
// File LocToGlobalDGMap.cpp
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
// Description: Local to Global Base Class mapping routines
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/LocalToGlobalDGMap.h>

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>


namespace Nektar
{
    namespace MultiRegions
    {
        LocalToGlobalDGMap::LocalToGlobalDGMap():
            m_numDirichletBndPhys(0)
        {
        }

        LocalToGlobalDGMap::~LocalToGlobalDGMap()
        {
        }


        /**
         *
         */
        LocalToGlobalDGMap::LocalToGlobalDGMap( const LibUtilities::CommSharedPtr &pComm,
                                                const SpatialDomains::MeshGraph1D &graph1D,
                                                const ExpList &locExp,
                                                const GlobalSysSolnType solnType,
                                                const Array<OneD, const MultiRegions::ExpListSharedPtr> &bndCondExp,
                                                const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndCond):
                LocalToGlobalBaseMap(pComm)
        {
            int i,j;
            int cnt, vid, gid;
            int nbnd = bndCondExp.num_elements();

            // set up Local to Continuous mapping
            Array<OneD,unsigned int> vmap;
            LocalRegions::SegExpSharedPtr locSegExp;

            const boost::shared_ptr<StdRegions::StdExpansionVector> exp1D = locExp.GetExp();

            m_numGlobalBndCoeffs  = exp1D->size()+1;
            m_numGlobalCoeffs = m_numGlobalBndCoeffs;
            m_numLocalBndCoeffs = 2*exp1D->size();
            m_numLocalCoeffs = m_numLocalBndCoeffs;
            m_localToGlobalBndMap   = Array<OneD, int>(m_numLocalBndCoeffs,-1);
            m_localToGlobalBndSign  = Array<OneD, NekDouble>(m_numLocalBndCoeffs,1.0);
            m_signChange = true;
            m_solnType = solnType;
            m_staticCondLevel = 0;
            m_numPatches = exp1D->size();
            m_numLocalBndCoeffsPerPatch =  Array<OneD, unsigned int>(m_numPatches);
            m_numLocalIntCoeffsPerPatch =  Array<OneD, unsigned int>(m_numPatches);
            for(i = 0; i < m_numPatches; ++i)
            {
                m_numLocalBndCoeffsPerPatch[i] = (unsigned int) (*exp1D)[i]->NumDGBndryCoeffs();
                m_numLocalIntCoeffsPerPatch[i] = (unsigned int) 0;
            }

            map<int, int> MeshVertToLocalVert;

            // Order the Dirichlet vertices first.
            gid = 0;
            for(i = 0; i < nbnd; i++)
            {
                if(bndCond[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                {
                    m_numDirichletBndPhys++;
                    vid = ((bndCondExp[i])->GetVertex())->GetVid();

                    MeshVertToLocalVert[vid] = gid++;
                }
            }

            // set up simple map based on vertex and edge id's
            cnt = 0;
            for(i = 0; i < exp1D->size(); ++i)
            {
                if(locSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>((*exp1D)[i]))
                {
                    locSegExp->GetBoundaryMap(vmap);

                    for(j = 0; j < locSegExp->GetNverts(); ++j)
                    {
                        vid = (locSegExp->GetGeom1D())->GetVid(j);

                        if(MeshVertToLocalVert.count(vid) == 0)
                        {
                            MeshVertToLocalVert[vid] = gid++;
                        }

                        m_localToGlobalBndMap[cnt+j] =
                            MeshVertToLocalVert.find(vid)->second;
                    }
                    cnt += locSegExp->NumBndryCoeffs();
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a segment expansion failed");
                }
            }

            // Set up boundary mapping
            m_bndCondCoeffsToGlobalCoeffsMap = Array<OneD, int >(nbnd);
            m_bndCondCoeffsToGlobalCoeffsSign = Array<OneD, NekDouble >(nbnd,1.0);
            m_numLocalDirBndCoeffs = 0;
            m_numDirichletBndPhys = 0;

            for(i = 0; i < nbnd; ++i)
            {
                vid = ((bndCondExp[i])->GetVertex())->GetVid();
                m_bndCondCoeffsToGlobalCoeffsMap[i] = MeshVertToLocalVert.find(vid)->second;

                if(bndCond[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                {
                    m_numLocalDirBndCoeffs += 1;
                    m_numDirichletBndPhys  += 1;
                }
            }

            m_numGlobalDirBndCoeffs = m_numLocalDirBndCoeffs;
            CalculateBndSystemBandWidth();
        }


        /**
         *
         */
        LocalToGlobalDGMap::LocalToGlobalDGMap(const LibUtilities::CommSharedPtr &pComm,
                                               SpatialDomains::MeshGraph2D &graph2D,
                                               const ExpList1DSharedPtr &trace,
                                               const ExpList &locExp,
                                               const GlobalSysSolnType solnType,
                                               const Array<OneD, MultiRegions::ExpListSharedPtr> &bndCondExp,
                                               const Array<OneD, SpatialDomains::BoundaryConditionShPtr> &bndCond,
                                               const map<int,int> &periodicEdges) :
                LocalToGlobalBaseMap(pComm)
        {


            int i,j,k,cnt,eid, id, id1, order_e,gid;
            int ntrace_exp = trace->GetExpSize();
            int nbnd = bndCondExp.num_elements();
            LocalRegions::SegExpSharedPtr  locSegExp,locSegExp1;
            LocalRegions::QuadExpSharedPtr locQuadExp;
            LocalRegions::TriExpSharedPtr  locTriExp;
            SpatialDomains::Geometry1DSharedPtr SegGeom;

            const boost::shared_ptr<StdRegions::StdExpansionVector> exp2D = locExp.GetExp();
            int nel        = exp2D->size();

            map<int, int> MeshEdgeId;

            m_signChange = true;

            // determine mapping from geometry edges to trace
            for(i = 0; i < ntrace_exp; ++i)
            {
                if(locSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(trace->GetExp(i)))
                {
                    id = (locSegExp->GetGeom1D())->GetEid();

                    if(periodicEdges.count(id) > 0)
                    {
                        if(MeshEdgeId.count(id) == 0)
                        {
                            id1 = periodicEdges.find(id)->second;
                            MeshEdgeId[id] = i;
                            MeshEdgeId[id1] = i;
                        }
                    }
                    else
                    {
                        MeshEdgeId[id] = i;
                    }
                }
                else
                {
                    ASSERTL0(false,"Dynamics cast to segment expansion failed");
                }
            }

            // Count total number of edges
            cnt = 0;
            for(i = 0; i < nel; ++i)
            {
                cnt += (*exp2D)[i]->GetNedges();
            }

            Array<OneD, StdRegions::StdExpansion1DSharedPtr> edgemap(cnt);
            m_elmtToTrace = Array<OneD, Array<OneD,StdRegions::StdExpansion1DSharedPtr> >(nel);

            // set up edge expansions links;
            cnt = 0;
            for(i = 0; i < nel; ++i)
            {
                m_elmtToTrace[i] = edgemap + cnt;

                if(locQuadExp = boost::dynamic_pointer_cast<LocalRegions::QuadExp>((*exp2D)[i]))
                {
                    for(j = 0; j < locQuadExp->GetNedges(); ++j)
                    {
                        SegGeom = (locQuadExp->GetGeom2D())->GetEdge(j);

                        id = SegGeom->GetEid();

                        if(MeshEdgeId.count(id) > 0)
                        {
                            m_elmtToTrace[i][j] = boost::dynamic_pointer_cast< LocalRegions::SegExp> ((*trace).GetExp(MeshEdgeId.find(id)->second));

                        }
                        else
                        {
                            ASSERTL0(false,"Failed to find edge map");
                        }
                    }
                }
                else if(locTriExp = boost::dynamic_pointer_cast<LocalRegions::TriExp>((*exp2D)[i]))
                {
                    for(j = 0; j < locTriExp->GetNedges(); ++j)
                    {
                        SegGeom = (locTriExp->GetGeom2D())->GetEdge(j);

                        id = SegGeom->GetEid();

                        if(MeshEdgeId.count(id) > 0)
                        {
                            m_elmtToTrace[i][j] = boost::dynamic_pointer_cast< LocalRegions::SegExp> ((*trace).GetExp((MeshEdgeId.find(id))->second));

                        }
                        else
                        {
                            ASSERTL0(false,"Failed to find edge map");
                        }
                    }

                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a local 2D expansion failed");
                }
                cnt += (*exp2D)[i]->GetNedges();
            }

            // Set up boundary mapping
            cnt = 0;
            for(i = 0; i < nbnd; ++i)
            {
                cnt += bndCondExp[i]->GetExpSize();
            }

#if OLDMAP
            m_bndCondCoeffsToGlobalCoeffsMap = Array<OneD,int >(cnt);
#endif
            m_bndExpAdjacentOrient = Array<OneD, AdjacentTraceOrientation > (cnt);
            m_numLocalDirBndCoeffs = 0;
            m_numDirichletBndPhys  = 0;

            cnt = 0;
            for(i = 0; i < bndCondExp.num_elements(); ++i)
            {
                for(j = 0; j < bndCondExp[i]->GetExpSize(); ++j)
                {
                    if(locSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(bndCondExp[i]->GetExp(j)))
                    {
                        SegGeom = locSegExp->GetGeom1D();
                        id = SegGeom->GetEid();

#if OLDMAP
                        id = SegGeom->GetEid();
                        if(MeshEdgeId.count(id) > 0)
                        {
                            m_bndCondCoeffsToGlobalCoeffsMap[cnt+j] = MeshEdgeId.find(id)->second;
                        }
                        else
                        {
                            ASSERTL0(false,"Failed to find edge map");
                        }
#endif
                        // Check to see which way boundary edge is
                        // orientated with respect to connecting
                        // element counter-clockwise convention.

                        SpatialDomains::ElementEdgeVectorSharedPtr con_elmt
                            = graph2D.GetElementsFromEdge(SegGeom);

                        if((boost::dynamic_pointer_cast<SpatialDomains::Geometry2D>((*con_elmt)[0]->m_Element))->GetEorient((*con_elmt)[0]->m_EdgeIndx) == StdRegions::eForwards)
                        {
                            m_bndExpAdjacentOrient[cnt+j] = eAdjacentEdgeIsForwards;
                        }
                        else
                        {
                            m_bndExpAdjacentOrient[cnt+j] = eAdjacentEdgeIsBackwards;
                        }

                    }
                    else
                    {
                        ASSERTL0(false,"dynamic cast to a local Segment expansion failed");
                    }

                    if(bndCond[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                    {
                        m_numLocalDirBndCoeffs  += locSegExp->GetNcoeffs();
                        m_numDirichletBndPhys   += locSegExp->GetTotPoints();
                    }

                }
                cnt += j;
            }

            // Set up integer mapping array and sign change for each
            // degree of freedom + initialise some more data members
            m_solnType = solnType;
            m_staticCondLevel = 0;
            m_numPatches = nel;
            m_numLocalBndCoeffsPerPatch =  Array<OneD, unsigned int>(nel);
            m_numLocalIntCoeffsPerPatch =  Array<OneD, unsigned int>(nel);
            int nbndry = 0;
            for(i = 0; i < nel; ++i) // count number of elements in array
            {
                eid = locExp.GetOffset_Elmt_Id(i);
                nbndry += (*exp2D)[eid]->NumDGBndryCoeffs();
                m_numLocalBndCoeffsPerPatch[i] = (unsigned int) (*exp2D)[eid]->NumDGBndryCoeffs();
                m_numLocalIntCoeffsPerPatch[i] = (unsigned int) 0;
            }

            m_numGlobalDirBndCoeffs = m_numLocalDirBndCoeffs;

            m_numLocalBndCoeffs = nbndry;
            m_numLocalCoeffs = nbndry;
            m_localToGlobalBndMap  = Array<OneD, int > (nbndry);
            m_localToGlobalBndSign = Array<OneD, NekDouble > (nbndry,1);

            // Set up array for potential mesh optimsation
            Array<OneD,int> TraceElmtGid(ntrace_exp,-1);
            int nDir = 0;
            cnt = 0;

            // We are now going to construct a graph of the mesh
            // which can be reordered depending on the type of solver we would
            // like to use.
            typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> BoostGraph;
            typedef boost::graph_traits<BoostGraph>::vertex_descriptor BoostVertex;

            BoostGraph boostGraphObj;
            int trace_id,trace_id1;

            // make trace edge renumbering map where first solved
            // edge starts at 0 so we can set up graph.
            for(i = 0; i < ntrace_exp; ++i)
            {
                if(trace->GetCoeff_Offset(i) >= m_numLocalDirBndCoeffs)
                {
                    // Initial put in element ordering (starting
                    // from zero) into TraceElmtGid
                    boost::add_vertex(boostGraphObj);
                    TraceElmtGid[i] = cnt++;
                }
                else
                {
                    // Use existing offset for Dirichlet edges
                    TraceElmtGid[i] = trace->GetCoeff_Offset(i);
                    nDir++;
                }
            }

            // Set up boost Graph
            for(i = 0; i < nel; ++i)
            {
                eid = locExp.GetOffset_Elmt_Id(i);
                nbndry += (*exp2D)[eid]->NumDGBndryCoeffs();

                for(j = 0; j < (*exp2D)[eid]->GetNedges(); ++j)
                {
                    locSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(m_elmtToTrace[eid][j]);
                    SegGeom = locSegExp->GetGeom1D();

                    // Add edge to boost graph for non-Dirichlet Boundary
                    id  = SegGeom->GetEid();
                    trace_id = MeshEdgeId.find(id)->second;
                    if(trace->GetCoeff_Offset(trace_id) >= m_numLocalDirBndCoeffs)
                    {
                        for(k = j+1; k < (*exp2D)[eid]->GetNedges(); ++k)
                        {
                            locSegExp1 = boost::dynamic_pointer_cast<LocalRegions::SegExp>(m_elmtToTrace[eid][k]);
                            SegGeom = locSegExp1->GetGeom1D();

                            id1  = SegGeom->GetEid();
                            trace_id1 = MeshEdgeId.find(id1)->second;
                            if(trace->GetCoeff_Offset(trace_id1)
                               >= m_numLocalDirBndCoeffs)
                            {
                                boost::add_edge( (size_t) TraceElmtGid[trace_id], (size_t) TraceElmtGid[trace_id1], boostGraphObj);
                            }
                        }
                    }
                }
            }


            int nGraphVerts = ntrace_exp-nDir;
            Array<OneD, int> perm(nGraphVerts);
            Array<OneD, int> iperm(nGraphVerts);
            BottomUpSubStructuredGraphSharedPtr bottomUpGraph;
            Array<OneD, int> vwgts(nGraphVerts);
            for(i = 0; i < nGraphVerts; ++i)
            {
                vwgts[i] = trace->GetExp(i+nDir)->GetNcoeffs();
            }

            if(nGraphVerts)
            {
                switch(solnType)
                {
                case eDirectFullMatrix:
                case eIterativeFull:
                    {
                        NoReordering(boostGraphObj,perm,iperm);
                    }
                    break;
                case eDirectStaticCond:
                    {
                        CuthillMckeeReordering(boostGraphObj,perm,iperm);
                    }
                    break;
                case eDirectMultiLevelStaticCond:
                    {
                        MultiLevelBisectionReordering(boostGraphObj,vwgts,perm,iperm,bottomUpGraph);
                    }
                    break;
                default:
                    {
                        ASSERTL0(false,"Unrecognised solution type");
                    }
                }
            }

            // Recast the permutation so that it can be
            // used as a map Form old trace edge ID to new trace
            // edge ID
            cnt = m_numLocalDirBndCoeffs;
            for(i = 0; i < ntrace_exp-nDir; ++i)
            {
                TraceElmtGid[perm[i]+nDir]=cnt;
                cnt += trace->GetExp(perm[i]+nDir)->GetNcoeffs();
            }

            // Now have trace edges Gid position
            nbndry = cnt = 0;
            for(i = 0; i < nel; ++i)
            {
                // order list according to m_offset_elmt_id details in
                // Exp2D so that triangules are listed first and then
                // quads
                eid = locExp.GetOffset_Elmt_Id(i);
                nbndry += (*exp2D)[eid]->NumDGBndryCoeffs();

                for(j = 0; j < (*exp2D)[eid]->GetNedges(); ++j)
                {
                    locSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(m_elmtToTrace[eid][j]);
                    SegGeom = locSegExp->GetGeom1D();

                    id  = SegGeom->GetEid();
                    gid = TraceElmtGid[MeshEdgeId.find(id)->second];

                    //Peter order_e = locSegExp->GetNcoeffs();
                    order_e = (*exp2D)[eid]->GetEdgeNcoeffs(j);

                    if((*exp2D)[eid]->GetEorient(j) == StdRegions::eForwards)
                    {
                        for(k = 0; k < order_e; ++k)
                        {
                            m_localToGlobalBndMap[k+cnt] = gid + k;
                        }
                    }
                    else // backwards orientated
                    {
                        switch(locSegExp->GetBasisType(0))
                        {
                        case LibUtilities::eModified_A:
                            // reverse vertex order
                            m_localToGlobalBndMap[cnt] = gid + 1;
                            m_localToGlobalBndMap[cnt+1] = gid;
                            for(k = 2; k < order_e; ++k)
                            {
                                m_localToGlobalBndMap[k+cnt] = gid + k;
                            }

                            // negate odd modes
                            for(k = 3; k < order_e; k+=2)
                            {
                                m_localToGlobalBndSign[cnt+k] = -1.0;
                            }


                            break;
                        case LibUtilities::eGLL_Lagrange:
                            // reverse  order
                            for(k = 0; k < order_e; ++k)
                            {
                                m_localToGlobalBndMap[cnt+order_e-k-1] = gid + k;
                            }
                            break;
                        default:
                            ASSERTL0(false,"Boundary type not permitted");

                        }
                    }
                    cnt += order_e;
                }
            }

            // set up m_bndCondCoeffsToGlobalCoeffsMap to align with map
            cnt = 0;
            for(i = 0; i < nbnd; ++i)
            {
                cnt += bndCondExp[i]->GetNcoeffs();
            }

            m_bndCondCoeffsToGlobalCoeffsMap = Array<OneD,int >(cnt);

            // Number of boundary expansions
            int nbndexp = 0;
            for(cnt = i = 0; i < nbnd; ++i)
            {
                for(j = 0; j < bndCondExp[i]->GetExpSize(); ++j)
                {
                    if(locSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(bndCondExp[i]->GetExp(j)))
                    {
                        nbndexp++;
                        SegGeom = locSegExp->GetGeom1D();
                        id      = SegGeom->GetEid();
                        gid     = TraceElmtGid[MeshEdgeId.find(id)->second];

                        order_e = locSegExp->GetNcoeffs();

                        // Since boundary information is defined to be
                        // aligned with the geometry just use forward
                        // defintiion for gid's
                        for(k = 0; k < order_e; ++k)
                        {
                            m_bndCondCoeffsToGlobalCoeffsMap[cnt++] = gid + k;
                        }
                    }
                }
            }

            m_numGlobalBndCoeffs = trace->GetNcoeffs();
            m_numGlobalCoeffs = m_numGlobalBndCoeffs;

            CalculateBndSystemBandWidth();

            if( (solnType == eDirectMultiLevelStaticCond) && nGraphVerts )
            {
                if(m_staticCondLevel < (bottomUpGraph->GetNlevels()-1))
                {

                    Array<OneD, int> vwgts_perm(nGraphVerts);
                    
                    for(int i = 0; i < nGraphVerts; i++)
                    {
                        vwgts_perm[i] = vwgts[perm[i]];
                    }

                    bottomUpGraph->ExpandGraphWithVertexWeights(vwgts_perm);

                    m_nextLevelLocalToGlobalMap = MemoryManager<LocalToGlobalBaseMap>::
                        AllocateSharedPtr(this,bottomUpGraph);
                }
            }

            cnt = 0;
            m_bndCondTraceToGlobalTraceMap = Array<OneD, int >(nbndexp);
            for(i = 0; i < bndCondExp.num_elements(); ++i)
            {
                for(j = 0; j < bndCondExp[i]->GetExpSize(); ++j)
                {
                    if(locSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(bndCondExp[i]->GetExp(j)))
                    {
                        SegGeom = locSegExp->GetGeom1D();
                        id = SegGeom->GetEid();

                        m_bndCondTraceToGlobalTraceMap[cnt++] = MeshEdgeId.find(id)->second;
                    }
                }
            }

            // Now set up mapping from global coefficients to universal.
            SetUpUniversalDGMap(locExp);

            // Initialise GSlib and populate the unique map.
            Nektar::Array<OneD, long> tmp(m_globalToUniversalBndMap.num_elements());
            for (unsigned int i = 0; i < m_globalToUniversalBndMap.num_elements(); ++i)
            {
                tmp[i] = m_globalToUniversalBndMap[i];
            }
            m_gsh = Gs::Init(tmp, pComm);
            Gs::Unique(tmp, pComm);
            for (unsigned int i = 0; i < m_globalToUniversalBndMap.num_elements(); ++i)
            {
                m_globalToUniversalBndMapUnique[i] = (tmp[i] >= 0 ? 1 : 0);
            }

        }


        /**
         * Constructs a mapping between the process-local global numbering and
         * a universal numbering of the trace space expansion. The universal
         * numbering is defined by the mesh edge IDs to enforce consistency
         * across processes.
         * @param       locExp  List of local elemental expansions.
         * @todo        Update to support 1D and 3D DG expansions.
         */
        void LocalToGlobalDGMap::SetUpUniversalDGMap(const ExpList &locExp)
        {
            LocalRegions::SegExpSharedPtr locSegExp;

            int eid = 0;
            int cnt = 0;
            int i,j,k;
            int id = 0;
            int order_e = 0;
            int vGlobalId = 0;
            int maxEdgeDof = 0;
            int dof = 0;
            const StdRegions::StdExpansionVector &locExpVector = *(locExp.GetExp());

            // Initialise the global to universal maps.
            m_globalToUniversalBndMap = Nektar::Array<OneD, int>(m_numGlobalBndCoeffs, -1);
            m_globalToUniversalBndMapUnique = Nektar::Array<OneD, int>(m_numGlobalBndCoeffs, -1);

            // Loop over all the elements in the domain and compute max edge
            // DOF. Reduce across all processes to get universal maximum.
            for(i = 0; i < locExpVector.size(); ++i)
            {
                // Loop over all edges of element i
                for(j = 0; j < locExpVector[i]->GetNedges(); ++j)
                {
                    dof = locExpVector[i]->GetEdgeNcoeffs(j);
                    maxEdgeDof = (dof > maxEdgeDof ? dof : maxEdgeDof);
                }
            }
            m_comm->AllReduce(maxEdgeDof, LibUtilities::ReduceMax);

            // Now have trace edges Gid position
            cnt = 0;
            for(i = 0; i < locExpVector.size(); ++i)
            {
                // order list according to m_offset_elmt_id details in
                // Exp2D so that triangules are listed first and then
                // quads
                eid = locExp.GetOffset_Elmt_Id(i);

                // Populate mapping for each edge of the element.
                for(j = 0; j < locExpVector[eid]->GetNedges(); ++j)
                {
                    locSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(m_elmtToTrace[eid][j]);

                    id  = locSegExp->GetGeom1D()->GetEid();
                    order_e = locExpVector[eid]->GetEdgeNcoeffs(j);

                    for(k = 0; k < order_e; ++k)
                    {
                        vGlobalId = m_localToGlobalBndMap[k+cnt];
                        m_globalToUniversalBndMap[vGlobalId]
                            = id * maxEdgeDof + k + 1;
                    }
                    cnt += order_e;
                }
            }
        }

        int LocalToGlobalDGMap::v_GetLocalToGlobalMap(const int i) const
        {
            return m_localToGlobalBndMap[i];
        }

        int LocalToGlobalDGMap::v_GetGlobalToUniversalMap(const int i) const
        {
            return m_globalToUniversalBndMap[i];
        }

        int LocalToGlobalDGMap::v_GetGlobalToUniversalMapUnique(const int i) const
        {
            return m_globalToUniversalBndMapUnique[i];
        }

        const Array<OneD,const int>& LocalToGlobalDGMap::v_GetLocalToGlobalMap()
        {
            return m_localToGlobalBndMap;
        }

        const Array<OneD,const int>& LocalToGlobalDGMap::v_GetGlobalToUniversalMap()
        {
            return m_globalToUniversalBndMap;
        }

        const Array<OneD,const int>& LocalToGlobalDGMap::v_GetGlobalToUniversalMapUnique()
        {
            return m_globalToUniversalBndMapUnique;
        }

        NekDouble LocalToGlobalDGMap::v_GetLocalToGlobalSign(
                    const int i) const
        {
            return GetLocalToGlobalBndSign(i);
        }

        const void LocalToGlobalDGMap::v_LocalToGlobal(
                    const Array<OneD, const NekDouble>& loc,
                          Array<OneD,       NekDouble>& global) const
        {
            AssembleBnd(loc,global);
        }

        const void LocalToGlobalDGMap::v_LocalToGlobal(
                    const NekVector<const NekDouble>& loc,
                          NekVector<      NekDouble>& global) const
        {
            AssembleBnd(loc,global);
        }

        const void LocalToGlobalDGMap::v_GlobalToLocal(
                    const Array<OneD, const NekDouble>& global,
                          Array<OneD,       NekDouble>& loc) const
        {
            GlobalToLocalBnd(global,loc);
        }

        const void LocalToGlobalDGMap::v_GlobalToLocal(
                    const NekVector<const NekDouble>& global,
                          NekVector<      NekDouble>& loc) const
        {
            GlobalToLocalBnd(global,loc);
        }

        const void LocalToGlobalDGMap::v_Assemble(
                    const Array<OneD, const NekDouble> &loc,
                          Array<OneD,       NekDouble> &global) const
        {
            AssembleBnd(loc,global);
        }

        const void LocalToGlobalDGMap::v_Assemble(
                    const NekVector<const NekDouble>& loc,
                          NekVector<      NekDouble>& global) const
        {
            AssembleBnd(loc,global);
        }

        const void LocalToGlobalDGMap::v_UniversalAssemble(
                      Array<OneD,     NekDouble>& pGlobal) const
        {
            Gs::Gather(pGlobal, Gs::gs_add, m_gsh);
        }

        const void LocalToGlobalDGMap::v_UniversalAssemble(
                      NekVector<      NekDouble>& pGlobal) const
        {
            UniversalAssemble(pGlobal.GetPtr());
        }

        const int LocalToGlobalDGMap::v_GetFullSystemBandWidth() const
        {
            return GetBndSystemBandWidth();
        }
    }


}
