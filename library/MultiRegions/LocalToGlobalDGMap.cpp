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
#include <LocalRegions/HexExp.h>
#include <LocalRegions/TetExp.h>
#include <LocalRegions/PrismExp.h>
#include <LocalRegions/PyrExp.h>

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
        LocalToGlobalDGMap::LocalToGlobalDGMap( const LibUtilities::SessionReaderSharedPtr &pSession,
                                                const SpatialDomains::MeshGraphSharedPtr &graph1D,
                                                const ExpList &locExp,
                                                const Array<OneD, const MultiRegions::ExpListSharedPtr> &bndCondExp,
                                                const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndCond):
                LocalToGlobalBaseMap(pSession)
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
			
			
            // Check to see which way boundary point is
            // orientated with respect to convention (croth)
            m_bndExpAdjacentOrient = Array<OneD, AdjacentTraceOrientation > (nbnd);
            
            for (int i=0; i<nbnd; i++)
            {
                vid = ((bndCondExp[i])->GetVertex())->GetVid();
                //cout << "VID = "<<vid<<endl;
		
                if(vid == 0)
                {
                    m_bndExpAdjacentOrient[i] = eAdjacentEdgeIsBackwards;
                }
                else
                {
                    m_bndExpAdjacentOrient[i] = eAdjacentEdgeIsForwards;
                }
            }
        }


        /**
         *
         */
        LocalToGlobalDGMap::LocalToGlobalDGMap(const LibUtilities::SessionReaderSharedPtr &pSession,
                                               const SpatialDomains::MeshGraphSharedPtr &graph2D,
                                               const ExpList1DSharedPtr &trace,
                                               const ExpList &locExp,
                                               const Array<OneD, MultiRegions::ExpListSharedPtr> &bndCondExp,
                                               const Array<OneD, SpatialDomains::BoundaryConditionShPtr> &bndCond,
                                               const map<int,int> &periodicEdges) :
                LocalToGlobalBaseMap(pSession)
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
                            id1 = abs(periodicEdges.find(id)->second);
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
                            = boost::dynamic_pointer_cast<SpatialDomains::MeshGraph2D>(graph2D)->GetElementsFromEdge(SegGeom);

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
                switch(m_solnType)
                {
                case eDirectFullMatrix:
                case eIterativeFull:
                    {
                        NoReordering(boostGraphObj,perm,iperm);
                    }
                    break;
                case eDirectStaticCond:
                case eIterativeStaticCond:
                    {
                        CuthillMckeeReordering(boostGraphObj,perm,iperm);
                    }
                    break;
                case eDirectMultiLevelStaticCond:
                case eIterativeMultiLevelStaticCond:
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

            if( m_solnType == eDirectMultiLevelStaticCond && nGraphVerts )
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
            m_gsh = Gs::Init(tmp, m_comm);
            Gs::Unique(tmp, m_comm);
            for (unsigned int i = 0; i < m_globalToUniversalBndMap.num_elements(); ++i)
            {
                m_globalToUniversalBndMapUnique[i] = (tmp[i] >= 0 ? 1 : 0);
            }

        }

        /**
         * Constructor for trace map for three-dimensional expansion.
         */
        LocalToGlobalDGMap::LocalToGlobalDGMap(const LibUtilities::SessionReaderSharedPtr &pSession,
                const SpatialDomains::MeshGraphSharedPtr &graph3D,
                const ExpList2DSharedPtr &trace,
                const ExpList &locExp,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &bndCondExp,
                const Array<OneD, SpatialDomains::BoundaryConditionShPtr> &bndCond,
                const map<int,int> &periodicFaces):
                LocalToGlobalBaseMap(pSession)
        {


            int i,j,k,cnt,eid, id, id1, order_e,gid;
            int ntrace_exp = trace->GetExpSize();
            int nbnd = bndCondExp.num_elements();
            LocalRegions::QuadExpSharedPtr locQuadExp, locQuadExp1;
            LocalRegions::TriExpSharedPtr  locTriExp, locTriExp1;
            LocalRegions::HexExpSharedPtr locHexExp;
            LocalRegions::PrismExpSharedPtr locPrismExp;
            LocalRegions::PyrExpSharedPtr locPyrExp;
            LocalRegions::TetExpSharedPtr locTetExp;
            SpatialDomains::Geometry2DSharedPtr FaceGeom;

            const boost::shared_ptr<StdRegions::StdExpansionVector> exp3D = locExp.GetExp();
            int nel = exp3D->size();

            map<int, int> MeshFaceId;

            m_signChange = true;

            // determine mapping from geometry edges to trace
            for(i = 0; i < ntrace_exp; ++i)
            {
                //quad face
                if(locQuadExp = boost::dynamic_pointer_cast<LocalRegions::QuadExp>(trace->GetExp(i)))
                {
                    id = (locQuadExp->GetGeom2D())->GetFid();
                }
                //tri face
                else if(locTriExp = boost::dynamic_pointer_cast<LocalRegions::TriExp>(trace->GetExp(i)))
                {
                    id = (locTriExp->GetGeom2D())->GetFid();
                }
                else
                {
                    ASSERTL0(false,"Dynamic cast to face expansion failed");
                }
				
                if(periodicFaces.count(id) > 0)
                {
                    if(MeshFaceId.count(id) == 0)
                    {
                        id1 = periodicFaces.find(id)->second;
                        MeshFaceId[id] = i;
                        MeshFaceId[id1] = i;
                    }
                }
                else
                {
                    MeshFaceId[id] = i;
                }
            }

            // Count total number of faces
            cnt = 0;
            for(i = 0; i < nel; ++i)
            {
                cnt += (*exp3D)[i]->GetNfaces();
            }

            Array<OneD, StdRegions::StdExpansion2DSharedPtr> facemap(cnt);
            m_elmtToFace = Array<OneD, Array<OneD,StdRegions::StdExpansion2DSharedPtr> >(nel);

            // set up face expansions links;
            cnt = 0;
            for(i = 0; i < nel; ++i)
            {
                m_elmtToFace[i] = facemap + cnt;
				//if Hex expansion
                if(locHexExp = boost::dynamic_pointer_cast<LocalRegions::HexExp>((*exp3D)[i]))
                {
                    for(j = 0; j < locHexExp->GetNfaces(); ++j)
                    {
                        FaceGeom = (locHexExp->GetGeom3D())->GetFace(j);

                        id = FaceGeom->GetFid();

                        if(MeshFaceId.count(id) > 0)
                        {
                            if(FaceGeom->GetGeomShapeType() == SpatialDomains::eQuadrilateral)
                            {
                                m_elmtToFace[i][j] = boost::dynamic_pointer_cast< LocalRegions::QuadExp> ((*trace).GetExp(MeshFaceId.find(id)->second));
                            }
                            else if(FaceGeom->GetGeomShapeType() == SpatialDomains::eTriangle)
                            {
                                m_elmtToFace[i][j] = boost::dynamic_pointer_cast< LocalRegions::TriExp> ((*trace).GetExp(MeshFaceId.find(id)->second));
                            }
                            else
                            {
                                ASSERTL0(false,"Unknown face geometry shape type");
                            }
                        }
                        else
                        {
                            ASSERTL0(false,"Failed to find face map");
                        }
                    }
                }
                //else if Tet expansion
                else if(locTetExp = boost::dynamic_pointer_cast<LocalRegions::TetExp>((*exp3D)[i]))
                {
                    for(j = 0; j < locTetExp->GetNfaces(); ++j)
                    {
                        FaceGeom = (locTetExp->GetGeom3D())->GetFace(j);
                        
                        id = FaceGeom->GetFid();
                        
                        if(MeshFaceId.count(id) > 0)
                        {
                            if(FaceGeom->GetGeomShapeType() == SpatialDomains::eQuadrilateral)
                            {
                                m_elmtToFace[i][j] = boost::dynamic_pointer_cast< LocalRegions::QuadExp> ((*trace).GetExp(MeshFaceId.find(id)->second));
                            }
                            else if(FaceGeom->GetGeomShapeType() == SpatialDomains::eTriangle)
                            {
                                m_elmtToFace[i][j] = boost::dynamic_pointer_cast< LocalRegions::TriExp> ((*trace).GetExp(MeshFaceId.find(id)->second));
                            }
                            else
                            {
                                ASSERTL0(false,"Unknown face geometry shape type");
                            }
                        }
                        else
                        {
                            ASSERTL0(false,"Failed to find face map");
                        }
                    }
                }
                //else if Pyramid expansion
                else if(locPyrExp = boost::dynamic_pointer_cast<LocalRegions::PyrExp>((*exp3D)[i]))
                {
                    for(j = 0; j < locPyrExp->GetNfaces(); ++j)
                    {
                        FaceGeom = (locPyrExp->GetGeom3D())->GetFace(j);
                        
                        id = FaceGeom->GetFid();
                        
                        if(MeshFaceId.count(id) > 0)
                        {
                            if(FaceGeom->GetGeomShapeType() == SpatialDomains::eQuadrilateral)
                            {
                                m_elmtToFace[i][j] = boost::dynamic_pointer_cast< LocalRegions::QuadExp> ((*trace).GetExp(MeshFaceId.find(id)->second));
                            }
                            else if(FaceGeom->GetGeomShapeType() == SpatialDomains::eTriangle)
                            {
                                m_elmtToFace[i][j] = boost::dynamic_pointer_cast< LocalRegions::TriExp> ((*trace).GetExp(MeshFaceId.find(id)->second));
                            }
                            else
                            {
                                ASSERTL0(false,"Unknown face geometry shape type");
                            }
                        }
                        else
                        {
                            ASSERTL0(false,"Failed to find face map");
                        }
                    }
                }
                //else if Prism expansion
                else if(locPrismExp = boost::dynamic_pointer_cast<LocalRegions::PrismExp>((*exp3D)[i]))
                {
                    for(j = 0; j < locPrismExp->GetNfaces(); ++j)
                    {
                        FaceGeom = (locPrismExp->GetGeom3D())->GetFace(j);
                        
                        id = FaceGeom->GetFid();
                        
                        if(MeshFaceId.count(id) > 0)
                        {
                            if(FaceGeom->GetGeomShapeType() == SpatialDomains::eQuadrilateral)
                            {
                                m_elmtToFace[i][j] = boost::dynamic_pointer_cast< LocalRegions::QuadExp> ((*trace).GetExp(MeshFaceId.find(id)->second));
                            }
                            else if(FaceGeom->GetGeomShapeType() == SpatialDomains::eTriangle)
                            {
                                m_elmtToFace[i][j] = boost::dynamic_pointer_cast< LocalRegions::TriExp> ((*trace).GetExp(MeshFaceId.find(id)->second));
                            }
                            else
                            {
                                ASSERTL0(false,"Unknown face geometry shape type");
                            }
                        }
                        else
                        {
                            ASSERTL0(false,"Failed to find face map");
                        }
                    }
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a local 3D expansion failed");
                }
                cnt += (*exp3D)[i]->GetNfaces();
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
            m_bndExpAdjacentFaceOrient = Array<OneD, AdjacentFaceOrientation > (cnt);
            m_numLocalDirBndCoeffs = 0;
            m_numDirichletBndPhys  = 0;

            cnt = 0;
            for(i = 0; i < bndCondExp.num_elements(); ++i)
            {
                for(j = 0; j < bndCondExp[i]->GetExpSize(); ++j)
                {
                    //if face is quad
                    if(locQuadExp = boost::dynamic_pointer_cast<LocalRegions::QuadExp>(bndCondExp[i]->GetExp(j)))
                    {
                        FaceGeom = locQuadExp->GetGeom2D();
                        id = FaceGeom->GetFid();

#if OLDMAP
                        id = FaceGeom->GetFid();
                        if(MeshFaceId.count(id) > 0)
                        {
                            m_bndCondCoeffsToGlobalCoeffsMap[cnt+j] = MeshFaceId.find(id)->second;
                        }
                        else
                        {
                            ASSERTL0(false,"Failed to find face map");
                        }
#endif
                        // Check to see which way boundary face is
                        // orientated with respect to connecting
                        // element.

                        SpatialDomains::ElementFaceVectorSharedPtr con_elmt
                            = boost::dynamic_pointer_cast<SpatialDomains::MeshGraph3D>(graph3D)->GetElementsFromFace(FaceGeom);

                        StdRegions::FaceOrientation cur_face_orientation
                            = (boost::dynamic_pointer_cast<SpatialDomains::Geometry3D>((*con_elmt)[0]->m_Element))->GetFaceOrient((*con_elmt)[0]->m_FaceIndx);	
                        
                        switch(cur_face_orientation)
                        {
                            case StdRegions::eDir1FwdDir1_Dir2FwdDir2:
                                m_bndExpAdjacentFaceOrient[cnt+j] = eAdjacentFaceDir1FwdDir1_Dir2FwdDir2;
                                break;
                            case StdRegions::eDir1FwdDir1_Dir2BwdDir2:
                                m_bndExpAdjacentFaceOrient[cnt+j] = eAdjacentFaceDir1FwdDir1_Dir2BwdDir2;
                                break;
                            case StdRegions::eDir1BwdDir1_Dir2FwdDir2:
                                m_bndExpAdjacentFaceOrient[cnt+j] = eAdjacentFaceDir1BwdDir1_Dir2FwdDir2;
                                break;
                            case StdRegions::eDir1BwdDir1_Dir2BwdDir2:
                                m_bndExpAdjacentFaceOrient[cnt+j] = eAdjacentFaceDir1BwdDir1_Dir2BwdDir2;
                                break;
                            case StdRegions::eDir1FwdDir2_Dir2FwdDir1:
                                m_bndExpAdjacentFaceOrient[cnt+j] = eAdjacentFaceDir1FwdDir2_Dir2FwdDir1;
                                break;
                            case StdRegions::eDir1FwdDir2_Dir2BwdDir1:
                                m_bndExpAdjacentFaceOrient[cnt+j] = eAdjacentFaceDir1FwdDir2_Dir2BwdDir1;
                                break;
                            case StdRegions::eDir1BwdDir2_Dir2FwdDir1:
                                m_bndExpAdjacentFaceOrient[cnt+j] = eAdjacentFaceDir1BwdDir2_Dir2FwdDir1;
                                break;
                            case StdRegions::eDir1BwdDir2_Dir2BwdDir1:
                                m_bndExpAdjacentFaceOrient[cnt+j] = eAdjacentFaceDir1BwdDir2_Dir2BwdDir1;
                                break;
                            default:
                                ASSERTL0(false, "Unknown adjacent face orientation");
                        };
                        
                        if(bndCond[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                        {
                            m_numLocalDirBndCoeffs  += locQuadExp->GetNcoeffs();
                            m_numDirichletBndPhys   += locQuadExp->GetTotPoints();
                        }
                    }
                    //else if face is triangle
                    else if(locTriExp = boost::dynamic_pointer_cast<LocalRegions::TriExp>(bndCondExp[i]->GetExp(j)))
                    {
                        FaceGeom = locTriExp->GetGeom2D();
                        id = FaceGeom->GetFid();
                        
#if OLDMAP
                        id = FaceGeom->GetFid();
                        if(MeshFaceId.count(id) > 0)
                        {
                            m_bndCondCoeffsToGlobalCoeffsMap[cnt+j] = MeshFaceId.find(id)->second;
                        }
                        else
                        {
                            ASSERTL0(false,"Failed to find face map");
                        }
#endif
                        // Check to see which way boundary face is
                        // orientated with respect to connecting
                        // element.
                        
                        SpatialDomains::ElementFaceVectorSharedPtr con_elmt
                            = boost::dynamic_pointer_cast<SpatialDomains::MeshGraph3D>(graph3D)->GetElementsFromFace(FaceGeom);
                        StdRegions::FaceOrientation cur_face_orientation
                            = (boost::dynamic_pointer_cast<SpatialDomains::Geometry3D>((*con_elmt)[0]->m_Element))->GetFaceOrient((*con_elmt)[0]->m_FaceIndx);	
                        
                        switch(cur_face_orientation)
                        {
                            case StdRegions::eDir1FwdDir1_Dir2FwdDir2:
                                m_bndExpAdjacentFaceOrient[cnt+j] = eAdjacentFaceDir1FwdDir1_Dir2FwdDir2;
                                break;
                            case StdRegions::eDir1FwdDir1_Dir2BwdDir2:
                                m_bndExpAdjacentFaceOrient[cnt+j] = eAdjacentFaceDir1FwdDir1_Dir2BwdDir2;
                                break;
                            case StdRegions::eDir1BwdDir1_Dir2FwdDir2:
                                m_bndExpAdjacentFaceOrient[cnt+j] = eAdjacentFaceDir1BwdDir1_Dir2FwdDir2;
                                break;
                            case StdRegions::eDir1BwdDir1_Dir2BwdDir2:
                                m_bndExpAdjacentFaceOrient[cnt+j] = eAdjacentFaceDir1BwdDir1_Dir2BwdDir2;
                                break;
                            case StdRegions::eDir1FwdDir2_Dir2FwdDir1:
                                m_bndExpAdjacentFaceOrient[cnt+j] = eAdjacentFaceDir1FwdDir2_Dir2FwdDir1;
                                break;
                            case StdRegions::eDir1FwdDir2_Dir2BwdDir1:
                                m_bndExpAdjacentFaceOrient[cnt+j] = eAdjacentFaceDir1FwdDir2_Dir2BwdDir1;
                                break;
                            case StdRegions::eDir1BwdDir2_Dir2FwdDir1:
                                m_bndExpAdjacentFaceOrient[cnt+j] = eAdjacentFaceDir1BwdDir2_Dir2FwdDir1;
                                break;
                            case StdRegions::eDir1BwdDir2_Dir2BwdDir1:
                                m_bndExpAdjacentFaceOrient[cnt+j] = eAdjacentFaceDir1BwdDir2_Dir2BwdDir1;
                                break;
                            default:
                                ASSERTL0(false, "Unknown adjacent face orientation");
                        };
			
                        if(bndCond[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                        {
                            m_numLocalDirBndCoeffs  += locTriExp->GetNcoeffs();
                            m_numDirichletBndPhys   += locTriExp->GetTotPoints();
                        }
                    }
                    else
                    {
                        ASSERTL0(false,"dynamic cast to a local face expansion failed");
                    }
                    
                    
                }
                cnt += j;
            }

//--------------> NumDGBndryCoeffs is currently only implemented for Hex elements

            // Set up integer mapping array and sign change for each
            // degree of freedom + initialise some more data members
            m_staticCondLevel = 0;
            m_numPatches = nel;
            m_numLocalBndCoeffsPerPatch =  Array<OneD, unsigned int>(nel);
            m_numLocalIntCoeffsPerPatch =  Array<OneD, unsigned int>(nel);
            int nbndry = 0;
            for(i = 0; i < nel; ++i) // count number of elements in array
            {
                eid = locExp.GetOffset_Elmt_Id(i);
                nbndry += (*exp3D)[eid]->NumDGBndryCoeffs();
                m_numLocalBndCoeffsPerPatch[i] = (unsigned int) (*exp3D)[eid]->NumDGBndryCoeffs();
                m_numLocalIntCoeffsPerPatch[i] = (unsigned int) 0;
            }

            m_numGlobalDirBndCoeffs = m_numLocalDirBndCoeffs;

            m_numLocalBndCoeffs = nbndry;
            m_numLocalCoeffs = nbndry;
            m_localToGlobalBndMap  = Array<OneD, int > (nbndry);
            m_localToGlobalBndSign = Array<OneD, NekDouble > (nbndry,1);

            // Set up array for potential mesh optimsation
            Array<OneD,int> FaceElmtGid(ntrace_exp,-1);
            int nDir = 0;
            cnt = 0;

            // We are now going to construct a graph of the mesh
            // which can be reordered depending on the type of solver we would
            // like to use.
            typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> BoostGraph;
            typedef boost::graph_traits<BoostGraph>::vertex_descriptor BoostVertex;

            BoostGraph boostGraphObj;
            int face_id,face_id1;

            // make trace face renumbering map where first solved
            // face starts at 0 so we can set up graph.
            for(i = 0; i < ntrace_exp; ++i)
            {
                if(trace->GetCoeff_Offset(i) >= m_numLocalDirBndCoeffs)
                {
                    // Initial put in element ordering (starting
                    // from zero) into FaceElmtGid
                    boost::add_vertex(boostGraphObj);
                    FaceElmtGid[i] = cnt++;
                }
                else
                {
                    // Use existing offset for Dirichlet edges
                    FaceElmtGid[i] = trace->GetCoeff_Offset(i);
                    nDir++;
                }
            }

            // Set up boost Graph
            for(i = 0; i < nel; ++i)
            {
                eid = locExp.GetOffset_Elmt_Id(i);

                for(j = 0; j < (*exp3D)[eid]->GetNfaces(); ++j)
                {
                    //if face is quad
                    if(locQuadExp = boost::dynamic_pointer_cast<LocalRegions::QuadExp>(m_elmtToFace[eid][j]))
                    {
                        FaceGeom = locQuadExp->GetGeom2D();
                    }
                    //else if face is triangle
                    else if(locTriExp = boost::dynamic_pointer_cast<LocalRegions::TriExp>(m_elmtToFace[eid][j]))
                    {
                        FaceGeom = locTriExp->GetGeom2D();
                    }
                    else
                    {
                        ASSERTL0(false,"dynamic cast to a local face expansion failed");
                    }
                    
                    // Add face to boost graph for non-Dirichlet Boundary
                    id = FaceGeom->GetFid();
                    face_id = MeshFaceId.find(id)->second;
                    if(trace->GetCoeff_Offset(face_id) >= m_numLocalDirBndCoeffs)
                    {
                        for(k = j+1; k < (*exp3D)[eid]->GetNfaces(); ++k)
                        {
                            //if face is quad
                            if(locQuadExp1 = boost::dynamic_pointer_cast<LocalRegions::QuadExp>(m_elmtToFace[eid][k]))
                            {
                                FaceGeom = locQuadExp1->GetGeom2D();
                            }
                            //else if face is triangle
                            else if(locTriExp1 = boost::dynamic_pointer_cast<LocalRegions::TriExp>(m_elmtToFace[eid][k]))
                            {
                                FaceGeom = locTriExp1->GetGeom2D();
                            }
                            else
                            {
                                ASSERTL0(false,"dynamic cast to a local face expansion failed");
                            }
                            id1  = FaceGeom->GetFid();
                            face_id1 = MeshFaceId.find(id1)->second;
                            if(trace->GetCoeff_Offset(face_id1)
                               >= m_numLocalDirBndCoeffs)
                            {
                                boost::add_edge( (size_t) FaceElmtGid[face_id], (size_t) FaceElmtGid[face_id1], boostGraphObj);
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
                switch(m_solnType)
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
            // used as a map Form old trace face ID to new trace
            // face ID
            cnt = m_numLocalDirBndCoeffs;
            for(i = 0; i < ntrace_exp-nDir; ++i)
            {
                FaceElmtGid[perm[i]+nDir]=cnt;
                cnt += trace->GetExp(perm[i]+nDir)->GetNcoeffs();
            }

            // Now have trace edges Gid position
            cnt = 0;
            for(i = 0; i < nel; ++i)
            {
                // order list according to m_offset_elmt_id details in Exp3D
                eid = locExp.GetOffset_Elmt_Id(i);
                
                for(j = 0; j < (*exp3D)[eid]->GetNfaces(); ++j)
                {
                    //if face is quad
                    if(locQuadExp = boost::dynamic_pointer_cast<LocalRegions::QuadExp>(m_elmtToFace[eid][j]))
                    {
                        FaceGeom = locQuadExp->GetGeom2D();
                    }
                    //else if face is triangle
                    else if(locTriExp = boost::dynamic_pointer_cast<LocalRegions::TriExp>(m_elmtToFace[eid][j]))
                    {
                        FaceGeom = locTriExp->GetGeom2D();
                    }
                    else
                    {
                        ASSERTL0(false,"dynamic cast to a local face expansion failed");
                    }

                    id  = FaceGeom->GetFid();
                    gid = FaceElmtGid[MeshFaceId.find(id)->second];
                    order_e = (*exp3D)[eid]->GetFaceNcoeffs(j);
                    
                    Array<OneD, unsigned int> elmMap1 (order_e);
                    Array<OneD,          int> elmSign1(order_e);
                    Array<OneD, unsigned int> elmMap2 (order_e);
                    Array<OneD,          int> elmSign2(order_e);
                    StdRegions::FaceOrientation fo = (*exp3D)[eid]->GetFaceOrient(j);
                    
                    // Construct mapping which will permute global IDs
                    // according to face orientations. 
                    (*exp3D)[eid]->GetFaceToElementMap(j,fo,elmMap1,elmSign2);
                    (*exp3D)[eid]->GetFaceToElementMap(j,(StdRegions::FaceOrientation)0,elmMap2,elmSign2);
                    
                    std::map<int,int> orientMap;
                    
                    for (k = 0; k < elmMap1.num_elements(); ++k)
                    {
                        // Find the elemental co-efficient in the original
                        // mapping.
                        int idx = -1;
                        for (int l = 0; l < elmMap2.num_elements(); ++l)
                        {
                            if (elmMap1[k] == elmMap2[l])
                            {
                                idx = l;
                                break;
                            }
                        }
                        
                        ASSERTL2(idx != -1, "Problem with face to element map!");
                        orientMap[k] = idx;
                    }
                    
                    for(k = 0; k < order_e; ++k)
                    {
                        m_localToGlobalBndMap [k+cnt] = gid + orientMap[k];
                        m_localToGlobalBndSign[k+cnt] = elmSign1[orientMap[k]];
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
                    //if face is quad
                    if(locQuadExp = boost::dynamic_pointer_cast<LocalRegions::QuadExp>(bndCondExp[i]->GetExp(j)))
                    {
                        nbndexp++;
                        FaceGeom = locQuadExp->GetGeom2D();
                        id      = FaceGeom->GetFid();
                        gid     = FaceElmtGid[MeshFaceId.find(id)->second];
                        
                        order_e = locQuadExp->GetNcoeffs();
                        
                        // Since boundary information is defined to be
                        // aligned with the geometry just use forward/forward
                        // (both coordinate directions) defintiion for gid's
                        for(k = 0; k < order_e; ++k)
                        {
                            m_bndCondCoeffsToGlobalCoeffsMap[cnt++] = gid + k;
                        }
                    }
                    //else if face is triangle
                    else if(locTriExp = boost::dynamic_pointer_cast<LocalRegions::TriExp>(bndCondExp[i]->GetExp(j)))
                    {
                        nbndexp++;
                        FaceGeom = locTriExp->GetGeom2D();
                        id      = FaceGeom->GetFid();
                        gid     = FaceElmtGid[MeshFaceId.find(id)->second];
                        
                        order_e = locTriExp->GetNcoeffs();
                        
                        // Since boundary information is defined to be
                        // aligned with the geometry just use forward/forward
                        // (both coordinate directions) defintiion for gid's
                        for(k = 0; k < order_e; ++k)
                        {
                            m_bndCondCoeffsToGlobalCoeffsMap[cnt++] = gid + k;
                        }
                    }
                    else
                    {
                        ASSERTL0(false,"dynamic cast to a local face expansion failed");
                    }
                }
            }
            
            m_numGlobalBndCoeffs = trace->GetNcoeffs();
            m_numGlobalCoeffs = m_numGlobalBndCoeffs;

            CalculateBndSystemBandWidth();

            if( m_solnType == eDirectMultiLevelStaticCond && nGraphVerts )
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
                    //if face is quad
                    if(locQuadExp = boost::dynamic_pointer_cast<LocalRegions::QuadExp>(bndCondExp[i]->GetExp(j)))
                    {
                        FaceGeom = locQuadExp->GetGeom2D();
                        id      = FaceGeom->GetFid();
                        m_bndCondTraceToGlobalTraceMap[cnt++] = MeshFaceId.find(id)->second;
                    }
                    //else if face is triangle
                    else if(locTriExp = boost::dynamic_pointer_cast<LocalRegions::TriExp>(bndCondExp[i]->GetExp(j)))
                    {
                        FaceGeom = locTriExp->GetGeom2D();
                        id      = FaceGeom->GetFid();
                        m_bndCondTraceToGlobalTraceMap[cnt++] = MeshFaceId.find(id)->second;
                    }
                    else
                    {
                        ASSERTL0(false,"dynamic cast to a local face expansion failed");
                    }
                }
            }
            
            // Now set up mapping from global coefficients to universal.
            SetUpUniversalDGMap3D(locExp);

            // Initialise GSlib and populate the unique map.
            Nektar::Array<OneD, long> tmp(m_globalToUniversalBndMap.num_elements());
            for (unsigned int i = 0; i < m_globalToUniversalBndMap.num_elements(); ++i)
            {
                tmp[i] = m_globalToUniversalBndMap[i];
            }
            m_gsh = Gs::Init(tmp, m_comm);
            Gs::Unique(tmp, m_comm);
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

        /**
		 * Temporary implementation for 3D trace, will be merged with the
		 * SetUpUniversalDGMap
         */
        void LocalToGlobalDGMap::SetUpUniversalDGMap3D(const ExpList &locExp)
        {
            LocalRegions::QuadExpSharedPtr locQuadExp;
            LocalRegions::TriExpSharedPtr locTriExp;

            int fid = 0;
            int cnt = 0;
            int i,j,k;
            int id = 0;
            int order_f = 0;
            int vGlobalId = 0;
            int maxFaceDof = 0;
            int dof = 0;
            const StdRegions::StdExpansionVector &locExpVector = *(locExp.GetExp());

            // Initialise the global to universal maps.
            m_globalToUniversalBndMap = Nektar::Array<OneD, int>(m_numGlobalBndCoeffs, -1);
            m_globalToUniversalBndMapUnique = Nektar::Array<OneD, int>(m_numGlobalBndCoeffs, -1);

            // Loop over all the elements in the domain and compute max face
            // DOF. Reduce across all processes to get universal maximum.
            for(i = 0; i < locExpVector.size(); ++i)
            {
                // Loop over all edges of element i
                for(j = 0; j < locExpVector[i]->GetNfaces(); ++j)
                {
                    dof = locExpVector[i]->GetFaceNcoeffs(j);
                    maxFaceDof = (dof > maxFaceDof ? dof : maxFaceDof);
                }
            }
            m_comm->AllReduce(maxFaceDof, LibUtilities::ReduceMax);

            // Now have trace faces Gid position
            cnt = 0;
            for(i = 0; i < locExpVector.size(); ++i)
            {
                // order list according to m_offset_elmt_id details in
                // Exp2D so that triangules are listed first and then
                // quads
                fid = locExp.GetOffset_Elmt_Id(i);

                // Populate mapping for each edge of the element.
                for(j = 0; j < locExpVector[fid]->GetNfaces(); ++j)
                {
                    //if face is a quad
                    if(locQuadExp = boost::dynamic_pointer_cast<LocalRegions::QuadExp>(m_elmtToFace[fid][j]))
                    {
                        id  = locQuadExp->GetGeom2D()->GetFid();
                    }
                    //else if face is a triangle
                    else if(locTriExp = boost::dynamic_pointer_cast<LocalRegions::TriExp>(m_elmtToFace[fid][j]))
                    {
                        id  = locTriExp->GetGeom2D()->GetFid();
                    }
                    else
                    {
                        ASSERTL0(false,"dynamic cast to a local face expansion failed");
                    }
                    
                    order_f = locExpVector[fid]->GetFaceNcoeffs(j);
                    
                    for(k = 0; k < order_f; ++k)
                    {
                        vGlobalId = m_localToGlobalBndMap[k+cnt];
                        m_globalToUniversalBndMap[vGlobalId]
                            = id * maxFaceDof + k + 1;
                    }
                    cnt += order_f;
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
                    const NekVector<NekDouble>& loc,
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
                    const NekVector<NekDouble>& global,
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
                    const NekVector<NekDouble>& loc,
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
