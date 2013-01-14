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

#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <MultiRegions/ExpList.h>
#include <LocalRegions/SegExp.h>
#include <LocalRegions/TriExp.h>
#include <LocalRegions/QuadExp.h>
#include <LocalRegions/HexExp.h>
#include <LocalRegions/TetExp.h>
#include <LocalRegions/PrismExp.h>
#include <LocalRegions/PyrExp.h>
#include <LocalRegions/PointExp.h>

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>


namespace Nektar
{
    namespace MultiRegions
    {
        AssemblyMapDG::AssemblyMapDG():
            m_numDirichletBndPhys(0)
        {
        }

        AssemblyMapDG::~AssemblyMapDG()
        {
        }


        /**
         *
         */
        AssemblyMapDG::AssemblyMapDG(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const SpatialDomains::MeshGraphSharedPtr &graph1D,
            const ExpList0DSharedPtr &trace,
            const ExpList &locExp,
            const Array<OneD, const MultiRegions::ExpListSharedPtr> &bndCondExp,
            const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndCond,
            const map<int,int> &periodicVertices)
        : AssemblyMap(pSession)
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
                if((locSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>((*exp1D)[i])))
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

            m_hash = boost::hash_range(m_localToGlobalBndMap.begin(),
                                       m_localToGlobalBndMap.end());
        }


        /**
         *
         */
        AssemblyMapDG::AssemblyMapDG(const LibUtilities::SessionReaderSharedPtr &pSession,
                                               const SpatialDomains::MeshGraphSharedPtr &graph2D,
                                               const ExpList1DSharedPtr &trace,
                                               const ExpList &locExp,
                                               const Array<OneD, MultiRegions::ExpListSharedPtr> &bndCondExp,
                                               const Array<OneD, SpatialDomains::BoundaryConditionShPtr> &bndCond,
                                               const map<int,int> &periodicEdges) :
                AssemblyMap(pSession)
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
                if((locSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(trace->GetExp(i))))
                {
                    id = (locSegExp->GetGeom1D())->GetEid();

                    /*
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
                    */
                    MeshEdgeId[id] = i;
                    /*
                    }
                    */
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

            Array<OneD, StdRegions::StdExpansionSharedPtr> edgemap(cnt);
            m_elmtToTrace = Array<OneD, Array<OneD,StdRegions::StdExpansionSharedPtr> >(nel);

            // set up edge expansions links;
            cnt = 0;
            for(i = 0; i < nel; ++i)
            {
                m_elmtToTrace[i] = edgemap + cnt;

                if((locQuadExp = boost::dynamic_pointer_cast<LocalRegions::QuadExp>((*exp2D)[i])))
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
                else if((locTriExp = boost::dynamic_pointer_cast<LocalRegions::TriExp>((*exp2D)[i])))
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
            m_numLocalDirBndCoeffs = 0;
            m_numDirichletBndPhys  = 0;

            cnt = 0;
            for(i = 0; i < bndCondExp.num_elements(); ++i)
            {
                for(j = 0; j < bndCondExp[i]->GetExpSize(); ++j)
                {
                    if((locSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(bndCondExp[i]->GetExp(j))))
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
                        break;
                    }
                    case eDirectStaticCond:
                    case eIterativeStaticCond:
                    {
                        CuthillMckeeReordering(boostGraphObj,perm,iperm);
                        break;
                    }
                    case eDirectMultiLevelStaticCond:
                    {
                        MultiLevelBisectionReordering(boostGraphObj,perm,iperm,bottomUpGraph);
                        break;
                    }
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
                            {
                                // reverse vertex order
                                m_localToGlobalBndMap[cnt]   = gid + 1;
                                m_localToGlobalBndMap[cnt+1] = gid;
                                for (k = 2; k < order_e; ++k)
                                {
                                    m_localToGlobalBndMap[k+cnt] = gid + k;
                                }

                                // negate odd modes
                                for(k = 3; k < order_e; k+=2)
                                {
                                    m_localToGlobalBndSign[cnt+k] = -1.0;
                                }
                                break;
                            }
                            case LibUtilities::eGLL_Lagrange:
                            {
                                // reverse  order
                                for(k = 0; k < order_e; ++k)
                                {
                                    m_localToGlobalBndMap[cnt+order_e-k-1] = gid + k;
                                }
                                break;
                            }   
                            case LibUtilities::eGauss_Lagrange:
                            {
                                // reverse  order
                                for(k = 0; k < order_e; ++k)
                                {
                                    m_localToGlobalBndMap[cnt+order_e-k-1] = gid + k;
                                }
                                break;
                            }
                            default:
                            {
                                ASSERTL0(false,"Boundary type not permitted");
                            }
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
                    if((locSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(bndCondExp[i]->GetExp(j))))
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

                    m_nextLevelLocalToGlobalMap = MemoryManager<AssemblyMap>::
                        AllocateSharedPtr(this,bottomUpGraph);
                }
            }

            cnt = 0;
            m_bndCondTraceToGlobalTraceMap = Array<OneD, int >(nbndexp);
            for(i = 0; i < bndCondExp.num_elements(); ++i)
            {
                for(j = 0; j < bndCondExp[i]->GetExpSize(); ++j)
                {
                    if((locSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(bndCondExp[i]->GetExp(j))))
                    {
                        SegGeom = locSegExp->GetGeom1D();
                        id = SegGeom->GetEid();

                        m_bndCondTraceToGlobalTraceMap[cnt++] = MeshEdgeId.find(id)->second;
                    }
                }
            }

            // Now set up mapping from global coefficients to universal.
            ExpListSharedPtr tr = boost::dynamic_pointer_cast<ExpList>(trace);
            SetUpUniversalDGMap   (locExp);
            SetUpUniversalTraceMap(locExp, tr);

            m_hash = boost::hash_range(m_localToGlobalBndMap.begin(),
                                       m_localToGlobalBndMap.end());
        }

        /**
         * Constructor for trace map for three-dimensional expansion.
         */
        AssemblyMapDG::AssemblyMapDG(
            const LibUtilities::SessionReaderSharedPtr                &pSession,
            const SpatialDomains::MeshGraphSharedPtr                  &graph3D,
            const ExpList2DSharedPtr                                  &trace,
            const ExpList                                             &locExp,
            const Array<OneD, MultiRegions::ExpListSharedPtr>         &bndCondExp,
            const Array<OneD, SpatialDomains::BoundaryConditionShPtr> &bndCond,
            const map<int,PeriodicFace>                               &periodicFaces):
            AssemblyMap(pSession)
        {
            int i,j,k,cnt,eid, id, id1, order_e,gid;
            int ntrace_exp = trace->GetExpSize();
            int nbnd = bndCondExp.num_elements();
            LocalRegions::QuadExpSharedPtr      locQuadExp, locQuadExp1;
            LocalRegions::TriExpSharedPtr       locTriExp, locTriExp1;
            LocalRegions::HexExpSharedPtr       locHexExp;
            LocalRegions::PrismExpSharedPtr     locPrismExp;
            LocalRegions::PyrExpSharedPtr       locPyrExp;
            LocalRegions::TetExpSharedPtr       locTetExp;
            SpatialDomains::Geometry2DSharedPtr FaceGeom;
            StdRegions::StdExpansionSharedPtr   locBndExp;

            const boost::shared_ptr<StdRegions::StdExpansionVector> exp3D = 
                locExp.GetExp();
            int nel = exp3D->size();

            map<int, int> MeshFaceId;

            m_signChange = true;

            // determine mapping from geometry edges to trace
            for(i = 0; i < ntrace_exp; ++i)
            {
                id = trace->GetExp(i)->GetGeom2D()->GetFid();
                /*
                if(periodicFaces.count(id) > 0)
                {
                    if(MeshFaceId.count(id) == 0)
                    {
                        id1 = periodicFaces.find(id)->second.first;
                        MeshFaceId[id] = i;
                        MeshFaceId[id1] = i;
                    }
                }
                else
                {
                */
                MeshFaceId[id] = i;
                /*
                }
                */
            }

            // Count total number of faces
            cnt = 0;
            for(i = 0; i < nel; ++i)
            {
                cnt += (*exp3D)[i]->GetNfaces();
            }

            Array<OneD, StdRegions::StdExpansionSharedPtr> facemap(cnt);
            m_elmtToTrace = Array<OneD, 
                Array<OneD, StdRegions::StdExpansionSharedPtr> >(nel);

            // set up face expansions links;
            cnt = 0;
            for(i = 0; i < nel; ++i)
            {
                m_elmtToTrace[i] = facemap + cnt;
                
                for(j = 0; j < (*exp3D)[i]->GetNfaces(); ++j)
                {
                    id = (*exp3D)[i]->GetGeom3D()->GetFid(j);
                    
                    if(MeshFaceId.count(id) > 0)
                    {
                        m_elmtToTrace[i][j] = 
                            trace->GetExp(MeshFaceId.find(id)->second);
                    }
                    else
                    {
                        ASSERTL0(false,"Failed to find face map");
                    }
                }

                cnt += (*exp3D)[i]->GetNfaces();
            }

            // Set up boundary mapping
            cnt = 0;
            for(i = 0; i < nbnd; ++i)
            {
                cnt += bndCondExp[i]->GetExpSize();
            }

            set<int> dirFaces;

            m_numLocalDirBndCoeffs = 0;
            m_numDirichletBndPhys  = 0;

            cnt = 0;
            for(i = 0; i < bndCondExp.num_elements(); ++i)
            {
                for(j = 0; j < bndCondExp[i]->GetExpSize(); ++j)
                {
                    locBndExp = bndCondExp[i]->GetExp(j);
                    FaceGeom  = locBndExp->GetGeom2D();
                    id        = FaceGeom->GetFid();
                    
                    if(bndCond[i]->GetBoundaryConditionType() == 
                           SpatialDomains::eDirichlet)
                    {
                        m_numLocalDirBndCoeffs += locBndExp->GetNcoeffs();
                        m_numDirichletBndPhys  += locBndExp->GetTotPoints();
                        dirFaces.insert(id);
                    }
                }

                cnt += j;
            }

            // Set up integer mapping array and sign change for each
            // degree of freedom + initialise some more data members
            m_staticCondLevel           = 0;
            m_numPatches                = nel;
            m_numLocalBndCoeffsPerPatch = Array<OneD, unsigned int>(nel);
            m_numLocalIntCoeffsPerPatch = Array<OneD, unsigned int>(nel);

            int nbndry = 0;
            for(i = 0; i < nel; ++i) // count number of elements in array
            {
                eid     = locExp.GetOffset_Elmt_Id(i);
                nbndry += (*exp3D)[eid]->NumDGBndryCoeffs();
                m_numLocalBndCoeffsPerPatch[i] = (unsigned int) (*exp3D)[eid]->NumDGBndryCoeffs();
                m_numLocalIntCoeffsPerPatch[i] = (unsigned int) 0;
            }

            m_numGlobalDirBndCoeffs = m_numLocalDirBndCoeffs;
            m_numLocalBndCoeffs     = nbndry;
            m_numLocalCoeffs        = nbndry;
            m_localToGlobalBndMap   = Array<OneD, int>       (nbndry);
            m_localToGlobalBndSign  = Array<OneD, NekDouble> (nbndry,1);

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
            int face_id, face_id1;
            int dirOffset = 0;

            // make trace face renumbering map where first solved
            // face starts at 0 so we can set up graph.
            for(i = 0; i < ntrace_exp; ++i)
            {
                id = trace->GetExp(i)->GetGeom2D()->GetFid();
                
                if (dirFaces.count(id) == 0)
                {
                    // Initial put in element ordering (starting
                    // from zero) into FaceElmtGid
                    boost::add_vertex(boostGraphObj);
                    FaceElmtGid[i] = cnt++;
                }
                else
                {
                    // Use existing offset for Dirichlet edges
                    FaceElmtGid[i] = dirOffset;
                    dirOffset     += trace->GetExp(i)->GetNcoeffs();
                    nDir++;
                }
            }

            // Set up boost Graph
            for(i = 0; i < nel; ++i)
            {
                eid = locExp.GetOffset_Elmt_Id(i);

                for(j = 0; j < (*exp3D)[eid]->GetNfaces(); ++j)
                {
                    // Add face to boost graph for non-Dirichlet Boundary
                    FaceGeom = m_elmtToTrace[eid][j]->GetGeom2D();
                    id       = FaceGeom->GetFid();
                    face_id  = MeshFaceId.find(id)->second;

                    if(dirFaces.count(id) == 0)
                    {
                        for(k = j+1; k < (*exp3D)[eid]->GetNfaces(); ++k)
                        {
                            FaceGeom = m_elmtToTrace[eid][k]->GetGeom2D();
                            id1      = FaceGeom->GetFid();
                            face_id1 = MeshFaceId.find(id1)->second;
                            
                            if(dirFaces.count(id1) == 0)
                            {
                                boost::add_edge((size_t) FaceElmtGid[face_id], 
                                                (size_t) FaceElmtGid[face_id1], 
                                                boostGraphObj);
                            }
                        }
                    }
                }
            }

            int                                 nGraphVerts = ntrace_exp - nDir;
            Array<OneD, int>                    perm (nGraphVerts);
            Array<OneD, int>                    iperm(nGraphVerts);
            Array<OneD, int>                    vwgts(nGraphVerts);
            BottomUpSubStructuredGraphSharedPtr bottomUpGraph;
            
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
                    case eIterativeStaticCond:
                    {
                        NoReordering(boostGraphObj,perm,iperm);
                        break;
                    }
                    case eDirectStaticCond:
                    {
                        CuthillMckeeReordering(boostGraphObj,perm,iperm);
                        break;
                    }
                    case eDirectMultiLevelStaticCond:
                    {
                        MultiLevelBisectionReordering(boostGraphObj,perm,iperm,
                                                      bottomUpGraph);
                        break;
                    }
                    default:
                    {
                        ASSERTL0(false,"Unrecognised solution type");
                    }
                }
            }

            // Recast the permutation so that it can be used as a map Form old
            // trace face ID to new trace face ID
            cnt = m_numLocalDirBndCoeffs;
            for(i = 0; i < ntrace_exp - nDir; ++i)
            {
                FaceElmtGid[perm[i]+nDir] = cnt;
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
                    FaceGeom = m_elmtToTrace[eid][j]->GetGeom2D();
                    id       = FaceGeom->GetFid();
                    gid      = FaceElmtGid[MeshFaceId.find(id)->second];
                    order_e  = (*exp3D)[eid]->GetFaceNcoeffs(j);
                    
                    std::map<int,int> orientMap;
                    
                    Array<OneD, unsigned int> elmMap1 (order_e);
                    Array<OneD,          int> elmSign1(order_e);
                    Array<OneD, unsigned int> elmMap2 (order_e);
                    Array<OneD,          int> elmSign2(order_e);
                    
                    StdRegions::Orientation fo = (*exp3D)[eid]->GetFaceOrient(j);
                    
                    // Construct mapping which will permute global IDs
                    // according to face orientations. 
                    (*exp3D)[eid]->GetFaceToElementMap(j,fo,elmMap1,elmSign1);
                    (*exp3D)[eid]->GetFaceToElementMap(
                        j,StdRegions::eDir1FwdDir1_Dir2FwdDir2,elmMap2,elmSign2);
                    
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
                        m_localToGlobalBndSign[k+cnt] = elmSign2[orientMap[k]];
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

            m_bndCondCoeffsToGlobalCoeffsMap = Array<OneD,int>(cnt);

            // Number of boundary expansions
            int nbndexp = 0, bndOffset, bndTotal = 0;
            for(cnt = i = 0; i < nbnd; ++i)
            {
                for(j = 0; j < bndCondExp[i]->GetExpSize(); ++j)
                {
                    locBndExp = bndCondExp[i]->GetExp(j);
                    id        = locBndExp->GetGeom2D()->GetFid();
                    gid       = FaceElmtGid[MeshFaceId.find(id)->second];
                    bndOffset = bndCondExp[i]->GetCoeff_Offset(j) + bndTotal;
                        
                    // Since boundary information is defined to be aligned with
                    // the geometry just use forward/forward (both coordinate
                    // directions) defintiion for gid's
                    for(k = 0; k < locBndExp->GetNcoeffs(); ++k)
                    {
                        m_bndCondCoeffsToGlobalCoeffsMap[bndOffset+k] = gid + k;
                    }
                }

                nbndexp  += bndCondExp[i]->GetExpSize();
                bndTotal += bndCondExp[i]->GetNcoeffs();
            }
            
            m_numGlobalBndCoeffs = trace->GetNcoeffs();
            m_numGlobalCoeffs    = m_numGlobalBndCoeffs;

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

                    m_nextLevelLocalToGlobalMap = MemoryManager<AssemblyMap>::
                        AllocateSharedPtr(this,bottomUpGraph);
                }
            }

            cnt = 0;
            m_bndCondTraceToGlobalTraceMap = Array<OneD, int>(nbndexp);
            for(i = 0; i < bndCondExp.num_elements(); ++i)
            {
                for(j = 0; j < bndCondExp[i]->GetExpSize(); ++j)
                {
                    locBndExp = bndCondExp[i]->GetExp(j);
                    FaceGeom  = locBndExp->GetGeom2D();
                    id        = FaceGeom->GetFid();
                    m_bndCondTraceToGlobalTraceMap[cnt++] = 
                        MeshFaceId.find(id)->second;
                }
            }

            // Now set up mapping from global coefficients to universal.
            ExpListSharedPtr tr = boost::dynamic_pointer_cast<ExpList>(trace);
            SetUpUniversalDGMap   (locExp);
            SetUpUniversalTraceMap(locExp, tr);

            m_hash = boost::hash_range(m_localToGlobalBndMap.begin(),
                                       m_localToGlobalBndMap.end());
        }

        /**
         * Constructs a mapping between the process-local global numbering and
         * a universal numbering of the trace space expansion. The universal
         * numbering is defined by the mesh edge IDs to enforce consistency
         * across processes.
         * 
         * @param       locExp  List of local elemental expansions.
         */
        void AssemblyMapDG::SetUpUniversalDGMap(const ExpList &locExp)
        {
            StdRegions::StdExpansionSharedPtr locExpansion;
            int eid       = 0;
            int cnt       = 0;
            int id        = 0;
            int order_e   = 0;
            int vGlobalId = 0;
            int maxDof    = 0;
            int dof       = 0;
            int nDim      = 0;
            int i,j,k;

            const StdRegions::StdExpansionVector &locExpVector = *(locExp.GetExp());

            // Initialise the global to universal maps.
            m_globalToUniversalBndMap = Nektar::Array<OneD, int>(m_numGlobalBndCoeffs, -1);
            m_globalToUniversalBndMapUnique = Nektar::Array<OneD, int>(m_numGlobalBndCoeffs, -1);

            // Loop over all the elements in the domain and compute max edge
            // DOF. Reduce across all processes to get universal maximum.
            for(i = 0; i < locExpVector.size(); ++i)
            {
                locExpansion = boost::dynamic_pointer_cast<
                    StdRegions::StdExpansion>(locExpVector[i]);
                nDim = locExpansion->GetShapeDimension();

                // Loop over all edges of element i
                if (nDim == 1)
                {
                    maxDof = (1 > maxDof ? 1 : maxDof);
                }
                else if (nDim == 2)
                {
                    for (j = 0; j < locExpansion->GetNedges(); ++j)
                    {
                        dof    = locExpansion->GetEdgeNcoeffs(j);
                        maxDof = (dof > maxDof ? dof : maxDof);
                    }
                }
                else if (nDim == 3)
                {
                    for (j = 0; j < locExpansion->GetNfaces(); ++j)
                    {
                        dof    = locExpansion->GetFaceNcoeffs(j);
                        maxDof = (dof > maxDof ? dof : maxDof);
                    }
                }
            }
            m_comm->AllReduce(maxDof, LibUtilities::ReduceMax);

            // Now have trace edges Gid position
            cnt = 0;
            for(i = 0; i < locExpVector.size(); ++i)
            {
                locExpansion = boost::dynamic_pointer_cast<
                    StdRegions::StdExpansion>(locExpVector[i]);
                nDim = locExpansion->GetShapeDimension();

                // Order list according to m_offset_elmt_id details in Exp2D
                // so that triangules are listed first and then quads
                eid = locExp.GetOffset_Elmt_Id(i);

                // Populate mapping for each edge of the element.
                if (nDim == 1)
                {
                    for(j = 0; j < locExpansion->GetNverts(); ++j, ++cnt)
                    {
                        LocalRegions::PointExpSharedPtr locPointExp = 
                            boost::dynamic_pointer_cast<
                                LocalRegions::PointExp>(m_elmtToTrace[eid][j]);
                        id = locPointExp->GetGeom()->GetEid();
                        vGlobalId = m_localToGlobalBndMap[cnt+j];
                        m_globalToUniversalBndMap[vGlobalId]
                            = id * maxDof + j + 1;
                    }
                } 
                else if (nDim == 2)
                {
                    for(j = 0; j < locExpansion->GetNedges(); ++j)
                    {
                        LocalRegions::SegExpSharedPtr locSegExp = 
                            boost::dynamic_pointer_cast<
                                LocalRegions::SegExp>(m_elmtToTrace[eid][j]);

                        id  = locSegExp->GetGeom1D()->GetEid();
                        order_e = locExpVector[eid]->GetEdgeNcoeffs(j);
                        
                        map<int,int> orientMap;
                        Array<OneD, unsigned int> map1(order_e), map2(order_e);
                        Array<OneD, int> sign1(order_e), sign2(order_e);
                        
                        locExpVector[eid]->GetEdgeToElementMap(j, StdRegions::eForwards, map1, sign1);
                        locExpVector[eid]->GetEdgeToElementMap(j, locExpVector[eid]->GetEorient(j), map2, sign2);
                        
                        for (k = 0; k < map1.num_elements(); ++k)
                        {
                            // Find the elemental co-efficient in the original
                            // mapping.
                            int idx = -1;
                            for (int l = 0; l < map2.num_elements(); ++l)
                            {
                                if (map1[k] == map2[l])
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
                            vGlobalId = m_localToGlobalBndMap[k+cnt];
                            m_globalToUniversalBndMap[vGlobalId]
                                = id * maxDof + orientMap[k] + 1;
                        }
                        cnt += order_e;
                    }
                }
                else if (nDim == 3)
                {
                    for(j = 0; j < locExpansion->GetNfaces(); ++j)
                    {
                        LocalRegions::Expansion2DSharedPtr locFaceExp = 
                            boost::dynamic_pointer_cast<
                                LocalRegions::Expansion2D>(m_elmtToTrace[eid][j]);

                        id  = locFaceExp->GetGeom2D()->GetFid();
                        order_e = locExpVector[eid]->GetFaceNcoeffs(j);

                        map<int,int> orientMap;
                        Array<OneD, unsigned int> map1(order_e), map2(order_e);
                        Array<OneD, int> sign1(order_e), sign2(order_e);
                        
                        locExpVector[eid]->GetFaceToElementMap(j, StdRegions::eDir1FwdDir1_Dir2FwdDir2, map1, sign1);
                        locExpVector[eid]->GetFaceToElementMap(j, locExpVector[eid]->GetFaceOrient(j), map2, sign2);
                        
                        for (k = 0; k < map1.num_elements(); ++k)
                        {
                            // Find the elemental co-efficient in the original
                            // mapping.
                            int idx = -1;
                            for (int l = 0; l < map2.num_elements(); ++l)
                            {
                                if (map1[k] == map2[l])
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
                            vGlobalId = m_localToGlobalBndMap[k+cnt];
                            m_globalToUniversalBndMap[vGlobalId]
                                = id * maxDof + orientMap[k] + 1;
                        }
                        cnt += order_e;
                    }
                }
            }
            
            // Initialise GSlib and populate the unique map.
            Array<OneD, long> tmp(m_globalToUniversalBndMap.num_elements());
            for (i = 0; i < m_globalToUniversalBndMap.num_elements(); ++i)
            {
                tmp[i] = m_globalToUniversalBndMap[i];
            }
            m_bndGsh = m_gsh = Gs::Init(tmp, m_comm);
            Gs::Unique(tmp, m_comm);
            for (i = 0; i < m_globalToUniversalBndMap.num_elements(); ++i)
            {
                m_globalToUniversalBndMapUnique[i] = (tmp[i] >= 0 ? 1 : 0);
            }
        }

        void AssemblyMapDG::SetUpUniversalTraceMap(const ExpList         &locExp,
                                                        const ExpListSharedPtr trace)
        {
            StdRegions::StdExpansionSharedPtr locExpansion;
            int i,j,k;
            int maxQuad = 0, quad = 0, nDim = 0, eid = 0, offset = 0;

            const StdRegions::StdExpansionVector &locExpVector = *(locExp.GetExp());

            int nTracePhys = trace->GetTotPoints();

            // Initialise the trace to universal maps.
            m_traceToUniversalMap       = 
                Nektar::Array<OneD, int>(nTracePhys, -1);
            m_traceToUniversalMapUnique = 
                Nektar::Array<OneD, int>(nTracePhys, -1);

            // Assume that each element of the expansion is of the same
            // dimension.
            nDim = locExpVector[0]->GetShapeDimension();

            if (nDim == 1)
            {
                maxQuad = (1 > maxQuad ? 1 : maxQuad);
            }
            else
            {
                for (i = 0; i < trace->GetExpSize(); ++i)
                {
                    quad = trace->GetExp(i)->GetTotPoints();
                    if (quad > maxQuad)
                    {
                        maxQuad = quad;
                    }
                }
            }
            m_comm->AllReduce(maxQuad, LibUtilities::ReduceMax);

            if (nDim == 1)
            {
                for (int i = 0; i < trace->GetExpSize(); ++i)
                {
                    eid = trace->GetExp(i)->GetGeom()->GetGlobalID();
                    offset = trace->GetPhys_Offset(i);
                    m_traceToUniversalMap[offset] = eid*maxQuad+1;
                }
            }
            else 
            {
                for (int i = 0; i < trace->GetExpSize(); ++i)
                {
                    eid    = trace->GetExp(i)->GetGeom()->GetGlobalID();
                    offset = trace->GetPhys_Offset(i);
                    quad   = trace->GetExp(i)->GetTotPoints();

                    for(int j = 0; j < quad; ++j)
                    {
                        m_traceToUniversalMap[j+offset] = eid*maxQuad+j+1;
                    }
                }
            }

            Array<OneD, long> tmp(nTracePhys);
            for (int i = 0; i < nTracePhys; ++i)
            {
                tmp[i] = m_traceToUniversalMap[i];
            }
            m_traceGsh = Gs::Init(tmp, m_comm);
            Gs::Unique(tmp, m_comm);
            for (int i = 0; i < nTracePhys; ++i)
            {
                m_traceToUniversalMapUnique[i] = tmp[i];
            }
        }

        void AssemblyMapDG::UniversalTraceAssemble(
            Array<OneD, NekDouble> &pGlobal) const
        {
            Gs::Gather(pGlobal, Gs::gs_add, m_traceGsh);
        }

        int AssemblyMapDG::v_GetLocalToGlobalMap(const int i) const
        {
            return m_localToGlobalBndMap[i];
        }

        int AssemblyMapDG::v_GetGlobalToUniversalMap(const int i) const
        {
            return m_globalToUniversalBndMap[i];
        }

        int AssemblyMapDG::v_GetGlobalToUniversalMapUnique(const int i) const
        {
            return m_globalToUniversalBndMapUnique[i];
        }

        const Array<OneD,const int>& AssemblyMapDG::v_GetLocalToGlobalMap()
        {
            return m_localToGlobalBndMap;
        }

        const Array<OneD,const int>& AssemblyMapDG::v_GetGlobalToUniversalMap()
        {
            return m_globalToUniversalBndMap;
        }

        const Array<OneD,const int>& AssemblyMapDG::v_GetGlobalToUniversalMapUnique()
        {
            return m_globalToUniversalBndMapUnique;
        }

        NekDouble AssemblyMapDG::v_GetLocalToGlobalSign(
                    const int i) const
        {
            return GetLocalToGlobalBndSign(i);
        }

        const void AssemblyMapDG::v_LocalToGlobal(
                    const Array<OneD, const NekDouble>& loc,
                          Array<OneD,       NekDouble>& global) const
        {
            AssembleBnd(loc,global);
        }

        const void AssemblyMapDG::v_LocalToGlobal(
                    const NekVector<NekDouble>& loc,
                          NekVector<      NekDouble>& global) const
        {
            AssembleBnd(loc,global);
        }

        const void AssemblyMapDG::v_GlobalToLocal(
                    const Array<OneD, const NekDouble>& global,
                          Array<OneD,       NekDouble>& loc) const
        {
            GlobalToLocalBnd(global,loc);
        }

        const void AssemblyMapDG::v_GlobalToLocal(
                    const NekVector<NekDouble>& global,
                          NekVector<      NekDouble>& loc) const
        {
            GlobalToLocalBnd(global,loc);
        }

        const void AssemblyMapDG::v_Assemble(
                    const Array<OneD, const NekDouble> &loc,
                          Array<OneD,       NekDouble> &global) const
        {
            AssembleBnd(loc,global);
        }

        const void AssemblyMapDG::v_Assemble(
                    const NekVector<NekDouble>& loc,
                          NekVector<      NekDouble>& global) const
        {
            AssembleBnd(loc,global);
        }

        const void AssemblyMapDG::v_UniversalAssemble(
                      Array<OneD,     NekDouble>& pGlobal) const
        {
            Gs::Gather(pGlobal, Gs::gs_add, m_gsh);
        }

        const void AssemblyMapDG::v_UniversalAssemble(
                      NekVector<      NekDouble>& pGlobal) const
        {
            UniversalAssemble(pGlobal.GetPtr());
        }

        const int AssemblyMapDG::v_GetFullSystemBandWidth() const
        {
            return GetBndSystemBandWidth();
        }

        int AssemblyMapDG::GetTraceToUniversalMap(int i)
        {
            return m_traceToUniversalMap[i];
        }

        int AssemblyMapDG::GetTraceToUniversalMapUnique(int i)
        {
            return m_traceToUniversalMapUnique[i];
        }

        int AssemblyMapDG::GetNumDirichletBndPhys()
        {
            return m_numDirichletBndPhys;
        }

        Array<OneD, StdRegions::StdExpansionSharedPtr>&
                    AssemblyMapDG::GetElmtToTrace(const int i)
        {
            ASSERTL1(i >= 0 && i < m_elmtToTrace.num_elements(),
                     "i is out of range");
            return m_elmtToTrace[i];
        }

        Array<OneD, Array< OneD, StdRegions::StdExpansionSharedPtr> >&
                    AssemblyMapDG::GetElmtToTrace()
        {
            return m_elmtToTrace;
        }


    } //namespace
} // namespace
