///////////////////////////////////////////////////////////////////////////////
//
// File ExpListConstructors.cpp
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
// Description: Constructors for Expansion list
//
///////////////////////////////////////////////////////////////////////////////

#include <iomanip>
#include <MultiRegions/ExpList.h>
#include <LocalRegions/PointExp.h>
#include <LocalRegions/SegExp.h>
#include <LocalRegions/TriExp.h>
#include <LocalRegions/QuadExp.h>
#include <LocalRegions/NodalTriExp.h>
#include <LocalRegions/HexExp.h>
#include <LocalRegions/PrismExp.h>
#include <LocalRegions/PyrExp.h>
#include <LocalRegions/TetExp.h>
#include <LocalRegions/Expansion0D.h>
#include <LocalRegions/Expansion3D.h>

using namespace std;

namespace Nektar
{
namespace MultiRegions
{
    //----------------------------------------------------------------------
    //                        0D Expansion Constructors 
    //----------------------------------------------------------------------
    ExpList::ExpList(const SpatialDomains::PointGeomSharedPtr &geom):
        m_expType(e0D),
        m_ncoeffs(1),
        m_npoints(1),
        m_physState(false),
        m_exp(MemoryManager<LocalRegions::ExpansionVector>
              ::AllocateSharedPtr()),
        m_blockMat(MemoryManager<BlockMatrixMap>::AllocateSharedPtr()),
        m_WaveSpace(false)
    {
        LocalRegions::PointExpSharedPtr Point =
            MemoryManager<LocalRegions::PointExp>::AllocateSharedPtr(geom);
        (*m_exp).push_back(Point);

        SetupCoeffPhys();
    }


    /**
     * Store expansions for the trace space expansions used in
     * DisContField2D
     *
     * @param  pSession      A session within information about expansion
     * @param  bndConstraint Array of ExpList1D objects each containing a
     *                       1D spectral/hp element expansion on a single
     *                       boundary region.
     * @param  bndCond       Array of BoundaryCondition objects which contain
     *                       information about the boundary conditions on the
     *                       different boundary regions.
     * @param  locexp        Complete domain expansion list.
     * @param  graph         mesh corresponding to the expansion list.
     * @param  DeclareCoeffPhysArrays Declare the coefficient and
     *                               phys space arrays
     * @param  variable      The variable name associated with the expansion
     * @param  ImpType       Detail about the implementation type to use 
     *                       in operators
     */
    ExpList::ExpList(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const Array<OneD,const ExpListSharedPtr>   &bndConstraint,
            const Array<OneD, const SpatialDomains::BoundaryConditionShPtr>
                                                       &bndCond,
            const LocalRegions::ExpansionVector        &locexp,
            const SpatialDomains::MeshGraphSharedPtr   &graph,
            const bool                                  DeclareCoeffPhysArrays,
            const std::string                           variable,
            const Collections::ImplementationType       ImpType):
        m_comm(pSession->GetComm()),
        m_session(pSession),
        m_graph(graph),
        m_physState(false),
        m_exp(MemoryManager<LocalRegions::ExpansionVector>
              ::AllocateSharedPtr()),
        m_blockMat(MemoryManager<BlockMatrixMap>::AllocateSharedPtr()),
        m_WaveSpace(false)
    {
        boost::ignore_unused(variable,ImpType);
        int i, j, id, elmtid = 0;
        set<int> tracesDone;
        
        SpatialDomains::PointGeomSharedPtr  PointGeom;
        SpatialDomains::Geometry1DSharedPtr segGeom;
        SpatialDomains::Geometry2DSharedPtr ElGeom;
        SpatialDomains::Geometry2DSharedPtr FaceGeom;
        SpatialDomains::QuadGeomSharedPtr   QuadGeom;
        SpatialDomains::TriGeomSharedPtr    TriGeom;
        
        LocalRegions::ExpansionSharedPtr    exp;
        LocalRegions::Expansion0DSharedPtr  exp0D;
        LocalRegions::Expansion1DSharedPtr  exp1D;
        LocalRegions::Expansion2DSharedPtr  exp2D;
        LocalRegions::Expansion3DSharedPtr  exp3D;
        
        // First loop over boundary conditions to reorder
        // Dirichlet boundaries
        for(i = 0; i < bndCond.size(); ++i)
        {
            if(bndCond[i]->GetBoundaryConditionType() == 
                                     SpatialDomains::eDirichlet)
            {
                for(j = 0; j < bndConstraint[i]->GetExpSize(); ++j)
                {
                    if((exp0D = std::dynamic_pointer_cast<
                        LocalRegions::Expansion0D>(
                                         bndConstraint[i]->GetExp(j))))
                    {
                        m_expType = e0D;
                        
                        PointGeom = exp0D->GetGeom()->GetVertex(0);
                        exp= MemoryManager<LocalRegions::PointExp>::
                            AllocateSharedPtr(PointGeom);
                        tracesDone.insert(PointGeom->GetVid());
                    }
                    else if((exp1D = std::dynamic_pointer_cast<
                             LocalRegions::Expansion1D>(
                                          bndConstraint[i]->GetExp(j))))
                    {
                        m_expType = e1D;

                        LibUtilities::BasisKey bkey = exp1D->
                            GetBasis(0)->GetBasisKey();
                        segGeom = exp1D->GetGeom1D();
                        exp = MemoryManager<LocalRegions::SegExp>
                            ::AllocateSharedPtr(bkey, segGeom);
                        tracesDone.insert(segGeom->GetGlobalID());
                    
                    }
                    else if ((exp2D = std::dynamic_pointer_cast
                              <LocalRegions::Expansion2D>(bndConstraint[i]->
                                                          GetExp(j))))
                    {
                        m_expType = e2D;

                        LibUtilities::BasisKey bkey0 = exp2D
                            ->GetBasis(0)->GetBasisKey();
                        LibUtilities::BasisKey bkey1 = exp2D
                            ->GetBasis(1)->GetBasisKey();
                        FaceGeom = exp2D->GetGeom2D();

                        //if face is a quad
                        if((QuadGeom = std::dynamic_pointer_cast<
                            SpatialDomains::QuadGeom>(FaceGeom)))
                        {
                            exp = MemoryManager<LocalRegions::QuadExp>
                                ::AllocateSharedPtr(bkey0, bkey1, QuadGeom);
                            tracesDone.insert(QuadGeom->GetGlobalID());
                        }
                        //if face is a triangle
                        else if((TriGeom = std::dynamic_pointer_cast<
                                 SpatialDomains::TriGeom>(FaceGeom)))
                        {
                            exp = MemoryManager<LocalRegions::TriExp>
                                ::AllocateSharedPtr(bkey0, bkey1, TriGeom);
                            tracesDone.insert(TriGeom->GetGlobalID());
                        }
                        else
                        {
                            ASSERTL0(false,"dynamic cast to a proper "
                                     "face geometry failed");
                        }
                    }
                    // Assign next id
                    exp->SetElmtId(elmtid++);
                
                    // Add the expansion
                    (*m_exp).push_back(exp);
                }
            }
        }

        map<int, pair<SpatialDomains::Geometry1DSharedPtr,
                      LibUtilities::BasisKey> > edgeOrders;

        map<int, pair<SpatialDomains::Geometry2DSharedPtr,
                      pair<LibUtilities::BasisKey,
                           LibUtilities::BasisKey> > > faceOrders;
        
        for(i = 0; i < locexp.size(); ++i)
        {
            if((exp1D =
                std::dynamic_pointer_cast<
                LocalRegions::Expansion1D>(locexp[i])))
            {
                m_expType = e0D;

                for(j = 0; j < 2; ++j)
                {
                    PointGeom = (exp1D->GetGeom1D())->GetVertex(j);
                    id = PointGeom->GetVid();
		
                    // Ignore Dirichlet edges
                    if (tracesDone.count(id) != 0)
                    {
                        continue;
                    }

                    exp = MemoryManager<LocalRegions::PointExp>::
                        AllocateSharedPtr(PointGeom);
                    tracesDone.insert(id);
                    exp->SetElmtId(elmtid++);
                    (*m_exp).push_back(exp);
                }
            }
            else if((exp2D =
                std::dynamic_pointer_cast<
                LocalRegions::Expansion2D>(locexp[i])))
            {
                m_expType = e1D;
                for(j = 0; j < locexp[i]->GetNtraces(); ++j)
                {
                    segGeom = exp2D->GetGeom2D()->GetEdge(j);
                    id      = segGeom->GetGlobalID();
                    // Ignore Dirichlet edges
                    if (tracesDone.count(id) != 0)
                    {
                        continue;
                    }
                    
                    auto it = edgeOrders.find(id);
                    
                    if (it == edgeOrders.end())
                    {
                        edgeOrders.insert(std::make_pair(id, std::make_pair(
                                                                            segGeom, locexp[i]->GetTraceBasisKey(j))));
                    }
                    else // variable modes/points
                    {
                        LibUtilities::BasisKey edge
                            = locexp[i]->GetTraceBasisKey(j);
                        LibUtilities::BasisKey existing
                            = it->second.second;
                        
                        int np1 = edge    .GetNumPoints();
                        int np2 = existing.GetNumPoints();
                        int nm1 = edge    .GetNumModes ();
                        int nm2 = existing.GetNumModes ();
                        
                        if (np2 >= np1 && nm2 >= nm1)
                        {
                            continue;
                        }
                        else if (np2 < np1 && nm2 < nm1)
                        {
                            it->second.second = edge;
                        }
                        else
                        {
                            ASSERTL0(false,
                                    "inappropriate number of points/modes (max"
                                    "num of points is not set with max order)");
                        }
                    }
                }
            }
            else if((exp3D =
                     dynamic_pointer_cast<
                     LocalRegions::Expansion3D>(locexp[i])))
            {
                m_expType = e2D;
                for (j = 0; j < exp3D->GetNtraces(); ++j)
                {
                    FaceGeom = exp3D->GetGeom3D()->GetFace(j);
                    id       = FaceGeom->GetGlobalID();

                    if(tracesDone.count(id) != 0)
                    {
                        continue;
                    }
                    auto it = faceOrders.find(id);

                    if (it == faceOrders.end())
                    {
                        LibUtilities::BasisKey face_dir0
                            = locexp[i]->GetTraceBasisKey(j,0);
                        LibUtilities::BasisKey face_dir1
                            = locexp[i]->GetTraceBasisKey(j,1);

                        faceOrders.insert(
                            std::make_pair(
                                id, std::make_pair(
                                    FaceGeom,
                                    std::make_pair(face_dir0, face_dir1))));
                    }
                    else // variable modes/points
                    {
                        LibUtilities::BasisKey face0     =
                            locexp[i]->GetTraceBasisKey(j,0);
                        LibUtilities::BasisKey face1     =
                            locexp[i]->GetTraceBasisKey(j,1);
                        LibUtilities::BasisKey existing0 =
                            it->second.second.first;
                        LibUtilities::BasisKey existing1 =
                            it->second.second.second;

                        int np11 = face0    .GetNumPoints();
                        int np12 = face1    .GetNumPoints();
                        int np21 = existing0.GetNumPoints();
                        int np22 = existing1.GetNumPoints();
                        int nm11 = face0    .GetNumModes ();
                        int nm12 = face1    .GetNumModes ();
                        int nm21 = existing0.GetNumModes ();
                        int nm22 = existing1.GetNumModes ();

                        if ((np22 >= np12 || np21 >= np11) &&
                            (nm22 >= nm12 || nm21 >= nm11))
                        {
                            continue;
                        }
                        else if((np22 < np12 || np21 < np11) &&
                                (nm22 < nm12 || nm21 < nm11))
                        {
                            it->second.second.first  = face0;
                            it->second.second.second = face1;
                        }
                        else
                        {
                            ASSERTL0(false,
                                    "inappqropriate number of points/modes (max"
                                    "num of points is not set with max order)");
                        }
                    }
                }
            }
        }
        
        int nproc   = m_comm->GetSize(); // number of processors
        int tracepr = m_comm->GetRank(); // ID processor
        
        if (nproc > 1)
        {
            int tCnt = 0;
            
            // Count the number of traces on each partition
            for(i = 0; i < locexp.size(); ++i)
            {
                tCnt += locexp[i]->GetNtraces();
            }
            
            // Set up the offset and the array that will contain the list of
            // edge IDs, then reduce this across processors.
            Array<OneD, int> tracesCnt(nproc, 0);
            tracesCnt[tracepr] = tCnt;
            m_comm->AllReduce(tracesCnt, LibUtilities::ReduceSum);
            
            // Set up offset array.
            int totTraceCnt = Vmath::Vsum(nproc, tracesCnt, 1);
            Array<OneD, int> tTotOffsets(nproc,0);

            for (i = 1; i < nproc; ++i)
            {
                tTotOffsets[i] = tTotOffsets[i-1] + tracesCnt[i-1];
            }
            
            // Local list of the edges per element
            Array<OneD, int> TracesTotID(totTraceCnt, 0);
            Array<OneD, int> TracesTotNm0(totTraceCnt, 0);
            Array<OneD, int> TracesTotNm1(totTraceCnt, 0);
            Array<OneD, int> TracesTotPnts0(totTraceCnt, 0);
            Array<OneD, int> TracesTotPnts1(totTraceCnt, 0);
            
            int cntr = tTotOffsets[tracepr];
            
            for(i = 0; i < locexp.size(); ++i)
            {
                if((exp2D = locexp[i]->as<LocalRegions::Expansion2D>()))
                {
                    
                    int nedges = locexp[i]->GetNtraces();
                
                    for(j = 0; j < nedges; ++j, ++cntr)
                    {
                        LibUtilities::BasisKey bkeyEdge =
                            locexp[i]->GetTraceBasisKey(j);
                        TracesTotID   [cntr] = exp2D->GetGeom2D()->GetEid(j);
                        TracesTotNm0  [cntr] = bkeyEdge.GetNumModes();
                        TracesTotPnts0[cntr] = bkeyEdge.GetNumPoints();
                    }
                }
                else if((exp3D = locexp[i]->as<LocalRegions::Expansion3D>()))
                {
                    int nfaces = locexp[i]->GetNtraces();
                    
                    for(j = 0; j < nfaces; ++j, ++cntr)
                    {
                        LibUtilities::BasisKey face_dir0
                            = locexp[i]->GetTraceBasisKey(j,0);
                        LibUtilities::BasisKey face_dir1
                            = locexp[i]->GetTraceBasisKey(j,1);

                        TracesTotID[cntr]    = exp3D->GetGeom3D()->GetFid(j);
                        TracesTotNm0[cntr]   = face_dir0.GetNumModes ();
                        TracesTotNm1[cntr]   = face_dir1.GetNumModes ();
                        TracesTotPnts0[cntr] = face_dir0.GetNumPoints();
                        TracesTotPnts1[cntr] = face_dir1.GetNumPoints();
                    }
                }
            }

            m_comm->AllReduce(TracesTotID,    LibUtilities::ReduceSum);
            m_comm->AllReduce(TracesTotNm0,   LibUtilities::ReduceSum);
            m_comm->AllReduce(TracesTotPnts0, LibUtilities::ReduceSum);
            if(m_expType == e2D)
            {
                m_comm->AllReduce(TracesTotNm1,   LibUtilities::ReduceSum);
                m_comm->AllReduce(TracesTotPnts1, LibUtilities::ReduceSum);
            }
            
            if(edgeOrders.size())
            {
                for (i = 0; i < totTraceCnt; ++i)
                {
                    auto it = edgeOrders.find(TracesTotID[i]);
                    
                    if (it == edgeOrders.end())
                    {
                        continue;
                    }
                    
                    LibUtilities::BasisKey existing
                        = it->second.second;
                    LibUtilities::BasisKey edge(existing.GetBasisType(),
                                                TracesTotNm0[i],
                                                LibUtilities::PointsKey(
                                                               TracesTotPnts0[i],
                                                existing.GetPointsType()));
                    
                    int np1 = edge    .GetNumPoints();
                    int np2 = existing.GetNumPoints();
                    int nm1 = edge    .GetNumModes ();
                    int nm2 = existing.GetNumModes ();
                    
                    if (np2 >= np1 && nm2 >= nm1)
                    {
                        continue;
                    }
                    else if (np2 < np1 && nm2 < nm1)
                    {
                        it->second.second = edge;
                    }
                    else
                    {
                        ASSERTL0(false,
                                 "inappropriate number of points/modes (max "
                                 "num of points is not set with max order)");
                    }
                }
            }
            else if(faceOrders.size())
            {
                for (i = 0; i < totTraceCnt; ++i)
                {
                    auto it = faceOrders.find(TracesTotID[i]);

                    if (it == faceOrders.end())
                    {
                        continue;
                    }

                    LibUtilities::BasisKey existing0 =
                        it->second.second.first;
                    LibUtilities::BasisKey existing1 =
                        it->second.second.second;
                    LibUtilities::BasisKey face0(
                        existing0.GetBasisType(), TracesTotNm0[i],
                        LibUtilities::PointsKey(TracesTotPnts0[i],
                                                existing0.GetPointsType()));
                    LibUtilities::BasisKey face1(
                        existing1.GetBasisType(), TracesTotNm1[i],
                        LibUtilities::PointsKey(TracesTotPnts1[i],
                                                existing1.GetPointsType()));

                    int np11 = face0    .GetNumPoints();
                    int np12 = face1    .GetNumPoints();
                    int np21 = existing0.GetNumPoints();
                    int np22 = existing1.GetNumPoints();
                    int nm11 = face0    .GetNumModes ();
                    int nm12 = face1    .GetNumModes ();
                    int nm21 = existing0.GetNumModes ();
                    int nm22 = existing1.GetNumModes ();

                    if ((np22 >= np12 || np21 >= np11) &&
                        (nm22 >= nm12 || nm21 >= nm11))
                    {
                        continue;
                    }
                    else if((np22 < np12 || np21 < np11) &&
                            (nm22 < nm12 || nm21 < nm11))
                    {
                        it->second.second.first  = face0;
                        it->second.second.second = face1;
                    }
                    else
                    {
                        ASSERTL0(false,
                                 "inappropriate number of points/modes (max "
                                 "num of points is not set with max order)");
                    }
                }
            }
        }

        if(edgeOrders.size())
        {
            for (auto &it : edgeOrders)
            {
                exp = MemoryManager<LocalRegions::SegExp>
                    ::AllocateSharedPtr(it.second.second, it.second.first);
                exp->SetElmtId(elmtid++);
                (*m_exp).push_back(exp);
            }        
        }
        else
        {
            for (auto &it : faceOrders)
            {
                FaceGeom = it.second.first;
            
                if ((QuadGeom = std::dynamic_pointer_cast<
                     SpatialDomains::QuadGeom>(FaceGeom)))
                {
                    exp = MemoryManager<LocalRegions::QuadExp>
                        ::AllocateSharedPtr(it.second.second.first,
                                            it.second.second.second,
                                            QuadGeom);
                }
                else if ((TriGeom = std::dynamic_pointer_cast<
                          SpatialDomains::TriGeom>(FaceGeom)))
                {
                    exp = MemoryManager<LocalRegions::TriExp>
                        ::AllocateSharedPtr(it.second.second.first,
                                            it.second.second.second,
                                            TriGeom);
                }
                exp->SetElmtId(elmtid++);
                (*m_exp).push_back(exp);
            }
        }
        
        // Set up m_coeffs, m_phys and offset arrays.
        SetupCoeffPhys(DeclareCoeffPhysArrays);
        
        
        // Set up collections
        if(m_expType != e0D)
        {
            CreateCollections(ImpType);
        }
    }

    /**
     * Fills the list of local expansions with the trace from the mesh
     * specified by \a domain. This CompositeMap contains a list of
     * Composites which define the boundary. It is also used to set up
     * expansion domains in the 1D Pulse Wave solver. 
     *
     * @param  pSession     A session within information about expansion
     * @param  domain       A domain, comprising of one or more composite
     *                      regions,
     * @param  graph        A mesh, containing information about the
     *                      domain and the spectral/hp element expansion.
     * @param DeclareCoeffPhysArrays Declare the coefficient and
     *                               phys space arrays. Default is true. 
     * @param  variable     The variable name associated with the expansion
     * @param  SetToOneSpaceDimension Reduce to one space dimension expansion
     * @param  comm         An optional communicator that can be used with the
     *                      boundary expansion in case of more global
     *                      parallel operations. Default to a Null Communicator
     * @param  ImpType      Detail about the implementation type to use 
     *                      in operators. Default is eNoImpType. 
     *
     */
    ExpList::ExpList(const LibUtilities::SessionReaderSharedPtr &pSession,
                     const SpatialDomains::CompositeMap &domain,
                     const SpatialDomains::MeshGraphSharedPtr &graph,
                     const bool DeclareCoeffPhysArrays,
                     const std::string variable,
                     bool SetToOneSpaceDimension,
                     const LibUtilities::CommSharedPtr comm,
                     const Collections::ImplementationType ImpType):
        m_comm(comm), 
        m_session(pSession),
        m_graph(graph),
        m_physState(false),
        m_exp(MemoryManager<LocalRegions::ExpansionVector>
              ::AllocateSharedPtr()),
        m_blockMat(MemoryManager<BlockMatrixMap>::AllocateSharedPtr()),
        m_WaveSpace(false)
    {
        int j, elmtid=0;
        SpatialDomains::PointGeomSharedPtr PtGeom;
        SpatialDomains::SegGeomSharedPtr   SegGeom;
        SpatialDomains::TriGeomSharedPtr   TriGeom;
        SpatialDomains::QuadGeomSharedPtr  QuadGeom;

        LocalRegions::ExpansionSharedPtr  exp;

        LibUtilities::PointsType TriNb;

        int meshdim = graph->GetMeshDimension();

        // Retrieve the list of expansions (needed of meshdim == 1
        const SpatialDomains::ExpansionInfoMap &expansions
            = graph->GetExpansionInfos(variable);
        
        // Retrieve the list of expansions
        // Process each composite region.
        for(auto &compIt : domain)
        {
            // Process each expansion in the region.
            for(j = 0; j < compIt.second->m_geomVec.size(); ++j)
            {
                if((PtGeom = std::dynamic_pointer_cast <
                    SpatialDomains::PointGeom>(compIt.second->m_geomVec[j])))
                {
                    m_expType = e0D;

                    exp = MemoryManager<LocalRegions::PointExp>
                        ::AllocateSharedPtr(PtGeom);
                }
                else  if((SegGeom = std::dynamic_pointer_cast<
                    SpatialDomains::SegGeom>(compIt.second->m_geomVec[j])))
                {
                    m_expType = e1D;
                    
                    // Retrieve the basis key from the expansion.
                    LibUtilities::BasisKey bkey = LibUtilities::NullBasisKey;

                    if(meshdim == 1)
                    {
                        auto expIt = expansions.find(SegGeom->GetGlobalID());
                        ASSERTL0(expIt != expansions.end(),
                                 "Failed to find basis key");
                        bkey = expIt->second->m_basisKeyVector[0];
                    }
                    else
                    {
                        bkey = graph->GetEdgeBasisKey(SegGeom, variable);
                    }

                    if(SetToOneSpaceDimension)
                    {
                        SpatialDomains::SegGeomSharedPtr OneDSegmentGeom = 
                            SegGeom->GenerateOneSpaceDimGeom();
                        
                        exp = MemoryManager<LocalRegions::SegExp>
                            ::AllocateSharedPtr(bkey, OneDSegmentGeom);
                    }
                    else
                    {
                        
                        exp = MemoryManager<LocalRegions::SegExp>
                            ::AllocateSharedPtr(bkey, SegGeom);
                    }
                }
                else if ((TriGeom = std::dynamic_pointer_cast<
                         SpatialDomains::TriGeom>(compIt.second->m_geomVec[j])))
                {
                    m_expType = e2D;
                    
                    LibUtilities::BasisKey TriBa
                        = graph->GetFaceBasisKey(TriGeom,0,variable);
                    LibUtilities::BasisKey TriBb
                        = graph->GetFaceBasisKey(TriGeom,1,variable);
                    
                    if (graph->GetExpansionInfos().begin()->second->
                        m_basisKeyVector[0].GetBasisType() == 
                        LibUtilities::eGLL_Lagrange)
                    {
                        ASSERTL0(false,"This method needs sorting");
                        TriNb = LibUtilities::eNodalTriElec;
                        
                        exp = MemoryManager<LocalRegions::NodalTriExp>
                            ::AllocateSharedPtr(TriBa,TriBb,TriNb,
                                                TriGeom);
                    }
                    else
                    {
                        exp = MemoryManager<LocalRegions::TriExp>
                            ::AllocateSharedPtr(TriBa, TriBb, TriGeom);
                    }
                }
                else if ((QuadGeom = std::dynamic_pointer_cast<
                        SpatialDomains::QuadGeom>(compIt.second->m_geomVec[j])))
                {
                    m_expType = e2D;

                    LibUtilities::BasisKey QuadBa
                        = graph->GetFaceBasisKey(QuadGeom, 0, variable);
                    LibUtilities::BasisKey QuadBb
                        = graph->GetFaceBasisKey(QuadGeom, 1, variable);
                    
                    exp = MemoryManager<LocalRegions::QuadExp>
                        ::AllocateSharedPtr(QuadBa, QuadBb, QuadGeom);
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a Geom (possibly 3D) failed");
                }

                exp->SetElmtId(elmtid++);
                (*m_exp).push_back(exp);
            }
        }
        
        // Set up m_coeffs, m_phys and offset arrays.
        SetupCoeffPhys(DeclareCoeffPhysArrays);

        if(m_expType != e0D)
        {
            CreateCollections(ImpType);
        }
    }

}
}
