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
#include <LocalRegions/SegExp.h>
#include <LocalRegions/QuadExp.h>

using namespace std;

namespace Nektar
{
namespace MultiRegions
{
    
    /**
     * Given a mesh \a graph1D, containing information about the domain and
     * the spectral/hp element expansion, this constructor fills the list
     * of local expansions \texttt{m_exp} with the proper expansions,
     * calculates the total number of quadrature points \f$x_i\f$ and local
     * expansion coefficients \f$\hat{u}^e_n\f$ and allocates memory for
     * the arrays #m_coeffs and #m_phys.
     *
     * @param  pSession    A session within information about expansion
     *
     * @param  graph       A mesh, containing information about the
     *                      domain and the spectral/hp element expansion.
     *
     * @param DeclareCoeffPhysArrays Declare the coefficient and
     *                               phys space arrays
     *
     * @param  ImpType     Detail about the implementation type to use 
     *                     in operators
     */
    ExpList::ExpList(const LibUtilities::SessionReaderSharedPtr &pSession,
                     const SpatialDomains::MeshGraphSharedPtr &graph,
                     const bool DeclareCoeffPhysArrays,
                     const Collections::ImplementationType ImpType):
        m_comm(pSession->GetComm()),
        m_session(pSession),
        m_graph(graph),
        m_ncoeffs(0),
        m_npoints(0),
        m_physState(false),
        m_exp(MemoryManager<LocalRegions::ExpansionVector>
              ::AllocateSharedPtr()),
        m_blockMat(MemoryManager<BlockMatrixMap>::AllocateSharedPtr()),
        m_WaveSpace(false)
    {
        int id=0;
        LocalRegions::SegExpSharedPtr seg;
        SpatialDomains::SegGeomSharedPtr SegmentGeom;
        
        // Retrieve the list of expansions
        const SpatialDomains::ExpansionMap &expansions
            = graph->GetExpansions();
        
        // Process each expansion in the graph
        for (auto &expIt : expansions)
        {
            // Retrieve basis key from expansion
            LibUtilities::BasisKey bkey = expIt.second->m_basisKeyVector[0];
            
            if ((SegmentGeom = std::dynamic_pointer_cast<
                 SpatialDomains::SegGeom>(expIt.second->m_geomShPtr)))
            {
                m_expType = e1D;
                
                seg = MemoryManager<LocalRegions::SegExp>
                    ::AllocateSharedPtr(bkey, SegmentGeom);
                
                // Assign next ID
                seg->SetElmtId(id++);
                
                // Add the expansion
                (*m_exp).push_back(seg);
            }
            else
            {
                ASSERTL0(false,"dynamic cast to a SegGeom failed");
            }
        }
        
        // set up offset arrays.
        SetCoeffPhysOffsets();
        
        if(DeclareCoeffPhysArrays)
        {
            // Set up m_coeffs, m_phys.
            m_coeffs = Array<OneD, NekDouble>(m_ncoeffs,0.0);
            m_phys   = Array<OneD, NekDouble>(m_npoints,0.0);
        }
        
        CreateCollections(ImpType);
    }
    

    /**
     * Given a meshgraph \a graph, and the spectral/hp element
     * expansion as well as and separate information about a \a
     * domain, this constructor fills the list of local
     * expansions \texttt{m_exp} with the proper expansions,
     * calculates the total number of quadrature points \f$x_i\f$
     * and local expansion coefficients \f$\hat{u}^e_n\f$ and
     * allocates memory for the arrays #m_coeffs and #m_phys.
     *
     *
     * @param  pSession     A session within information about expansion
     * @param   graph       A mesh, containing information about the
     *                      domain and the spectral/hp element expansion.
     * @param  domain       A Composite list describing the domain of the 
     *                      expansion 
     * @param DeclareCoeffPhysArrays Declare the coefficient and
     *                               phys space arrays
     * @param  var          The variable name associated with the expansion
     * @param  SetToOneSpaceDimnesion Reduce to one space dimension expansion
     * @param  ImpType      Detail about the implementation type to use 
     *                      in operators
     */
    ExpList::ExpList(const LibUtilities::SessionReaderSharedPtr &pSession,
                         const SpatialDomains::MeshGraphSharedPtr &graph1D,
                         const SpatialDomains::CompositeMap &domain,
                         const bool DeclareCoeffPhysArrays,
                         const std::string var,
                         bool SetToOneSpaceDimension,
                         const Collections::ImplementationType ImpType):
        m_comm(pSession->GetComm()),
        m_session(pSession),
        m_graph(graph1D),
        m_ncoeffs(0),
        m_npoints(0),
        m_physState(false),
        m_exp(MemoryManager<LocalRegions::ExpansionVector>
              ::AllocateSharedPtr()),
        m_blockMat(MemoryManager<BlockMatrixMap>::AllocateSharedPtr()),
        m_WaveSpace(false)
    {
        int j,id=0;
        LocalRegions::SegExpSharedPtr seg;
        SpatialDomains::SegGeomSharedPtr SegmentGeom;
        SpatialDomains::Composite comp;

        // Retrieve the list of expansions
        const SpatialDomains::ExpansionMap &expansions
            = graph1D->GetExpansions(var);
        
        // Process each composite region in domain
        for(auto &compIt : domain)
        {
            // Process each expansion in the graph
            for(j = 0; j < compIt.second->m_geomVec.size(); ++j)
            {
                if((SegmentGeom = std::dynamic_pointer_cast<
                    SpatialDomains::SegGeom>(compIt.second->m_geomVec[j])))
                {
                    m_expType = e1D;
                    // Retrieve basis key from expansion and
                    // define expansion
                    auto expIt = expansions.find(SegmentGeom->GetGlobalID());
                    if(expIt != expansions.end())
                    {
                        LibUtilities::BasisKey bkey =
                            expIt->second->m_basisKeyVector[0];
                        
                        if(SetToOneSpaceDimension)
                        {
                            SpatialDomains::SegGeomSharedPtr OneDSegmentGeom = 
                                SegmentGeom->GenerateOneSpaceDimGeom();
                            
                            seg = MemoryManager<LocalRegions::SegExp>
                                ::AllocateSharedPtr(bkey, OneDSegmentGeom);
                        }
                        else
                        {
                            seg = MemoryManager<LocalRegions::SegExp>
                                ::AllocateSharedPtr(bkey, SegmentGeom);
                        }
                    }
                    else
                    {
                        ASSERTL0(false,"Failed to find basis key");
                    }
                }
                else
                {
                    ASSERTL0(false,"Failed to dynamic cast geometry to SegGeom");
                }
                
                // Assign next ID
                seg->SetElmtId(id++);
                
                // Add the expansion
                (*m_exp).push_back(seg);
            }
        }
        
        // set up offset arrays.
        SetCoeffPhysOffsets();
        
        if(DeclareCoeffPhysArrays)
        {
            // Set up m_coeffs, m_phys.
            m_coeffs = Array<OneD, NekDouble>(m_ncoeffs,0.0);
            m_phys   = Array<OneD, NekDouble>(m_npoints,0.0);
        }
        
        CreateCollections(ImpType);
    }
    

    /**
     * After initialising the basic data populate the expansion list
     * from the segments defined in the supplied
     * SpatialDomains#MeshGraph1D. All expansions in the graph are
     * defined using the same LibUtilities#BasisKey which overrides
     * that specified in \a graph1D.
     *
     * @param  pSession     A session within information about expansion
     * @param  Ba           BasisKey describing quadrature points and
     *                      number of modes to impose on this ExpList. 
     * @param  graph1D      Domain and expansion definitions.
     * @param  ImpType      Detail about the implementation type to use 
     *                      in operators
     */
    ExpList::ExpList(const LibUtilities::SessionReaderSharedPtr &pSession,
                         const LibUtilities::BasisKey &Ba,
                         const SpatialDomains::MeshGraphSharedPtr &graph1D,
                         const Collections::ImplementationType ImpType):
        m_comm(pSession->GetComm()),
        m_session(pSession),
        m_graph(graph1D),
        m_ncoeffs(0),
        m_npoints(0),
        m_physState(false),
        m_exp(MemoryManager<LocalRegions::ExpansionVector>
              ::AllocateSharedPtr()),
        m_blockMat(MemoryManager<BlockMatrixMap>::AllocateSharedPtr()),
        m_WaveSpace(false)
    {
        m_expType = e1D; 
        
        int id=0;
        LocalRegions::SegExpSharedPtr seg;
        SpatialDomains::SegGeomSharedPtr SegmentGeom;
        
        const SpatialDomains::ExpansionMap &expansions
            = graph1D->GetExpansions();
        
        // For each element in the mesh, create a segment expansion using
        // the supplied BasisKey and segment geometry.
        for (auto &expIt : expansions)
        {
            if ((SegmentGeom = std::dynamic_pointer_cast<
                 SpatialDomains::SegGeom>(expIt.second->m_geomShPtr)))
            {
                seg = MemoryManager<LocalRegions::SegExp>
                    ::AllocateSharedPtr(Ba,SegmentGeom);
                seg->SetElmtId(id++);
                (*m_exp).push_back(seg);
            }
            else
            {
                ASSERTL0(false,"dynamic cast to a SegGeom failed");
            }
        }
        
        // Allocate storage for data and populate element offset lists.
        SetCoeffPhysOffsets();
        
        m_coeffs = Array<OneD, NekDouble>(m_ncoeffs);
        m_phys   = Array<OneD, NekDouble>(m_npoints);
        
        CreateCollections(ImpType);
    }
    
    /**
     * Fills the list of local expansions with the segments from the 2D
     * mesh specified by \a domain. This CompositeMap contains a list of
     * Composites which define the Neumann boundary.
     *
     * @param  pSession     A session within information about expansion
     * @param  domain       A domain, comprising of one or more composite
     *                      regions,
     * @param  graph2D      A mesh, containing information about the
     *                      domain and the spectral/hp element expansion.
     * @param DeclareCoeffPhysArrays Declare the coefficient and
     *                               phys space arrays
     * @param  variable     The variable name associated with the expansion
     * @param  comm         A communicator for this expansion
     * @param  ImpType      Detail about the implementation type to use 
     *                      in operators
     *
     * Note: comm is available thorugh the pSession !!
     */
    ExpList::ExpList(const LibUtilities::SessionReaderSharedPtr &pSession,
                     const SpatialDomains::CompositeMap &domain,
                     const SpatialDomains::MeshGraphSharedPtr &graph2D,
                     const bool DeclareCoeffPhysArrays,
                     const std::string variable,
                     const LibUtilities::CommSharedPtr comm,
                     const Collections::ImplementationType ImpType):
        m_comm(pSession->GetComm()), 
        m_session(pSession),
        m_graph(graph2D),
        m_ncoeffs(0),
        m_npoints(0),
        m_physState(false),
        m_exp(MemoryManager<LocalRegions::ExpansionVector>
              ::AllocateSharedPtr()),
        m_blockMat(MemoryManager<BlockMatrixMap>::AllocateSharedPtr()),
        m_WaveSpace(false)
    {
        if (comm)
        {
            m_comm = comm;
        }

        int j, id=0;
        SpatialDomains::Composite comp;
        SpatialDomains::SegGeomSharedPtr SegmentGeom;
        LocalRegions::SegExpSharedPtr seg;
        
        // Process each composite region.
        for(auto &compIt : domain)
        {
            // Process each expansion in the region.
            for(j = 0; j < compIt.second->m_geomVec.size(); ++j)
            {
                if((SegmentGeom = std::dynamic_pointer_cast<
                    SpatialDomains::SegGeom>(compIt.second->m_geomVec[j])))
                {
                    m_expType = e1D;
                    
                    // Retrieve the basis key from the expansion.
                    LibUtilities::BasisKey bkey = graph2D->GetEdgeBasisKey(
                                               SegmentGeom, variable);
                    
                    seg = MemoryManager<LocalRegions::SegExp>
                        ::AllocateSharedPtr(bkey, SegmentGeom);
                    
                    // Add the segment to the expansion list.
                    seg->SetElmtId(id++);
                    (*m_exp).push_back(seg);
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a SegGeom failed");
                }
            }
        }
        
        // Allocate storage for data and populate element offset lists.
        SetCoeffPhysOffsets();

            // Set up m_coeffs, m_phys.
        if(DeclareCoeffPhysArrays)
        {
            m_coeffs = Array<OneD, NekDouble>(m_ncoeffs,0.0);
            m_phys   = Array<OneD, NekDouble>(m_npoints,0.0);
        }
        
        CreateCollections(ImpType);
    }
    
    /**
     * Store expansions for the trace space expansions used in
     * DisContField2D.
     *
     * @param  pSession      A session within information about expansion
     * @param  bndConstraint Array of ExpList1D objects each containing a
     *                       1D spectral/hp element expansion on a single
     *                       boundary region.
     * @param  bndCond       Array of BoundaryCondition objects which contain
     *                       information about the boundary conditions on the
     *                       different boundary regions.
     * @param  locexp        Complete domain expansion list.
     * @param  graph2D       2D mesh corresponding to the expansion list.
     * @param  periodicEdges List of periodic edges.
     * @param DeclareCoeffPhysArrays Declare the coefficient and
     *                               phys space arrays
     * @param  variable      The variable name associated with the expansion
     * @param  ImpType       Detail about the implementation type to use 
     *                       in operators
     */
    ExpList::ExpList(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const Array<OneD,const ExpListSharedPtr>   &bndConstraint,
            const Array<OneD, const SpatialDomains::BoundaryConditionShPtr>  &bndCond,
            const LocalRegions::ExpansionVector &locexp,
            const SpatialDomains::MeshGraphSharedPtr &graph2D,
            const PeriodicMap &periodicEdges,
            const bool DeclareCoeffPhysArrays,
            const std::string variable,
            const Collections::ImplementationType ImpType):
        m_comm(pSession->GetComm()), 
        m_session(pSession),
        m_graph(graph2D),
        m_ncoeffs(0),
        m_npoints(0),
        m_physState(false),
        m_exp(MemoryManager<LocalRegions::ExpansionVector>
              ::AllocateSharedPtr()),
        m_blockMat(MemoryManager<BlockMatrixMap>::AllocateSharedPtr()),
        m_WaveSpace(false)
    {
        int i, j, id, elmtid = 0;
        set<int> edgesDone;
        
        SpatialDomains::Geometry1DSharedPtr segGeom;
        SpatialDomains::Geometry2DSharedPtr ElGeom;
        LocalRegions::SegExpSharedPtr       seg;
        LocalRegions::SegExpSharedPtr       seg_tmp;
        LocalRegions::Expansion1DSharedPtr  exp1D;
        LocalRegions::Expansion2DSharedPtr  exp2D;
        
        map<int,int> EdgeDone;
        map<int,int> NormalSet;
        
        LocalRegions::SegExpSharedPtr Seg;
        
        m_expType = e1D;

        // First loop over boundary conditions to renumber
        // Dirichlet boundaries
        for(i = 0; i < bndCond.num_elements(); ++i)
        {
            if(bndCond[i]->GetBoundaryConditionType()
               == SpatialDomains::eDirichlet)
            {
                for(j = 0; j < bndConstraint[i]->GetExpSize(); ++j)
                {
                    LibUtilities::BasisKey bkey = bndConstraint[i]
                        ->GetExp(j)->GetBasis(0)->GetBasisKey();
                    exp1D = bndConstraint[i]->GetExp(j)->
                        as<LocalRegions::Expansion1D>();
                    segGeom = exp1D->GetGeom1D();
                    
                    seg = MemoryManager<LocalRegions::SegExp>
                        ::AllocateSharedPtr(bkey, segGeom);
                    edgesDone.insert(segGeom->GetGlobalID());
                    
                    seg->SetElmtId(elmtid++);
                    (*m_exp).push_back(seg);
                }
            }
        }

        map<int, pair<SpatialDomains::Geometry1DSharedPtr,
                      LibUtilities::BasisKey> > edgeOrders;
        
        for(i = 0; i < locexp.size(); ++i)
        {
            exp2D = locexp[i]->as<LocalRegions::Expansion2D>();
            
            for(j = 0; j < locexp[i]->GetNedges(); ++j)
            {
                segGeom = exp2D->GetGeom2D()->GetEdge(j);
                id      = segGeom->GetGlobalID();
                // Ignore Dirichlet edges
                if (edgesDone.count(id) != 0)
                {
                    continue;
                }
                
                auto it = edgeOrders.find(id);
                
                if (it == edgeOrders.end())
                {
                    edgeOrders.insert(std::make_pair(id, std::make_pair(
                               segGeom, locexp[i]->DetEdgeBasisKey(j))));
                }
                else // variable modes/points
                {
                    LibUtilities::BasisKey edge
                        = locexp[i]->DetEdgeBasisKey(j);
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
                                 "inappropriate number of points/modes (max "
                                 "num of points is not set with max order)");
                    }
                }
            }
        }
        
        LibUtilities::CommSharedPtr vComm =
            pSession->GetComm()->GetRowComm();
        int nproc = vComm->GetSize(); // number of processors
        int edgepr = vComm->GetRank(); // ID processor
        
        if (nproc > 1)
        {
            int eCnt = 0;
            
            // Count the number of edges on each partition
            for(i = 0; i < locexp.size(); ++i)
            {
                eCnt += locexp[i]->GetNedges();
            }
            
            // Set up the offset and the array that will contain the list of
            // edge IDs, then reduce this across processors.
            Array<OneD, int> edgesCnt(nproc, 0);
            edgesCnt[edgepr] = eCnt;
            vComm->AllReduce(edgesCnt, LibUtilities::ReduceSum);
            
            // Set up offset array.
            int totEdgeCnt = Vmath::Vsum(nproc, edgesCnt, 1);
            Array<OneD, int> eTotOffsets(nproc,0);
            for (i = 1; i < nproc; ++i)
            {
                eTotOffsets[i] = eTotOffsets[i-1] + edgesCnt[i-1];
            }
            
            // Local list of the edges per element
            Array<OneD, int> EdgesTotID(totEdgeCnt, 0);
            Array<OneD, int> EdgesTotNm(totEdgeCnt, 0);
            Array<OneD, int> EdgesTotPnts(totEdgeCnt, 0);
            
            int cntr = eTotOffsets[edgepr];
            
            for(i = 0; i < locexp.size(); ++i)
            {
                exp2D = locexp[i]->as<LocalRegions::Expansion2D>();
                
                int nedges = locexp[i]->GetNedges();
                
                for(j = 0; j < nedges; ++j, ++cntr)
                {
                    LibUtilities::BasisKey bkeyEdge =
                        locexp[i]->DetEdgeBasisKey(j);
                    EdgesTotID  [cntr] = exp2D->GetGeom2D()->GetEid(j);
                    EdgesTotNm  [cntr] = bkeyEdge.GetNumModes();
                    EdgesTotPnts[cntr] = bkeyEdge.GetNumPoints();
                }
            }
            
            vComm->AllReduce(EdgesTotID, LibUtilities::ReduceSum);
            vComm->AllReduce(EdgesTotNm, LibUtilities::ReduceSum);
            vComm->AllReduce(EdgesTotPnts, LibUtilities::ReduceSum);
            
            for (i = 0; i < totEdgeCnt; ++i)
            {
                auto it = edgeOrders.find(EdgesTotID[i]);
                
                if (it == edgeOrders.end())
                {
                    continue;
                }
                
                LibUtilities::BasisKey existing
                    = it->second.second;
                LibUtilities::BasisKey edge(
                          existing.GetBasisType(), EdgesTotNm[i],
                          LibUtilities::PointsKey(EdgesTotPnts[i],
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
        
        for (auto &it : edgeOrders)
        {
            seg = MemoryManager<LocalRegions::SegExp>
                ::AllocateSharedPtr(it.second.second, it.second.first);
            seg->SetElmtId(elmtid++);
            (*m_exp).push_back(seg);
        }
        
        // Set up offset information and array sizes
        SetCoeffPhysOffsets();
        
        // Set up m_coeffs, m_phys.
        if(DeclareCoeffPhysArrays)
        {
            m_coeffs = Array<OneD, NekDouble>(m_ncoeffs,0.0);
            m_phys   = Array<OneD, NekDouble>(m_npoints,0.0);
        }
        
        CreateCollections(ImpType);
    }
}
}
