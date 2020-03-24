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
        
        SetupCoeffPhys(true);
    }

    /**
     * Store expansions for the trace space expansions used in
     * DisContField1D.
     *
     * @param   bndConstraint   Array of ExpList1D objects each containing a
     *                      1D spectral/hp element expansion on a single
     *                      boundary region.
     * @param   bndCond     Array of BoundaryCondition objects which contain
     *                      information about the boundary conditions on the
     *                      different boundary regions.
     * @param   locexp      Complete domain expansion list.
     * @param   graph1D     1D mesh corresponding to the expansion list.
     * @param   UseGenSegExp If true, create general segment expansions
     *                      instead of just normal segment expansions.
     */
    ExpList::ExpList(
            const Array<OneD, const ExpListSharedPtr> &bndConstraint,
            const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> 
                                                      &bndCond,
            const LocalRegions::ExpansionVector       &locexp,
            const SpatialDomains::MeshGraphSharedPtr  &graph1D,
            const bool                                 DeclareCoeffPhysArrays):
        m_expType(e0D),
        m_ncoeffs(1),
        m_npoints(1),
        m_physState(false),
        m_exp(MemoryManager<LocalRegions::ExpansionVector>
              ::AllocateSharedPtr()),
        m_blockMat(MemoryManager<BlockMatrixMap>::AllocateSharedPtr()),
        m_WaveSpace(false)
    {
        int i, j, id, elmtid=0;
        map<int,int> EdgeDone;
        map<int,int> NormalSet;
        
        SpatialDomains::PointGeomSharedPtr PointGeom;
        LocalRegions::PointExpSharedPtr Point;
        LocalRegions::Expansion1DSharedPtr exp;
	
        // First loop over boundary conditions to renumber Dirichlet boundaries
        for(i = 0; i < bndCond.num_elements(); ++i)
        {
            if(bndCond[i]->GetBoundaryConditionType() ==
               SpatialDomains::eDirichlet)
            {
                for(j = 0; j < bndConstraint[i]->GetExpSize(); ++j)
                {
                    PointGeom = bndConstraint[i]->GetExp(0)->GetGeom()->
                        GetVertex(0);
                    Point = MemoryManager<LocalRegions::PointExp>::
                        AllocateSharedPtr(PointGeom);
                    
                    EdgeDone[PointGeom->GetVid()] = elmtid;
                    
                    Point->SetElmtId(elmtid++);
                    (*m_exp).push_back(Point);
                }
            }
        }
	
        // loop over all other edges and fill out other connectivities
        for(i = 0; i < locexp.size(); ++i)
        {
            for(j = 0; j < 2; ++j)
            {
                exp = locexp[i]->as<LocalRegions::Expansion1D>();
                PointGeom = (exp->GetGeom1D())->GetVertex(j);
                id = PointGeom->GetVid();
		
                if(EdgeDone.count(id)==0)
                {						
                    Point = MemoryManager<LocalRegions::PointExp>::
                        AllocateSharedPtr(PointGeom);
                    EdgeDone[id] = elmtid;
                    
                    Point->SetElmtId(elmtid++);
                    (*m_exp).push_back(Point);
                }
            }
        }
		 
        // Set up offset information and array sizes
        SetupCoeffPhys(DeclareCoeffPhysArrays);
    }

    //----------------------------------------------------------------------
    //                        1D Expansion Constructors 
    //----------------------------------------------------------------------
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
                     const std::string &var,
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
        // Retrieve the list of expansions
        const SpatialDomains::ExpansionMap &expansions
            = graph->GetExpansions(var);

        int id=0;
        LocalRegions::SegExpSharedPtr      seg;
        LocalRegions::TriExpSharedPtr      tri;
        LocalRegions::NodalTriExpSharedPtr Ntri;
        LibUtilities::PointsType           TriNb;
        LocalRegions::QuadExpSharedPtr     quad;
        SpatialDomains::Composite          comp;
        
        SpatialDomains::SegGeomSharedPtr  SegmentGeom;
        SpatialDomains::TriGeomSharedPtr  TriangleGeom;
        SpatialDomains::QuadGeomSharedPtr QuadrilateralGeom;

        
        // Process each expansion in the graph
        for (auto &expIt : expansions)
        {
            
            if ((SegmentGeom = std::dynamic_pointer_cast<
                 SpatialDomains::SegGeom>(expIt.second->m_geomShPtr)))
            {
                m_expType = e1D;
                
                // Retrieve basis key from expansion
                LibUtilities::BasisKey bkey = expIt.second->m_basisKeyVector[0];

                seg = MemoryManager<LocalRegions::SegExp>
                    ::AllocateSharedPtr(bkey, SegmentGeom);
                
                // Assign next ID
                seg->SetElmtId(id++);
                
                // Add the expansion
                (*m_exp).push_back(seg);
            }
            else if ((TriangleGeom = std::dynamic_pointer_cast<SpatialDomains
                        ::TriGeom>(expIt.second->m_geomShPtr)))
            {
                m_expType = e2D;

                LibUtilities::BasisKey TriBa
                    = expIt.second->m_basisKeyVector[0];
                LibUtilities::BasisKey TriBb
                    = expIt.second->m_basisKeyVector[1];

                // This is not elegantly implemented needs re-thinking.
                if (TriBa.GetBasisType() == LibUtilities::eGLL_Lagrange)
                {
                    LibUtilities::BasisKey newBa(LibUtilities::eOrtho_A,
                                                 TriBa.GetNumModes(),
                                                 TriBa.GetPointsKey());
                    
                    TriNb = LibUtilities::eNodalTriElec;
                    Ntri = MemoryManager<LocalRegions::NodalTriExp>
                        ::AllocateSharedPtr(newBa,TriBb,TriNb,
                                            TriangleGeom);
                    Ntri->SetElmtId(id++);
                    (*m_exp).push_back(Ntri);
                }
                else
                {
                    tri = MemoryManager<LocalRegions::TriExp>
                        ::AllocateSharedPtr(TriBa,TriBb,
                                            TriangleGeom);
                    tri->SetElmtId(id++);
                    (*m_exp).push_back(tri);
                }
            }
            else if ((QuadrilateralGeom = std::dynamic_pointer_cast<
                      SpatialDomains::QuadGeom>(expIt.second->m_geomShPtr)))
            {
                m_expType = e2D;
                LibUtilities::BasisKey QuadBa
                    = expIt.second->m_basisKeyVector[0];
                LibUtilities::BasisKey QuadBb
                    = expIt.second->m_basisKeyVector[1];
                
                quad = MemoryManager<LocalRegions::QuadExp>
                    ::AllocateSharedPtr(QuadBa,QuadBb,
                                        QuadrilateralGeom);
                quad->SetElmtId(id++);
                (*m_exp).push_back(quad);

            }
            else
            {
                ASSERTL0(false,"dynamic cast to a 1D or 2D Geom failed");
            }
        }
        
        // set up offset arrays.
        SetupCoeffPhys(DeclareCoeffPhysArrays);
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
        SetupCoeffPhys(DeclareCoeffPhysArrays);
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
        SetupCoeffPhys(true);

        CreateCollections(ImpType);
    }
    
    //-----------------------------------------------------------------------------
    //                         1D and 2D Constructors
    //-----------------------------------------------------------------------------

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
            const Array<OneD, const SpatialDomains::BoundaryConditionShPtr>  &bndCond,
            const LocalRegions::ExpansionVector &locexp,
            const SpatialDomains::MeshGraphSharedPtr &graph,
            const bool DeclareCoeffPhysArrays,
            const std::string variable,
            const Collections::ImplementationType ImpType):
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
        int i, j, id, elmtid = 0;
        set<int> tracesDone;
        
        SpatialDomains::Geometry1DSharedPtr segGeom;
        SpatialDomains::Geometry2DSharedPtr ElGeom;
        SpatialDomains::Geometry2DSharedPtr FaceGeom;
        SpatialDomains::QuadGeomSharedPtr   FaceQuadGeom;
        SpatialDomains::TriGeomSharedPtr    FaceTriGeom;
        LocalRegions::SegExpSharedPtr       seg;
        LocalRegions::SegExpSharedPtr       seg_tmp;
        LocalRegions::QuadExpSharedPtr      FaceQuadExp;
        LocalRegions::TriExpSharedPtr       FaceTriExp;
        LocalRegions::Expansion1DSharedPtr  exp1D;
        LocalRegions::Expansion2DSharedPtr  exp2D;
        LocalRegions::Expansion3DSharedPtr  exp3D;
        
        m_expType = e1D;

        // First loop over boundary conditions to reorder
        // Dirichlet boundaries
        for(i = 0; i < bndCond.num_elements(); ++i)
        {
            if(bndCond[i]->GetBoundaryConditionType()
               == SpatialDomains::eDirichlet)
            {
                for(j = 0; j < bndConstraint[i]->GetExpSize(); ++j)
                {
                    if((exp1D = std::dynamic_pointer_cast<LocalRegions::Expansion1D>
                        (bndConstraint[i]->GetExp(j))))
                    {
                        LibUtilities::BasisKey bkey = bndConstraint[i]
                            ->GetExp(j)->GetBasis(0)->GetBasisKey();
                        segGeom = exp1D->GetGeom1D();
                    
                        seg = MemoryManager<LocalRegions::SegExp>
                            ::AllocateSharedPtr(bkey, segGeom);
                        tracesDone.insert(segGeom->GetGlobalID());
                    
                        seg->SetElmtId(elmtid++);
                        (*m_exp).push_back(seg);
                    }
                    else if ((exp2D = std::dynamic_pointer_cast
                              <LocalRegions::Expansion2D>(bndConstraint[i]->GetExp(j))))
                    {
                        LibUtilities::BasisKey bkey0 = bndConstraint[i]
                            ->GetExp(j)->GetBasis(0)->GetBasisKey();
                        LibUtilities::BasisKey bkey1 = bndConstraint[i]
                            ->GetExp(j)->GetBasis(1)->GetBasisKey();
                        FaceGeom = exp2D->GetGeom2D();

                        //if face is a quad
                        if((FaceQuadGeom = std::dynamic_pointer_cast<
                            SpatialDomains::QuadGeom>(FaceGeom)))
                        {
                            FaceQuadExp = MemoryManager<LocalRegions::QuadExp>
                                ::AllocateSharedPtr(bkey0, bkey1, FaceQuadGeom);
                            tracesDone.insert(FaceQuadGeom->GetGlobalID());
                            FaceQuadExp->SetElmtId(elmtid++);
                            (*m_exp).push_back(FaceQuadExp);
                        }
                        //if face is a triangle
                        else if((FaceTriGeom = std::dynamic_pointer_cast<
                                 SpatialDomains::TriGeom>(FaceGeom)))
                        {
                            FaceTriExp = MemoryManager<LocalRegions::TriExp>
                                ::AllocateSharedPtr(bkey0, bkey1, FaceTriGeom);
                            tracesDone.insert(FaceTriGeom->GetGlobalID());
                            FaceTriExp->SetElmtId(elmtid++);
                            (*m_exp).push_back(FaceTriExp);
                        }
                        else
                        {
                            ASSERTL0(false,"dynamic cast to a proper face geometry failed");
                        }

                    }
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
            if((exp2D =
                std::dynamic_pointer_cast<LocalRegions::Expansion2D>(locexp[i])))
            {
            
                for(j = 0; j < locexp[i]->GetNedges(); ++j)
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
            else if((exp3D =
                     dynamic_pointer_cast<LocalRegions::Expansion3D>(locexp[i])))
            {
               for (j = 0; j < exp3D->GetNfaces(); ++j)
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
                            = locexp[i]->DetFaceBasisKey(j,0);
                        LibUtilities::BasisKey face_dir1
                            = locexp[i]->DetFaceBasisKey(j,1);

                        faceOrders.insert(
                            std::make_pair(
                                id, std::make_pair(
                                    FaceGeom,
                                    std::make_pair(face_dir0, face_dir1))));
                    }
                    else // variable modes/points
                    {
                        LibUtilities::BasisKey face0     =
                            locexp[i]->DetFaceBasisKey(j,0);
                        LibUtilities::BasisKey face1     =
                            locexp[i]->DetFaceBasisKey(j,1);
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
                                     "inappropriate number of points/modes (max "
                                     "num of points is not set with max order)");
                        }
                    }
                }
            }
        }

        
        LibUtilities::CommSharedPtr vComm =
            pSession->GetComm()->GetRowComm();
        int nproc   = vComm->GetSize(); // number of processors
        int tracepr = vComm->GetRank(); // ID processor
        
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
            vComm->AllReduce(tracesCnt, LibUtilities::ReduceSum);
            
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
                    
                    int nedges = locexp[i]->GetNedges();
                
                    for(j = 0; j < nedges; ++j, ++cntr)
                    {
                        LibUtilities::BasisKey bkeyEdge =
                            locexp[i]->DetEdgeBasisKey(j);
                        TracesTotID   [cntr] = exp2D->GetGeom2D()->GetEid(j);
                        TracesTotNm0  [cntr] = bkeyEdge.GetNumModes();
                        TracesTotPnts0[cntr] = bkeyEdge.GetNumPoints();
                    }
                }
                else if((exp3D = locexp[i]->as<LocalRegions::Expansion3D>()))
                {
                    int nfaces = locexp[i]->GetNfaces();
                    
                    for(j = 0; j < nfaces; ++j, ++cntr)
                    {
                        LibUtilities::BasisKey face_dir0
                            = locexp[i]->DetFaceBasisKey(j,0);
                        LibUtilities::BasisKey face_dir1
                            = locexp[i]->DetFaceBasisKey(j,1);

                        TracesTotID[cntr]    = exp3D->GetGeom3D()->GetFid(j);
                        TracesTotNm0[cntr]   = face_dir0.GetNumModes ();
                        TracesTotNm1[cntr]   = face_dir1.GetNumModes ();
                        TracesTotPnts0[cntr] = face_dir0.GetNumPoints();
                        TracesTotPnts1[cntr] = face_dir1.GetNumPoints();
                    }
                }
            }

            vComm->AllReduce(TracesTotID,    LibUtilities::ReduceSum);
            vComm->AllReduce(TracesTotNm0,   LibUtilities::ReduceSum);
            vComm->AllReduce(TracesTotPnts0, LibUtilities::ReduceSum);
            if(m_expType == e2D)
            {
                vComm->AllReduce(TracesTotNm1,   LibUtilities::ReduceSum);
                vComm->AllReduce(TracesTotPnts1, LibUtilities::ReduceSum);
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
                seg = MemoryManager<LocalRegions::SegExp>
                    ::AllocateSharedPtr(it.second.second, it.second.first);
                seg->SetElmtId(elmtid++);
                (*m_exp).push_back(seg);
            }        
        }
        else
        {
            for (auto &it : faceOrders)
            {
                FaceGeom = it.second.first;
            
                if ((FaceQuadGeom = std::dynamic_pointer_cast<
                     SpatialDomains::QuadGeom>(FaceGeom)))
                {
                    FaceQuadExp = MemoryManager<LocalRegions::QuadExp>
                        ::AllocateSharedPtr(it.second.second.first,
                                            it.second.second.second,
                                            FaceQuadGeom);
                    FaceQuadExp->SetElmtId(elmtid++);
                    (*m_exp).push_back(FaceQuadExp);
                }
                else if ((FaceTriGeom = std::dynamic_pointer_cast<
                          SpatialDomains::TriGeom>(FaceGeom)))
                {
                    FaceTriExp = MemoryManager<LocalRegions::TriExp>
                        ::AllocateSharedPtr(it.second.second.first,
                                            it.second.second.second,
                                            FaceTriGeom);
                    FaceTriExp->SetElmtId(elmtid++);
                    (*m_exp).push_back(FaceTriExp);
                }
            }
        }
        
        SetupCoeffPhys(DeclareCoeffPhysArrays);

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
                     const SpatialDomains::MeshGraphSharedPtr &graph,
                     const bool DeclareCoeffPhysArrays,
                     const std::string variable,
                     const LibUtilities::CommSharedPtr comm,
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
        if (comm)
        {
            m_comm = comm;
        }

        int j, elmtid=0;
        SpatialDomains::Composite comp;
        SpatialDomains::SegGeomSharedPtr SegmentGeom;
        LocalRegions::SegExpSharedPtr seg;
        SpatialDomains::TriGeomSharedPtr TriangleGeom;
        SpatialDomains::QuadGeomSharedPtr QuadrilateralGeom;
        
        LocalRegions::TriExpSharedPtr tri;
        LocalRegions::NodalTriExpSharedPtr Ntri;
        LibUtilities::PointsType TriNb;
        LocalRegions::QuadExpSharedPtr quad;
        
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
                    LibUtilities::BasisKey bkey = graph->GetEdgeBasisKey(
                                               SegmentGeom, variable);
                    
                    seg = MemoryManager<LocalRegions::SegExp>
                        ::AllocateSharedPtr(bkey, SegmentGeom);
                    
                    // Add the segment to the expansion list.
                    seg->SetElmtId(elmtid++);
                    (*m_exp).push_back(seg);
                }
                else if ((TriangleGeom = std::dynamic_pointer_cast<
                          SpatialDomains::TriGeom>(compIt.second->m_geomVec[j])))
                {
                    m_expType = e2D;
                    
                    LibUtilities::BasisKey TriBa
                        = graph->GetFaceBasisKey(TriangleGeom, 0, variable);
                    LibUtilities::BasisKey TriBb
                        = graph->GetFaceBasisKey(TriangleGeom,1,variable);
                    
                    if (graph->GetExpansions().begin()->second->
                        m_basisKeyVector[0].GetBasisType() == 
                        LibUtilities::eGLL_Lagrange)
                    {
                        ASSERTL0(false,"This method needs sorting");
                        TriNb = LibUtilities::eNodalTriElec;
                        
                        Ntri = MemoryManager<LocalRegions::NodalTriExp>
                            ::AllocateSharedPtr(TriBa,TriBb,TriNb,
                                                TriangleGeom);
                        Ntri->SetElmtId(elmtid++);
                        (*m_exp).push_back(Ntri);
                    }
                    else
                    {
                        tri = MemoryManager<LocalRegions::TriExp>
                            ::AllocateSharedPtr(TriBa, TriBb,
                                                TriangleGeom);
                        tri->SetElmtId(elmtid++);
                        (*m_exp).push_back(tri);
                    }
                }
                else if ((QuadrilateralGeom = std::dynamic_pointer_cast<
                          SpatialDomains::QuadGeom>(compIt.second->m_geomVec[j])))
                {
                    m_expType = e2D;

                    LibUtilities::BasisKey QuadBa
                        = graph->GetFaceBasisKey(QuadrilateralGeom, 0,
                                                   variable);
                    LibUtilities::BasisKey QuadBb
                        = graph->GetFaceBasisKey(QuadrilateralGeom, 1,
                                                   variable);
                    
                    quad = MemoryManager<LocalRegions::QuadExp>
                        ::AllocateSharedPtr(QuadBa, QuadBb,
                                            QuadrilateralGeom);
                    quad->SetElmtId(elmtid++);
                    (*m_exp).push_back(quad);
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a SegGeom failed");
                }
            }
        }
        
        SetupCoeffPhys(DeclareCoeffPhysArrays);
        
        CreateCollections(ImpType);
    }

    //-----------------------------------------------------------------------------
    //                                2D Constructors
    //-----------------------------------------------------------------------------


    /**
     * Given an expansion vector \a expansions, containing
     * information about the domain and the spectral/hp element
     * expansion, this constructor fills the list of local
     * expansions \texttt{m_exp} with the proper expansions,
     * calculates the total number of quadrature points
     * \f$\boldsymbol{x}_i\f$ and local expansion coefficients
     * \f$\hat{u}^e_n\f$ and allocates memory for the arrays
     * #m_coeffs and #m_phys.
     *
     * @param  pSession      A session within information about expansion
     * @param expansions     A vector containing information about the
     *                       domain and the spectral/hp element
     *                       expansion.
     * @param DeclareCoeffPhysArrays Declare the coefficient and
     *                               phys space arrays
     * @param  ImpType       Detail about the implementation type to use 
     *                       in operators
     */
    ExpList::ExpList(const LibUtilities::SessionReaderSharedPtr &pSession,
                     const SpatialDomains::ExpansionMap &expansions,
                     const bool DeclareCoeffPhysArrays,
                     const Collections::ImplementationType ImpType):
        m_expType(e2D),
        m_comm(pSession->GetComm()), 
        m_session(pSession),
        m_ncoeffs(0),
        m_npoints(0),
        m_physState(false),
        m_exp(MemoryManager<LocalRegions::ExpansionVector>
              ::AllocateSharedPtr()),
        m_blockMat(MemoryManager<BlockMatrixMap>::AllocateSharedPtr()),
        m_WaveSpace(false)
    {
        int elmtid=0;
        LocalRegions::TriExpSharedPtr      tri;
        LocalRegions::NodalTriExpSharedPtr Ntri;
        LibUtilities::PointsType           TriNb;
        LocalRegions::QuadExpSharedPtr     quad;
        SpatialDomains::Composite          comp;
        
        for (auto &expIt : expansions)
        {
            SpatialDomains::TriGeomSharedPtr  TriangleGeom;
            SpatialDomains::QuadGeomSharedPtr QuadrilateralGeom;
            
            if ((TriangleGeom = std::dynamic_pointer_cast<SpatialDomains
                 ::TriGeom>(expIt.second->m_geomShPtr)))
            {
                LibUtilities::BasisKey TriBa
                    = expIt.second->m_basisKeyVector[0];
                LibUtilities::BasisKey TriBb
                    = expIt.second->m_basisKeyVector[1];
                
                // This is not elegantly implemented needs re-thinking.
                if (TriBa.GetBasisType() == LibUtilities::eGLL_Lagrange)
                {
                    LibUtilities::BasisKey newBa(LibUtilities::eOrtho_A,
                                                 TriBa.GetNumModes(),
                                                 TriBa.GetPointsKey());
                    
                    TriNb = LibUtilities::eNodalTriElec;
                    Ntri = MemoryManager<LocalRegions::NodalTriExp>
                        ::AllocateSharedPtr(newBa,TriBb,TriNb,
                                            TriangleGeom);
                    Ntri->SetElmtId(elmtid++);
                    (*m_exp).push_back(Ntri);
                }
                else
                {
                    tri = MemoryManager<LocalRegions::TriExp>
                        ::AllocateSharedPtr(TriBa,TriBb,
                                            TriangleGeom);
                    tri->SetElmtId(elmtid++);
                    (*m_exp).push_back(tri);
                }
            }
            else if ((QuadrilateralGeom = std::dynamic_pointer_cast<
                      SpatialDomains::QuadGeom>(expIt.second->m_geomShPtr)))
            {
                LibUtilities::BasisKey QuadBa
                    = expIt.second->m_basisKeyVector[0];
                LibUtilities::BasisKey QuadBb
                    = expIt.second->m_basisKeyVector[1];
                
                quad = MemoryManager<LocalRegions::QuadExp>
                    ::AllocateSharedPtr(QuadBa,QuadBb,
                                        QuadrilateralGeom);
                quad->SetElmtId(elmtid++);
                (*m_exp).push_back(quad);

            }
            else
            {
                ASSERTL0(false, "dynamic cast to a proper Geometry2D "
                         "failed");
            }
            
        }
        
        SetupCoeffPhys(DeclareCoeffPhysArrays);
        
        CreateCollections(ImpType);
    }

    /**
     * Given a mesh \a graph2D, containing information about the domain and
     * the a list of basiskeys, this constructor fills the list
     * of local expansions \texttt{m_exp} with the proper expansions,
     * calculates the total number of quadrature points
     * \f$\boldsymbol{x}_i\f$ and local expansion coefficients
     * \f$\hat{u}^e_n\f$ and allocates memory for the arrays #m_coeffs
     * and #m_phys.
     *
     * @param   pSession    A session within information about expansion
     * @param   TriBa       A BasisKey, containing the definition of the
     *                      basis in the first coordinate direction for
     *                      triangular elements
     * @param   TriBb       A BasisKey, containing the definition of the
     *                      basis in the second coordinate direction for
     *                      triangular elements
     * @param   QuadBa      A BasisKey, containing the definition of the
     *                      basis in the first coordinate direction for
     *                      quadrilateral elements
     * @param   QuadBb      A BasisKey, containing the definition of the
     *                      basis in the second coordinate direction for
     *                      quadrilateral elements
     * @param   graph2D     A mesh, containing information about the domain
     *                      and the spectral/hp element expansion.
     * @param   TriNb       The PointsType of possible nodal points
     * @param   ImpType     Detail about the implementation type to use 
     *                      in operators
     */
    ExpList::ExpList(const LibUtilities::SessionReaderSharedPtr &pSession,
                     const LibUtilities::BasisKey &TriBa,
                     const LibUtilities::BasisKey &TriBb,
                     const LibUtilities::BasisKey &QuadBa,
                     const LibUtilities::BasisKey &QuadBb,
                     const SpatialDomains::MeshGraphSharedPtr &graph2D,
                     const LibUtilities::PointsType TriNb,
                     const Collections::ImplementationType ImpType):
        m_expType(e2D),
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
        
        int elmtid=0;
        LocalRegions::TriExpSharedPtr tri;
        LocalRegions::NodalTriExpSharedPtr Ntri;
        LocalRegions::QuadExpSharedPtr quad;
        SpatialDomains::Composite comp;
        
        const SpatialDomains::ExpansionMap &expansions = 
            graph2D->GetExpansions();
        m_ncoeffs = 0;
        m_npoints = 0;
        
        m_physState = false;
        
        for (auto &expIt : expansions)
        {
            SpatialDomains::TriGeomSharedPtr TriangleGeom;
            SpatialDomains::QuadGeomSharedPtr QuadrilateralGeom;
            
            if ((TriangleGeom = std::dynamic_pointer_cast<SpatialDomains::
                 TriGeom>(expIt.second->m_geomShPtr)))
            {
                if (TriNb < LibUtilities::SIZE_PointsType)
                {
                    Ntri = MemoryManager<LocalRegions::NodalTriExp>::
                        AllocateSharedPtr(TriBa, TriBb, TriNb, 
                                          TriangleGeom);
                        Ntri->SetElmtId(elmtid++);
                        (*m_exp).push_back(Ntri);
                }
                else
                {
                    tri = MemoryManager<LocalRegions::TriExp>::
                        AllocateSharedPtr(TriBa, TriBb, TriangleGeom);
                    tri->SetElmtId(elmtid++);
                    (*m_exp).push_back(tri);
                }
                
                m_ncoeffs += (TriBa.GetNumModes()*(TriBa.GetNumModes()+1))/2
                    + TriBa.GetNumModes() * (TriBb.GetNumModes() - 
                                             TriBa.GetNumModes());
                m_npoints += TriBa.GetNumPoints()*TriBb.GetNumPoints();
            }
            else if ((QuadrilateralGeom = std::dynamic_pointer_cast<
                      SpatialDomains::QuadGeom>(expIt.second->m_geomShPtr)))
            {
                quad = MemoryManager<LocalRegions::QuadExp>::
                    AllocateSharedPtr(QuadBa, QuadBb, QuadrilateralGeom);
                quad->SetElmtId(elmtid++);
                (*m_exp).push_back(quad);
                
                m_ncoeffs += QuadBa.GetNumModes()*QuadBb.GetNumModes();
                m_npoints += QuadBa.GetNumPoints()*QuadBb.GetNumPoints();
            }
            else
            {
                ASSERTL0(false,
                         "dynamic cast to a proper Geometry2D failed");
            }
            
        }
        
        SetupCoeffPhys(true);
        
        CreateCollections(ImpType);
    }
}
}
