///////////////////////////////////////////////////////////////////////////////
//
// File ExpList1D.cpp
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
// Description: Expansion list 1D definition
//
///////////////////////////////////////////////////////////////////////////////

#include <iomanip>

#include <boost/core/ignore_unused.hpp>

#include <MultiRegions/ExpList1D.h>
#include <LibUtilities/Polylib/Polylib.h>
#include <LocalRegions/SegExp.h>
#include <LocalRegions/QuadExp.h>
#include <LocalRegions/Expansion2D.h>
#include <LibUtilities/Foundations/ManagerAccess.h>  // for PointsManager, etc
#include <LibUtilities/Foundations/Interp.h>

using namespace std;

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class ExpList1D
         * This multi-elemental expansion, which does not exhibit any coupling
         * between the expansion on the separate elements, can be formulated
         * as, \f[u^{\delta}(x_i)=\sum_{e=1}^{{N_{\mathrm{el}}}}
         * \sum_{n=0}^{N^{e}_m-1}\hat{u}_n^e\phi_n^e(x_i).\f]
         * where \f${N_{\mathrm{el}}}\f$ is the number of elements and
         * \f$N^{e}_m\f$ is the local elemental number of expansion modes.
         *
         * Instances of this class may be optionally constructed which use
         * generalised segment expansions (LocalRegions#GenSegExp), rather than
         * the standard segment expansions (LocalRegions#SegExp).
         * LocalRegions#GenSegExp provides additional spatial
         * data including segment normals and is enabled using the \a
         * UseGenSegExp flag.
         *
         * This class inherits all its variables and member functions from the
         * base class MultiRegions#ExpList.
         */

        /**
         * Assumes use of standard segment expansions only. All data storage
         * areas are initialised to empty arrays by the default ExpList
         * constructor.
         */
        ExpList1D::ExpList1D():
            ExpList()
        {
            SetExpType(e1D);
        }


        /**
         * Creates an identical copy of another ExpList1D object.
         */
        ExpList1D::ExpList1D(const ExpList1D &In, const bool DeclareCoeffPhysArrays):
            ExpList(In,DeclareCoeffPhysArrays)
        {
            SetExpType(e1D);
        }

        /**
         *
         */
        ExpList1D::ExpList1D(const ExpList1D &In,
                             const std::vector<unsigned int> &eIDs,
                             const bool DeclareCoeffPhysArrays,
                             const Collections::ImplementationType ImpType):
            ExpList(In, eIDs, DeclareCoeffPhysArrays)
        {
            SetExpType(e1D);

            // Allocate storage for data and populate element offset lists.
            SetCoeffPhysOffsets();

            CreateCollections(ImpType);
        }


        /**
         * After initialising the data inherited through MultiRegions#ExpList,
         * populate the expansion list from the segments defined in the supplied
         * SpatialDomains#MeshGraph1D. All expansions in the graph are defined
         * using the same LibUtilities#BasisKey which overrides that specified
         * in \a graph1D.
         *
         * @see     ExpList1D#ExpList1D(SpatialDomains::MeshGraph1D&, bool)
         *          for details.
         * @deprecated          This constructor is no longer required since
         *                      the basis key is now available from the graph.
         * @param   Ba          BasisKey describing quadrature points and
         *                      number of modes.
         * @param   graph1D     Domain and expansion definitions.
         */
        ExpList1D::ExpList1D(const LibUtilities::SessionReaderSharedPtr &pSession,
                             const LibUtilities::BasisKey &Ba,
                             const SpatialDomains::MeshGraphSharedPtr &graph1D,
                             const Collections::ImplementationType ImpType):
            ExpList(pSession,graph1D)
        {
            SetExpType(e1D);

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
                         SpatialDomains::SegGeom>(
                             expIt.second->m_geomShPtr)))
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

            m_coeffs = Array<OneD, NekDouble>{size_t(m_ncoeffs), 0.0};
            m_phys   = Array<OneD, NekDouble>{size_t(m_npoints), 0.0};

            CreateCollections(ImpType);
        }


        /**
         * Given a mesh \a graph1D, containing information about the domain and
         * the spectral/hp element expansion, this constructor fills the list
         * of local expansions \texttt{m_exp} with the proper expansions,
         * calculates the total number of quadrature points \f$x_i\f$ and local
         * expansion coefficients \f$\hat{u}^e_n\f$ and allocates memory for
         * the arrays #m_coeffs and #m_phys.
         *
         * For each element its corresponding LibUtilities#BasisKey is
         * retrieved and this is used to construct either a standard segment
         * (LocalRegions#SegExp) or a generalised segment
         * (LocalRegions#GenSegExp) which is stored in the list #m_exp.
         * Finally, ExpList#SetCoeffPhys is called to initialise the data
         * storage areas and set up the offset arrays.
         *
         * @param   graph1D     A mesh, containing information about the
         *                      domain and the spectral/hp element expansion.
         * @param   UseGenSegExp If true, create general segment expansions
         *                      instead of just normal segment expansions.
         */
        ExpList1D::ExpList1D(const LibUtilities::SessionReaderSharedPtr &pSession,
                const SpatialDomains::MeshGraphSharedPtr &graph1D,
                             const bool DeclareCoeffPhysArrays,
                             const Collections::ImplementationType ImpType):
            ExpList(pSession,graph1D)
        {
            SetExpType(e1D);

            int id=0;
            LocalRegions::SegExpSharedPtr seg;
            SpatialDomains::SegGeomSharedPtr SegmentGeom;

            // Retrieve the list of expansions
            const SpatialDomains::ExpansionMap &expansions
                                                    = graph1D->GetExpansions();

            // Process each expansion in the graph
            for (auto &expIt : expansions)
            {
                // Retrieve basis key from expansion
                LibUtilities::BasisKey bkey = expIt.second->m_basisKeyVector[0];

                if ((SegmentGeom = std::dynamic_pointer_cast<
                         SpatialDomains::SegGeom>(
                             expIt.second->m_geomShPtr)))
                {
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
                m_coeffs = Array<OneD, NekDouble> {size_t(m_ncoeffs), 0.0};
                m_phys   = Array<OneD, NekDouble> {size_t(m_npoints), 0.0};
            }

            CreateCollections(ImpType);
        }



        /**
         * Given a mesh \a graph1D, and the spectral/hp element
         * expansion as well as and separate information about a \a
         * domain, this constructor fills the list of local
         * expansions \texttt{m_exp} with the proper expansions,
         * calculates the total number of quadrature points \f$x_i\f$
         * and local expansion coefficients \f$\hat{u}^e_n\f$ and
         * allocates memory for the arrays #m_coeffs and #m_phys.
         *
         * For each element its corresponding LibUtilities#BasisKey is
         * retrieved and this is used to construct either a standard segment
         * (LocalRegions#SegExp) or a generalised segment
         * (LocalRegions#GenSegExp) which is stored in the list #m_exp.
         * Finally, ExpList#SetCoeffPhys is called to initialise the data
         * storage areas and set up the offset arrays.
         *
         * @param   graph1D     A mesh, containing information about the
         *                      domain and the spectral/hp element expansion.
         * @param   UseGenSegExp If true, create general segment expansions
         *                      instead of just normal segment expansions.
         */
        ExpList1D::ExpList1D(const LibUtilities::SessionReaderSharedPtr &pSession,
                             const SpatialDomains::MeshGraphSharedPtr &graph1D,
                             const SpatialDomains::CompositeMap &domain,
                             const bool DeclareCoeffPhysArrays,
                             const std::string var,
                             bool SetToOneSpaceDimension,
                             const Collections::ImplementationType ImpType):
            ExpList(pSession,graph1D)
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
                            SpatialDomains::SegGeom>(
                                compIt.second->m_geomVec[j])))
                    {
                        // Retrieve basis key from expansion and define expansion
                        auto expIt = expansions.find(SegmentGeom->GetGlobalID());
                        if(expIt != expansions.end())
                        {
                            LibUtilities::BasisKey bkey = expIt->second->m_basisKeyVector[0];

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
                m_coeffs = Array<OneD, NekDouble> {size_t(m_ncoeffs), 0.0};
                m_phys   = Array<OneD, NekDouble> {size_t(m_npoints), 0.0};
            }

            CreateCollections(ImpType);
        }


        /**
         * Fills the list of local expansions with the segments from the 2D
         * mesh specified by \a domain. This CompositeMap contains a list of
         * Composites which define the Neumann boundary.
         * @see     ExpList1D#ExpList1D(SpatialDomains::MeshGraph1D&, bool)
         *          for details.
         * @param   domain      A domain, comprising of one or more composite
         *                      regions,
         * @param   graph2D     A mesh, containing information about the
         *                      domain and the spectral/hp element expansion.
         * @param   UseGenSegExp If true, create general segment expansions
         *                      instead of just normal segment expansions.
         */
        ExpList1D::ExpList1D(const LibUtilities::SessionReaderSharedPtr &pSession,
                             const SpatialDomains::CompositeMap &domain,
                             const SpatialDomains::MeshGraphSharedPtr &graph2D,
                             const bool DeclareCoeffPhysArrays,
                             const std::string variable,
                             const LibUtilities::CommSharedPtr comm,
                             const Collections::ImplementationType ImpType)
            : ExpList(pSession,graph2D)
        {
            SetExpType(e1D);

            m_graph = graph2D;

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
                            SpatialDomains::SegGeom>(
                                compIt.second->m_geomVec[j])))
                    {
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
                m_coeffs = Array<OneD, NekDouble> {size_t(m_ncoeffs), 0.0};
                m_phys   = Array<OneD, NekDouble> {size_t(m_npoints), 0.0};
            }

            CreateCollections(ImpType);
        }

        /**
         * Store expansions for the trace space expansions used in
         * DisContField2D.
         *
         * @param   bndConstraint   Array of ExpList1D objects each containing a
         *                      1D spectral/hp element expansion on a single
         *                      boundary region.
         * @param   bndCond     Array of BoundaryCondition objects which contain
         *                      information about the boundary conditions on the
         *                      different boundary regions.
         * @param   locexp      Complete domain expansion list.
         * @param   graph2D     2D mesh corresponding to the expansion list.
         * @param   periodicEdges   List of periodic edges.
         * @param   UseGenSegExp If true, create general segment expansions
         *                      instead of just normal segment expansions.
         */
        ExpList1D::ExpList1D(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const Array<OneD,const ExpListSharedPtr>  &bndConstraint,
            const Array<OneD, const SpatialDomains::BoundaryConditionShPtr>  &bndCond,
            const LocalRegions::ExpansionVector &locexp,
            const SpatialDomains::MeshGraphSharedPtr &graph2D,
            const PeriodicMap &periodicEdges,
            const bool DeclareCoeffPhysArrays,
            const std::string variable,
            const Collections::ImplementationType ImpType):
            ExpList(pSession,graph2D)
        {
            boost::ignore_unused(periodicEdges, variable);

            int i, j, id, elmtid = 0;
            set<int> edgesDone;

            SpatialDomains::Geometry1DSharedPtr segGeom;
            SpatialDomains::Geometry2DSharedPtr ElGeom;
            LocalRegions::SegExpSharedPtr       seg;
            LocalRegions::SegExpSharedPtr       seg_tmp;
            LocalRegions::Expansion1DSharedPtr  exp1D;
            LocalRegions::Expansion2DSharedPtr  exp2D;

            SetExpType(e1D);

            map<int,int> EdgeDone;
            map<int,int> NormalSet;

            LocalRegions::SegExpSharedPtr Seg;

            // First loop over boundary conditions to renumber
            // Dirichlet boundaries
            for(i = 0; i < bndCond.size(); ++i)
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
                m_coeffs = Array<OneD, NekDouble> {size_t(m_ncoeffs), 0.0};
                m_phys   = Array<OneD, NekDouble> {size_t(m_npoints), 0.0};
            }

            CreateCollections(ImpType);
        }

        /**
         *
         */
        ExpList1D::~ExpList1D()
        {
        }

        /**
         * To perform post-processing on the entire domain use \a elmtId = 0.
         * @param   kernel      The post-processing kernel.
         * @param   inarray     The set of evaluation points.
         * @param   outarray    Contains the resulting post-processed
         *                      solution for element \a elmId.
         * @param   h           The mesh spacing.
         * @param   elmId       Optionally specifies which element to perform
         *                      the post-processing on (0=whole domain).
         */
        void ExpList1D::PostProcess(LibUtilities::KernelSharedPtr kernel,
                                    Array<OneD,NekDouble> &inarray,
                                    Array<OneD,NekDouble> &outarray,
                                    NekDouble h,
                                    int elmId)

        {
            int i,j,r;

            // get the local element expansion of the elmId element
            LocalRegions::ExpansionSharedPtr elmExp = GetExp( elmId );

            // Get the quadrature points and weights required for integration
            int quad_npoints = elmExp->GetTotPoints();
            LibUtilities::PointsKey quadPointsKey(quad_npoints,
                                                    elmExp->GetPointsType(0));
            Array<OneD,NekDouble> quad_points
                        = LibUtilities::PointsManager()[quadPointsKey]->GetZ();
            Array<OneD,NekDouble> quad_weights
                        = LibUtilities::PointsManager()[quadPointsKey]->GetW();

            // Declare variable for the local kernel breaks
            int kernel_width = kernel->GetKernelWidth();
            Array<OneD,NekDouble> local_kernel_breaks(kernel_width+1);

            // Declare variable for the transformed quadrature points
            Array<OneD,NekDouble> mapped_quad_points(quad_npoints);

            // For each evaluation point
            for(i = 0; i < inarray.size(); i++)
            {
                // Move the center of the kernel to the current point
                kernel->MoveKernelCenter(inarray[i],local_kernel_breaks);

                // Find the mesh breaks under the kernel support
                Array<OneD,NekDouble> mesh_breaks;
                kernel->FindMeshUnderKernel(local_kernel_breaks,h,mesh_breaks);

                // Sort the total breaks for integration purposes
                int total_nbreaks = local_kernel_breaks.size() +
                                    mesh_breaks.size();
                                    // number of the total breaks
                Array<OneD,NekDouble> total_breaks(total_nbreaks);
                kernel->Sort(local_kernel_breaks,mesh_breaks,total_breaks);

                // Integrate the product of kernel and function over the total
                // breaks
                NekDouble integral_value = 0.0;
                for(j = 0; j < total_breaks.size()-1; j++)
                {
                    NekDouble a = total_breaks[j];
                    NekDouble b = total_breaks[j+1];

                    // Map the quadrature points to the appropriate interval
                    for(r = 0; r < quad_points.size(); r++)
                    {
                        mapped_quad_points[r]
                                = (quad_points[r] + 1.0) * 0.5 * (b - a) + a;
                    }

                    // Evaluate the function at the transformed quadrature
                    // points
                    Array<OneD,NekDouble> u_value(quad_npoints);
                    Array<OneD,NekDouble> coeffs = GetCoeffs();

                    PeriodicEval(coeffs,mapped_quad_points,h,
                                 elmExp->GetBasisNumModes(0),u_value);

                    // Evaluate the kernel at the transformed quadrature points
                    Array<OneD,NekDouble> k_value(quad_npoints);
                    kernel->EvaluateKernel(mapped_quad_points,h,k_value);

                    // Integrate
                    for(r = 0; r < quad_npoints; r++)
                    {
                        integral_value += (b - a) * 0.5 * k_value[r]
                                                * u_value[r] * quad_weights[r];
                    }
                }
                outarray[i] = integral_value/h;
            }
        }


        /**
         * Given the elemental coefficients \f$\hat{u}_n^e\f$ of an expansion,
         * periodically evaluate the spectral/hp expansion
         * \f$u^{\delta}(\boldsymbol{x})\f$ at arbitrary points.
         * @param   inarray1    An array of size \f$N_{\mathrm{eof}}\f$
         *                      containing the local coefficients
         *                      \f$\hat{u}_n^e\f$.
         * @param   inarray2    Contains the set of evaluation points.
         * @param   h           The mesh spacing.
         * @param   nmodes      The number of polynomial modes for each element
         *                      (we consider that each element has the same
         *                      number of polynomial modes).
         * @param   outarray    Contains the resulting values at the
         *                      evaluation points
         */
        void ExpList1D::PeriodicEval(Array<OneD,NekDouble> &inarray1,
                                     Array<OneD,NekDouble> &inarray2,
                                     NekDouble h, int nmodes,
                                     Array<OneD,NekDouble> &outarray)
        {
            int i,j,r;

            // Get the number of elements in the domain
            int num_elm = GetExpSize();

            // initializing the outarray
            for(i = 0; i < outarray.size(); i++)
            {
                outarray[i] = 0.0;
            }

            // Make a copy for further modification
            int x_size = inarray2.size();
            Array<OneD,NekDouble> x_values_cp(x_size);

            // Determining the element to which the x belongs
            Array<OneD,int> x_elm(x_size);
            for(i = 0; i < x_size; i++ )
            {
                x_elm[i] = (int)floor(inarray2[i]/h);
            }

            // Clamp indices periodically
            for(i = 0; i < x_size; i++)
            {
                while(x_elm[i] < 0)
                {
                    x_elm[i] += num_elm;
                }
                while(x_elm[i] >= num_elm)
                {
                    x_elm[i] -= num_elm ;
                }
            }

            // Map the values of x to [-1 1] on its interval
            for(i = 0; i < x_size; i++)
            {
                x_values_cp[i] = (inarray2[i]/h - floor(inarray2[i]/h))*2 - 1.0;
            }

            // Evaluate the jocobi polynomials
            // (Evaluating the base at some points other than the quadrature
            // points). Should it be added to the base class????
            Array<TwoD,NekDouble> jacobi_poly(nmodes,x_size);
            for(i = 0; i < nmodes; i++)
            {
                Polylib::jacobfd(x_size,x_values_cp.get(),
                                    jacobi_poly.get()+i*x_size,NULL,i,0.0,0.0);
            }

            // Evaluate the function values
            for(r = 0; r < nmodes; r++)
            {
                for(j = 0; j < x_size; j++)
                {
                    int index = ((x_elm[j])*nmodes)+r;
                    outarray[j] += inarray1[index]*jacobi_poly[r][j];
                }
            }

        }


        /**
         * Sets up the normals on all edges of expansions in the domain.
         * @param   locexp      Complete list of domain expansions.
         */
//        void ExpList1D::SetUpPhysNormals(
//                                const StdRegions::StdExpansionVector &locexp)
//        {
//            map<int, int> EdgeGID;
//            int i,cnt,n,id;
//
//            // setup map of all global ids along boundary
//            for(cnt = i = 0; i < (*m_exp).size(); ++i)
//            {
//                id =  (*m_exp)[i]->GetGeom1D()->GetEid();
//                EdgeGID[id] = cnt++;
//            }
//
//            // Loop over elements and find edges that match;
//            for(cnt = n = 0; n < locexp.size(); ++n)
//            {
//                for(i = 0; i < locexp[n]->GetNedges(); ++i)
//                {
//                    id = locexp[n]->GetGeom2D()->GetEid(i);
//
//                    if(EdgeGID.count(id) > 0)
//                    {
//                        (*m_exp)[EdgeGID.find(id)->second]
//                                            ->SetUpPhysNormals(locexp[n],i);
//                    }
//                }
//            }
//        }

		//croth
		void ExpList1D::v_SetUpPhysNormals()
        {
            int i, j;
            for (i = 0; i < m_exp->size(); ++i)
            {
                for (j = 0; j < (*m_exp)[i]->GetNverts(); ++j)
                {
                    (*m_exp)[i]->ComputeVertexNormal(j);
                }
            }
        }

        /**
         * Upwind the left and right states given by the Arrays Fwd and Bwd
         * using the vector quantity Vec and ouput the upwinded value in the
         * array upwind.
         *
         * @param   Vec         Velocity field.
         * @param   Fwd         Left state.
         * @param   Bwd         Right state.
         * @param   Upwind      Output vector.
         */
        void ExpList1D::v_Upwind(
            const Array<OneD, const Array<OneD,       NekDouble> > &Vec,
            const Array<OneD,                   const NekDouble>   &Fwd,
            const Array<OneD,                   const NekDouble>   &Bwd,
                  Array<OneD,                         NekDouble>   &Upwind)
        {
            int i,j,k,e_npoints,offset;
            Array<OneD,NekDouble> normals;
            NekDouble Vn;

            // Assume whole array is of same coordimate dimension
            int coordim = GetCoordim(0);

            ASSERTL1(Vec.size() >= coordim,
                    "Input vector does not have sufficient dimensions to "
                    "match coordim");

            // Process each expansion
            for(i = 0; i < m_exp->size(); ++i)
            {
                // Get the number of points in the expansion and the normals.
                e_npoints = (*m_exp)[i]->GetNumPoints(0);
                normals   = (*m_exp)[i]->GetPhysNormals();

                // Get the physical data offset of the expansion in m_phys.
                offset = m_phys_offset[i];

                // Compute each data point.
                for(j = 0; j < e_npoints; ++j)
                {
                    // Calculate normal velocity.
                    Vn = 0.0;
                    for(k = 0; k < coordim; ++k)
                    {
                        Vn += Vec[k][offset+j]*normals[k*e_npoints + j];
                    }

                    // Upwind based on direction of normal velocity.
                    if(Vn > 0.0)
                    {
                        Upwind[offset + j] = Fwd[offset + j];
                    }
                    else
                    {
                        Upwind[offset + j] = Bwd[offset + j];
                    }
                }
            }
        }

        /**
         * One-dimensional upwind.
         * \see    ExpList1D::Upwind(
         *           const Array<OneD, const Array<OneD, NekDouble> >,
         *           const Array<OneD, const NekDouble>,
         *           const Array<OneD, const NekDouble>,
         *                 Array<OneD, NekDouble>, int)
         *
         * @param   Vn          Velocity field.
         * @param   Fwd         Left state.
         * @param   Bwd         Right state.
         * @param   Upwind      Output vector.
         */
        void ExpList1D::v_Upwind(
            const Array<OneD, const NekDouble> &Vn,
            const Array<OneD, const NekDouble> &Fwd,
            const Array<OneD, const NekDouble> &Bwd,
                  Array<OneD,       NekDouble> &Upwind)
        {
            int i,j,e_npoints,offset;
            Array<OneD,NekDouble> normals;

            // Process each expansion.
            for(i = 0; i < m_exp->size(); ++i)
            {
                // Get the number of points and the data offset.
                e_npoints = (*m_exp)[i]->GetNumPoints(0);
                offset = m_phys_offset[i];

                // Process each point in the expansion.
                for(j = 0; j < e_npoints; ++j)
                {
                    // Upwind based on one-dimensional velocity.
                    if(Vn[offset + j] > 0.0)
                    {
                        Upwind[offset + j] = Fwd[offset + j];
                    }
                    else
                    {
                        Upwind[offset + j] = Bwd[offset + j];
                    }
                }
            }
        }


        /**
         * For each local element, copy the normals stored in the element list
         * into the array \a normals.
         * @param   normals     Multidimensional array in which to copy normals
         *                      to. Must have dimension equal to or larger than
         *                      the spatial dimension of the elements.
         */
        void ExpList1D::v_GetNormals(
            Array<OneD, Array<OneD, NekDouble> > &normals)
        {
            int i,j,k,e_npoints,offset;
            SpatialDomains::Geometry1DSharedPtr segGeom;
            Array<OneD,Array<OneD,NekDouble> > locnormals;
            Array<OneD,Array<OneD,NekDouble> > locnormals2;
            Array<OneD,Array<OneD,NekDouble> > Norms;
            // Assume whole array is of same coordinate dimension
            int coordim = GetCoordim(0);

            ASSERTL1(normals.size() >= coordim,
                     "Output vector does not have sufficient dimensions to "
                     "match coordim");

            for (i = 0; i < m_exp->size(); ++i)
            {
                LocalRegions::Expansion1DSharedPtr loc_exp = (*m_exp)[i]->as<LocalRegions::Expansion1D>();

                LocalRegions::Expansion2DSharedPtr loc_elmt =
                    loc_exp->GetLeftAdjacentElementExp();

                int edgeNumber = loc_exp->GetLeftAdjacentElementEdge();

                // Get the number of points and normals for this expansion.
                e_npoints  = (*m_exp)[i]->GetNumPoints(0);

                locnormals = loc_elmt->GetEdgeNormal(edgeNumber);
                int e_nmodes   = loc_exp->GetBasis(0)->GetNumModes();
                int loc_nmodes = loc_elmt->GetBasis(0)->GetNumModes();

                if (e_nmodes != loc_nmodes)
                {
                    if (loc_exp->GetRightAdjacentElementEdge() >= 0)
                    {
                        LocalRegions::Expansion2DSharedPtr loc_elmt =
                                       loc_exp->GetRightAdjacentElementExp();

                        int EdgeNumber = loc_exp->GetRightAdjacentElementEdge();
                        // Serial case: right element is connected so we can
                        // just grab that normal.
                        locnormals = loc_elmt->GetEdgeNormal(EdgeNumber);

                        offset = m_phys_offset[i];

                        // Process each point in the expansion.
                        for (j = 0; j < e_npoints; ++j)
                        {
                            // Process each spatial dimension and copy the values
                            // into the output array.
                            for (k = 0; k < coordim; ++k)
                            {
                                normals[k][offset + j] = -locnormals[k][j];
                            }
                        }
                    }
                    else
                    {
                        // Parallel case: need to interpolate normal.
                        Array<OneD, Array<OneD, NekDouble> > normal(coordim);

                        for (int p = 0; p < coordim; ++p)
                        {
                            normal[p] = Array<OneD, NekDouble>(e_npoints,0.0);
                            LibUtilities::PointsKey to_key =
                                loc_exp->GetBasis(0)->GetPointsKey();
                            LibUtilities::PointsKey from_key =
                                loc_elmt->GetBasis(0)->GetPointsKey();
                            LibUtilities::Interp1D(from_key,
                                                   locnormals[p],
                                                   to_key,
                                                   normal[p]);
                        }

                        offset = m_phys_offset[i];

                        // Process each point in the expansion.
                        for (j = 0; j < e_npoints; ++j)
                        {
                            // Process each spatial dimension and copy the values
                            // into the output array.
                            for (k = 0; k < coordim; ++k)
                            {
                                normals[k][offset + j] = normal[k][j];
                            }
                        }
                    }
                }
                else
                {
                    // Get the physical data offset for this expansion.
                    offset = m_phys_offset[i];

                    // Process each point in the expansion.
                    for (j = 0; j < e_npoints; ++j)
                    {
                        // Process each spatial dimension and copy the values
                        // into the output array.
                        for (k = 0; k < coordim; ++k)
                        {
                            normals[k][offset + j] = locnormals[k][j];
                        }
                    }
                }
            }
        }

        void ExpList1D::v_GetElmtNormalLength(
            Array<OneD, NekDouble>  &lengthsFwd,
            Array<OneD, NekDouble>  &lengthsBwd)
        {
            int e_npoints;

            Array<OneD, NekDouble> locLeng;
            Array<OneD, NekDouble> lengintp;
            Array<OneD, NekDouble> lengAdd;
            Array<OneD, int      > LRbndnumbs(2);
            Array<OneD, Array<OneD,NekDouble> > lengLR(2);
            lengLR[0]   =   lengthsFwd;
            lengLR[1]   =   lengthsBwd;
            Array<OneD, LocalRegions::Expansion2DSharedPtr> LRelmts(2);
            LocalRegions::Expansion2DSharedPtr loc_elmt;
            LocalRegions::Expansion1DSharedPtr loc_exp;
            int e_npoints0  =   -1; 
            for (int i = 0; i < m_exp->size(); ++i)
            {
                loc_exp = (*m_exp)[i]->as<LocalRegions::Expansion1D>();
                int offset = m_phys_offset[i];
                
                int e_nmodes   = loc_exp->GetBasis(0)->GetNumModes();
                e_npoints  = (*m_exp)[i]->GetNumPoints(0);
                if ( e_npoints0 < e_npoints)
                {
                    lengintp = Array<OneD, NekDouble>{size_t(e_npoints), 0.0};
                    e_npoints0 = e_npoints;
                }
                
                LRelmts[0] = loc_exp->GetLeftAdjacentElementExp();
                LRelmts[1] = loc_exp->GetRightAdjacentElementExp();

                LRbndnumbs[0]   =   loc_exp->GetLeftAdjacentElementEdge();
                LRbndnumbs[1]   =   loc_exp->GetRightAdjacentElementEdge();
                for (int nlr = 0; nlr < 2; ++nlr)
                {
                    Vmath::Zero(e_npoints0, lengintp, 1);
                    lengAdd     =   lengintp;
                    int bndNumber = LRbndnumbs[nlr];
                    loc_elmt = LRelmts[nlr];
                    if (bndNumber >= 0)
                    {
                        locLeng  = loc_elmt->GetElmtBndNormDirElmtLen(
                                                bndNumber);
                        lengAdd  =   locLeng;

                        int loc_nmodes  = loc_elmt->GetBasis(0)->GetNumModes();

                        if (e_nmodes != loc_nmodes)
                        {
                            // Parallel case: need to interpolate.
                            LibUtilities::PointsKey to_key =
                                loc_exp->GetBasis(0)->GetPointsKey();
                            LibUtilities::PointsKey from_key =
                                loc_elmt->GetBasis(0)->GetPointsKey();
                            LibUtilities::Interp1D(from_key,
                                                locLeng,
                                                to_key,
                                                lengintp);
                            lengAdd     =   lengintp;
                        }
                    }
                    for (int j = 0; j < e_npoints; ++j)
                    {
                        lengLR[nlr][offset + j] = lengAdd[j];
                    }
                }
            }
        }

        /**
         * Writes out the header for a <PIECE> VTK XML segment describing the
         * geometric information which comprises this element. This includes
         * vertex coordinates for each quadrature point, vertex connectivity
         * information, cell types and cell offset data.
         *
         * @param   outfile     Output stream to write data to.
         */
        void ExpList1D::v_WriteVtkPieceHeader(std::ostream &outfile, int expansion, int)
        {
            int i,j;
            int nquad0 = (*m_exp)[expansion]->GetNumPoints(0);
            int ntot = nquad0;
            int ntotminus = (nquad0-1);

            Array<OneD,NekDouble> coords[3];
            coords[0] = Array<OneD,NekDouble>(ntot, 0.0);
            coords[1] = Array<OneD,NekDouble>(ntot, 0.0);
            coords[2] = Array<OneD,NekDouble>(ntot, 0.0);
            (*m_exp)[expansion]->GetCoords(coords[0],coords[1],coords[2]);

            outfile << "    <Piece NumberOfPoints=\""
                    << ntot << "\" NumberOfCells=\""
                    << ntotminus << "\">" << endl;
            outfile << "      <Points>" << endl;
            outfile << "        <DataArray type=\"Float64\" "
                    << "NumberOfComponents=\"3\" format=\"ascii\">" << endl;
            outfile << "          ";
            for (i = 0; i < ntot; ++i)
            {
                for (j = 0; j < 3; ++j)
                {
                    outfile << setprecision(8) << scientific
                            << (float)coords[j][i] << " ";
                }
                outfile << endl;
            }
            outfile << endl;
            outfile << "        </DataArray>" << endl;
            outfile << "      </Points>" << endl;
            outfile << "      <Cells>" << endl;
            outfile << "        <DataArray type=\"Int32\" "
                    << "Name=\"connectivity\" format=\"ascii\">" << endl;
            for (i = 0; i < nquad0-1; ++i)
            {
                outfile << i << " " << i+1 << endl;
            }
            outfile << endl;
            outfile << "        </DataArray>" << endl;
            outfile << "        <DataArray type=\"Int32\" "
                    << "Name=\"offsets\" format=\"ascii\">" << endl;
            for (i = 0; i < ntotminus; ++i)
            {
                outfile << i*2+2 << " ";
            }
            outfile << endl;
            outfile << "        </DataArray>" << endl;
            outfile << "        <DataArray type=\"UInt8\" "
                    << "Name=\"types\" format=\"ascii\">" << endl;
            for (i = 0; i < ntotminus; ++i)
            {
                outfile << "3 ";
            }
            outfile << endl;
            outfile << "        </DataArray>" << endl;
            outfile << "      </Cells>" << endl;
            outfile << "      <PointData>" << endl;
        }

    } //end of namespace
} //end of namespace
