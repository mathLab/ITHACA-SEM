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
// Description: Expansion list 1D definition
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ExpList1D.h>
#include <LibUtilities/Polylib/Polylib.h>

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
            ExpList(),
            m_UseGenSegExp(false)
        {
        }


        /**
         * Creates an identical copy of another ExpList1D object.
         */
        ExpList1D::ExpList1D(const ExpList1D &In):
            ExpList(In),
            m_UseGenSegExp(In.m_UseGenSegExp)
        {
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
         * @param   UseGenSegExp If true, create general segment expansions
         *                      instead of just normal segment expansions.
         */
        ExpList1D::ExpList1D(const LibUtilities::BasisKey &Ba,
                             const SpatialDomains::MeshGraph1D &graph1D,
                             bool UseGenSegExp):
            ExpList(),
            m_UseGenSegExp(UseGenSegExp)
        {
            int i,j, id=0;
            LocalRegions::SegExpSharedPtr seg;
            SpatialDomains::SegGeomSharedPtr SegmentGeom;

            const SpatialDomains::ExpansionVector &expansions
                                                    = graph1D.GetExpansions();

            // For each element in the mesh, create a segment expansion using
            // the supplied BasisKey and segment geometry.
            for(i = 0; i < expansions.size(); ++i)
            {
                if(SegmentGeom = boost
                            ::dynamic_pointer_cast<SpatialDomains::SegGeom>(
                                                expansions[i]->m_GeomShPtr))
                {
                    // Use a general segment expansion with normal and binormal?
                    if (UseGenSegExp)
                    {
                        seg = MemoryManager<LocalRegions::GenSegExp>
                                            ::AllocateSharedPtr(Ba,SegmentGeom);
                    }
                    else {
                        seg = MemoryManager<LocalRegions::SegExp>
                                            ::AllocateSharedPtr(Ba,SegmentGeom);
                    }
                    seg->SetElmtId(id++);
                    (*m_exp).push_back(seg);
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a SegGeom failed");
                }
            }

            // Allocate storage for data and populate element offset lists.
            ExpList::SetCoeffPhys();
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
        ExpList1D::ExpList1D(SpatialDomains::MeshGraph1D &graph1D,
                             bool UseGenSegExp):
            ExpList(),
            m_UseGenSegExp(UseGenSegExp)
        {
            int i,id=0;
            LocalRegions::SegExpSharedPtr seg;
            SpatialDomains::SegGeomSharedPtr SegmentGeom;

            // Retrieve the list of expansions
            const SpatialDomains::ExpansionVector &expansions
                                                    = graph1D.GetExpansions();

            // Process each expansion in the graph
            for(i = 0; i < expansions.size(); ++i)
            {
                // Retrieve basis key from expansion
                LibUtilities::BasisKey bkey
                                        = expansions[i]->m_BasisKeyVector[0];

                if(SegmentGeom = boost
                            ::dynamic_pointer_cast<SpatialDomains::SegGeom>(
                                                expansions[i]->m_GeomShPtr))
                {
                    // Switch depending on if using general segment expansions.
                    if (UseGenSegExp)
                    {
                        seg = MemoryManager<LocalRegions::GenSegExp>
                                        ::AllocateSharedPtr(bkey, SegmentGeom);
                    }
                    else {
                        seg = MemoryManager<LocalRegions::SegExp>
                                        ::AllocateSharedPtr(bkey, SegmentGeom);
                    }
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

            // Allocate storage for data and populate element offset lists.
            ExpList::SetCoeffPhys();
        }


        /**
         * Fills the list of local expansions with the segments from the 2D
         * mesh specified by \a domain. This CompositeVector contains a list of
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
        ExpList1D::ExpList1D(const SpatialDomains::CompositeVector &domain,
                             SpatialDomains::MeshGraph2D &graph2D,
                             bool UseGenSegExp):
            ExpList(),
            m_UseGenSegExp(UseGenSegExp)
        {
            int i,j,cnt,id=0;
            SpatialDomains::Composite comp;
            SpatialDomains::SegGeomSharedPtr SegmentGeom;
            LocalRegions::SegExpSharedPtr seg;

            // Process each composite region.
            cnt = 0;
            for(i = 0; i < domain.size(); ++i)
            {
                comp = domain[i];

                // Process each expansion in the region.
                for(j = 0; j < comp->size(); ++j)
                {
                    if(SegmentGeom = boost
                            ::dynamic_pointer_cast<SpatialDomains::SegGeom>(
                                                                    (*comp)[j]))
                    {
                        // Retrieve the basis key from the expansion.
                        LibUtilities::BasisKey bkey
                                        = graph2D.GetEdgeBasisKey(SegmentGeom);

                        // Use general expansions?
                        if (UseGenSegExp)
                        {
                            seg = MemoryManager<LocalRegions::GenSegExp>
                                        ::AllocateSharedPtr(bkey, SegmentGeom);
                        }
                        else
                        {
                            seg = MemoryManager<LocalRegions::SegExp>
                                        ::AllocateSharedPtr(bkey, SegmentGeom);
                        }

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
            ExpList::SetCoeffPhys();
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
                    const Array<OneD,const ExpList1DSharedPtr>  &bndConstraint,
                    const Array<OneD, const SpatialDomains
                                           ::BoundaryConditionShPtr>  &bndCond,
                    const StdRegions::StdExpansionVector &locexp,
                    SpatialDomains::MeshGraph2D &graph2D,
                    const map<int,int> &periodicEdges,
                    bool UseGenSegExp):
            ExpList(),
            m_UseGenSegExp(UseGenSegExp)
        {
            int i,j,k,cnt,id, elmtid=0;
            map<int,int> EdgeDone;
            map<int,int> NormalSet;

            SpatialDomains::Geometry1DSharedPtr SegGeom;
            LocalRegions::SegExpSharedPtr Seg;

            // First loop over boundary conditions to renumber
            // Dirichlet boundaries
            cnt = 0;
            for(i = 0; i < bndCond.num_elements(); ++i)
            {
                if(bndCond[i]->GetBoundaryConditionType()
                                            == SpatialDomains::eDirichlet)
                {
                    for(j = 0; j < bndConstraint[i]->GetExpSize(); ++j)
                    {
                        LibUtilities::BasisKey bkey = bndConstraint[i]
                                    ->GetExp(j)->GetBasis(0)->GetBasisKey();
                        SegGeom = bndConstraint[i]->GetExp(j)->GetGeom1D();

                        if (UseGenSegExp)
                        {
                            Seg = MemoryManager<LocalRegions::GenSegExp>
                                            ::AllocateSharedPtr(bkey, SegGeom);
                            EdgeDone[SegGeom->GetEid()] = elmtid;
                        }
                        else
                        {
                            Seg = MemoryManager<LocalRegions::SegExp>
                                            ::AllocateSharedPtr(bkey, SegGeom);
                            EdgeDone[SegGeom->GetEid()] = 1;
                        }

                        Seg->SetElmtId(elmtid++);
                        (*m_exp).push_back(Seg);
                    }
                }
            }

            // loop over all other edges and fill out other connectivities
            for(i = 0; i < locexp.size(); ++i)
            {
                for(j = 0; j < locexp[i]->GetNedges(); ++j)
                {
                    SegGeom = (locexp[i]->GetGeom2D())->GetEdge(j);

                    id = SegGeom->GetEid();

                    if(EdgeDone.count(id)==0)
                    {
                        LibUtilities::BasisKey EdgeBkey
                                    = locexp[i]->DetEdgeBasisKey(j);

                        if (UseGenSegExp)
                        {
                            Seg = MemoryManager<LocalRegions::GenSegExp>
                                        ::AllocateSharedPtr(EdgeBkey, SegGeom);
                            EdgeDone[id] = elmtid;

                            if (periodicEdges.count(id) > 0)
                            {
                                EdgeDone[periodicEdges.find(id)->second]
                                        = elmtid;
                            }
                        }
                        else
                        {
                            Seg = MemoryManager<LocalRegions::SegExp>
                                        ::AllocateSharedPtr(EdgeBkey, SegGeom);
                            EdgeDone[id] = 1;
                        }

                        Seg->SetElmtId(elmtid++);
                        (*m_exp).push_back(Seg);
                    }
                    else // variable modes/points
                    {
                        if (UseGenSegExp) {
                            LibUtilities::BasisKey EdgeBkey
                                    = locexp[i]->DetEdgeBasisKey(j);

                            if((*m_exp)[EdgeDone[id]]->GetNumPoints(0)
                                    >= EdgeBkey.GetNumPoints()
                                && (*m_exp)[EdgeDone[id]]->GetBasisNumModes(0)
                                    >= EdgeBkey.GetNumModes())
                            {
                            }
                            else if((*m_exp)[EdgeDone[id]]->GetNumPoints(0)
                                    <= EdgeBkey.GetNumPoints()
                                && (*m_exp)[EdgeDone[id]]->GetBasisNumModes(0)
                                    <= EdgeBkey.GetNumModes())
                            {
                                Seg = MemoryManager<LocalRegions::GenSegExp>
                                        ::AllocateSharedPtr(EdgeBkey, SegGeom);
                                Seg->SetElmtId(EdgeDone[id]);
                                (*m_exp)[EdgeDone[id]] = Seg;
                                NormalSet.erase(id);
                            }
                            else
                            {
                                ASSERTL0(false,
                                    "inappropriate number of points/modes (max "
                                    "num of points is not set with max order)")
                            }
                        }
                    }

                    if (UseGenSegExp && NormalSet.count(id) == 0)
                    {
                        Seg = boost::dynamic_pointer_cast
                                    <LocalRegions::GenSegExp>(
                                        (*m_exp)[EdgeDone.find(id)->second]);

                        // Set up normals at all Segment Quadrature points
                        Seg->SetUpPhysNormals(locexp[i],j);
                        NormalSet[id] = 1;
                    }
                }
            }

            // Set up offset information and array sizes
            ExpList::SetCoeffPhys();

        }


        /**
         *
         */
        ExpList1D::~ExpList1D()
        {
        }


        /**
         * @param   graph1D     A mesh containing information about the domain
         *                      and the Spectral/hp element expansion.
         * @param   bcs         Information about the boundary conditions.
         * @param   variable    Specifies the field.
         * @param   bndCondExpansions   Array of ExpList1D objects each
         *                      containing a 1D spectral/hp element expansion
         *                      on a single boundary region.
         * @param   bncConditions   Array of BoundaryCondition objects which
         *                      contain information about the boundary
         *                      conditions on the different boundary regions.
         */
        void ExpList1D::SetBoundaryConditionExpansion(
                                const SpatialDomains::MeshGraph1D &graph1D,
                                      SpatialDomains::BoundaryConditions &bcs,
                                const std::string variable,
                                Array<OneD, LocalRegions::PointExpSharedPtr>
                                                            &bndCondExpansions,
                                Array<OneD, SpatialDomains
                                    ::BoundaryConditionShPtr> &bndConditions)
        {
            int i,j,k;
            int cnt  = 0;

            SpatialDomains::BoundaryRegionCollection &bregions
                                                = bcs.GetBoundaryRegions();
            SpatialDomains::BoundaryConditionCollection &bconditions
                                                = bcs.GetBoundaryConditions();

            LocalRegions::PointExpSharedPtr          locPointExp;
            SpatialDomains::BoundaryConditionShPtr   locBCond;
            SpatialDomains::VertexComponentSharedPtr vert;

            int nbnd = bregions.size();

            cnt=0;
            // list Dirichlet boundaries first
            for(i = 0; i < nbnd; ++i)
            {
                locBCond = (*(bconditions[i]))[variable];
                if(locBCond->GetBoundaryConditionType()
                        == SpatialDomains::eDirichlet)
                {
                    for(j = 0; j < bregions[i]->size(); j++)
                    {
                        for(k = 0; k < ((*bregions[i])[j])->size(); k++)
                        {
                            if(vert = boost::dynamic_pointer_cast
                                    <SpatialDomains::VertexComponent>(
                                        (*(*bregions[i])[j])[k]))
                            {
                                locPointExp
                                    = MemoryManager<LocalRegions::PointExp>
                                                ::AllocateSharedPtr(vert);
                                bndCondExpansions[cnt]  = locPointExp;
                                bndConditions[cnt++]    = locBCond;
                            }
                            else
                            {
                                ASSERTL0(false,
                                         "dynamic cast to a vertex failed");
                            }
                        }
                    }
                }
            } // end if Dirichlet

            // then, list the other (non-periodic) boundaries
            for(i = 0; i < nbnd; ++i)
            {
                locBCond = (*(bconditions[i]))[variable];
                if(locBCond->GetBoundaryConditionType()
                        == SpatialDomains::eNeumann)
                {
                    for(j = 0; j < bregions[i]->size(); j++)
                    {
                        for(k = 0; k < ((*bregions[i])[j])->size(); k++)
                        {
                            if(vert = boost::dynamic_pointer_cast
                                    <SpatialDomains::VertexComponent>(
                                        (*(*bregions[i])[j])[k]))
                            {
                                locPointExp
                                    = MemoryManager<LocalRegions::PointExp>
                                                ::AllocateSharedPtr(vert);
                                bndCondExpansions[cnt]  = locPointExp;
                                bndConditions[cnt++]    = locBCond;
                            }
                            else
                            {
                                ASSERTL0(false,
                                         "dynamic cast to a vertex failed");
                            }
                        }
                    }
                }
                else if((locBCond->GetBoundaryConditionType()
                            != SpatialDomains::eDirichlet) &&
                        (locBCond->GetBoundaryConditionType()
                            != SpatialDomains::ePeriodic))
                {
                    ASSERTL0(false,"This type of BC not implemented yet");
                }
            }
        }


        /**
         * @param   graph1D     A mesh containing information about the domain
         *                      and the spectral/hp element expansion.
         * @param   bcs         Information about the boundary conditions.
         * @param   variable    Specified the field.
         * @param   periodicVertices    Map into which the list of periodic
         *                      vertices is placed.
         */
        void ExpList1D::GetPeriodicVertices(
                                const SpatialDomains::MeshGraph1D &graph1D,
                                      SpatialDomains::BoundaryConditions &bcs,
                                const std::string variable,
                                      map<int,int>& periodicVertices)
        {

            int i,j,k;

            SpatialDomains::BoundaryRegionCollection &bregions
                    = bcs.GetBoundaryRegions();
            SpatialDomains::BoundaryConditionCollection &bconditions
                    = bcs.GetBoundaryConditions();

            int region1ID;
            int region2ID;

            SpatialDomains::Composite comp1;
            SpatialDomains::Composite comp2;

            SpatialDomains::VertexComponentSharedPtr vert1;
            SpatialDomains::VertexComponentSharedPtr vert2;

            SpatialDomains::BoundaryConditionShPtr locBCond;

            // This std::map is a check so that the periodic pairs
            // are not treated twice
            map<int, int> doneBndRegions;

            int nbnd = bregions.size();

            for(i = 0; i < nbnd; ++i)
            {
                locBCond = (*(bconditions[i]))[variable];
                if(locBCond->GetBoundaryConditionType()
                        == SpatialDomains::ePeriodic)
                {
                    region1ID = i;
                    region2ID = (boost::static_pointer_cast<SpatialDomains::
                                    PeriodicBoundaryCondition>(locBCond))
                                        ->m_ConnectedBoundaryRegion;

                    if(doneBndRegions.count(region1ID)==0)
                    {
                        ASSERTL0(bregions[region1ID]->size()
                                    == bregions[region2ID]->size(),
                                 "Size of the 2 periodic boundary regions "
                                 "should be equal");

                        for(j = 0; j < bregions[region1ID]->size(); j++)
                        {
                            comp1 = (*(bregions[region1ID]))[j];
                            comp2 = (*(bregions[region2ID]))[j];

                            ASSERTL0(comp1->size() == comp2->size(),
                                     "Size of the 2 periodic composites should "
                                     "be equal");

                            for(k = 0; k < comp1->size(); k++)
                            {
                                if(!(vert1 = boost::dynamic_pointer_cast
                                        <SpatialDomains::VertexComponent>(
                                            (*comp1)[k]))||
                                   !(vert2 = boost::dynamic_pointer_cast
                                        <SpatialDomains::VertexComponent>(
                                            (*comp2)[k])))
                                {
                                    ASSERTL0(false,"dynamic cast to a "
                                                   "VertexComponent failed");
                                }

                                // Extract the periodic vertices
                                periodicVertices[vert1->GetVid()]
                                    = vert2->GetVid();
                            }
                        }
                    }
                    else
                    {
                        ASSERTL0(doneBndRegions[region1ID]==region2ID,
                                 "Boundary regions are not mutually periodic");
                    }
                    doneBndRegions[region2ID] = region1ID;
                }
            }
        }


        /**
         * Evaluates the boundary condition expansions, \a bndCondExpansions,
         * given the information provided by \a bndConditions.
         * @param   time        The time at which the boundary conditions
         *                      should be evaluated.
         * @param   bndCondExpansions   List of boundary expansions.
         * @param   bndConditions   Information about the boundary conditions.
         */
        void ExpList1D::EvaluateBoundaryConditions(
                                const NekDouble time,
                                Array<OneD, LocalRegions::PointExpSharedPtr>
                                                            &bndCondExpansions,
                                Array<OneD, SpatialDomains
                                    ::BoundaryConditionShPtr> &bndConditions)
        {
            int i;

            NekDouble x0;
            NekDouble x1;
            NekDouble x2;

            for(i = 0; i < bndCondExpansions.num_elements(); ++i)
            {
                bndCondExpansions[i]->GetCoords(x0,x1,x2);

                if(bndConditions[i]->GetBoundaryConditionType()
                        == SpatialDomains::eDirichlet)
                {
                    bndCondExpansions[i]->SetValue(
                            (boost::static_pointer_cast<SpatialDomains
                             ::DirichletBoundaryCondition>(bndConditions[i])
                             ->m_DirichletCondition).Evaluate(x0,x1,x2,time));
                }
                else if(bndConditions[i]->GetBoundaryConditionType()
                        == SpatialDomains::eNeumann)
                {
                    bndCondExpansions[i]->SetValue(
                            (boost::static_pointer_cast<SpatialDomains
                             ::NeumannBoundaryCondition>(bndConditions[i])
                             ->m_NeumannCondition).Evaluate(x0,x1,x2,time));
                }
                else
                {
                    ASSERTL0(false,"This type of BC not implemented yet");
                }
            }
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
            StdRegions::StdExpansionSharedPtr elmExp = GetExp(elmId);

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
            for(i = 0; i < inarray.num_elements(); i++)
            {
                // Move the center of the kernel to the current point
                kernel->MoveKernelCenter(inarray[i],local_kernel_breaks);

                // Find the mesh breaks under the kernel support
                Array<OneD,NekDouble> mesh_breaks;
                kernel->FindMeshUnderKernel(local_kernel_breaks,h,mesh_breaks);

                // Sort the total breaks for integration purposes
                int total_nbreaks = local_kernel_breaks.num_elements() +
                                    mesh_breaks.num_elements();
                                    // number of the total breaks
                Array<OneD,NekDouble> total_breaks(total_nbreaks);
                kernel->Sort(local_kernel_breaks,mesh_breaks,total_breaks);

                // Integrate the product of kernel and function over the total
                // breaks
                NekDouble integral_value = 0.0;
                for(j = 0; j < total_breaks.num_elements()-1; j++)
                {
                    double a = total_breaks[j];
                    double b = total_breaks[j+1];

                    // Map the quadrature points to the appropriate interval
                    for(r = 0; r < quad_points.num_elements(); r++)
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
            for(i = 0; i < outarray.num_elements(); i++)
            {
                outarray[i] = 0.0;
            }

            // Make a copy for further modification
            int x_size = inarray2.num_elements();
            Array<OneD,NekDouble> x_values_cp(x_size);

            // Determining the element to which the x belongs
            Array<OneD,int> x_elm(x_size);
            for(i = 0; i < x_size; i++ )
            {
                x_elm[i] = floor(inarray2[i]/h);
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
        void ExpList1D::SetUpPhysNormals(
                                const StdRegions::StdExpansionVector &locexp)
        {
            ASSERTL0(m_UseGenSegExp, "Must use GenSegExp to use normals.");

            map<int, int> EdgeGID;
            int i,j,cnt,n,id;

            // setup map of all global ids along boundary
            for(cnt = i = 0; i < (*m_exp).size(); ++i)
            {
                id =  (*m_exp)[i]->GetGeom1D()->GetEid();
                EdgeGID[id] = cnt++;
            }

            // Loop over elements and find edges that match;
            for(cnt = n = 0; n < locexp.size(); ++n)
            {
                for(i = 0; i < locexp[n]->GetNedges(); ++i)
                {
                    id = locexp[n]->GetGeom2D()->GetEid(i);

                    if(EdgeGID.count(id) > 0)
                    {
                        (*m_exp)[EdgeGID.find(id)->second]
                                            ->SetUpPhysNormals(locexp[n],i);
                    }
                }
            }
        }


        /**
         * Upwind the left and right states given by the Arrays Fwd
         * and Bwd using the vector quantity Vec and ouput the
         * upwinded value in the array upwind.
         * @param   Vec         Velocity field.
         * @param   Fwd         Left state.
         * @param   Bwd         Right state.
         * @param   Upwind      Output vector.
         * @param   direction   (Unused)
         */
        void ExpList1D::Upwind(
                        const Array<OneD, const Array<OneD, NekDouble> > &Vec,
                        const Array<OneD, const NekDouble> &Fwd,
                        const Array<OneD, const NekDouble> &Bwd,
                              Array<OneD, NekDouble> &Upwind,
                        int direction)
        {
            ASSERTL0(m_UseGenSegExp, "Must use GenSegExp to use normals.");

            int i,j,k,e_npoints,offset;
            Array<OneD,NekDouble> normals;
            NekDouble Vn;

            // Assume whole array is of same coordimate dimension
            int coordim = (*m_exp)[0]->GetGeom1D()->GetCoordim();

            ASSERTL1(Vec.num_elements() >= coordim,
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
         * @param   Vn          Velocity field.
         * @param   Fwd         Left state.
         * @param   Bwd         Right state.
         * @param   Upwind      Output vector.
         * @param   direction   (Unused).
         */
        void ExpList1D::Upwind(   const Array<OneD, const NekDouble> &Vn,
                                  const Array<OneD, const NekDouble> &Fwd,
                                  const Array<OneD, const NekDouble> &Bwd,
                                  Array<OneD, NekDouble> &Upwind,
                                  int direction)
        {
            ASSERTL0(m_UseGenSegExp, "Must use GenSegExp to use normals.");

            int i,j,k,e_npoints,offset;
            Array<OneD,NekDouble> normals;

            // Assume whole array is of same coordimate dimention
            int coordim = (*m_exp)[0]->GetGeom1D()->GetCoordim();

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
        void ExpList1D::GetNormals(
                                Array<OneD, Array<OneD, NekDouble> > &normals)
        {
            ASSERTL0(m_UseGenSegExp, "Must use GenSegExp to use normals.");

            int i,j,k,e_npoints,offset;
            Array<OneD,NekDouble> locnormals;

            // Assume whole array is of same coordinate dimension
            int coordim = (*m_exp)[0]->GetGeom1D()->GetCoordim();

            ASSERTL1(normals.num_elements() >= coordim,
                     "Output vector does not have sufficient dimensions to "
                     "match coordim");

            // Process each expansion.
            for(i = 0; i < m_exp->size(); ++i)
            {
                // Get the number of points and normals for this expansion.
                e_npoints  = (*m_exp)[i]->GetNumPoints(0);
                locnormals = (*m_exp)[i]->GetPhysNormals();

                // Get the physical data offset for this expansion.
                offset = m_phys_offset[i];

                // Process each point in the expansion.
                for(j = 0; j < e_npoints; ++j)
                {
                    // Process each spatial dimension and copy the values into
                    // the output array.
                    for(k = 0; k < coordim; ++k)
                    {
                        normals[k][offset+j] = locnormals[k*e_npoints + j];
                    }
                }
            }
        }


        /**
         *
         */
        void ExpList1D::v_SetUpPhysNormals(
                                const StdRegions::StdExpansionVector &locexp)
        {
            SetUpPhysNormals(locexp);
        }

    } //end of namespace
} //end of namespace

/**
* $Log: ExpList1D.cpp,v $
* Revision 1.40  2009/11/04 20:30:15  cantwell
* Added documentation to ExpList and ExpList1D and tidied up code.
*
* Revision 1.39  2009/11/04 12:33:38  cantwell
* Fix for HDGHelmholtz2D solver.
*
* Revision 1.38  2009/11/02 19:15:43  cantwell
* Moved ContField1D to inherit from DisContField1D.
* Moved ContField3D to inherit from DisContField3D.
* Incorporated GenExpList1D functionality into ExpList1D.
* Tidied up and added documentation to various classes.
* Moved Namespace documentation and introductions to separate files along with
* doxygen configuration.
* Added option to use system ZLIB library instead of libboost_zlib on UNIX.
* Added extra search paths to FindMetis.cmake and FindNektar++.cmake.
* Updated Linux compiling instructions.
* Updated regDemo to use Helmholtz2D-g when built as debug.
*
* Revision 1.37  2009/09/06 22:28:45  sherwin
* Updates for Navier-Stokes solver
*
* Revision 1.36  2009/04/20 16:14:06  sherwin
* Updates for optimising bandwidth of DG solver and allowing write import on explist
*
* Revision 1.35  2009/02/08 09:11:49  sherwin
* General updates to introduce multiple matrix definitions based on different boundary types
*
* Revision 1.34  2009/01/13 02:50:10  mirzaee
* Added definitions for the PostProcessing functions and PeriodicEval
*
* Revision 1.33  2008/09/09 15:06:03  sherwin
* Modifications related to curved elements.
*
* Revision 1.32  2008/08/14 22:15:51  sherwin
* Added LocalToglobalMap and DGMap and depracted LocalToGlobalBndryMap1D,2D. Made DisContField classes compatible with updated ContField formats
*
* Revision 1.31  2008/07/31 11:17:13  sherwin
* Changed GetEdgeBasis with DetEdgeBasisKey
*
* Revision 1.30  2008/07/29 22:27:33  sherwin
* Updates for DG solvers, including using GenSegExp, fixed forcing function on UDG HelmSolve and started to tidy up the mapping arrays to be 1D rather than 2D
*
* Revision 1.29  2008/07/12 17:31:39  sherwin
* Added m_phys_offset and rename m_exp_offset to m_coeff_offset
*
* Revision 1.28  2008/06/23 14:21:01  pvos
* updates for 1D ExpLists
*
* Revision 1.27  2008/05/14 18:06:50  sherwin
* mods to fix Seggeom to Geometry1D casting
*
* Revision 1.26  2008/05/13 22:06:58  sherwin
* Changed SegGeom to Geometry1D
*
* Revision 1.25  2008/05/10 18:27:33  sherwin
* Modifications necessary for QuadExp Unified DG Solver
*
* Revision 1.24  2008/03/18 14:14:13  pvos
* Update for nodal triangular helmholtz solver
*
* Revision 1.23  2008/03/12 15:25:45  pvos
* Clean up of the code
*
* Revision 1.22  2007/12/06 22:52:30  pvos
* 2D Helmholtz solver updates
*
* Revision 1.21  2007/11/07 20:29:53  jfrazier
* Modified to use new expansion list contained in meshgraph.
*
* Revision 1.20  2007/09/25 14:25:29  pvos
* Update for helmholtz1D with different expansion orders
*
* Revision 1.19  2007/07/22 23:04:20  bnelson
* Backed out Nektar::ptr.
*
* Revision 1.18  2007/07/20 02:04:12  bnelson
* Replaced boost::shared_ptr with Nektar::ptr
*
* Revision 1.17  2007/07/13 16:48:47  pvos
* Another HelmHoltz update (homogeneous dir BC multi-elemental solver does work)
*
* Revision 1.16  2007/07/10 08:54:29  pvos
* Updated ContField1D constructor
*
* Revision 1.15  2007/07/06 18:39:34  pvos
* ContField1D constructor updates
*
* Revision 1.14  2007/06/05 16:36:55  pvos
* Updated Explist2D ContExpList2D and corresponding demo-codes
*
**/
