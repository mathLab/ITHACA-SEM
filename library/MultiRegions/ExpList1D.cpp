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

#include <iomanip>
#include <MultiRegions/ExpList1D.h>
#include <LibUtilities/Polylib/Polylib.h>
#include <LocalRegions/SegExp.h>
#include <LocalRegions/Expansion2D.h>
#include <SpatialDomains/MeshGraph2D.h>
#include <LibUtilities/Foundations/ManagerAccess.h>  // for PointsManager, etc

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
        }


        /**
         * Creates an identical copy of another ExpList1D object.
         */
        ExpList1D::ExpList1D(const ExpList1D &In, const bool DeclareCoeffPhysArrays):
            ExpList(In,DeclareCoeffPhysArrays)
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
         */
        ExpList1D::ExpList1D(const LibUtilities::SessionReaderSharedPtr &pSession,
                             const LibUtilities::BasisKey &Ba,
                             const SpatialDomains::MeshGraphSharedPtr &graph1D):
            ExpList(pSession,graph1D)
        {
            int id=0;
            LocalRegions::SegExpSharedPtr seg;
            SpatialDomains::SegGeomSharedPtr SegmentGeom;

            const SpatialDomains::ExpansionMap &expansions
                                                    = graph1D->GetExpansions();

            // For each element in the mesh, create a segment expansion using
            // the supplied BasisKey and segment geometry.
            SpatialDomains::ExpansionMap::const_iterator expIt;
            for (expIt = expansions.begin(); expIt != expansions.end(); ++expIt)
            {
                if ((SegmentGeom = boost::dynamic_pointer_cast<
                         SpatialDomains::SegGeom>(
                             expIt->second->m_geomShPtr)))
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

            // Setup Default optimisation information.
            int nel = GetExpSize();
            m_globalOptParam = MemoryManager<NekOptimize::GlobalOptParam>
                ::AllocateSharedPtr(nel);

            // Allocate storage for data and populate element offset lists.
            SetCoeffPhysOffsets();

            m_coeffs = Array<OneD, NekDouble>(m_ncoeffs);
            m_phys   = Array<OneD, NekDouble>(m_npoints);

            ReadGlobalOptimizationParameters();
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
                const bool DeclareCoeffPhysArrays):
            ExpList(pSession,graph1D)
        {
            int id=0;
            LocalRegions::SegExpSharedPtr seg;
            SpatialDomains::SegGeomSharedPtr SegmentGeom;

            // Retrieve the list of expansions
            const SpatialDomains::ExpansionMap &expansions
                                                    = graph1D->GetExpansions();

            // Process each expansion in the graph
            SpatialDomains::ExpansionMap::const_iterator expIt;
            for (expIt = expansions.begin(); expIt != expansions.end(); ++expIt)
            {
                // Retrieve basis key from expansion
                LibUtilities::BasisKey bkey = expIt->second->m_basisKeyVector[0];

                if ((SegmentGeom = boost::dynamic_pointer_cast<
                         SpatialDomains::SegGeom>(
                             expIt->second->m_geomShPtr)))
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

            // Setup Default optimisation information.
            int nel = GetExpSize();
            m_globalOptParam = MemoryManager<NekOptimize::GlobalOptParam>
                ::AllocateSharedPtr(nel);

            // set up offset arrays.
            SetCoeffPhysOffsets();

            if(DeclareCoeffPhysArrays)
            {
                // Set up m_coeffs, m_phys.
                m_coeffs = Array<OneD, NekDouble>(m_ncoeffs);
                m_phys   = Array<OneD, NekDouble>(m_npoints);
            }

            ReadGlobalOptimizationParameters();
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
        ExpList1D::ExpList1D(const SpatialDomains::CompositeMap &domain,
                             const SpatialDomains::MeshGraphSharedPtr &graph2D,
                             const bool DeclareCoeffPhysArrays,
                             const std::string variable):
            ExpList()
        {
            int j, id=0;
            SpatialDomains::Composite comp;
            SpatialDomains::CompositeMap::const_iterator compIt;
            SpatialDomains::SegGeomSharedPtr SegmentGeom;
            LocalRegions::SegExpSharedPtr seg;

            // Process each composite region.
            for(compIt = domain.begin(); compIt != domain.end(); ++compIt)
            {
                comp = compIt->second;

                // Process each expansion in the region.
                for(j = 0; j < compIt->second->size(); ++j)
                {
                    if((SegmentGeom = boost::dynamic_pointer_cast<
                            SpatialDomains::SegGeom>(
                                (*compIt->second)[j])))
                    {
                        // Retrieve the basis key from the expansion.
                        LibUtilities::BasisKey bkey
                            = boost::dynamic_pointer_cast<SpatialDomains::MeshGraph2D>(graph2D)->GetEdgeBasisKey(SegmentGeom, variable);

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

            // Setup Default optimisation information.
            int nel = GetExpSize();
            m_globalOptParam = MemoryManager<NekOptimize::GlobalOptParam>
                ::AllocateSharedPtr(nel);

            // Allocate storage for data and populate element offset lists.
            SetCoeffPhysOffsets();

            // Set up m_coeffs, m_phys.
            if(DeclareCoeffPhysArrays)
            {
                m_coeffs = Array<OneD, NekDouble>(m_ncoeffs);
                m_phys   = Array<OneD, NekDouble>(m_npoints);
            }
        }

        /**
         * Fills the list of local expansions with the segments in one
         * subdomain specified in an inputfile by \a domain. This 
         * CompositeMap contains a list of Composites which define the 
         * subdomains.
         * @param   domain      A domain, comprising of one or more composite
         *                      regions.
         * @param   i           Index of currently processed subdomain
         * @param   graph1D     A mesh, containing information about the
         *                      domain and the spectral/hp element expansion.
         * @param   DeclareCoeffPhysArrays If true, create general segment expansions
         *                      instead of just normal segment expansions.
         */
        ExpList1D::ExpList1D(const LibUtilities::SessionReaderSharedPtr &pSession,
							 const SpatialDomains::CompositeMap &domain,
                             const SpatialDomains::MeshGraphSharedPtr &graph1D,
							 int i,
                             const bool DeclareCoeffPhysArrays):
		ExpList(pSession)
        {
            int id=0;
            SpatialDomains::Composite comp;
            SpatialDomains::CompositeMap::const_iterator compIt;
            SpatialDomains::SegGeomSharedPtr SegmentGeom;
            LocalRegions::SegExpSharedPtr seg;
			
            int offset = 0;
            const SpatialDomains::ExpansionMap &expansions = graph1D->GetExpansions();
            SpatialDomains::ExpansionMap::const_iterator expIt;
			
			
            // Find the correct composite region to process
            compIt = domain.begin();
            for(int k = 0; k < i; ++k)
            {
                offset += compIt->second->size();
                ++compIt;
            }	
            comp = compIt->second;
			
            //Find the correct expansion start point for the current composite
            expIt = expansions.begin();
            for(int k = 0; k < offset; ++k)
            {
                ++expIt;
            }	
			
            // Process each expansion in the region.
            for(int j = 0; j < compIt->second->size(); ++j)
            {
                if ((SegmentGeom = boost::dynamic_pointer_cast<
                         SpatialDomains::SegGeom>(
                             (*compIt->second)[j])))
                {					
                    // Retrieve the basis key from the expansion.
                    LibUtilities::BasisKey bkey = expIt->second->m_basisKeyVector[0];
																				
                    seg = MemoryManager<LocalRegions::SegExp>
                        ::AllocateSharedPtr(bkey, SegmentGeom);
					
                    // Add the segment to the expansion list.
                    seg->SetElmtId(id++);
                    (*m_exp).push_back(seg);
					
                    expIt++;
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a SegGeom failed");
                }
            }
			
			
            // Setup Default optimisation information.
            int nel = GetExpSize();
            m_globalOptParam = MemoryManager<NekOptimize::GlobalOptParam>
			::AllocateSharedPtr(nel);
			
            // Allocate storage for data and populate element offset lists.
            SetCoeffPhysOffsets();
			
            // Set up m_coeffs, m_phys.
            if(DeclareCoeffPhysArrays)
            {
                m_coeffs = Array<OneD, NekDouble>(m_ncoeffs);
                m_phys   = Array<OneD, NekDouble>(m_npoints);
            }
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
                    const Array<OneD,const ExpListSharedPtr>  &bndConstraint,
                    const Array<OneD, const SpatialDomains
                                           ::BoundaryConditionShPtr>  &bndCond,
                    const StdRegions::StdExpansionVector &locexp,
                    const SpatialDomains::MeshGraphSharedPtr &graph2D,
                    const map<int,int> &periodicEdges,
                    const bool DeclareCoeffPhysArrays,
                    const std::string variable):
            ExpList()
        {
            int i, j, id, elmtid=0;
            map<int,int> EdgeDone;
            map<int,int> NormalSet;

            SpatialDomains::Geometry1DSharedPtr SegGeom;
            LocalRegions::SegExpSharedPtr Seg;

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
                        SegGeom = bndConstraint[i]->GetExp(j)->GetGeom1D();

                        Seg = MemoryManager<LocalRegions::SegExp>
                                            ::AllocateSharedPtr(bkey, SegGeom);
                        EdgeDone[SegGeom->GetEid()] = elmtid;

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

                        Seg = MemoryManager<LocalRegions::SegExp>
                                        ::AllocateSharedPtr(EdgeBkey, SegGeom);
                        EdgeDone[id] = elmtid;

                        /*
                        if (periodicEdges.count(id) > 0)
                        {
                            EdgeDone[abs(periodicEdges.find(id)->second)] = elmtid;
                        }
                        */

                        Seg->SetElmtId(elmtid++);
                        (*m_exp).push_back(Seg);
                    }
                    else // variable modes/points
                    {
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
                            Seg = MemoryManager<LocalRegions::SegExp>
                                    ::AllocateSharedPtr(EdgeBkey, SegGeom);
                            Seg->SetElmtId(EdgeDone[id]);
                            (*m_exp)[EdgeDone[id]] = Seg;
                            NormalSet.erase(id);
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

            // Setup Default optimisation information.
            int nel = GetExpSize();
            m_globalOptParam = MemoryManager<NekOptimize::GlobalOptParam>
                ::AllocateSharedPtr(nel);

            // Set up offset information and array sizes
            SetCoeffPhysOffsets();

            // Set up m_coeffs, m_phys.
            if(DeclareCoeffPhysArrays)
            {
                m_coeffs = Array<OneD, NekDouble>(m_ncoeffs);
                m_phys   = Array<OneD, NekDouble>(m_npoints);
            }
        }

        /**
         * Each expansion (local element) is processed in turn to
         * determine the number of coefficients and physical data
         * points it contributes to the domain. Three arrays,
         * #m_coeff_offset, #m_phys_offset and #m_offset_elmt_id, are
         * also initialised and updated to store the data offsets of
         * each element in the #m_coeffs and #m_phys arrays, and the
         * element id that each consecutive block is associated
         * respectively.
         */
        void ExpList1D::SetCoeffPhysOffsets()
        {
            int i;

            // Set up offset information and array sizes
            m_coeff_offset   = Array<OneD,int>(m_exp->size());
            m_phys_offset    = Array<OneD,int>(m_exp->size());
            m_offset_elmt_id = Array<OneD,int>(m_exp->size());

            m_ncoeffs = m_npoints = 0;

            for(i = 0; i < m_exp->size(); ++i)
            {
                m_coeff_offset[i]   = m_ncoeffs;
                m_phys_offset [i]   = m_npoints;
                m_offset_elmt_id[i] = i;
                m_ncoeffs += (*m_exp)[i]->GetNcoeffs();
                m_npoints += (*m_exp)[i]->GetTotPoints();
            }
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
                    NekDouble a = total_breaks[j];
                    NekDouble b = total_breaks[j+1];

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
		

	void ExpList1D::SetUpPhysTangents(
		const StdRegions::StdExpansionVector &locexp)
	{
	    map<int, int> EdgeGID;
	    int i,cnt,n,id;
	    
	    //setup map of all global ids along booundary
	    for(cnt = i=0; i< (*m_exp).size(); ++i)
	    {
	        id = (*m_exp)[i]->GetGeom1D()->GetEid();
	        EdgeGID[id] = cnt++;
	    }
	    
	    //loop over elements and find edges that match
	    for(cnt = n =0; n< locexp.size(); ++n)
	    {
	       for(i=0; i < locexp[n]->GetNedges(); ++i)
	       {
	       	  id = locexp[n]->GetGeom2D()->GetEid(i);
	       	  if(EdgeGID.count(id)> 0)
	       	  {
	       	      (*m_exp)[EdgeGID.find(id)->second]
	       	      			->SetUpPhysTangents(locexp[n],i);
	       	  }
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
            Array<OneD,Array<OneD,NekDouble> > locnormals;

            // Assume whole array is of same coordinate dimension
            int coordim = (*m_exp)[0]->GetGeom1D()->GetCoordim();

            ASSERTL1(normals.num_elements() >= coordim,
                     "Output vector does not have sufficient dimensions to "
                     "match coordim");

            // Process each expansion.
            for(i = 0; i < m_exp->size(); ++i)
            {
                LocalRegions::Expansion1DSharedPtr loc_exp = 
                    boost::dynamic_pointer_cast<
                        LocalRegions::Expansion1D>((*m_exp)[i]);
                LocalRegions::Expansion2DSharedPtr loc_elmt = 
                    loc_exp->GetLeftAdjacentElementExp();

                // Get the number of points and normals for this expansion.
                e_npoints  = (*m_exp)[i]->GetNumPoints(0);
                locnormals = loc_elmt->GetEdgeNormal(
                    loc_exp->GetLeftAdjacentElementEdge());

                // Get the physical data offset for this expansion.
                offset = m_phys_offset[i];
                
                // Process each point in the expansion.
                for(j = 0; j < e_npoints; ++j)
                {
                    // Process each spatial dimension and copy the values into
                    // the output array.
                    for(k = 0; k < coordim; ++k)
                    {
                        //normals[k][offset+j] = locnormals[k*e_npoints + j];
                        normals[k][offset+j] = locnormals[k][j];
                    }
                }
            }
        }

/*        void ExpList1D::v_GetTangents(
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &tangents)
        {
            int i,j,k,e_npoints,offset;
            Array<OneD,Array<OneD, NekDouble> > loctangent;

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
                for (j = 0; j < 2; ++j)
                {
                    loctangent = (*m_exp)[i]->GetMetricInfo()->GetTangent(j);

                    // Get the physical data offset for this expansion.
                    offset = m_phys_offset[i];
                    for (k = 0; k < coordim; ++k)
                    {
                        Vmath::Vcopy(e_npoints, &(loctangent[k][0]), 1,
                                                &(tangents[j][k][offset]), 1);
                    }
                }
            }
*/                // Process each point in the expansion.
/*                for(j = 0; j < e_npoints; ++j)
                {
                    // Process each spatial dimension and copy the values into
                    // the output array.
                    for(k = 0; k < coordim; ++k)
                    {
                        //normals[k][offset+j] = locnormals[k*e_npoints + j];
                        normals[k][offset+j] = locnormals[k][j];
                    }
                }*/

//        }

        /**
         *
         */
//        void ExpList1D::v_SetUpPhysNormals(
//                                const StdRegions::StdExpansionVector &locexp)
//        {
//            SetUpPhysNormals(locexp);
//        }

        void ExpList1D::v_SetUpPhysTangents(
                    const StdRegions::StdExpansionVector &locexp)
        {
            SetUpPhysTangents(locexp);
        }

        void ExpList1D::v_ReadGlobalOptimizationParameters()
        {
//            Array<OneD, int> NumShape(1,0);
//            NumShape[0] = GetExpSize();
//
//            int one = 1;
//            m_globalOptParam = MemoryManager<NekOptimize::GlobalOptParam>
//                ::AllocateSharedPtr(m_session,one,NumShape);
        }


        /**
         * Writes out the header for a <PIECE> VTK XML segment describing the
         * geometric information which comprises this element. This includes
         * vertex coordinates for each quadrature point, vertex connectivity
         * information, cell types and cell offset data.
         *
         * @param   outfile     Output stream to write data to.
         */
        void ExpList1D::v_WriteVtkPieceHeader(std::ofstream &outfile, int expansion)
        {
            int i,j;
            int nquad0 = (*m_exp)[expansion]->GetNumPoints(0);
            int ntot = nquad0;
            int ntotminus = (nquad0-1);

            Array<OneD,NekDouble> coords[3];
            coords[0] = Array<OneD,NekDouble>(ntot);
            coords[1] = Array<OneD,NekDouble>(ntot);
            coords[2] = Array<OneD,NekDouble>(ntot);
            (*m_exp)[expansion]->GetCoords(coords[0],coords[1],coords[2]);

            outfile << "    <Piece NumberOfPoints=\""
                    << ntot << "\" NumberOfCells=\""
                    << ntotminus << "\">" << endl;
            outfile << "      <Points>" << endl;
            outfile << "        <DataArray type=\"Float32\" "
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

/**
* $Log: ExpList1D.cpp,v $
* Revision 1.44  2009/12/18 18:53:14  bnelson
* Fixed windows compiler warnings.
*
* Revision 1.43  2009/12/15 18:09:02  cantwell
* Split GeomFactors into 1D, 2D and 3D
* Added generation of tangential basis into GeomFactors
* Updated ADR2DManifold solver to use GeomFactors for tangents
* Added <GEOMINFO> XML session section support in MeshGraph
* Fixed const-correctness in VmathArray
* Cleaned up LocalRegions code to generate GeomFactors
* Removed GenSegExp
* Temporary fix to SubStructuredGraph
* Documentation for GlobalLinSys and GlobalMatrix classes
*
* Revision 1.42  2009/11/19 23:30:36  cantwell
* Documentation for ExpList2D and GlobalMatrixKey
* Updated doxygen pages.
*
* Revision 1.41  2009/11/18 17:12:29  cantwell
* Added documentation to ExpList1D.
*
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
