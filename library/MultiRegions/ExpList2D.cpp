///////////////////////////////////////////////////////////////////////////////
//
// File ExpList2D.cpp
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
// Description: Expansion list 2D definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LocalRegions/TriExp.h>
//#include <LocalRegions/GenTriExp.h>
#include <LocalRegions/QuadExp.h>
//#include <LocalRegions/GenQuadExp.h>
#include <LocalRegions/NodalTriExp.h>
//#include <LocalRegions/GenNodalTriExp.h>

#include <MultiRegions/ExpList2D.h>

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class ExpList2D
         *
         * This multi-elemental expansion, which does not exhibit any coupling
         * between the expansion on the separate elements, can be formulated
         * as,
         * \f[u^{\delta}(\boldsymbol{x}_i)=\sum_{e=1}^{{N_{\mathrm{el}}}}
         * \sum_{n=0}^{N^{e}_m-1}\hat{u}_n^e\phi_n^e(\boldsymbol{x}_i).\f]
         * where \f${N_{\mathrm{el}}}\f$ is the number of elements and
         * \f$N^{e}_m\f$ is the local elemental number of expansion modes.
         * This class inherits all its variables and member functions from the
         * base class #ExpList.
         */

        /**
         *
         */
        ExpList2D::ExpList2D():
            ExpList()
        {
        }


        /**
         *
         */
        ExpList2D::~ExpList2D()
        {
        }


        /**
         * @param   In          ExpList2D object to copy.
         */
        ExpList2D::ExpList2D(const ExpList2D &In):
            ExpList(In)
        {
        }


        /**
         * Given a mesh \a graph2D, containing information about the domain and
         * the spectral/hp element expansion, this constructor fills the list
         * of local expansions \texttt{m_exp} with the proper expansions,
         * calculates the total number of quadrature points
         * \f$\boldsymbol{x}_i\f$ and local expansion coefficients
         * \f$\hat{u}^e_n\f$ and allocates memory for the arrays #m_coeffs
         * and #m_phys.
         *
         * @param   graph2D     A mesh, containing information about the domain
         *                      and the spectral/hp element expansion.
         */
        ExpList2D::ExpList2D(SpatialDomains::MeshGraph2D &graph2D,
                                bool UseGenExp):
            ExpList()
        {
            int i,elmtid=0;
            LocalRegions::TriExpSharedPtr      tri;
            LocalRegions::NodalTriExpSharedPtr Ntri;
            LibUtilities::PointsType           TriNb;
            LocalRegions::QuadExpSharedPtr     quad;
            SpatialDomains::Composite          comp;

            const SpatialDomains::ExpansionVector &expansions
                                        = graph2D.GetExpansions();

            m_coeff_offset = Array<OneD,int>(expansions.size());
            m_phys_offset = Array<OneD,int>(expansions.size());

            for(i = 0; i < expansions.size(); ++i)
            {
                SpatialDomains::TriGeomSharedPtr  TriangleGeom;
                SpatialDomains::QuadGeomSharedPtr QuadrilateralGeom;

                if(TriangleGeom = boost::dynamic_pointer_cast<SpatialDomains
                                        ::TriGeom>(expansions[i]->m_GeomShPtr))
                {
                    LibUtilities::BasisKey TriBa
                                        = expansions[i]->m_BasisKeyVector[0];
                    LibUtilities::BasisKey TriBb
                                        = expansions[i]->m_BasisKeyVector[1];

                    // This is not elegantly implemented needs re-thinking.
                    if(TriBa.GetBasisType() == LibUtilities::eGLL_Lagrange)
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

                    m_coeff_offset[i] = m_ncoeffs;
                    m_phys_offset[i] = m_npoints;
                    m_ncoeffs += (TriBa.GetNumModes()*(TriBa.GetNumModes()+1))/2
                                    + TriBa.GetNumModes()*(TriBb.GetNumModes()
                                    -TriBa.GetNumModes());
                    m_npoints += TriBa.GetNumPoints()*TriBb.GetNumPoints();
                }
                else if(QuadrilateralGeom = boost::dynamic_pointer_cast<
                        SpatialDomains::QuadGeom>(expansions[i]->m_GeomShPtr))
                {
                    LibUtilities::BasisKey QuadBa 
                                        = expansions[i]->m_BasisKeyVector[0];
                    LibUtilities::BasisKey QuadBb 
                                        = expansions[i]->m_BasisKeyVector[1];

                    quad = MemoryManager<LocalRegions::QuadExp>
                                        ::AllocateSharedPtr(QuadBa,QuadBb,
                                                            QuadrilateralGeom);
                    quad->SetElmtId(elmtid++);
                    (*m_exp).push_back(quad);

                    m_coeff_offset[i] = m_ncoeffs;
                    m_phys_offset[i] = m_npoints;
                    m_ncoeffs += QuadBa.GetNumModes()*QuadBb.GetNumModes();
                    m_npoints += QuadBa.GetNumPoints()*QuadBb.GetNumPoints();
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a proper Geometry2D "
                                   "failed");
                }

            }
            m_coeffs = Array<OneD, NekDouble>(m_ncoeffs);
            m_phys   = Array<OneD, NekDouble>(m_npoints);
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
         */
         ExpList2D::ExpList2D(const LibUtilities::BasisKey &TriBa, 
                              const LibUtilities::BasisKey &TriBb, 
                              const LibUtilities::BasisKey &QuadBa, 
                              const LibUtilities::BasisKey &QuadBb, 
                              const SpatialDomains::MeshGraph2D &graph2D,
                              const LibUtilities::PointsType TriNb):
             ExpList()
         {
             int i,elmtid=0;
             LocalRegions::TriExpSharedPtr tri;
             LocalRegions::NodalTriExpSharedPtr Ntri;
             LocalRegions::QuadExpSharedPtr quad;
             SpatialDomains::Composite comp;
 
             const SpatialDomains::ExpansionVector &expansions = graph2D.GetExpansions();
             m_ncoeffs = 0;
             m_npoints = 0;
             
             m_transState = eNotSet; 
             m_physState  = false;
             
             m_coeff_offset = Array<OneD,int>(expansions.size());
             m_phys_offset = Array<OneD,int>(expansions.size());
 
             for(i = 0; i < expansions.size(); ++i)
             {
                 SpatialDomains::TriGeomSharedPtr TriangleGeom;
                 SpatialDomains::QuadGeomSharedPtr QuadrilateralGeom;
                 
                 if(TriangleGeom = boost::dynamic_pointer_cast<SpatialDomains::TriGeom>(expansions[i]->m_GeomShPtr))
                 {
                     if(TriNb < LibUtilities::SIZE_PointsType)
                     {
                         Ntri = MemoryManager<LocalRegions::NodalTriExp>::AllocateSharedPtr(TriBa,TriBb,TriNb,TriangleGeom);
                         Ntri->SetElmtId(elmtid++);
                         (*m_exp).push_back(Ntri);
                     }
                     else
                     {
                         tri = MemoryManager<LocalRegions::TriExp>::AllocateSharedPtr(TriBa,TriBb,TriangleGeom);
                         tri->SetElmtId(elmtid++);
                         (*m_exp).push_back(tri);
                     }
                      
                     m_coeff_offset[i] = m_ncoeffs;
                     m_phys_offset[i] = m_npoints;
                     m_ncoeffs += (TriBa.GetNumModes()*(TriBa.GetNumModes()+1))/2 
                         + TriBa.GetNumModes()*(TriBb.GetNumModes()-TriBa.GetNumModes());
                     m_npoints += TriBa.GetNumPoints()*TriBb.GetNumPoints();
                 }
                 else if(QuadrilateralGeom = boost::dynamic_pointer_cast<SpatialDomains::QuadGeom>(expansions[i]->m_GeomShPtr))
                 {
                     quad = MemoryManager<LocalRegions::QuadExp>::AllocateSharedPtr(QuadBa,QuadBb,QuadrilateralGeom);
                     quad->SetElmtId(elmtid++);
                     (*m_exp).push_back(quad);
 
                     m_coeff_offset[i] = m_ncoeffs;
                     m_phys_offset[i] = m_npoints;
                     m_ncoeffs += QuadBa.GetNumModes()*QuadBb.GetNumModes();
                     m_npoints += QuadBa.GetNumPoints()*QuadBb.GetNumPoints();
                 }
                 else
                 {
                     ASSERTL0(false,"dynamic cast to a proper Geometry2D failed");
                 }  
                 
             }
             
             m_coeffs = Array<OneD, NekDouble>(m_ncoeffs);
             m_phys   = Array<OneD, NekDouble>(m_npoints);
         }


        /**
         * Fills the list of local expansions with the segments from the 3D
         * mesh specified by \a domain. This CompositeVector contains a list of
         * Composites which define the Neumann boundary.
         * @see     ExpList2D#ExpList2D(SpatialDomains::MeshGraph2D&)
         *          for details.
         * @param   domain      A domain, comprising of one or more composite
         *                      regions.
         * @param   graph3D     A mesh, containing information about the domain
         *                      and the spectral/hp element expansions.
         */
        ExpList2D::ExpList2D(   const SpatialDomains::CompositeVector &domain,
                                SpatialDomains::MeshGraph3D &graph3D,
                                bool UseGenExp):
            ExpList()
        {
            int i,j,elmtid=0;
            int nel=0;
            int cnt = 0;

            SpatialDomains::Composite comp;
            SpatialDomains::TriGeomSharedPtr TriangleGeom;
            SpatialDomains::QuadGeomSharedPtr QuadrilateralGeom;

            LocalRegions::TriExpSharedPtr tri;
            LocalRegions::NodalTriExpSharedPtr Ntri;
            LibUtilities::PointsType TriNb;
            LocalRegions::QuadExpSharedPtr quad;

            for(i = 0; i < domain.size(); ++i)
            {
                nel += (domain[i])->size();
            }

            m_coeff_offset = Array<OneD,int>(nel,0);
            m_phys_offset  = Array<OneD,int>(nel,0);

            for(i = 0; i < domain.size(); ++i)
            {
                comp = domain[i];

                for(j = 0; j < comp->size(); ++j)
                {
                    if(TriangleGeom = boost::dynamic_pointer_cast<
                                        SpatialDomains::TriGeom>((*comp)[j]))
                    {
                        LibUtilities::BasisKey TriBa 
                                    = graph3D.GetFaceBasisKey(TriangleGeom,0);
                        LibUtilities::BasisKey TriBb 
                                    = graph3D.GetFaceBasisKey(TriangleGeom,1);

                        if((graph3D.GetExpansions())[0]->m_BasisKeyVector[0]
                                .GetBasisType() == LibUtilities::eGLL_Lagrange)
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
                                        ::AllocateSharedPtr(TriBa,TriBb,
                                                            TriangleGeom);
                            tri->SetElmtId(elmtid++);
                            (*m_exp).push_back(tri);
                        }

                        m_coeff_offset[cnt] = m_ncoeffs;
                        m_phys_offset[cnt++] = m_npoints;
                        m_ncoeffs 
                            += (TriBa.GetNumModes()*(TriBa.GetNumModes()+1))/2
                                + TriBa.GetNumModes()*(TriBb.GetNumModes()
                                -TriBa.GetNumModes());
                        m_npoints += TriBa.GetNumPoints()*TriBb.GetNumPoints();
                    }
                    else if(QuadrilateralGeom = boost::dynamic_pointer_cast<
                                        SpatialDomains::QuadGeom>((*comp)[j]))
                    {
                        LibUtilities::BasisKey QuadBa 
                                = graph3D.GetFaceBasisKey(QuadrilateralGeom,0);
                        LibUtilities::BasisKey QuadBb 
                                = graph3D.GetFaceBasisKey(QuadrilateralGeom,0);

                        quad = MemoryManager<LocalRegions::QuadExp>
                                        ::AllocateSharedPtr(QuadBa,QuadBb,
                                                            QuadrilateralGeom);
                        quad->SetElmtId(elmtid++);
                        (*m_exp).push_back(quad);

                        m_coeff_offset[cnt] = m_ncoeffs;
                        m_phys_offset[cnt++] = m_npoints;
                        m_ncoeffs += QuadBa.GetNumModes()*QuadBb.GetNumModes();
                        m_npoints += QuadBa.GetNumPoints()
                                        * QuadBb.GetNumPoints();
                    }
                    else
                    {
                        ASSERTL0(false,"dynamic cast to a proper Geometry2D "
                                       "failed");
                    }
                }

            }
            m_coeffs = Array<OneD, NekDouble>(m_ncoeffs);
            m_phys   = Array<OneD, NekDouble>(m_npoints);

        }


        /**
         * @param   graph2D     A mesh containing information about the domain
         *                      and the spectral/hp element expansions.
         * @param   bcs         Information about the boundary conditions.
         * @param   variable    Specifies the field.
         * @param   bndCondExpansions   Array of ExpList1D objects each
         *                      containing a 1D spectral/hp element expansion
         *                      on a single boundary region.
         * @param   bncCondition    Array of BoundaryCondition objects which
         *                      contain information about the boundary
         *                      conditions on the different boundary regions.
         */
        void ExpList2D::SetBoundaryConditionExpansion(
                        SpatialDomains::MeshGraph2D &graph2D,
                        SpatialDomains::BoundaryConditions &bcs,
                        const std::string variable,
                        Array<OneD, ExpList1DSharedPtr> &bndCondExpansions,
                        Array<OneD, SpatialDomains::BoundaryConditionShPtr>
                                                                &bndConditions)
        {
            int i;
            int cnt  = 0;

            SpatialDomains::BoundaryRegionCollection &bregions 
                                        = bcs.GetBoundaryRegions();
            SpatialDomains::BoundaryConditionCollection &bconditions 
                                        = bcs.GetBoundaryConditions();

            MultiRegions::ExpList1DSharedPtr       locExpList;
            SpatialDomains::BoundaryConditionShPtr locBCond;

            int nbnd = bregions.size();

            cnt=0;
            // list Dirichlet boundaries first
            for(i = 0; i < nbnd; ++i)
            {
                locBCond = (*(bconditions[i]))[variable];
                if(locBCond->GetBoundaryConditionType() 
                                        == SpatialDomains::eDirichlet)
                {
                    locExpList = MemoryManager<MultiRegions::ExpList1D>
                                        ::AllocateSharedPtr(*(bregions[i]),
                                                            graph2D);
                    bndCondExpansions[cnt]  = locExpList;
                    bndConditions[cnt++]    = locBCond;
                } // end if Dirichlet
            }
            // then, list the other (non-periodic) boundaries
            for(i = 0; i < nbnd; ++i)
            {
                locBCond = (*(bconditions[i]))[variable];
                if(locBCond->GetBoundaryConditionType() 
                                        == SpatialDomains::eNeumann)
                {
                    bool UseGenSegExp = true;
                    locExpList = MemoryManager<MultiRegions::ExpList1D>
                                        ::AllocateSharedPtr(*(bregions[i]),
                                                            graph2D);
                    bndCondExpansions[cnt]  = locExpList;
                    bndConditions[cnt++]    = locBCond;
                }
                else if((locBCond->GetBoundaryConditionType() 
                                        != SpatialDomains::eDirichlet) 
                        && (locBCond->GetBoundaryConditionType() 
                                        != SpatialDomains::ePeriodic))
                {
                    ASSERTL0(false,"This type of BC not implemented yet");
                }
            }
        }


        /**
         * Evaluates the boundary condition expansions, \a bndCondExpansions,
         * given the information provided by \a bndConditions.
         * @param   time        The time at which the boundary conditions
         *                      should be evaluated.
         * @param   bndCondExpansions   List of boundary conditions.
         * @param   bndConditions   Information about the boundary conditions.
         */
        void ExpList2D::EvaluateBoundaryConditions(
                        const NekDouble time,
                        Array<OneD, ExpList1DSharedPtr> &bndCondExpansions,
                        Array<OneD, SpatialDomains::BoundaryConditionShPtr>
                                                                &bndConditions)
        {
            int i,j;
            int npoints;
            int nbnd = bndCondExpansions.num_elements();
			
            MultiRegions::ExpList1DSharedPtr locExpList;

            for(i = 0; i < nbnd; ++i)
            {
			  if(time == 0.0 || bndConditions[i]->GetUserDefined().GetEquation()=="TimeDependent")
			  {
                locExpList = bndCondExpansions[i];
                npoints = locExpList->GetNpoints();
				Array<OneD,NekDouble> x0(npoints,0.0);
                Array<OneD,NekDouble> x1(npoints,0.0);
                Array<OneD,NekDouble> x2(npoints,0.0);

                locExpList->GetCoords(x0,x1,x2);

                if(bndConditions[i]->GetBoundaryConditionType() 
                                        == SpatialDomains::eDirichlet)
                {
                    for(j = 0; j < npoints; j++)
                    {
                        (locExpList->UpdatePhys())[j] 
                            = (boost::static_pointer_cast<
                                    SpatialDomains::DirichletBoundaryCondition
                                >(bndConditions[i])->m_DirichletCondition
                              ).Evaluate(x0[j],x1[j],x2[j],time);
                    }

                    locExpList->FwdTrans_BndConstrained(locExpList->GetPhys(),
                                        locExpList->UpdateCoeffs());
                }
                else if(bndConditions[i]->GetBoundaryConditionType() 
                                        == SpatialDomains::eNeumann)
                {
                    for(j = 0; j < npoints; j++)
                    {
                        (locExpList->UpdatePhys())[j] 
                                = (boost::static_pointer_cast<
                                        SpatialDomains::NeumannBoundaryCondition
                                    >(bndConditions[i])->m_NeumannCondition
                                  ).Evaluate(x0[j],x1[j],x2[j],time);
                    }

                    locExpList->IProductWRTBase(locExpList->GetPhys(),
                                        locExpList->UpdateCoeffs());
                }
                else
                {
                    ASSERTL0(false,"This type of BC not implemented yet");
                }
			  }
            }
        }


        /**
         * @param   graph2D     A mesh containing information about the domain
         *                      and the spectral/hp element expansion.
         * @param   bcs         Information about the boundary conditions.
         * @param   variable    Specifies the field.
         * @param   periodicVerts   Vector of Maps into which the list of
         *                      periodic vertices is placed, one map for each
         *                      boundary region.
         * @param   periodicEdges   Map into which the list of periodic
         *                      edges is placed.
         */
        void ExpList2D::GetPeriodicEdges(
                        SpatialDomains::MeshGraph2D &graph2D,
                        SpatialDomains::BoundaryConditions &bcs,
                        const std::string variable,
                        vector<map<int,int> >& periodicVerts,
                        map<int,int>& periodicEdges)
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

            SpatialDomains::SegGeomSharedPtr segmentGeom1;
            SpatialDomains::SegGeomSharedPtr segmentGeom2;

            SpatialDomains::ElementEdgeVectorSharedPtr element1;
            SpatialDomains::ElementEdgeVectorSharedPtr element2;

            StdRegions::EdgeOrientation orient1;
            StdRegions::EdgeOrientation orient2;

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
                    region2ID = (boost::static_pointer_cast<
                                    SpatialDomains::PeriodicBoundaryCondition
                                    >(locBCond))->m_ConnectedBoundaryRegion;

                    if(doneBndRegions.count(region1ID)==0)
                    {
                        ASSERTL0(bregions[region1ID]->size() 
                                        == bregions[region2ID]->size(),
                                 "Size of the 2 periodic boundary regions "
                                 "should be equal");


                        map<int,int> periodicVertices;

                        for(j = 0; j < bregions[region1ID]->size(); j++)
                        {
                            comp1 = (*(bregions[region1ID]))[j];
                            comp2 = (*(bregions[region2ID]))[j];

                            ASSERTL0(comp1->size() == comp2->size(),
                                     "Size of the 2 periodic composites should "
                                     "be equal");

                            for(k = 0; k < comp1->size(); k++)
                            {
                                if(!(segmentGeom1 
                                        = boost::dynamic_pointer_cast<
                                          SpatialDomains::SegGeom>((*comp1)[k]))
                                    || !(segmentGeom2 
                                        = boost::dynamic_pointer_cast<
                                          SpatialDomains::SegGeom>((*comp2)[k]))
                                    )
                                {
                                    ASSERTL0(false,"dynamic cast to a SegGeom "
                                                   "failed");
                                }

                                // Extract the periodic edges
                                periodicEdges[segmentGeom1->GetEid()] 
                                        = segmentGeom2->GetEid();
                                periodicEdges[segmentGeom2->GetEid()] 
                                        = segmentGeom1->GetEid();

                                // Extract the periodic vertices
                                element1 = graph2D
                                            .GetElementsFromEdge(segmentGeom1);
                                element2 = graph2D
                                            .GetElementsFromEdge(segmentGeom2);

                                ASSERTL0(element1->size()==1,
                                         "The periodic boundaries belong to "
                                         "more than one element of the mesh");
                                ASSERTL0(element2->size()==1,
                                         "The periodic boundaries belong to "
                                         "more than one element of the mesh");

                                orient1 = (boost::dynamic_pointer_cast<
                                            SpatialDomains::Geometry2D>(
                                                (*element1)[0]->m_Element)
                                          )->GetEorient((*element1)[0]
                                                        ->m_EdgeIndx);
                                orient2 = (boost::dynamic_pointer_cast<
                                            SpatialDomains::Geometry2D>(
                                                (*element2)[0]->m_Element)
                                          )->GetEorient((*element2)[0]
                                                        ->m_EdgeIndx);

                                if(orient1!=orient2)
                                {
                                    periodicVertices[segmentGeom1->GetVid(0)] 
                                        = segmentGeom2->GetVid(0);
                                    periodicVertices[segmentGeom1->GetVid(1)] 
                                        = segmentGeom2->GetVid(1);
                                }
                                else
                                {
                                    periodicVertices[segmentGeom1->GetVid(0)] 
                                        = segmentGeom2->GetVid(1);
                                    periodicVertices[segmentGeom1->GetVid(1)] 
                                        = segmentGeom2->GetVid(0);
                                }
                            }
                        }
                        periodicVerts.push_back(periodicVertices);
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
         * Sets up the normals on all edges of expansions in the domain.
         * @param   locexp      Complete list of domain expansions.
         */
        void ExpList2D::SetUpPhysNormals(
                                const StdRegions::StdExpansionVector &locexp)
        {

        }
        
        /**
         * For each local element, copy the normals stored in the element list
         * into the array \a normals.
         * @param   normals     Multidimensional array in which to copy normals
         *                      to. Must have dimension equal to or larger than
         *                      the spatial dimension of the elements.
         */
        void ExpList2D::GetNormals(
                                Array<OneD, Array<OneD, NekDouble> > &normals)
        {

        }
        
        /**
         *
         */
        void ExpList2D::v_SetUpPhysNormals(
                                const StdRegions::StdExpansionVector &locexp)
        {
            SetUpPhysNormals(locexp);
        }

    } //end of namespace
} //end of namespace

/**
* $Log: ExpList2D.cpp,v $
* Revision 1.36  2010/01/28 15:20:46  abolis
* Time-Dependent boundary conditions
*
* Revision 1.35  2009/12/18 18:53:14  bnelson
* Fixed windows compiler warnings.
*
* Revision 1.34  2009/12/15 18:09:03  cantwell
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
* Revision 1.33  2009/11/25 14:51:29  pvos
* Updates for added Timings directory
*
* Revision 1.32  2009/11/23 22:11:07  cantwell
* Documentation.
*
* Revision 1.31  2009/11/19 23:30:36  cantwell
* Documentation for ExpList2D and GlobalMatrixKey
* Updated doxygen pages.
*
* Revision 1.30  2009/11/02 19:15:43  cantwell
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
* Revision 1.29  2009/09/06 22:28:45  sherwin
* Updates for Navier-Stokes solver
*
* Revision 1.28  2009/05/10 23:17:12  sherwin
* Updated mainly to handle doubly periodic meshes which required modification to vertex handling from a numbering perspective
*
* Revision 1.27  2009/04/27 15:02:04  pvos
* From h-to-p efficiently updates
*
* Revision 1.26  2009/04/20 16:14:06  sherwin
* Updates for optimising bandwidth of DG solver and allowing write import on explist
*
* Revision 1.25  2009/03/04 14:17:38  pvos
* Removed all methods that take and Expansion as argument
*
* Revision 1.24  2009/01/12 10:26:35  pvos
* Added input tags for nodal expansions
*
* Revision 1.23  2009/01/06 21:05:57  sherwin
* Added virtual function calls for BwdTrans, FwdTrans and IProductWRTBase from the class ExpList. Introduced _IterPerExp versions of these methods in ExpList.cpp§
*
* Revision 1.22  2008/09/23 18:21:00  pvos
* Updates for working ProjectContField3D demo
*
* Revision 1.21  2008/08/14 22:15:51  sherwin
* Added LocalToglobalMap and DGMap and depracted LocalToGlobalBndryMap1D,2D. Made DisContField classes compatible with updated ContField formats
*
* Revision 1.20  2008/07/12 17:31:39  sherwin
* Added m_phys_offset and rename m_exp_offset to m_coeff_offset
*
* Revision 1.19  2008/05/10 18:27:33  sherwin
* Modifications necessary for QuadExp Unified DG Solver
*
* Revision 1.18  2008/04/02 22:19:54  pvos
* Update for 2D local to global mapping
*
* Revision 1.17  2008/03/18 14:14:13  pvos
* Update for nodal triangular helmholtz solver
*
* Revision 1.16  2008/03/12 15:25:45  pvos
* Clean up of the code
*
* Revision 1.15  2007/12/06 22:52:30  pvos
* 2D Helmholtz solver updates
*
* Revision 1.14  2007/07/20 02:04:12  bnelson
* Replaced boost::shared_ptr with Nektar::ptr
*
* Revision 1.13  2007/06/07 15:54:19  pvos
* Modificications to make Demos/MultiRegions/ProjectCont2D work correctly.
* Also made corrections to various ASSERTL2 calls
*
* Revision 1.12  2007/06/05 16:36:55  pvos
* Updated Explist2D ContExpList2D and corresponding demo-codes
*
**/
