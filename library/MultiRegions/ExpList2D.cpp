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

#include <iomanip>
#include <LocalRegions/TriExp.h>
#include <LocalRegions/QuadExp.h>
#include <LocalRegions/NodalTriExp.h>
#include <LocalRegions/Expansion3D.h>
#include <MultiRegions/ExpList2D.h>
#include <LibUtilities/Foundations/Interp.h>
#include <LibUtilities/Foundations/PhysGalerkinProject.h>
#include <SpatialDomains/MeshGraph3D.h>


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
        ExpList2D::ExpList2D(const ExpList2D &In, const bool DeclareCoeffPhysArrays):
            ExpList(In,DeclareCoeffPhysArrays)
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
        ExpList2D::ExpList2D(const LibUtilities::SessionReaderSharedPtr &pSession,
                const SpatialDomains::MeshGraphSharedPtr &graph2D,
                const bool DeclareCoeffPhysArrays,
                const std::string &var):
            ExpList(pSession,graph2D)
        {
            int elmtid=0;
            LocalRegions::TriExpSharedPtr      tri;
            LocalRegions::NodalTriExpSharedPtr Ntri;
            LibUtilities::PointsType           TriNb;
            LocalRegions::QuadExpSharedPtr     quad;
            SpatialDomains::Composite          comp;

            const SpatialDomains::ExpansionMap &expansions
                                        = graph2D->GetExpansions(var);

            SpatialDomains::ExpansionMap::const_iterator expIt;
            for (expIt = expansions.begin(); expIt != expansions.end(); ++expIt)
            {
                SpatialDomains::TriGeomSharedPtr  TriangleGeom;
                SpatialDomains::QuadGeomSharedPtr QuadrilateralGeom;

                if((TriangleGeom = boost::dynamic_pointer_cast<SpatialDomains
                        ::TriGeom>(expIt->second->m_geomShPtr)))
                {
                    LibUtilities::BasisKey TriBa
                                        = expIt->second->m_basisKeyVector[0];
                    LibUtilities::BasisKey TriBb
                                        = expIt->second->m_basisKeyVector[1];

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
                    m_ncoeffs += (TriBa.GetNumModes()*(TriBa.GetNumModes()+1))/2
                                    + TriBa.GetNumModes()*(TriBb.GetNumModes()
                                    -TriBa.GetNumModes());
                    m_npoints += TriBa.GetNumPoints()*TriBb.GetNumPoints();
                }
                else if((QuadrilateralGeom = boost::dynamic_pointer_cast<
                         SpatialDomains::QuadGeom>(expIt->second->m_geomShPtr)))
                {
                    LibUtilities::BasisKey QuadBa
                                        = expIt->second->m_basisKeyVector[0];
                    LibUtilities::BasisKey QuadBb
                                        = expIt->second->m_basisKeyVector[1];

                    quad = MemoryManager<LocalRegions::QuadExp>
                        ::AllocateSharedPtr(QuadBa,QuadBb,
                                            QuadrilateralGeom);
                    quad->SetElmtId(elmtid++);
                    (*m_exp).push_back(quad);

                    m_ncoeffs += QuadBa.GetNumModes()*QuadBb.GetNumModes();
                    m_npoints += QuadBa.GetNumPoints()*QuadBb.GetNumPoints();
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a proper Geometry2D "
                                   "failed");
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
         * Given an expansion vector \a expansions, containing
         * information about the domain and the spectral/hp element
         * expansion, this constructor fills the list of local
         * expansions \texttt{m_exp} with the proper expansions,
         * calculates the total number of quadrature points
         * \f$\boldsymbol{x}_i\f$ and local expansion coefficients
         * \f$\hat{u}^e_n\f$ and allocates memory for the arrays
         * #m_coeffs and #m_phys.
         *
         * @param expansions A vector containing information about the
         *                      domain and the spectral/hp element
         *                      expansion.
         */
        ExpList2D::ExpList2D(
                             const LibUtilities::SessionReaderSharedPtr &pSession,
                             const SpatialDomains::ExpansionMap &expansions,
                             const bool DeclareCoeffPhysArrays):
            ExpList(pSession)
        {
            int elmtid=0;
            LocalRegions::TriExpSharedPtr      tri;
            LocalRegions::NodalTriExpSharedPtr Ntri;
            LibUtilities::PointsType           TriNb;
            LocalRegions::QuadExpSharedPtr     quad;
            SpatialDomains::Composite          comp;

            SpatialDomains::ExpansionMapConstIter expIt;
            for(expIt = expansions.begin(); expIt != expansions.end(); ++expIt)
            {
                SpatialDomains::TriGeomSharedPtr  TriangleGeom;
                SpatialDomains::QuadGeomSharedPtr QuadrilateralGeom;

                if((TriangleGeom = boost::dynamic_pointer_cast<SpatialDomains
                        ::TriGeom>(expIt->second->m_geomShPtr)))
                {
                    LibUtilities::BasisKey TriBa
                                        = expIt->second->m_basisKeyVector[0];
                    LibUtilities::BasisKey TriBb
                                        = expIt->second->m_basisKeyVector[1];

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
                    m_ncoeffs += (TriBa.GetNumModes()*(TriBa.GetNumModes()+1))/2
                                    + TriBa.GetNumModes()*(TriBb.GetNumModes()
                                    -TriBa.GetNumModes());
                    m_npoints += TriBa.GetNumPoints()*TriBb.GetNumPoints();
                }
                else if((QuadrilateralGeom = boost::dynamic_pointer_cast<
                         SpatialDomains::QuadGeom>(expIt->second->m_geomShPtr)))
                {
                    LibUtilities::BasisKey QuadBa
                        = expIt->second->m_basisKeyVector[0];
                    LibUtilities::BasisKey QuadBb
                        = expIt->second->m_basisKeyVector[1];

                    quad = MemoryManager<LocalRegions::QuadExp>
                        ::AllocateSharedPtr(QuadBa,QuadBb,
                                            QuadrilateralGeom);
                    quad->SetElmtId(elmtid++);
                    (*m_exp).push_back(quad);

                    m_ncoeffs += QuadBa.GetNumModes()*QuadBb.GetNumModes();
                    m_npoints += QuadBa.GetNumPoints()*QuadBb.GetNumPoints();
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a proper Geometry2D "
                                   "failed");
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
          ExpList2D::ExpList2D(const LibUtilities::SessionReaderSharedPtr &pSession,
                               const LibUtilities::BasisKey &TriBa,
                               const LibUtilities::BasisKey &TriBb,
                               const LibUtilities::BasisKey &QuadBa,
                               const LibUtilities::BasisKey &QuadBb,
                               const SpatialDomains::MeshGraphSharedPtr &graph2D,
                               const LibUtilities::PointsType TriNb):
              ExpList(pSession,graph2D)
          {
              int elmtid=0;
              LocalRegions::TriExpSharedPtr tri;
              LocalRegions::NodalTriExpSharedPtr Ntri;
              LocalRegions::QuadExpSharedPtr quad;
              SpatialDomains::Composite comp;

              const SpatialDomains::ExpansionMap &expansions = graph2D->GetExpansions();
              m_ncoeffs = 0;
              m_npoints = 0;

              m_physState  = false;

              SpatialDomains::ExpansionMap::const_iterator expIt;
              for (expIt = expansions.begin(); expIt != expansions.end(); ++expIt)
              {
                  SpatialDomains::TriGeomSharedPtr TriangleGeom;
                  SpatialDomains::QuadGeomSharedPtr QuadrilateralGeom;

                  if((TriangleGeom = boost::dynamic_pointer_cast<SpatialDomains::TriGeom>(expIt->second->m_geomShPtr)))
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

                      m_ncoeffs += (TriBa.GetNumModes()*(TriBa.GetNumModes()+1))/2
                          + TriBa.GetNumModes()*(TriBb.GetNumModes()-TriBa.GetNumModes());
                      m_npoints += TriBa.GetNumPoints()*TriBb.GetNumPoints();
                  }
                  else if((QuadrilateralGeom = boost::dynamic_pointer_cast<SpatialDomains::QuadGeom>(expIt->second->m_geomShPtr)))
                  {
                      quad = MemoryManager<LocalRegions::QuadExp>::AllocateSharedPtr(QuadBa,QuadBb,QuadrilateralGeom);
                      quad->SetElmtId(elmtid++);
                      (*m_exp).push_back(quad);

                      m_ncoeffs += QuadBa.GetNumModes()*QuadBb.GetNumModes();
                      m_npoints += QuadBa.GetNumPoints()*QuadBb.GetNumPoints();
                  }
                  else
                  {
                      ASSERTL0(false,"dynamic cast to a proper Geometry2D failed");
                  }

              }

              // Setup Default optimisation information.
              int nel = GetExpSize();
              m_globalOptParam = MemoryManager<NekOptimize::GlobalOptParam>
                  ::AllocateSharedPtr(nel);

             // Set up m_coeffs, m_phys and offset arrays.
             SetCoeffPhysOffsets();
             m_coeffs = Array<OneD, NekDouble>(m_ncoeffs);
             m_phys   = Array<OneD, NekDouble>(m_npoints);

             ReadGlobalOptimizationParameters();
          }

        /**
         * Specialized constructor for trace expansions. Store 
         * expansions for the trace space used in DisContField3D
         *
         * @param   bndConstraint   Array of ExpList2D objects each containing a
         *                      2D spectral/hp element expansion on a single
         *                      boundary region.
         * @param   bndCond     Array of BoundaryCondition objects which contain
         *                      information about the boundary conditions on the
         *                      different boundary regions.
         * @param   locexp      Complete domain expansion list.
         * @param   graph3D     3D mesh corresponding to the expansion list.
         * @param   periodicFaces   List of periodic faces.
         * @param   DeclareCoeffPhysArrays If true, set up m_coeffs, m_phys arrays
         **/
        ExpList2D::ExpList2D(
            const Array<OneD,const ExpListSharedPtr> &bndConstraint,
            const Array<OneD,const SpatialDomains::BoundaryConditionShPtr>  &bndCond,
            const StdRegions::StdExpansionVector &locexp,
            const SpatialDomains::MeshGraphSharedPtr &graph3D,
            const map<int,PeriodicFace> &periodicFaces,
            const bool DeclareCoeffPhysArrays, 
            const std::string variable):
            ExpList()
        {
            int i, j, id, elmtid=0;
            map<int,int> FaceDone;
            map<int,int> NormalSet;
            SpatialDomains::Geometry2DSharedPtr FaceGeom;
            SpatialDomains::QuadGeomSharedPtr FaceQuadGeom;
            SpatialDomains::TriGeomSharedPtr FaceTriGeom;
            LocalRegions::QuadExpSharedPtr FaceQuadExp;
            LocalRegions::TriExpSharedPtr FaceTriExp;
            
            // First loop over boundary conditions to renumber
            // Dirichlet boundaries
            for(i = 0; i < bndCond.num_elements(); ++i)
            {
                if(bndCond[i]->GetBoundaryConditionType()
                                            == SpatialDomains::eDirichlet)
                {
                    for(j = 0; j < bndConstraint[i]->GetExpSize(); ++j)
                    {
                        LibUtilities::BasisKey bkey0 = bndConstraint[i]
                                    ->GetExp(j)->GetBasis(0)->GetBasisKey();
                        LibUtilities::BasisKey bkey1 = bndConstraint[i]
                                    ->GetExp(j)->GetBasis(1)->GetBasisKey();
                        FaceGeom = bndConstraint[i]->GetExp(j)->GetGeom2D();

                        //if face is a quad
                        if((FaceQuadGeom = boost::dynamic_pointer_cast<SpatialDomains::QuadGeom>(FaceGeom)))
                        {
                            FaceQuadExp = MemoryManager<LocalRegions::QuadExp>::AllocateSharedPtr(bkey0, bkey1, FaceQuadGeom);
                            FaceDone[FaceGeom->GetFid()] = elmtid;
                            FaceQuadExp->SetElmtId(elmtid++);
                            (*m_exp).push_back(FaceQuadExp);
                        }
                        //if face is a triangle
                        else if((FaceTriGeom = boost::dynamic_pointer_cast<SpatialDomains::TriGeom>(FaceGeom)))
                        {
                            FaceTriExp = MemoryManager<LocalRegions::TriExp>::AllocateSharedPtr(bkey0, bkey1, FaceTriGeom);
                            FaceDone[FaceGeom->GetFid()] = elmtid;
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
            
            // loop over all other faces and fill out other connectivities
            for(i = 0; i < locexp.size(); ++i)
            {
                for(j = 0; j < locexp[i]->GetNfaces(); ++j)
                {
                    FaceGeom = (locexp[i]->GetGeom3D())->GetFace(j);

                    id = FaceGeom->GetFid();

                    if(FaceDone.count(id)==0)
                    {
                        LibUtilities::BasisKey bkey0 = 
                            boost::dynamic_pointer_cast<SpatialDomains::MeshGraph3D>(graph3D)->GetFaceBasisKey(FaceGeom, 0, variable); 
                        LibUtilities::BasisKey bkey1 = 
                            boost::dynamic_pointer_cast<SpatialDomains::MeshGraph3D>(graph3D)->GetFaceBasisKey(FaceGeom, 1);
                        
                        //if face is a quad
                        if((FaceQuadGeom = boost::dynamic_pointer_cast<SpatialDomains::QuadGeom>(FaceGeom)))
                        {
                            FaceQuadExp = MemoryManager<LocalRegions::QuadExp>::AllocateSharedPtr(bkey0, bkey1, FaceQuadGeom);
                            
                            FaceDone[id] = elmtid;
                            /*
                            if (periodicFaces.count(id) > 0)
                            {
                                FaceDone[periodicFaces.find(id)->second.first] = elmtid;
                            }
                            */
                            FaceQuadExp->SetElmtId(elmtid++);
                            (*m_exp).push_back(FaceQuadExp);
                        }
                        //if face is a triangle
                        else if((FaceTriGeom = boost::dynamic_pointer_cast<SpatialDomains::TriGeom>(FaceGeom)))
                        {
                            FaceTriExp = MemoryManager<LocalRegions::TriExp>::AllocateSharedPtr(bkey0, bkey1, FaceTriGeom);
                            
                            FaceDone[id] = elmtid;
                            /*
                            if (periodicFaces.count(id) > 0)
                            {
                                FaceDone[periodicFaces.find(id)->second.first] = elmtid;
                            }
                            */
                            FaceTriExp->SetElmtId(elmtid++);
                            (*m_exp).push_back(FaceTriExp);
                        }
                        else
                        {
                            ASSERTL0(false,"dynamic cast to a proper face geometry failed"); 
                        }
                    }
                    //variable modes/points
                    //if for the current edge we have more modes/points at least in one direction
                    //we replace the old edge in the trace expansion with the current one
                    else
                    {
                        LibUtilities::BasisKey bkey0 = 
                            boost::dynamic_pointer_cast<SpatialDomains::MeshGraph3D>(graph3D)->GetFaceBasisKey(FaceGeom, 0); 
                        LibUtilities::BasisKey bkey1 = 
                            boost::dynamic_pointer_cast<SpatialDomains::MeshGraph3D>(graph3D)->GetFaceBasisKey(FaceGeom, 1);
                        
                        if( ((*m_exp)[FaceDone[id]]->GetNumPoints(0)
                                >= bkey0.GetNumPoints()
                             || (*m_exp)[FaceDone[id]]->GetNumPoints(1)
                                >= bkey1.GetNumPoints())
                            && ((*m_exp)[FaceDone[id]]->GetBasisNumModes(0)
                                >= bkey0.GetNumModes()
                                || (*m_exp)[FaceDone[id]]->GetBasisNumModes(1)
                                >= bkey1.GetNumModes()) )
                        {
                        }
                        else if( ((*m_exp)[FaceDone[id]]->GetNumPoints(0)
                                <= bkey0.GetNumPoints()
                                  || (*m_exp)[FaceDone[id]]->GetNumPoints(1)
                                <= bkey1.GetNumPoints())
                            && ((*m_exp)[FaceDone[id]]->GetBasisNumModes(0)
                                <= bkey0.GetNumModes()
                                || (*m_exp)[FaceDone[id]]->GetBasisNumModes(1)
                                <= bkey1.GetNumModes()) )
                        {
                            //if face is a quad
                            if((FaceQuadGeom = boost::dynamic_pointer_cast<SpatialDomains::QuadGeom>(FaceGeom)))
                            {
                                FaceQuadExp = MemoryManager<LocalRegions::QuadExp>::AllocateSharedPtr(bkey0, bkey1, FaceQuadGeom);
                                FaceQuadExp->SetElmtId(FaceDone[id]);
                                (*m_exp)[FaceDone[id]] = FaceQuadExp;
                            }
                            //if face is a triangle
                            else if((FaceTriGeom = boost::dynamic_pointer_cast<SpatialDomains::TriGeom>(FaceGeom)))
                            {
                                FaceTriExp = MemoryManager<LocalRegions::TriExp>::AllocateSharedPtr(bkey0, bkey1, FaceTriGeom);
                                FaceTriExp->SetElmtId(FaceDone[id]);
                                (*m_exp)[FaceDone[id]] = FaceTriExp;
                            }
                            else
                            {
                                ASSERTL0(false,"dynamic cast to a proper face geometry failed"); 
                            }
                            
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
          * Fills the list of local expansions with the segments from the 3D
          * mesh specified by \a domain. This CompositeMap contains a list of
          * Composites which define the Neumann boundary.
          * @see     ExpList2D#ExpList2D(SpatialDomains::MeshGraph2D&)
          *          for details.
          * @param   domain      A domain, comprising of one or more composite
          *                      regions.
          * @param   graph3D     A mesh, containing information about the domain
          *                      and the spectral/hp element expansions.
          */
         ExpList2D::ExpList2D(   const LibUtilities::SessionReaderSharedPtr &pSession,
                                 const SpatialDomains::CompositeMap &domain,
                                 const SpatialDomains::MeshGraphSharedPtr &graph3D,
                                 const std::string variable):
             ExpList(pSession,graph3D)
         {
             ASSERTL0(boost::dynamic_pointer_cast<SpatialDomains::MeshGraph3D>(graph3D),
                     "Expected a MeshGraph3D object.");

             int j,elmtid=0;
             int nel = 0;

             SpatialDomains::Composite comp;
             SpatialDomains::TriGeomSharedPtr TriangleGeom;
             SpatialDomains::QuadGeomSharedPtr QuadrilateralGeom;

             LocalRegions::TriExpSharedPtr tri;
             LocalRegions::NodalTriExpSharedPtr Ntri;
             LibUtilities::PointsType TriNb;
             LocalRegions::QuadExpSharedPtr quad;

             SpatialDomains::CompositeMap::const_iterator compIt;
             for (compIt = domain.begin(); compIt != domain.end(); ++compIt)
             {
                 nel += (compIt->second)->size();
             }

             for (compIt = domain.begin(); compIt != domain.end(); ++compIt)
             {
                 for(j = 0; j < compIt->second->size(); ++j)
                 {
                     if((TriangleGeom = boost::dynamic_pointer_cast<
                             SpatialDomains::TriGeom>((*compIt->second)[j])))
                     {
                         LibUtilities::BasisKey TriBa
                                     = boost::dynamic_pointer_cast<SpatialDomains::MeshGraph3D>(graph3D)->GetFaceBasisKey(TriangleGeom,0);
                         LibUtilities::BasisKey TriBb
                                     = boost::dynamic_pointer_cast<SpatialDomains::MeshGraph3D>(graph3D)->GetFaceBasisKey(TriangleGeom,1);

                         if(graph3D->GetExpansions().begin()->second->m_basisKeyVector[0]
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

                         m_ncoeffs
                             += (TriBa.GetNumModes()*(TriBa.GetNumModes()+1))/2
                                 + TriBa.GetNumModes()*(TriBb.GetNumModes()
                                 -TriBa.GetNumModes());
                         m_npoints += TriBa.GetNumPoints()*TriBb.GetNumPoints();
                     }
                     else if((QuadrilateralGeom = boost::dynamic_pointer_cast<
                              SpatialDomains::QuadGeom>((*compIt->second)[j])))
                     {
                         LibUtilities::BasisKey QuadBa
                                 = boost::dynamic_pointer_cast<SpatialDomains::MeshGraph3D>(graph3D)->GetFaceBasisKey(QuadrilateralGeom,0);
                         LibUtilities::BasisKey QuadBb
                                 = boost::dynamic_pointer_cast<SpatialDomains::MeshGraph3D>(graph3D)->GetFaceBasisKey(QuadrilateralGeom,0);

                         quad = MemoryManager<LocalRegions::QuadExp>
                             ::AllocateSharedPtr(QuadBa,QuadBb,
                                                 QuadrilateralGeom);
                         quad->SetElmtId(elmtid++);
                         (*m_exp).push_back(quad);

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

             // Setup Default optimisation information.
             nel = GetExpSize();
             m_globalOptParam = MemoryManager<NekOptimize::GlobalOptParam>
                 ::AllocateSharedPtr(nel);

             // Set up m_coeffs, m_phys and offset arrays.
            SetCoeffPhysOffsets();
            m_coeffs = Array<OneD, NekDouble>(m_ncoeffs);
            m_phys   = Array<OneD, NekDouble>(m_npoints);

            ReadGlobalOptimizationParameters();
            
         }
        
        /**
         * One-dimensional upwind.
         * @param   Vn          Velocity field.
         * @param   Fwd         Left state.
         * @param   Bwd         Right state.
         * @param   Upwind      Output vector.
         */
        void ExpList2D::v_Upwind(
            const Array<OneD, const NekDouble> &Vn,
            const Array<OneD, const NekDouble> &Fwd,
            const Array<OneD, const NekDouble> &Bwd,
                  Array<OneD,       NekDouble> &Upwind)
        {
            int i,j,f_npoints,offset;

            // Process each expansion.
            for(i = 0; i < m_exp->size(); ++i)
            {
                // Get the number of points and the data offset.
                f_npoints = (*m_exp)[i]->GetNumPoints(0)*
                            (*m_exp)[i]->GetNumPoints(1);
                offset    = m_phys_offset[i];
                
                // Process each point in the expansion.
                for(j = 0; j < f_npoints; ++j)
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
        void ExpList2D::v_GetNormals(
            Array<OneD, Array<OneD, NekDouble> > &normals)
        {
            int i,k,offset;
            Array<OneD,Array<OneD,NekDouble> > locnormals;
            Array<OneD, NekDouble> tmp;
            // Assume whole array is of same coordinate dimension
            int coordim = (*m_exp)[0]->GetGeom2D()->GetCoordim();
            
            ASSERTL1(normals.num_elements() >= coordim,
                     "Output vector does not have sufficient dimensions to "
                     "match coordim");
            
            // Process each expansion.
            for(i = 0; i < m_exp->size(); ++i)
            {
                LocalRegions::Expansion2DSharedPtr loc_exp = 
                    boost::dynamic_pointer_cast<
                        LocalRegions::Expansion2D>((*m_exp)[i]);
                LocalRegions::Expansion3DSharedPtr loc_elmt = 
                    loc_exp->GetLeftAdjacentElementExp();
                int faceNumber = loc_exp->GetLeftAdjacentElementFace();
                
                // Get the number of points and normals for this expansion.
                locnormals = loc_elmt->GetFaceNormal(faceNumber);
                
                // Get the physical data offset for this expansion.
                offset = m_phys_offset[i];
                
                for (k = 0; k < coordim; ++k)
                {
                    LibUtilities::Interp2D(
                        loc_elmt->GetFacePointsKey(faceNumber, 0),
                        loc_elmt->GetFacePointsKey(faceNumber, 1),
                        locnormals[k],
                        (*m_exp)[i]->GetBasis(0)->GetPointsKey(),
                        (*m_exp)[i]->GetBasis(1)->GetPointsKey(),
                        tmp = normals[k]+offset);
                }
            }
        }

        /**
         * Each expansion (local element) is processed in turn to
         * determine the number of coefficients and physical data
         * points it contributes to the domain. Three arrays,
         * #m_coeff_offset, #m_phys_offset and #m_offset_elmt_id, are
         * initialised and updated to store the data offsets of each
         * element in the #m_coeffs and #m_phys arrays, and the
         * element id that each consecutive block is associated
         * respectively.
         */
        void ExpList2D::SetCoeffPhysOffsets()
        {
            int i;

            // Set up offset information and array sizes
            m_coeff_offset   = Array<OneD,int>(m_exp->size());
            m_phys_offset    = Array<OneD,int>(m_exp->size());
            m_offset_elmt_id = Array<OneD,int>(m_exp->size());

            m_ncoeffs = m_npoints = 0;

            int cnt = 0;
            for(i = 0; i < m_exp->size(); ++i)
            {
                if((*m_exp)[i]->DetExpansionType() == StdRegions::eTriangle)
                {
                    m_coeff_offset[i]   = m_ncoeffs;
                    m_phys_offset [i]   = m_npoints;
                    m_offset_elmt_id[cnt++] = i;
                    m_ncoeffs += (*m_exp)[i]->GetNcoeffs();
                    m_npoints += (*m_exp)[i]->GetTotPoints();
                }
            }

            for(i = 0; i < m_exp->size(); ++i)
            {
                if((*m_exp)[i]->DetExpansionType() == StdRegions::eQuadrilateral)
                {
                    m_coeff_offset[i]   = m_ncoeffs;
                    m_phys_offset [i]   = m_npoints;
                    m_offset_elmt_id[cnt++] = i;
                    m_ncoeffs += (*m_exp)[i]->GetNcoeffs();
                    m_npoints += (*m_exp)[i]->GetTotPoints();
                }
            }
        }


        /**
         *
         */
        void ExpList2D::v_SetUpPhysNormals()
        {
            int i, j;
            for (i = 0; i < m_exp->size(); ++i)
            {
                for (j = 0; j < (*m_exp)[i]->GetNedges(); ++j)
                {
                    (*m_exp)[i]->ComputeEdgeNormal(j);
                }
            }
        }


        void ExpList2D::v_ReadGlobalOptimizationParameters()
        {
            Array<OneD, int> NumShape(2,0);

            for(int i = 0; i < GetExpSize(); ++i)
            {
                if((*m_exp)[i]->DetExpansionType() == StdRegions::eTriangle)
                {
                    NumShape[0] += 1;
                }
                else  // Quadrilateral element
                {
                    NumShape[1] += 1;
                }
            }

            m_globalOptParam = MemoryManager<NekOptimize::GlobalOptParam>
                ::AllocateSharedPtr(m_session,2,NumShape);
        }

        void ExpList2D::v_WriteVtkPieceHeader(std::ofstream &outfile, int expansion)
        {
            int i,j;
            int nquad0 = (*m_exp)[expansion]->GetNumPoints(0);
            int nquad1 = (*m_exp)[expansion]->GetNumPoints(1);
            int ntot = nquad0*nquad1;
            int ntotminus = (nquad0-1)*(nquad1-1);

            Array<OneD,NekDouble> coords[3];
            coords[0] = Array<OneD,NekDouble>(ntot,0.0);
            coords[1] = Array<OneD,NekDouble>(ntot,0.0);
            coords[2] = Array<OneD,NekDouble>(ntot,0.0);
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
                    outfile << setprecision(8)     << scientific 
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
                for (j = 0; j < nquad1-1; ++j)
                {
                    outfile << j*nquad0 + i << " "
                            << j*nquad0 + i + 1 << " "
                            << (j+1)*nquad0 + i + 1 << " "
                            << (j+1)*nquad0 + i << endl;
                }
            }
            outfile << endl;
            outfile << "        </DataArray>" << endl;
            outfile << "        <DataArray type=\"Int32\" "
                    << "Name=\"offsets\" format=\"ascii\">" << endl;
            for (i = 0; i < ntotminus; ++i)
            {
                outfile << i*4+4 << " ";
            }
            outfile << endl;
            outfile << "        </DataArray>" << endl;
            outfile << "        <DataArray type=\"UInt8\" "
                    << "Name=\"types\" format=\"ascii\">" << endl;
            for (i = 0; i < ntotminus; ++i)
            {
                outfile << "9 ";
            }
            outfile << endl;
            outfile << "        </DataArray>" << endl;
            outfile << "      </Cells>" << endl;
            outfile << "      <PointData>" << endl;
        }


        void ExpList2D::v_PhysInterp1DScaled(const NekDouble scale, const Array<OneD, NekDouble> &inarray, Array<OneD, NekDouble> &outarray)
        {
            int cnt,cnt1;

            cnt = cnt1 = 0;
            for(int i = 0; i < GetExpSize(); ++i)
            {
                // get new points key
                int pt0 = (*m_exp)[i]->GetNumPoints(0);
                int pt1 = (*m_exp)[i]->GetNumPoints(1);
                int npt0 = (int) pt0*scale;
                int npt1 = (int) pt1*scale;
                
                LibUtilities::PointsKey newPointsKey0(npt0,(*m_exp)[i]->GetPointsType(0));
                LibUtilities::PointsKey newPointsKey1(npt1,(*m_exp)[i]->GetPointsType(1));

                // Interpolate points; 
                LibUtilities::Interp2D((*m_exp)[i]->GetBasis(0)->GetPointsKey(),
                                       (*m_exp)[i]->GetBasis(1)->GetPointsKey(),
                                       &inarray[cnt],newPointsKey0,
                                       newPointsKey1,&outarray[cnt1]);

                cnt  += pt0*pt1;
                cnt1 += npt0*npt1;
            }
        }
        
        void ExpList2D::v_PhysGalerkinProjection1DScaled(const NekDouble scale, const Array<OneD, NekDouble> &inarray, Array<OneD, NekDouble> &outarray)
        {
            int cnt,cnt1;

            cnt = cnt1 = 0;
            for(int i = 0; i < GetExpSize(); ++i)
            {
                // get new points key
                int pt0 = (*m_exp)[i]->GetNumPoints(0);
                int pt1 = (*m_exp)[i]->GetNumPoints(1);
                int npt0 = (int) pt0*scale;
                int npt1 = (int) pt1*scale;
                
                LibUtilities::PointsKey newPointsKey0(npt0,(*m_exp)[i]->GetPointsType(0));
                LibUtilities::PointsKey newPointsKey1(npt1,(*m_exp)[i]->GetPointsType(1));

                // Project points; 
                LibUtilities::PhysGalerkinProject2D(newPointsKey0, 
                                                    newPointsKey1,
                                       &inarray[cnt],
                                       (*m_exp)[i]->GetBasis(0)->GetPointsKey(),
                                       (*m_exp)[i]->GetBasis(1)->GetPointsKey(),
                                       &outarray[cnt1]);
                
                cnt  += npt0*npt1;
                cnt1 += pt0*pt1;
            }

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
