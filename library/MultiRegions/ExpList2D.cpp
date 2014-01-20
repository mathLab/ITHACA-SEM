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
            SetExpType(e2D);
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
            SetExpType(e2D);
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
            SetExpType(e2D);

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
            SetExpType(e2D);

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
              SetExpType(e2D);

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
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const Array<OneD,const ExpListSharedPtr> &bndConstraint,
            const Array<OneD,const SpatialDomains::BoundaryConditionShPtr>  &bndCond,
            const LocalRegions::ExpansionVector &locexp,
            const SpatialDomains::MeshGraphSharedPtr &graph3D,
            const PeriodicMap &periodicFaces,
            const bool DeclareCoeffPhysArrays, 
            const std::string variable):
            ExpList()
        {
            SetExpType(e2D);

            int i, j, id, elmtid=0;
            set<int> facesDone;
            bool var_p = false;
            
            SpatialDomains::Geometry2DSharedPtr FaceGeom;
            SpatialDomains::QuadGeomSharedPtr   FaceQuadGeom;
            SpatialDomains::TriGeomSharedPtr    FaceTriGeom;
            LocalRegions::QuadExpSharedPtr      FaceQuadExp;
            LocalRegions::TriExpSharedPtr       FaceTriExp;
            LocalRegions::Expansion2DSharedPtr  exp2D;
            LocalRegions::Expansion3DSharedPtr  exp3D;
            
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
                        FaceGeom = LocalRegions::Expansion2D::FromStdExp(
                            bndConstraint[i]->GetExp(j))->GetGeom2D();

                        //if face is a quad
                        if((FaceQuadGeom = boost::dynamic_pointer_cast<
                            SpatialDomains::QuadGeom>(FaceGeom)))
                        {
                            FaceQuadExp = MemoryManager<LocalRegions::QuadExp>
                                ::AllocateSharedPtr(bkey0, bkey1, FaceQuadGeom);
                            facesDone.insert(FaceQuadGeom->GetFid());
                            FaceQuadExp->SetElmtId(elmtid++);
                            (*m_exp).push_back(FaceQuadExp);
                        }
                        //if face is a triangle
                        else if((FaceTriGeom = boost::dynamic_pointer_cast<
                                 SpatialDomains::TriGeom>(FaceGeom)))
                        {
                            FaceTriExp = MemoryManager<LocalRegions::TriExp>
                                ::AllocateSharedPtr(bkey0, bkey1, FaceTriGeom);
                            facesDone.insert(FaceTriGeom->GetFid());
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
            
            map<int, pair<SpatialDomains::Geometry2DSharedPtr,
                          pair<LibUtilities::BasisKey,
                               LibUtilities::BasisKey> > > faceOrders;
            map<int, pair<SpatialDomains::Geometry2DSharedPtr,
                          pair<LibUtilities::BasisKey,
                               LibUtilities::BasisKey> > >::iterator it;

            for(i = 0; i < locexp.size(); ++i)
            {
                exp3D = LocalRegions::Expansion3D::FromStdExp(locexp[i]);
                
                for(j = 0; j < exp3D->GetNfaces(); ++j)
                {
                    FaceGeom = exp3D->GetGeom3D()->GetFace(j);
                    id       = FaceGeom->GetFid();
                    
                    if(facesDone.count(id) != 0)
                    {
                        continue;
                    }
                    it = faceOrders.find(id);
                    
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
                        var_p = true;
                        
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
        
            LibUtilities::CommSharedPtr vComm = pSession->GetComm();
            int nproc = vComm->GetSize(); // number of processors
            int edgepr = vComm->GetRank(); // ID processor
            
            m_parallel = false;
            
            if (nproc > 1 && var_p == true)
            {
                m_parallel = true;
                
                int numModes_ref0,
                numModes_ref1,
                numModes_com0,
                numModes_com1,
                numMods_max0,
                numMods_max1,
                pnts_com0,
                pnts_com1,
                pnts_ref0,
                pnts_ref1,
                pnts_max0,
                pnts_max1;
                
                int fCnt = 0;
                
                //count the number of face on each partition
                
                for(i = 0; i < locexp.size(); ++i)
                {
                    for(j = 0; j < locexp[i]->GetNfaces(); ++j)
                    {
                        fCnt = fCnt + 1;
                    }
                }
                
                // set-up the offset and the array that will contain the list of
                // edge ID's
                
                Array<OneD, int> faceCnt(nproc,0);
                faceCnt[edgepr] = fCnt;
                
                vComm->AllReduce(faceCnt, LibUtilities::ReduceSum);
                
                int totFaceCnt = Vmath::Vsum(nproc, faceCnt, 1);
                
                Array<OneD, int> fTotOffsets(nproc,0);
                
                for (i = 1; i < nproc; ++i)
                {
                    fTotOffsets[i] = fTotOffsets[i-1] + faceCnt[i-1];
                }
                
                // Local list of the edges per element
                
                Array<OneD, int> faceID(fCnt,0);
                Array<OneD, int> faceNm0(fCnt,0);
                Array<OneD, int> faceNm1(fCnt,0);
                Array<OneD, int> facePnts0(fCnt,0);
                Array<OneD, int> facePnts1(fCnt,0);
                
                int cntr = 0;
                
                for(i = 0; i < locexp.size(); ++i)
                {
                    exp3D = LocalRegions::Expansion3D::FromStdExp(locexp[i]);
                    
                    int nfaces = locexp[i]->GetNfaces();
                    
                    for(j = 0; j < nfaces; ++j)
                    {
                        LibUtilities::BasisKey face_dir0
                        = locexp[i]->DetFaceBasisKey(j,0);
                        LibUtilities::BasisKey face_dir1
                        = locexp[i]->DetFaceBasisKey(j,1);
                        
                        int nm_ref0         = face_dir0.GetNumModes();
                        int nm_ref1         = face_dir1.GetNumModes();
                        
                        int pnts_ref0       = face_dir0.GetNumPoints();
                        int pnts_ref1       = face_dir1.GetNumPoints();
                        
                        FaceGeom            = exp3D->GetGeom3D()->GetFace(j);
                        
                        faceID[cntr]      = FaceGeom->GetFid();
                        faceNm0[cntr]      = nm_ref0;
                        faceNm1[cntr]      = nm_ref1;
                        facePnts0[cntr]    = pnts_ref0;
                        facePnts1[cntr]    = pnts_ref1;
                        
                        cntr = cntr + 1;
                    }
                }
                
                // Make a list of the edges per element. In this case double edge
                // ID's can be identified which are the trace edges
                
                Array<OneD, int> FacesTotID(totFaceCnt, 0);
                Array<OneD, int> FacesTotNm0(totFaceCnt, 0);
                Array<OneD, int> FacesTotNm1(totFaceCnt, 0);
                Array<OneD, int> FacesTotPnts0(totFaceCnt, 0);
                Array<OneD, int> FacesTotPnts1(totFaceCnt, 0);
                
                for (i = 0; i < fCnt; ++i)
                {
                    FacesTotID[fTotOffsets[edgepr] + i]   = faceID[i];
                    FacesTotNm0[fTotOffsets[edgepr] + i]   = faceNm0[i];
                    FacesTotNm1[fTotOffsets[edgepr] + i]   = faceNm1[i];
                    FacesTotPnts0[fTotOffsets[edgepr] + i] = facePnts0[i];
                    FacesTotPnts1[fTotOffsets[edgepr] + i] = facePnts1[i];
                }
                
                vComm->AllReduce(FacesTotID, LibUtilities::ReduceSum);
                vComm->AllReduce(FacesTotNm0, LibUtilities::ReduceSum);
                vComm->AllReduce(FacesTotNm1, LibUtilities::ReduceSum);
                vComm->AllReduce(FacesTotPnts0, LibUtilities::ReduceSum);
                vComm->AllReduce(FacesTotPnts1, LibUtilities::ReduceSum);
                
                Array<OneD, int> tmp_traceGlobalID(totFaceCnt, 0);
                Array<OneD, int> tmp_traceGlobalNm0(totFaceCnt, 0);
                Array<OneD, int> tmp_traceGlobalNm1(totFaceCnt, 0);
                Array<OneD, int> tmp_traceGlobalPnts0(totFaceCnt, 0);
                Array<OneD, int> tmp_traceGlobalPnts1(totFaceCnt, 0);
                
                int cnt = 0;
                
                for (int u = 0; u < totFaceCnt; ++u)
                {
                    for (int k = 0; k < totFaceCnt; ++k)
                    {
                        if (FacesTotID[k] == FacesTotID[u] && u != k)
                        {
                            numModes_ref0 = FacesTotNm0[u];
                            numModes_com0 = FacesTotNm0[k];
                            
                            numModes_ref1 = FacesTotNm1[u];
                            numModes_com1 = FacesTotNm1[k];
                            
                            pnts_ref0 = FacesTotPnts0[u];
                            pnts_com0 = FacesTotPnts0[k];
                            
                            pnts_ref1 = FacesTotPnts1[u];
                            pnts_com1 = FacesTotPnts1[k];
                            
                            // determine the maximum nummodes
                            
                            if (numModes_ref0 > numModes_com0
                                && numModes_ref1 > numModes_com1)
                            {
                                numMods_max0 = numModes_ref0;
                                numMods_max1 = numModes_ref1;
                                
                                pnts_max0    = pnts_ref0;
                                pnts_max1    = pnts_ref1;
                            }
                            else
                            {
                                numMods_max0 = numModes_com0;
                                numMods_max1 = numModes_com1;
                                
                                pnts_max0    = pnts_com0;
                                pnts_max1    = pnts_com1;
                            }
                            
                            tmp_traceGlobalID[cnt]       = FacesTotID[k];
                            tmp_traceGlobalNm0[cnt]      = numMods_max0;
                            tmp_traceGlobalNm1[cnt]      = numMods_max1;
                            tmp_traceGlobalPnts0[cnt]    = pnts_max0;
                            tmp_traceGlobalPnts1[cnt]    = pnts_max1;
                            
                            cnt++;
                        }
                    }
                }
                
                int nTraceGlobal = cnt;
                
                for (int i = 0; i < nTraceGlobal; ++i)
                {
                    for (int j = i+1; j < nTraceGlobal;)
                    {
                        if (tmp_traceGlobalID[j] == tmp_traceGlobalID[i])
                        {
                            for (int k = j; k < nTraceGlobal; ++k)
                            {
                                tmp_traceGlobalID[k]    = tmp_traceGlobalID[k+1];
                                tmp_traceGlobalNm0[k]   = tmp_traceGlobalNm0[k+1];
                                tmp_traceGlobalNm1[k]   = tmp_traceGlobalNm1[k+1];
                                tmp_traceGlobalPnts0[k] = tmp_traceGlobalPnts0[k+1];
                                tmp_traceGlobalPnts1[k] = tmp_traceGlobalPnts1[k+1];
                            }
                            nTraceGlobal--;
                        }
                        else
                        {
                            j++;
                        }
                    }
                }
                
                Array<OneD, int> traceGlobalID(nTraceGlobal, 0);
                Array<OneD, int> traceGlobalNm0(nTraceGlobal, 0);
                Array<OneD, int> traceGlobalNm1(nTraceGlobal, 0);
                Array<OneD, int> traceGlobalPnts0(nTraceGlobal, 0);
                Array<OneD, int> traceGlobalPnts1(nTraceGlobal, 0);
                
                Vmath::Vcopy(nTraceGlobal,
                             tmp_traceGlobalID, 1,
                             traceGlobalID, 1);
                
                Vmath::Vcopy(nTraceGlobal,
                             tmp_traceGlobalNm0, 1,
                             traceGlobalNm0, 1);
                
                Vmath::Vcopy(nTraceGlobal,
                             tmp_traceGlobalNm1, 1,
                             traceGlobalNm1, 1);
                
                Vmath::Vcopy(nTraceGlobal,
                             tmp_traceGlobalPnts0, 1,
                             traceGlobalPnts0, 1);
                
                Vmath::Vcopy(nTraceGlobal,
                             tmp_traceGlobalPnts1, 1,
                             traceGlobalPnts1, 1);
                
                // Loop to set up the correct trace expansions
                // locexp.size() gives the number of elements on this node
                
                for(i = 0; i < locexp.size(); ++i)
                {
                    exp3D = LocalRegions::Expansion3D::FromStdExp(locexp[i]);
                    
                    int nfaces = locexp[i]->GetNfaces();
                    
                    for(j = 0; j < nfaces; ++j)
                    {
                        FaceGeom = exp3D->GetGeom3D()->GetFace(j);
                        id       = FaceGeom->GetFid();
                        
                        LibUtilities::BasisKey face_dir0
                        = locexp[i]->DetFaceBasisKey(j,0);
                        LibUtilities::BasisKey face_dir1
                        = locexp[i]->DetFaceBasisKey(j,1);
                        
                        int nm_ref0 = face_dir0.GetNumModes();
                        int nm_ref1 = face_dir1.GetNumModes();
                        
                        if (facesDone.count(id) != 0)
                        {
                            continue;
                        }
                        
                        it = faceOrders.find(id);
                        
                        for (int u = 0; u < nTraceGlobal; ++u)
                        {
                            if(id == traceGlobalID[u]
                               && nm_ref0 != traceGlobalNm0[u]
                               && nm_ref1 != traceGlobalNm1[u])
                            {
                                // set-up the correct PointsKey for the newly
                                // defined basiskey for the incorrect defined trace edges
                                
                                const LibUtilities::PointsKey newPkey0(
                                        traceGlobalPnts0[u],
                                        face_dir0.GetPointsType());
                                
                                const LibUtilities::PointsKey newPkey1(
                                        traceGlobalPnts1[u],
                                        face_dir1.GetPointsType());
                                
                                // set-up the correct basiskey for the trace edges
                                LibUtilities::BasisKey traceGlobalBkey0(
                                        face_dir0.GetBasisType(),
                                        traceGlobalNm0[u],
                                        newPkey0);
                                
                                LibUtilities::BasisKey traceGlobalBkey1(
                                        face_dir1.GetBasisType(),
                                        traceGlobalNm1[u],
                                        newPkey1);
                                
                                faceOrders.insert(std::make_pair(id,
                                        std::make_pair(FaceGeom,std::make_pair(
                                        traceGlobalBkey0, traceGlobalBkey1))));
                                
                                it->second.second.first  = traceGlobalBkey0;
                                it->second.second.second = traceGlobalBkey1;
                            }
                        }
                    }
                }
            }

            for (it = faceOrders.begin(); it != faceOrders.end(); ++it)
            {
                FaceGeom = it->second.first;
                
                if ((FaceQuadGeom = boost::dynamic_pointer_cast<
                     SpatialDomains::QuadGeom>(FaceGeom)))
                {
                    FaceQuadExp = MemoryManager<LocalRegions::QuadExp>
                        ::AllocateSharedPtr(it->second.second.first,
                                            it->second.second.second,
                                            FaceQuadGeom);
                    FaceQuadExp->SetElmtId(elmtid++);
                    (*m_exp).push_back(FaceQuadExp);
                }
                else if ((FaceTriGeom = boost::dynamic_pointer_cast<
                          SpatialDomains::TriGeom>(FaceGeom)))
                {
                    FaceTriExp = MemoryManager<LocalRegions::TriExp>
                        ::AllocateSharedPtr(it->second.second.first,
                                            it->second.second.second,
                                            FaceTriGeom);
                    FaceTriExp->SetElmtId(elmtid++);
                    (*m_exp).push_back(FaceTriExp);
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

             SetExpType(e2D);

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
                             = boost::dynamic_pointer_cast<SpatialDomains::MeshGraph3D>(graph3D)->GetFaceBasisKey(TriangleGeom,0,variable);
                         LibUtilities::BasisKey TriBb
                             = boost::dynamic_pointer_cast<SpatialDomains::MeshGraph3D>(graph3D)->GetFaceBasisKey(TriangleGeom,1,variable);

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
                             = boost::dynamic_pointer_cast<SpatialDomains::MeshGraph3D>(graph3D)->GetFaceBasisKey(QuadrilateralGeom,0,variable);
                         LibUtilities::BasisKey QuadBb
                             = boost::dynamic_pointer_cast<SpatialDomains::MeshGraph3D>(graph3D)->GetFaceBasisKey(QuadrilateralGeom,1,variable);

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
            int i,j,k,offset;
            Array<OneD,Array<OneD,NekDouble> > locnormals;
            Array<OneD, NekDouble> tmp;
            // Assume whole array is of same coordinate dimension
            int coordim = GetCoordim(0);
            
            ASSERTL1(normals.num_elements() >= coordim,
                     "Output vector does not have sufficient dimensions to "
                     "match coordim");
            
            // Process each expansion.
            if(m_parallel == false)
            {
                for(i = 0; i < m_exp->size(); ++i)
                {
                    int e_npoints = (*m_exp)[i]->GetTotPoints();
                
                    LocalRegions::Expansion2DSharedPtr loc_exp =
                        boost::dynamic_pointer_cast<
                            LocalRegions::Expansion2D>((*m_exp)[i]);
                    LocalRegions::Expansion3DSharedPtr loc_elmt = 
                        loc_exp->GetLeftAdjacentElementExp();
                    int faceNumber = loc_exp->GetLeftAdjacentElementFace();
                
                    // Get the number of points and normals for this expansion.
                    locnormals = loc_elmt->GetFaceNormal(faceNumber);
                    int e_nmodes   = loc_exp->GetBasis(0)->GetNumModes();
                    int loc_nmodes = loc_elmt->GetBasis(0)->GetNumModes();
                    //if (e_npoints != locnormals[0].num_elements())
                    if (e_nmodes != loc_nmodes)
                    {
                        LocalRegions::Expansion3DSharedPtr loc_elmt =
                            loc_exp->GetRightAdjacentElementExp();
                        int faceNumber = loc_exp->GetRightAdjacentElementFace();
                        
                        // Get the number of points and normals for this expansion.
                        locnormals = loc_elmt->GetFaceNormal(faceNumber);

                        offset = m_phys_offset[i];
                        // Process each point in the expansion.
                        for(int j = 0; j < e_npoints; ++j)
                        {
                            for(k = 0; k < coordim; ++k)
                            {
                                normals[k][offset+j] = -locnormals[k][j];
                            }
                        }
                    }
                    else
                    {
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
            }
            if(m_parallel == true)
            {
                for(i = 0; i < m_exp->size(); ++i)
                {
                    int e_npoints = (*m_exp)[i]->GetTotPoints();
                    
                    LocalRegions::Expansion2DSharedPtr loc_exp =
                    boost::dynamic_pointer_cast<
                    LocalRegions::Expansion2D>((*m_exp)[i]);
                    
                    LocalRegions::Expansion3DSharedPtr loc_elmt =
                    loc_exp->GetLeftAdjacentElementExp();
                    
                    LocalRegions::Expansion2DSharedPtr to_exp =
                    boost::dynamic_pointer_cast<
                    LocalRegions::Expansion2D>(loc_elmt->GetFaceExp(
                                loc_exp->GetLeftAdjacentElementFace()));
                    
                    int faceNumber = loc_exp->GetLeftAdjacentElementFace();
                    
                    int e_nmodes   = loc_exp->GetBasis(0)->GetNumModes();
                    int loc_nmodes = loc_elmt->GetBasis(0)->GetNumModes();
                    
                    // Get the number of points and normals for this expansion.
                    locnormals = loc_elmt->GetFaceNormal(faceNumber);
                    
                    if (e_nmodes != loc_nmodes)
                    {
                        Array<OneD, Array<OneD, NekDouble> > normal(coordim);
                        
                        for (int p = 0; p < coordim; ++p)
                        {
                            normal[p] = Array<OneD, NekDouble>(e_npoints,0.0);
                            
                            LibUtilities::PointsKey to_key0 =
                            loc_exp->GetBasis(0)->GetPointsKey();
                            LibUtilities::PointsKey to_key1 =
                            loc_exp->GetBasis(1)->GetPointsKey();
                            
                            LibUtilities::PointsKey from_key0 =
                            loc_elmt->GetBasis(0)->GetPointsKey();
                            LibUtilities::PointsKey from_key1 =
                            loc_elmt->GetBasis(1)->GetPointsKey();
                            
                            LibUtilities::Interp2D(from_key0,
                                                   from_key1,
                                                   locnormals[p],
                                                   to_key0,
                                                   to_key1,
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
                    else
                    {
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
                if((*m_exp)[i]->DetShapeType() == LibUtilities::eTriangle)
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
                if((*m_exp)[i]->DetShapeType() == LibUtilities::eQuadrilateral)
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
                if((*m_exp)[i]->DetShapeType() == LibUtilities::eTriangle)
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
