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

#include <boost/core/ignore_unused.hpp>

#include <LocalRegions/TriExp.h>
#include <LocalRegions/QuadExp.h>
#include <LocalRegions/NodalTriExp.h>
#include <LocalRegions/Expansion3D.h>
#include <MultiRegions/ExpList2D.h>
#include <LibUtilities/Foundations/Interp.h>
#include <LibUtilities/Foundations/PhysGalerkinProject.h>

using namespace std;

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
         * @param   In   ExpList2D object to copy.
         */
        ExpList2D::ExpList2D(
            const ExpList2D &In,
            const bool DeclareCoeffPhysArrays):
            ExpList(In,DeclareCoeffPhysArrays)
        {
            SetExpType(e2D);
        }

        /**
         * @param   In   ExpList2D object to copy.
         * @param   eIDs Id of elements that should be copied.
         */
        ExpList2D::ExpList2D(
            const ExpList2D &In,
            const std::vector<unsigned int> &eIDs,
            const bool DeclareCoeffPhysArrays,
            const Collections::ImplementationType ImpType):
            ExpList(In,eIDs,DeclareCoeffPhysArrays)
        {
            SetExpType(e2D);

            // set up offset arrays.
            SetCoeffPhysOffsets();

            CreateCollections(ImpType);
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
        ExpList2D::ExpList2D(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const SpatialDomains::MeshGraphSharedPtr &graph2D,
                const bool DeclareCoeffPhysArrays,
                const std::string &var,
                const Collections::ImplementationType ImpType):
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
                    m_ncoeffs += (TriBa.GetNumModes()*(TriBa.GetNumModes()+1))/2
                                    + TriBa.GetNumModes()*(TriBb.GetNumModes()
                                    -TriBa.GetNumModes());
                    m_npoints += TriBa.GetNumPoints()*TriBb.GetNumPoints();
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

                    m_ncoeffs += QuadBa.GetNumModes()*QuadBb.GetNumModes();
                    m_npoints += QuadBa.GetNumPoints()*QuadBb.GetNumPoints();
                }
                else
                {
                    ASSERTL0(false, "dynamic cast to a proper Geometry2D "
                                    "failed");
                }

            }

            // set up element numbering
            for(int i = 0; i < (*m_exp).size(); ++i)
            {
                (*m_exp)[i]->SetElmtId(i);
            }

            // set up offset arrays.
            SetCoeffPhysOffsets();

            if (DeclareCoeffPhysArrays)
            {
                // Set up m_coeffs, m_phys.
                m_coeffs = Array<OneD, NekDouble> {size_t(m_ncoeffs), 0.0};
                m_phys   = Array<OneD, NekDouble> {size_t(m_npoints), 0.0};
             }

            CreateCollections(ImpType);
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
            const bool DeclareCoeffPhysArrays,
            const Collections::ImplementationType ImpType):
            ExpList(pSession)
        {
            SetExpType(e2D);

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
                    m_ncoeffs += (TriBa.GetNumModes()*(TriBa.GetNumModes()+1))/2
                                    + TriBa.GetNumModes()*(TriBb.GetNumModes()
                                    -TriBa.GetNumModes());
                    m_npoints += TriBa.GetNumPoints()*TriBb.GetNumPoints();
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

                    m_ncoeffs += QuadBa.GetNumModes()*QuadBb.GetNumModes();
                    m_npoints += QuadBa.GetNumPoints()*QuadBb.GetNumPoints();
                }
                else
                {
                    ASSERTL0(false, "dynamic cast to a proper Geometry2D "
                                    "failed");
                }

            }

            // set up offset arrays.
            SetCoeffPhysOffsets();

            if (DeclareCoeffPhysArrays)
            {
                // Set up m_coeffs, m_phys.
                m_coeffs = Array<OneD, NekDouble> {size_t(m_ncoeffs), 0.0};
                m_phys   = Array<OneD, NekDouble> {size_t(m_npoints), 0.0};
             }

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
        ExpList2D::ExpList2D(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const LibUtilities::BasisKey &TriBa,
            const LibUtilities::BasisKey &TriBb,
            const LibUtilities::BasisKey &QuadBa,
            const LibUtilities::BasisKey &QuadBb,
            const SpatialDomains::MeshGraphSharedPtr &graph2D,
            const LibUtilities::PointsType TriNb,
            const Collections::ImplementationType ImpType):
            ExpList(pSession, graph2D)
        {
            SetExpType(e2D);

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

            // Set up m_coeffs, m_phys and offset arrays.
            SetCoeffPhysOffsets();
            m_coeffs = Array<OneD, NekDouble> {size_t(m_ncoeffs), 0.0};
            m_phys   = Array<OneD, NekDouble> {size_t(m_npoints), 0.0};

            CreateCollections(ImpType);
        }

        /**
         * Specialized constructor for trace expansions. Store
         * expansions for the trace space used in DisContField3D
         *
         * @param   bndConstraint  Array of ExpList2D objects each containing a
         *                         2D spectral/hp element expansion on a single
         *                         boundary region.
         * @param   bndCond   Array of BoundaryCondition objects which contain
         *                    information about the boundary conditions on the
         *                    different boundary regions.
         * @param   locexp   Complete domain expansion list.
         * @param   graph3D   3D mesh corresponding to the expansion list.
         * @param   periodicFaces   List of periodic faces.
         * @param   DeclareCoeffPhysArrays   If true, set up m_coeffs,
         *                                   m_phys arrays
         **/
        ExpList2D::ExpList2D(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const Array<OneD,const ExpListSharedPtr> &bndConstraint,
            const Array<OneD,const SpatialDomains::BoundaryConditionShPtr>  &bndCond,
            const LocalRegions::ExpansionVector &locexp,
            const SpatialDomains::MeshGraphSharedPtr &graph3D,
            const PeriodicMap &periodicFaces,
            const bool DeclareCoeffPhysArrays,
            const std::string variable,
            const Collections::ImplementationType ImpType):
            ExpList(pSession, graph3D)
        {
            boost::ignore_unused(periodicFaces, variable);

            SetExpType(e2D);

            int i, j, id, elmtid=0;
            set<int> facesDone;

            SpatialDomains::Geometry2DSharedPtr FaceGeom;
            SpatialDomains::QuadGeomSharedPtr   FaceQuadGeom;
            SpatialDomains::TriGeomSharedPtr    FaceTriGeom;
            LocalRegions::QuadExpSharedPtr      FaceQuadExp;
            LocalRegions::TriExpSharedPtr       FaceTriExp;
            LocalRegions::Expansion2DSharedPtr  exp2D;
            LocalRegions::Expansion3DSharedPtr  exp3D;

            // First loop over boundary conditions to renumber
            // Dirichlet boundaries
            for (i = 0; i < bndCond.size(); ++i)
            {
                if (bndCond[i]->GetBoundaryConditionType()
                                            == SpatialDomains::eDirichlet)
                {
                    for (j = 0; j < bndConstraint[i]->GetExpSize(); ++j)
                    {
                        LibUtilities::BasisKey bkey0 = bndConstraint[i]
                            ->GetExp(j)->GetBasis(0)->GetBasisKey();
                        LibUtilities::BasisKey bkey1 = bndConstraint[i]
                            ->GetExp(j)->GetBasis(1)->GetBasisKey();
                        exp2D = bndConstraint[i]->GetExp(j)
                                    ->as<LocalRegions::Expansion2D>();
                        FaceGeom = exp2D->GetGeom2D();

                        //if face is a quad
                        if((FaceQuadGeom = std::dynamic_pointer_cast<
                            SpatialDomains::QuadGeom>(FaceGeom)))
                        {
                            FaceQuadExp = MemoryManager<LocalRegions::QuadExp>
                                ::AllocateSharedPtr(bkey0, bkey1, FaceQuadGeom);
                            facesDone.insert(FaceQuadGeom->GetGlobalID());
                            FaceQuadExp->SetElmtId(elmtid++);
                            (*m_exp).push_back(FaceQuadExp);
                        }
                        //if face is a triangle
                        else if((FaceTriGeom = std::dynamic_pointer_cast<
                                 SpatialDomains::TriGeom>(FaceGeom)))
                        {
                            FaceTriExp = MemoryManager<LocalRegions::TriExp>
                                ::AllocateSharedPtr(bkey0, bkey1, FaceTriGeom);
                            facesDone.insert(FaceTriGeom->GetGlobalID());
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

            for(i = 0; i < locexp.size(); ++i)
            {
                exp3D = locexp[i]->as<LocalRegions::Expansion3D>();
                for (j = 0; j < exp3D->GetNfaces(); ++j)
                {
                    FaceGeom = exp3D->GetGeom3D()->GetFace(j);
                    id       = FaceGeom->GetGlobalID();

                    if(facesDone.count(id) != 0)
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

            LibUtilities::CommSharedPtr vComm = pSession->GetComm();
            int nproc = vComm->GetSize(); // number of processors
            int facepr = vComm->GetRank(); // ID processor

            if (nproc > 1)
            {
                int fCnt = 0;

                // Count the number of faces on each partition
                for(i = 0; i < locexp.size(); ++i)
                {
                    fCnt += locexp[i]->GetNfaces();
                }

                // Set up the offset and the array that will contain the list of
                // face IDs, then reduce this across processors.
                Array<OneD, int> faceCnt(nproc,0);
                faceCnt[facepr] = fCnt;
                vComm->AllReduce(faceCnt, LibUtilities::ReduceSum);

                int totFaceCnt = Vmath::Vsum(nproc, faceCnt, 1);
                Array<OneD, int> fTotOffsets(nproc,0);

                for (i = 1; i < nproc; ++i)
                {
                    fTotOffsets[i] = fTotOffsets[i-1] + faceCnt[i-1];
                }

                // Local list of the edges per element

                Array<OneD, int> FacesTotID   (totFaceCnt, 0);
                Array<OneD, int> FacesTotNm0  (totFaceCnt, 0);
                Array<OneD, int> FacesTotNm1  (totFaceCnt, 0);
                Array<OneD, int> FacesTotPnts0(totFaceCnt, 0);
                Array<OneD, int> FacesTotPnts1(totFaceCnt, 0);

                int cntr = fTotOffsets[facepr];

                for(i = 0; i < locexp.size(); ++i)
                {
                    exp3D = locexp[i]->as<LocalRegions::Expansion3D>();

                    int nfaces = locexp[i]->GetNfaces();

                    for(j = 0; j < nfaces; ++j, ++cntr)
                    {
                        LibUtilities::BasisKey face_dir0
                            = locexp[i]->DetFaceBasisKey(j,0);
                        LibUtilities::BasisKey face_dir1
                            = locexp[i]->DetFaceBasisKey(j,1);

                        FacesTotID[cntr]    = exp3D->GetGeom3D()->GetFid(j);
                        FacesTotNm0[cntr]   = face_dir0.GetNumModes ();
                        FacesTotNm1[cntr]   = face_dir1.GetNumModes ();
                        FacesTotPnts0[cntr] = face_dir0.GetNumPoints();
                        FacesTotPnts1[cntr] = face_dir1.GetNumPoints();
                    }
                }

                vComm->AllReduce(FacesTotID,    LibUtilities::ReduceSum);
                vComm->AllReduce(FacesTotNm0,   LibUtilities::ReduceSum);
                vComm->AllReduce(FacesTotNm1,   LibUtilities::ReduceSum);
                vComm->AllReduce(FacesTotPnts0, LibUtilities::ReduceSum);
                vComm->AllReduce(FacesTotPnts1, LibUtilities::ReduceSum);

                for (i = 0; i < totFaceCnt; ++i)
                {
                    auto it = faceOrders.find(FacesTotID[i]);

                    if (it == faceOrders.end())
                    {
                        continue;
                    }

                    LibUtilities::BasisKey existing0 =
                        it->second.second.first;
                    LibUtilities::BasisKey existing1 =
                        it->second.second.second;
                    LibUtilities::BasisKey face0(
                        existing0.GetBasisType(), FacesTotNm0[i],
                        LibUtilities::PointsKey(FacesTotPnts0[i],
                                                existing0.GetPointsType()));
                    LibUtilities::BasisKey face1(
                        existing1.GetBasisType(), FacesTotNm1[i],
                        LibUtilities::PointsKey(FacesTotPnts1[i],
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
        ExpList2D::ExpList2D(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const SpatialDomains::CompositeMap &domain,
            const SpatialDomains::MeshGraphSharedPtr &graph3D,
            const std::string variable,
            const LibUtilities::CommSharedPtr comm,
            const Collections::ImplementationType ImpType)
            : ExpList(pSession, graph3D)
         {
             SetExpType(e2D);

             if (comm)
             {
                 m_comm = comm;
             }

             int j, elmtid=0;
             int nel = 0;

             SpatialDomains::Composite comp;
             SpatialDomains::TriGeomSharedPtr TriangleGeom;
             SpatialDomains::QuadGeomSharedPtr QuadrilateralGeom;

             LocalRegions::TriExpSharedPtr tri;
             LocalRegions::NodalTriExpSharedPtr Ntri;
             LibUtilities::PointsType TriNb;
             LocalRegions::QuadExpSharedPtr quad;

             for (auto &compIt : domain)
             {
                 nel += compIt.second->m_geomVec.size();
             }

             for (auto &compIt : domain)
             {
                 for (j = 0; j < compIt.second->m_geomVec.size(); ++j)
                 {
                     if ((TriangleGeom = std::dynamic_pointer_cast<
                             SpatialDomains::TriGeom>(compIt.second->m_geomVec[j])))
                     {
                         LibUtilities::BasisKey TriBa
                             = graph3D->GetFaceBasisKey(TriangleGeom, 0, variable);
                         LibUtilities::BasisKey TriBb
                             = graph3D->GetFaceBasisKey(TriangleGeom,1,variable);

                         if (graph3D->GetExpansions().begin()->second->
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

                         m_ncoeffs
                             += (TriBa.GetNumModes()*(TriBa.GetNumModes()+1))/2
                                 + TriBa.GetNumModes()*(TriBb.GetNumModes()
                                 -TriBa.GetNumModes());
                         m_npoints += TriBa.GetNumPoints()*TriBb.GetNumPoints();
                     }
                     else if ((QuadrilateralGeom = std::dynamic_pointer_cast<
                              SpatialDomains::QuadGeom>(compIt.second->m_geomVec[j])))
                     {
                         LibUtilities::BasisKey QuadBa
                             = graph3D->GetFaceBasisKey(QuadrilateralGeom, 0,
                                                        variable);
                         LibUtilities::BasisKey QuadBb
                             = graph3D->GetFaceBasisKey(QuadrilateralGeom, 1,
                                                        variable);

                         quad = MemoryManager<LocalRegions::QuadExp>
                             ::AllocateSharedPtr(QuadBa, QuadBb,
                                                 QuadrilateralGeom);
                         quad->SetElmtId(elmtid++);
                         (*m_exp).push_back(quad);

                         m_ncoeffs += QuadBa.GetNumModes()*QuadBb.GetNumModes();
                         m_npoints += QuadBa.GetNumPoints()
                                         * QuadBb.GetNumPoints();
                     }
                     else
                     {
                         ASSERTL0(false,
                                  "dynamic cast to a proper Geometry2D failed");
                     }
                 }

             }

             // Setup Default optimisation information.
             nel = GetExpSize();

             // Set up m_coeffs, m_phys and offset arrays.
            SetCoeffPhysOffsets();
            m_coeffs = Array<OneD, NekDouble> {size_t(m_ncoeffs), 0.0};
            m_phys   = Array<OneD, NekDouble> {size_t(m_npoints), 0.0};

            CreateCollections(ImpType);
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
         * @brief Helper function to re-align face to a given orientation.
         */
        void AlignFace(const StdRegions::Orientation       orient,
                       const int                           nquad1,
                       const int                           nquad2,
                       const Array<OneD, const NekDouble> &in,
                             Array<OneD,       NekDouble> &out)
        {
            // Copy transpose.
            if (orient == StdRegions::eDir1FwdDir2_Dir2FwdDir1 ||
                orient == StdRegions::eDir1BwdDir2_Dir2FwdDir1 ||
                orient == StdRegions::eDir1FwdDir2_Dir2BwdDir1 ||
                orient == StdRegions::eDir1BwdDir2_Dir2BwdDir1)
            {
                for (int i = 0; i < nquad2; ++i)
                {
                    for (int j = 0; j < nquad1; ++j)
                    {
                        out[i*nquad1 + j] = in[j*nquad2 + i];
                    }
                }
            }
            else
            {
                for (int i = 0; i < nquad2; ++i)
                {
                    for (int j = 0; j < nquad1; ++j)
                    {
                        out[i*nquad1 + j] = in[i*nquad1 + j];
                    }
                }
            }

            if (orient == StdRegions::eDir1BwdDir1_Dir2FwdDir2 ||
                orient == StdRegions::eDir1BwdDir1_Dir2BwdDir2 ||
                orient == StdRegions::eDir1BwdDir2_Dir2FwdDir1 ||
                orient == StdRegions::eDir1BwdDir2_Dir2BwdDir1)
            {
                // Reverse x direction
                for (int i = 0; i < nquad2; ++i)
                {
                    for (int j = 0; j < nquad1/2; ++j)
                    {
                        swap(out[i*nquad1 + j],
                             out[i*nquad1 + nquad1-j-1]);
                    }
                }
            }

            if (orient == StdRegions::eDir1FwdDir1_Dir2BwdDir2 ||
                orient == StdRegions::eDir1BwdDir1_Dir2BwdDir2 ||
                orient == StdRegions::eDir1FwdDir2_Dir2BwdDir1 ||
                orient == StdRegions::eDir1BwdDir2_Dir2BwdDir1)
            {
                // Reverse y direction
                for (int j = 0; j < nquad1; ++j)
                {
                    for (int i = 0; i < nquad2/2; ++i)
                    {
                        swap(out[i*nquad1 + j],
                             out[(nquad2-i-1)*nquad1 + j]);
                    }
                }
            }
        }

        /**
         * @brief For each local element, copy the normals stored in the element
         * list into the array \a normals.
         *
         * @param   normals     Multidimensional array in which to copy normals
         *                      to. Must have dimension equal to or larger than
         *                      the spatial dimension of the elements.
         */
        void ExpList2D::v_GetNormals(
            Array<OneD, Array<OneD, NekDouble> > &normals)
        {
            Array<OneD, NekDouble> tmp;
            int i, j;
            const int coordim = GetCoordim(0);

            ASSERTL1(normals.size() >= coordim,
                     "Output vector does not have sufficient dimensions to "
                     "match coordim");

            // Process each expansion.
            for (i = 0; i < m_exp->size(); ++i)
            {
                LocalRegions::Expansion2DSharedPtr traceExp = (*m_exp)[i]->as<
                    LocalRegions::Expansion2D>();
                LocalRegions::Expansion3DSharedPtr exp3D =
                    traceExp->GetLeftAdjacentElementExp();

                // Get the number of points and normals for this expansion.
                int faceNum = traceExp->GetLeftAdjacentElementFace();
                int offset  = m_phys_offset[i];

                const Array<OneD, const Array<OneD, NekDouble> > &locNormals
                    = exp3D->GetFaceNormal(faceNum);

                // Project normals from 3D element onto the same orientation as
                // the trace expansion.
                StdRegions::Orientation orient = exp3D->GetForient(faceNum);


                int fromid0,fromid1;

                if(orient < StdRegions::eDir1FwdDir2_Dir2FwdDir1)
                {
                    fromid0 = 0;
                    fromid1 = 1;
                }
                else
                {
                    fromid0 = 1;
                    fromid1 = 0;
                }

                LibUtilities::BasisKey faceBasis0
                    = exp3D->DetFaceBasisKey(faceNum, fromid0);
                LibUtilities::BasisKey faceBasis1
                    = exp3D->DetFaceBasisKey(faceNum, fromid1);
                LibUtilities::BasisKey traceBasis0
                    = traceExp->GetBasis(0)->GetBasisKey();
                LibUtilities::BasisKey traceBasis1
                    = traceExp->GetBasis(1)->GetBasisKey();

                const int faceNq0 = faceBasis0.GetNumPoints();
                const int faceNq1 = faceBasis1.GetNumPoints();

                for (j = 0; j < coordim; ++j)
                {
                    Array<OneD, NekDouble> traceNormals(faceNq0 * faceNq1);
                    AlignFace(orient, faceNq0, faceNq1,
                              locNormals[j], traceNormals);
                    LibUtilities::Interp2D(
                        faceBasis0.GetPointsKey(),
                        faceBasis1.GetPointsKey(),
                        traceNormals,
                        traceBasis0.GetPointsKey(),
                        traceBasis1.GetPointsKey(),
                        tmp = normals[j]+offset);
                }
            }
        }

        /**
         * For each local element, copy the element length in boundary normal 
         * direction stored in the element list into the array \a normals.
         */
        void ExpList2D::v_GetElmtNormalLength(
            Array<OneD, NekDouble>  &lengthsFwd,
            Array<OneD, NekDouble>  &lengthsBwd)
        {
            int e_npoints;

            Array<OneD, NekDouble> locLeng;
            Array<OneD, NekDouble> lengintp;
            Array<OneD, int      > LRbndnumbs(2);
            Array<OneD, Array<OneD,NekDouble> > lengLR(2);
            lengLR[0]   =   lengthsFwd;
            lengLR[1]   =   lengthsBwd;
            Array<OneD, LocalRegions::Expansion3DSharedPtr> LRelmts(2);
            LocalRegions::Expansion3DSharedPtr loc_elmt;
            LocalRegions::Expansion2DSharedPtr loc_exp;
            int e_npoints0  =   -1; 
            for (int i = 0; i < m_exp->size(); ++i)
            {
                loc_exp = (*m_exp)[i]->as<LocalRegions::Expansion2D>();
                int offset = m_phys_offset[i];

                LibUtilities::BasisKey traceBasis0
                    = loc_exp->GetBasis(0)->GetBasisKey();
                LibUtilities::BasisKey traceBasis1
                    = loc_exp->GetBasis(1)->GetBasisKey();
                const int TraceNq0 = traceBasis0.GetNumPoints();
                const int TraceNq1 = traceBasis1.GetNumPoints();
                e_npoints  =   TraceNq0*TraceNq1;
                if (e_npoints0 < e_npoints)
                {
                    lengintp = Array<OneD, NekDouble>{size_t(e_npoints), 0.0};
                    e_npoints0 = e_npoints;
                }
                
                LRelmts[0] = loc_exp->GetLeftAdjacentElementExp();
                LRelmts[1] = loc_exp->GetRightAdjacentElementExp();

                LRbndnumbs[0]   =   loc_exp->GetLeftAdjacentElementFace();
                LRbndnumbs[1]   =   loc_exp->GetRightAdjacentElementFace();
                for (int nlr = 0; nlr < 2; ++nlr)
                {
                    Vmath::Zero(e_npoints0, lengintp, 1);
                    int bndNumber = LRbndnumbs[nlr];
                    loc_elmt = LRelmts[nlr];
                    if (bndNumber >= 0)
                    {
                        locLeng = loc_elmt->GetElmtBndNormDirElmtLen(
                                                bndNumber);
                        // Project normals from 3D element onto the same orientation as
                        // the trace expansion.
                        StdRegions::Orientation orient = loc_elmt->GetForient(
                                                            bndNumber);

                        int fromid0,fromid1;
                        if (orient < StdRegions::eDir1FwdDir2_Dir2FwdDir1)
                        {
                            fromid0 = 0;
                            fromid1 = 1;
                        }
                        else
                        {
                            fromid0 = 1;
                            fromid1 = 0;
                        }

                        LibUtilities::BasisKey faceBasis0 
                            = loc_elmt->DetFaceBasisKey(bndNumber, fromid0);
                        LibUtilities::BasisKey faceBasis1 
                            = loc_elmt->DetFaceBasisKey(bndNumber, fromid1);
                        const int faceNq0 = faceBasis0.GetNumPoints();
                        const int faceNq1 = faceBasis1.GetNumPoints();
                        Array<OneD, NekDouble> alignedLeng(faceNq0 * faceNq1);
                        
                        AlignFace(orient, faceNq0, faceNq1,
                                locLeng, alignedLeng);
                        LibUtilities::Interp2D(
                            faceBasis0.GetPointsKey(),
                            faceBasis1.GetPointsKey(),
                            alignedLeng,
                            traceBasis0.GetPointsKey(),
                            traceBasis1.GetPointsKey(),
                            lengintp);
                    }
                    for (int j = 0; j < e_npoints; ++j)
                    {
                        lengLR[nlr][offset + j] = lengintp[j];
                    }
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


        void ExpList2D::v_WriteVtkPieceHeader(
            std::ostream &outfile,
            int expansion,
            int istrip)
        {
            boost::ignore_unused(istrip);

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


        void ExpList2D::v_PhysInterp1DScaled(
            const NekDouble scale,
            const Array<OneD, NekDouble> &inarray,
                  Array<OneD, NekDouble> &outarray)
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

                LibUtilities::PointsKey newPointsKey0(npt0,
                                                (*m_exp)[i]->GetPointsType(0));
                LibUtilities::PointsKey newPointsKey1(npt1,
                                                (*m_exp)[i]->GetPointsType(1));

                // Interpolate points;
                LibUtilities::Interp2D((*m_exp)[i]->GetBasis(0)->GetPointsKey(),
                                       (*m_exp)[i]->GetBasis(1)->GetPointsKey(),
                                       &inarray[cnt],newPointsKey0,
                                       newPointsKey1,&outarray[cnt1]);

                cnt  += pt0*pt1;
                cnt1 += npt0*npt1;
            }
        }

        void ExpList2D::v_PhysGalerkinProjection1DScaled(
            const NekDouble scale,
            const Array<OneD, NekDouble> &inarray,
                  Array<OneD, NekDouble> &outarray)
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

                LibUtilities::PointsKey newPointsKey0(npt0,
                                                (*m_exp)[i]->GetPointsType(0));
                LibUtilities::PointsKey newPointsKey1(npt1,
                                                (*m_exp)[i]->GetPointsType(1));

                // Project points;
                LibUtilities::PhysGalerkinProject2D(
                                        newPointsKey0,
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
