 //////////////////////////////////////////////////////////////////////////////
 //
 // File DisContField3D.cpp
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
 // Description: Field definition for 3D domain with boundary
 // conditions using LDG flux
 //
 ///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/DisContField3D.h>
#include <LocalRegions/Expansion3D.h>
#include <LocalRegions/Expansion2D.h>
#include <SpatialDomains/MeshGraph.h>
#include <LocalRegions/HexExp.h>
#include <LocalRegions/TetExp.h>
#include <LocalRegions/PrismExp.h>
#include <LibUtilities/Foundations/Interp.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <tuple>

using namespace std;

 namespace Nektar
 {
     namespace MultiRegions
     {
         /**
          * @class DisContField3D
          * Abstraction of a global discontinuous three-dimensional spectral/hp
          * element expansion which approximates the solution of a set of
          * partial differential equations.
          */

         /**
          * @brief Default constructor.
          */
         DisContField3D::DisContField3D() :
             DisContField          ()
         {
         }

         /**
          * @brief Constructs a global discontinuous field based on an input
          * mesh with boundary conditions.
          */
         DisContField3D::DisContField3D(
             const LibUtilities::SessionReaderSharedPtr &pSession,
             const SpatialDomains::MeshGraphSharedPtr   &graph3D,
             const std::string                          &variable,
             const bool                                  SetUpJustDG,
             const bool                                  DeclareCoeffPhysArrays, 
             const Collections::ImplementationType       ImpType):
             DisContField       (pSession, graph3D, variable, SetUpJustDG,
                                 DeclareCoeffPhysArrays, ImpType)
         {
         }

         /*
          * @brief Copy type constructor which declares new boundary conditions
          * and re-uses mapping info and trace space if possible
          */
         DisContField3D::DisContField3D( 
             const DisContField3D                     &In,
             const SpatialDomains::MeshGraphSharedPtr &graph3D,
             const std::string                        &variable,
             const bool                                SetUpJustDG, 
             const bool                                DeclareCoeffPhysArrays)
             : DisContField(In,graph3D,variable,SetUpJustDG,DeclareCoeffPhysArrays)
         {
         }

         /**
          *
          */
         DisContField3D::DisContField3D(const DisContField3D &In) :
             DisContField(In)
         {
         }

         /**
          * @brief Destructor.
          */
         DisContField3D::~DisContField3D()
         {
         }




        /**
         * @brief Add trace contributions into elemental coefficient spaces.
         * 
         * Given some quantity \f$ \vec{Fn} \f$, which conatins this
         * routine calculates the integral
         * 
         * \f[ 
         * \int_{\Omega^e} \vec{Fn}, \mathrm{d}S
         * \f] 
         * 
         * and adds this to the coefficient space provided by
         * outarray. The value of q is determined from the routine
         * IsLeftAdjacentFace() which if true we use Fwd else we use
         * Bwd
         *
         * @see Expansion3D::AddFaceNormBoundaryInt
         * 
         * @param Fwd       The trace quantities associated with left (fwd)
         *                  adjancent elmt.
         * @param Bwd       The trace quantities associated with right (bwd)
         *                  adjacent elet.
         * @param outarray  Resulting 3D coefficient space.
         */
        void DisContField3D::v_AddFwdBwdTraceIntegral(
            const Array<OneD, const NekDouble> &Fwd, 
            const Array<OneD, const NekDouble> &Bwd, 
                  Array<OneD,       NekDouble> &outarray)
        {
            Array<OneD, NekDouble> Coeffs(m_trace->GetNcoeffs());

            m_trace->IProductWRTBase(Fwd,Coeffs);
            m_locTraceToTraceMap->AddTraceCoeffsToFieldCoeffs(0,Coeffs,outarray);
            m_trace->IProductWRTBase(Bwd,Coeffs);
            m_locTraceToTraceMap->AddTraceCoeffsToFieldCoeffs(1,Coeffs,outarray);
        }

        /**
         * @brief Calculates the result of the multiplication of a global matrix
         * of type specified by @a mkey with a vector given by @a inarray.
         * 
         * @param mkey      Key representing desired matrix multiplication.
         * @param inarray   Input vector.
         * @param outarray  Resulting multiplication.
         */
        void DisContField3D::v_GeneralMatrixOp(
               const GlobalMatrixKey             &gkey,
               const Array<OneD,const NekDouble> &inarray,
               Array<OneD,      NekDouble> &outarray)
        {
            int     LocBndCoeffs = m_traceMap->GetNumLocalBndCoeffs();
            Array<OneD, NekDouble> loc_lambda(LocBndCoeffs);
            DNekVec LocLambda(LocBndCoeffs,loc_lambda,eWrapper);
            const DNekScalBlkMatSharedPtr& HDGHelm = GetBlockMatrix(gkey);

            m_traceMap->GlobalToLocalBnd(inarray, loc_lambda);
            LocLambda = (*HDGHelm) * LocLambda;
            m_traceMap->AssembleBnd(loc_lambda,outarray);
        }
        
        /**
         * @brief Evaluate HDG post-processing to increase polynomial order of
         * solution.
         * 
         * This function takes the solution (assumed to be one order lower) in
         * physical space, and postprocesses at the current polynomial order by
         * solving the system:
         * 
         * \f[
         * \begin{aligned}
         *   (\nabla w, \nabla u^*) &= (\nabla w, u), \\
         *   \langle \nabla u^*, 1 \rangle &= \langle \nabla u, 1 \rangle
         * \end{aligned}
         * \f]
         * 
         * where \f$ u \f$ corresponds with the current solution as stored
         * inside #m_coeffs.
         * 
         * @param outarray  The resulting field \f$ u^* \f$.
         */
        void  DisContField3D::EvaluateHDGPostProcessing(
            Array<OneD, NekDouble> &outarray)
        {
            int    i,cnt,f,ncoeff_face;
            Array<OneD, NekDouble> force, out_tmp,qrhs,qrhs1;
            Array<OneD, Array< OneD, LocalRegions::ExpansionSharedPtr> > 
                &elmtToTrace = m_traceMap->GetElmtToTrace();
            
            int     nq_elmt, nm_elmt;
            int     LocBndCoeffs = m_traceMap->GetNumLocalBndCoeffs();
            Array<OneD, NekDouble> loc_lambda(LocBndCoeffs), face_lambda;
            Array<OneD, NekDouble> tmp_coeffs;
            m_traceMap->GlobalToLocalBnd(m_trace->GetCoeffs(),loc_lambda);

            face_lambda = loc_lambda;

            // Calculate Q using standard DG formulation.
            for(i = cnt = 0; i < GetExpSize(); ++i)
            {
                LocalRegions::Expansion3DSharedPtr exp =
                        (*m_exp)[i]->as<LocalRegions::Expansion3D>();

                nq_elmt = (*m_exp)[i]->GetTotPoints();
                nm_elmt = (*m_exp)[i]->GetNcoeffs();
                qrhs    = Array<OneD, NekDouble>(nq_elmt);
                qrhs1   = Array<OneD, NekDouble>(nq_elmt);
                force   = Array<OneD, NekDouble>(2*nm_elmt);
                out_tmp = force + nm_elmt;
                LocalRegions::ExpansionSharedPtr ppExp;

                int num_points0 = (*m_exp)[i]->GetBasis(0)->GetNumPoints();
                int num_points1 = (*m_exp)[i]->GetBasis(1)->GetNumPoints();
                int num_points2 = (*m_exp)[i]->GetBasis(2)->GetNumPoints();
                int num_modes0 = (*m_exp)[i]->GetBasis(0)->GetNumModes();
                int num_modes1 = (*m_exp)[i]->GetBasis(1)->GetNumModes();
                int num_modes2 = (*m_exp)[i]->GetBasis(2)->GetNumModes();

                // Probably a better way of setting up lambda than this.  Note
                // cannot use PutCoeffsInToElmts since lambda space is mapped
                // during the solve.
                int nFaces = (*m_exp)[i]->GetNfaces();
                Array<OneD, Array<OneD, NekDouble> > faceCoeffs(nFaces);
                for(f = 0; f < nFaces; ++f)
                {
                    ncoeff_face = elmtToTrace[i][f]->GetNcoeffs();
                    faceCoeffs[f] = Array<OneD, NekDouble>(ncoeff_face);
                    Vmath::Vcopy(ncoeff_face, face_lambda, 1, faceCoeffs[f], 1);
                    exp->SetFaceToGeomOrientation(f, faceCoeffs[f]);
                    face_lambda = face_lambda + ncoeff_face;
                }

                //creating orthogonal expansion (checking if we have quads or triangles)
                LibUtilities::ShapeType shape = (*m_exp)[i]->DetShapeType();
                switch(shape)
                {
                    case LibUtilities::eHexahedron:
                    {
                        const LibUtilities::PointsKey PkeyH1(num_points0,LibUtilities::eGaussLobattoLegendre);
                        const LibUtilities::PointsKey PkeyH2(num_points1,LibUtilities::eGaussLobattoLegendre);
                        const LibUtilities::PointsKey PkeyH3(num_points2,LibUtilities::eGaussLobattoLegendre);
                        LibUtilities::BasisKey  BkeyH1(LibUtilities::eOrtho_A, num_modes0, PkeyH1);
                        LibUtilities::BasisKey  BkeyH2(LibUtilities::eOrtho_A, num_modes1, PkeyH2);
                        LibUtilities::BasisKey  BkeyH3(LibUtilities::eOrtho_A, num_modes2, PkeyH3);
                        SpatialDomains::HexGeomSharedPtr hGeom = std::dynamic_pointer_cast<SpatialDomains::HexGeom>((*m_exp)[i]->GetGeom());
                        ppExp = MemoryManager<LocalRegions::HexExp>::AllocateSharedPtr(BkeyH1, BkeyH2, BkeyH3, hGeom);
                    }
                    break;
                    case LibUtilities::eTetrahedron:
                    {
                        const LibUtilities::PointsKey PkeyT1(num_points0,LibUtilities::eGaussLobattoLegendre);
                        const LibUtilities::PointsKey PkeyT2(num_points1,LibUtilities::eGaussRadauMAlpha1Beta0);
                        const LibUtilities::PointsKey PkeyT3(num_points2,LibUtilities::eGaussRadauMAlpha2Beta0);
                        LibUtilities::BasisKey  BkeyT1(LibUtilities::eOrtho_A, num_modes0, PkeyT1);
                        LibUtilities::BasisKey  BkeyT2(LibUtilities::eOrtho_B, num_modes1, PkeyT2);
                        LibUtilities::BasisKey  BkeyT3(LibUtilities::eOrtho_C, num_modes2, PkeyT3);
                        SpatialDomains::TetGeomSharedPtr tGeom = std::dynamic_pointer_cast<SpatialDomains::TetGeom>((*m_exp)[i]->GetGeom());
                        ppExp = MemoryManager<LocalRegions::TetExp>::AllocateSharedPtr(BkeyT1, BkeyT2, BkeyT3, tGeom);
                    }
                    break;
                    case LibUtilities::ePrism:
                    {
                        const LibUtilities::PointsKey PkeyP1(num_points0,LibUtilities::eGaussLobattoLegendre);
                        const LibUtilities::PointsKey PkeyP2(num_points1,LibUtilities::eGaussLobattoLegendre);
                        const LibUtilities::PointsKey PkeyP3(num_points2,LibUtilities::eGaussRadauMAlpha1Beta0);
                        LibUtilities::BasisKey  BkeyP1(LibUtilities::eOrtho_A, num_modes0, PkeyP1);
                        LibUtilities::BasisKey  BkeyP2(LibUtilities::eOrtho_A, num_modes1, PkeyP2);
                        LibUtilities::BasisKey  BkeyP3(LibUtilities::eOrtho_B, num_modes2, PkeyP3);
                        SpatialDomains::PrismGeomSharedPtr pGeom = std::dynamic_pointer_cast<SpatialDomains::PrismGeom>((*m_exp)[i]->GetGeom());
                        ppExp = MemoryManager<LocalRegions::PrismExp>::AllocateSharedPtr(BkeyP1, BkeyP2, BkeyP3, pGeom);
                    }
                    break;
                    default:
                        ASSERTL0(false, "Wrong shape type, HDG postprocessing is not implemented");
                };


                //DGDeriv	
                // (d/dx w, q_0)
                (*m_exp)[i]->DGDeriv(
                    0,tmp_coeffs = m_coeffs + m_coeff_offset[i],
                    elmtToTrace[i], faceCoeffs, out_tmp);
                (*m_exp)[i]->BwdTrans(out_tmp,qrhs);
                ppExp->IProductWRTDerivBase(0,qrhs,force);


                // + (d/dy w, q_1)
                (*m_exp)[i]->DGDeriv(
                    1,tmp_coeffs = m_coeffs + m_coeff_offset[i],
                    elmtToTrace[i], faceCoeffs, out_tmp);
                (*m_exp)[i]->BwdTrans(out_tmp,qrhs);
                ppExp->IProductWRTDerivBase(1,qrhs,out_tmp);

                Vmath::Vadd(nm_elmt,force,1,out_tmp,1,force,1);

                // + (d/dz w, q_2)
                (*m_exp)[i]->DGDeriv(
                    2,tmp_coeffs = m_coeffs + m_coeff_offset[i],
                    elmtToTrace[i], faceCoeffs, out_tmp);
                (*m_exp)[i]->BwdTrans(out_tmp,qrhs);
                ppExp->IProductWRTDerivBase(2,qrhs,out_tmp);

                Vmath::Vadd(nm_elmt,force,1,out_tmp,1,force,1);
                // determine force[0] = (1,u)
                (*m_exp)[i]->BwdTrans(
                    tmp_coeffs = m_coeffs + m_coeff_offset[i],qrhs);
                force[0] = (*m_exp)[i]->Integral(qrhs);

                // multiply by inverse Laplacian matrix
                // get matrix inverse
                LocalRegions::MatrixKey  lapkey(StdRegions::eInvLaplacianWithUnityMean, ppExp->DetShapeType(), *ppExp);
                DNekScalMatSharedPtr lapsys = ppExp->GetLocMatrix(lapkey); 

                NekVector<NekDouble> in (nm_elmt, force, eWrapper);
                NekVector<NekDouble> out(nm_elmt);

                out = (*lapsys)*in;

                // Transforming back to modified basis
                Array<OneD, NekDouble> work(nq_elmt);
                ppExp->BwdTrans(out.GetPtr(), work);
                (*m_exp)[i]->FwdTrans(work,
                                tmp_coeffs = outarray + m_coeff_offset[i]);
            }
        }

    } // end of namespace
} // end of namespace
