//////////////////////////////////////////////////////////////////////////////
//
// File DisContField2D.cpp
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
// Description: Field definition for 2D domain with boundary conditions using
// LDG flux.
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/DisContField2D.h>
#include <LocalRegions/MatrixKey.h>
#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/Expansion.h> 
#include <LocalRegions/QuadExp.h>   
#include <LocalRegions/TriExp.h>    
#include <SpatialDomains/MeshGraph.h>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>

using namespace std;

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class DisContField2D
         * Abstraction of a global discontinuous two-dimensional spectral/hp
         * element expansion which approximates the solution of a set of
         * partial differential equations.
         */
        
        /**
         * @brief Default constructor.
         */
        DisContField2D::DisContField2D(void)
            : DisContField            ()
        {
        }

        DisContField2D::DisContField2D(
            const DisContField2D &In, 
            const bool            DeclareCoeffPhysArrays)
            : DisContField         (In,DeclareCoeffPhysArrays)
        {
        }

        /**
         * @brief Constructs a global discontinuous field based on an input
         * mesh with boundary conditions.
         */
        DisContField2D::DisContField2D(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const SpatialDomains::MeshGraphSharedPtr   &graph2D,
            const std::string                          &variable,
            const bool                                  SetUpJustDG,
            const bool                                  DeclareCoeffPhysArrays,
            const Collections::ImplementationType       ImpType)
            : DisContField(pSession, graph2D, variable, SetUpJustDG,
                           DeclareCoeffPhysArrays,ImpType)
        {
        }


        /*
         * @brief Copy type constructor which declares new boundary conditions
         * and re-uses mapping info and trace space if possible
         */
        DisContField2D::DisContField2D(
            const DisContField2D                     &In,
            const SpatialDomains::MeshGraphSharedPtr &graph2D,
            const std::string                        &variable,
            const bool                                SetUpJustDG,
            const bool                                DeclareCoeffPhysArrays):
            DisContField(In,graph2D, variable, SetUpJustDG, DeclareCoeffPhysArrays)
        {
        }

        /**
         *
         */
        DisContField2D::~DisContField2D()
        {
        }
        

        /**
         * @brief Add trace contributions into elemental coefficient spaces.
         * 
         * Given some quantity \f$ \vec{q} \f$, calculate the elemental integral
         * 
         * \f[ 
         * \int_{\Omega^e} \vec{q}, \mathrm{d}S
         * \f] 
         * 
         * and adds this to the coefficient space provided by
         * outarray. The value of q is determined from the routine
         * IsLeftAdjacentTrace() which if true we use Fwd else we use
         * Bwd
         * 
         * @see Expansion2D::AddEdgeNormBoundaryInt
         * 
         * @param Fwd       The trace quantities associated with left (fwd)
         *                  adjancent elmt.
         * @param Bwd       The trace quantities associated with right (bwd)
         *                  adjacent elet.
         * @param outarray  Resulting 2D coefficient space.
         */
        void DisContField2D::v_AddFwdBwdTraceIntegral(
            const Array<OneD, const NekDouble> &Fwd, 
            const Array<OneD, const NekDouble> &Bwd, 
                  Array<OneD,       NekDouble> &outarray)
        {
            int e,n,offset, t_offset;
            Array<OneD, NekDouble> e_outarray;
            Array<OneD, Array<OneD, LocalRegions::ExpansionSharedPtr> >
                &elmtToTrace = m_traceMap->GetElmtToTrace();

            for (n = 0; n < GetExpSize(); ++n)
            {
                offset = GetCoeff_Offset(n);
                for (e = 0; e < (*m_exp)[n]->GetNedges(); ++e)
                {
                    t_offset = GetTrace()->GetPhys_Offset(
                                            elmtToTrace[n][e]->GetElmtId());
                    
                    // Evaluate upwind flux less local edge 
                    if (IsLeftAdjacentTrace(n, e))
                    {
                        (*m_exp)[n]->AddEdgeNormBoundaryInt(
                        e, elmtToTrace[n][e], Fwd+t_offset,
                        e_outarray = outarray+offset);
                    }
                    else
                    {
                        (*m_exp)[n]->AddEdgeNormBoundaryInt(
                        e, elmtToTrace[n][e], Bwd+t_offset,
                        e_outarray = outarray+offset);
                    }

                }
            }
        }

        /**
         * @brief Set up a list of element IDs and edge IDs that link to the
         * boundary conditions.
         */
        void DisContField2D::v_GetBoundaryToElmtMap(
            Array<OneD, int> &ElmtID, 
            Array<OneD, int> &EdgeID)
        {
            if (m_BCtoElmMap.num_elements() == 0)
            {
                map<int, int> globalIdMap;
                int i,n;
                int cnt;
                int nbcs = 0;

                // Populate global ID map (takes global geometry ID to local
                // expansion list ID).
                for (i = 0; i < GetExpSize(); ++i)
                {
                    globalIdMap[(*m_exp)[i]->GetGeom()->GetGlobalID()] = i;
                }

                // Determine number of boundary condition expansions.
                for(i = 0; i < m_bndConditions.num_elements(); ++i)
                {
                    nbcs += m_bndCondExpansions[i]->GetExpSize();
                }

                // Initialize arrays
                m_BCtoElmMap = Array<OneD, int>(nbcs);
                m_BCtoEdgMap = Array<OneD, int>(nbcs);

                LocalRegions::Expansion1DSharedPtr exp1d;
                for (cnt = n = 0; n < m_bndCondExpansions.num_elements(); ++n)
                {
                    for (i = 0; i < m_bndCondExpansions[n]->GetExpSize(); 
                         ++i, ++cnt)
                    {
                        exp1d = m_bndCondExpansions[n]->GetExp(i)->
                                            as<LocalRegions::Expansion1D>();
                        // Use edge to element map from MeshGraph.
                        SpatialDomains::GeometryLinkSharedPtr tmp =
                            m_graph->GetElementsFromEdge(exp1d->GetGeom1D());

                        m_BCtoElmMap[cnt] = globalIdMap[
                            (*tmp)[0].first->GetGlobalID()];
                        m_BCtoEdgMap[cnt] = (*tmp)[0].second;
                    }
                }
            }
            ElmtID = m_BCtoElmMap;
            EdgeID = m_BCtoEdgMap;
        }
        
        void DisContField2D::v_GetBndElmtExpansion(int i,
                            std::shared_ptr<ExpList> &result,
                            const bool DeclareCoeffPhysArrays)
        {
            int n, cnt, nq;
            int offsetOld, offsetNew;
            std::vector<unsigned int> eIDs;
            
            Array<OneD, int> ElmtID,EdgeID;
            GetBoundaryToElmtMap(ElmtID,EdgeID);
            
            // Skip other boundary regions
            for (cnt = n = 0; n < i; ++n)
            {
                cnt += m_bndCondExpansions[n]->GetExpSize();
            }

            // Populate eIDs with information from BoundaryToElmtMap
            for (n = 0; n < m_bndCondExpansions[i]->GetExpSize(); ++n)
            {
                eIDs.push_back(ElmtID[cnt+n]);
            }
            
            // Create expansion list
            result = 
                MemoryManager<ExpList>::AllocateSharedPtr
                    (*this, eIDs, DeclareCoeffPhysArrays);
            
            // Copy phys and coeffs to new explist
            if( DeclareCoeffPhysArrays)
            {
                Array<OneD, NekDouble> tmp1, tmp2;
                for (n = 0; n < result->GetExpSize(); ++n)
                {
                    nq = GetExp(ElmtID[cnt+n])->GetTotPoints();
                    offsetOld = GetPhys_Offset(ElmtID[cnt+n]);
                    offsetNew = result->GetPhys_Offset(n);
                    Vmath::Vcopy(nq, tmp1 = GetPhys()+ offsetOld, 1,
                                tmp2 = result->UpdatePhys()+ offsetNew, 1);

                    nq = GetExp(ElmtID[cnt+n])->GetNcoeffs();
                    offsetOld = GetCoeff_Offset(ElmtID[cnt+n]);
                    offsetNew = result->GetCoeff_Offset(n);
                    Vmath::Vcopy(nq, tmp1 = GetCoeffs()+ offsetOld, 1,
                                tmp2 = result->UpdateCoeffs()+ offsetNew, 1);
                }
            }
        }

        /**
         * @brief Reset this field, so that geometry information can be updated.
         */
        void DisContField2D::v_Reset()
        {
            ExpList::v_Reset();

            // Reset boundary condition expansions.
            for (int n = 0; n < m_bndCondExpansions.num_elements(); ++n)
            {
                m_bndCondExpansions[n]->Reset();
            }
        }

        /** 
         * @brief Calculate the \f$ L^2 \f$ error of the \f$ Q_{\rm dir} \f$
         * derivative using the consistent DG evaluation of \f$ Q_{\rm dir} \f$.
         * 
         * The solution provided is of the primative variation at the quadrature
         * points and the derivative is compared to the discrete derivative at
         * these points, which is likely to be undesirable unless using a much
         * higher number of quadrature points than the polynomial order used to
         * evaluate \f$ Q_{\rm dir} \f$.
        */
        NekDouble DisContField2D::L2_DGDeriv(
            const int                           dir,
            const Array<OneD, const NekDouble> &soln)
        {
            int    i,e,ncoeff_edge;
            Array<OneD, const NekDouble> tmp_coeffs;
            Array<OneD, NekDouble> out_d(m_ncoeffs), out_tmp;

            Array<OneD, Array<OneD, LocalRegions::ExpansionSharedPtr> > 
                &elmtToTrace = m_traceMap->GetElmtToTrace();

            StdRegions::Orientation edgedir;

            int     cnt;
            int     LocBndCoeffs = m_traceMap->GetNumLocalBndCoeffs();
            Array<OneD, NekDouble> loc_lambda(LocBndCoeffs), edge_lambda;

            
            m_traceMap->GlobalToLocalBnd(m_trace->GetCoeffs(),loc_lambda);

            edge_lambda = loc_lambda;
            
            // Calculate Q using standard DG formulation.
            for(i = cnt = 0; i < GetExpSize(); ++i)
            {
                // Probably a better way of setting up lambda than this.
                // Note cannot use PutCoeffsInToElmts since lambda space
                // is mapped during the solve.
                int nEdges = (*m_exp)[i]->GetNedges();
                Array<OneD, Array<OneD, NekDouble> > edgeCoeffs(nEdges);

                for(e = 0; e < nEdges; ++e)
                {
                    edgedir = (*m_exp)[i]->GetEorient(e);
                    ncoeff_edge = elmtToTrace[i][e]->GetNcoeffs();
                    edgeCoeffs[e] = Array<OneD, NekDouble>(ncoeff_edge);
                    Vmath::Vcopy(ncoeff_edge, edge_lambda, 1, edgeCoeffs[e], 1);
                    elmtToTrace[i][e]->SetCoeffsToOrientation(
                        edgedir, edgeCoeffs[e], edgeCoeffs[e]);
                    edge_lambda = edge_lambda + ncoeff_edge;
                }

                (*m_exp)[i]->DGDeriv(dir,
                                       tmp_coeffs=m_coeffs+m_coeff_offset[i],
                                       elmtToTrace[i],
                                       edgeCoeffs,
                                       out_tmp = out_d+cnt);
                cnt  += (*m_exp)[i]->GetNcoeffs();
            }
            
            BwdTrans(out_d,m_phys);
            Vmath::Vsub(m_npoints,m_phys,1,soln,1,m_phys,1);
            return L2(m_phys);
        }


        /**
         * @brief Calculates the result of the multiplication of a global matrix
         * of type specified by @a mkey with a vector given by @a inarray.
         * 
         * @param mkey      Key representing desired matrix multiplication.
         * @param inarray   Input vector.
         * @param outarray  Resulting multiplication.
         */
        void DisContField2D::v_GeneralMatrixOp(
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
         * @brief Search through the edge expansions and identify which ones
         * have Robin/Mixed type boundary conditions.
         * 
         * If a Robin boundary is found then store the edge ID of the boundary
         * condition and the array of points of the physical space boundary
         * condition which are hold the boundary condition primitive variable
         * coefficient at the quatrature points
         *
         * @return A map containing the Robin boundary condition information
         *         using a key of the element ID.
         */
        map<int, RobinBCInfoSharedPtr> DisContField2D::v_GetRobinBCInfo(void)
        {
            int i,cnt;
            map<int, RobinBCInfoSharedPtr> returnval;
            Array<OneD, int> ElmtID,EdgeID;
            GetBoundaryToElmtMap(ElmtID,EdgeID);

            for(cnt = i = 0; i < m_bndCondExpansions.num_elements(); ++i)
            {
                MultiRegions::ExpListSharedPtr locExpList;

                if(m_bndConditions[i]->GetBoundaryConditionType() == 
                       SpatialDomains::eRobin)
                {
                    int e,elmtid;
                    Array<OneD, NekDouble> Array_tmp;

                    locExpList = m_bndCondExpansions[i];

                    int npoints    = locExpList->GetNpoints();
                    Array<OneD, NekDouble> x0(npoints, 0.0);
                    Array<OneD, NekDouble> x1(npoints, 0.0);
                    Array<OneD, NekDouble> x2(npoints, 0.0);
                    Array<OneD, NekDouble> coeffphys(npoints);

                    locExpList->GetCoords(x0, x1, x2);

                    LibUtilities::Equation coeffeqn =
                        std::static_pointer_cast<
                            SpatialDomains::RobinBoundaryCondition>
                        (m_bndConditions[i])->m_robinPrimitiveCoeff;

                    // evalaute coefficient 
                    coeffeqn.Evaluate(x0, x1, x2, 0.0, coeffphys);

                    for(e = 0; e < locExpList->GetExpSize(); ++e)
                    {
                        RobinBCInfoSharedPtr rInfo =
                            MemoryManager<RobinBCInfo>
                            ::AllocateSharedPtr(
                                EdgeID[cnt+e],
                                Array_tmp = coeffphys + 
                                locExpList->GetPhys_Offset(e));
                        
                        elmtid = ElmtID[cnt+e];
                        // make link list if necessary
                        if(returnval.count(elmtid) != 0)
                        {
                            rInfo->next = returnval.find(elmtid)->second;
                        }
                        returnval[elmtid] = rInfo;
                    }
                }
                cnt += m_bndCondExpansions[i]->GetExpSize();
            }

            return returnval;
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
        void  DisContField2D::EvaluateHDGPostProcessing(
            Array<OneD, NekDouble> &outarray)
        {
            int    i,cnt,e,ncoeff_edge;
            Array<OneD, NekDouble> force, out_tmp, qrhs, qrhs1;
            Array<OneD, Array< OneD, LocalRegions::ExpansionSharedPtr> > 
                &elmtToTrace = m_traceMap->GetElmtToTrace();

            StdRegions::Orientation edgedir;

            int     nq_elmt, nm_elmt;
            int     LocBndCoeffs = m_traceMap->GetNumLocalBndCoeffs();
            Array<OneD, NekDouble> loc_lambda(LocBndCoeffs), edge_lambda;
            Array<OneD, NekDouble> tmp_coeffs;
            m_traceMap->GlobalToLocalBnd(m_trace->GetCoeffs(),loc_lambda);

            edge_lambda = loc_lambda;

            // Calculate Q using standard DG formulation.
            for(i = cnt = 0; i < GetExpSize(); ++i)
            {
                nq_elmt = (*m_exp)[i]->GetTotPoints();
                nm_elmt = (*m_exp)[i]->GetNcoeffs();
                qrhs  = Array<OneD, NekDouble>(nq_elmt);
                qrhs1  = Array<OneD, NekDouble>(nq_elmt);
                force = Array<OneD, NekDouble>(2*nm_elmt);
                out_tmp = force + nm_elmt;
                LocalRegions::ExpansionSharedPtr ppExp;

                int num_points0 = (*m_exp)[i]->GetBasis(0)->GetNumPoints();
                int num_points1 = (*m_exp)[i]->GetBasis(1)->GetNumPoints();
                int num_modes0 = (*m_exp)[i]->GetBasis(0)->GetNumModes();
                int num_modes1 = (*m_exp)[i]->GetBasis(1)->GetNumModes();

                // Probably a better way of setting up lambda than this.  Note
                // cannot use PutCoeffsInToElmts since lambda space is mapped
                // during the solve.
                int nEdges = (*m_exp)[i]->GetNedges();
                Array<OneD, Array<OneD, NekDouble> > edgeCoeffs(nEdges);

                for(e = 0; e < (*m_exp)[i]->GetNedges(); ++e)
                {
                    edgedir = (*m_exp)[i]->GetEorient(e);
                    ncoeff_edge = elmtToTrace[i][e]->GetNcoeffs();
                    edgeCoeffs[e] = Array<OneD, NekDouble>(ncoeff_edge);
                    Vmath::Vcopy(ncoeff_edge, edge_lambda, 1, edgeCoeffs[e], 1);
                    elmtToTrace[i][e]->SetCoeffsToOrientation(
                        edgedir, edgeCoeffs[e], edgeCoeffs[e]);
                    edge_lambda = edge_lambda + ncoeff_edge;
                }

                //creating orthogonal expansion (checking if we have quads or triangles)
                LibUtilities::ShapeType shape = (*m_exp)[i]->DetShapeType();
                switch(shape)
                {
                    case LibUtilities::eQuadrilateral:
                    {
                        const LibUtilities::PointsKey PkeyQ1(num_points0,LibUtilities::eGaussLobattoLegendre);
                        const LibUtilities::PointsKey PkeyQ2(num_points1,LibUtilities::eGaussLobattoLegendre);
                        LibUtilities::BasisKey  BkeyQ1(LibUtilities::eOrtho_A, num_modes0, PkeyQ1);
                        LibUtilities::BasisKey  BkeyQ2(LibUtilities::eOrtho_A, num_modes1, PkeyQ2);
                        SpatialDomains::QuadGeomSharedPtr qGeom = std::dynamic_pointer_cast<SpatialDomains::QuadGeom>((*m_exp)[i]->GetGeom());
                        ppExp = MemoryManager<LocalRegions::QuadExp>::AllocateSharedPtr(BkeyQ1, BkeyQ2, qGeom);
                    }
                    break;
                    case LibUtilities::eTriangle:
                    {
                        const LibUtilities::PointsKey PkeyT1(num_points0,LibUtilities::eGaussLobattoLegendre);
                        const LibUtilities::PointsKey PkeyT2(num_points1,LibUtilities::eGaussRadauMAlpha1Beta0);
                        LibUtilities::BasisKey  BkeyT1(LibUtilities::eOrtho_A, num_modes0, PkeyT1);
                        LibUtilities::BasisKey  BkeyT2(LibUtilities::eOrtho_B, num_modes1, PkeyT2);
                        SpatialDomains::TriGeomSharedPtr tGeom = std::dynamic_pointer_cast<SpatialDomains::TriGeom>((*m_exp)[i]->GetGeom());
                        ppExp = MemoryManager<LocalRegions::TriExp>::AllocateSharedPtr(BkeyT1, BkeyT2, tGeom);
                    }
                    break;
                    default:
                        ASSERTL0(false, "Wrong shape type, HDG postprocessing is not implemented");
                };

               
                //DGDeriv    
                // (d/dx w, d/dx q_0)
                (*m_exp)[i]->DGDeriv(
                    0,tmp_coeffs = m_coeffs + m_coeff_offset[i],
                    elmtToTrace[i], edgeCoeffs, out_tmp);
                (*m_exp)[i]->BwdTrans(out_tmp,qrhs);
                //(*m_exp)[i]->IProductWRTDerivBase(0,qrhs,force);
                ppExp->IProductWRTDerivBase(0,qrhs,force);


                // + (d/dy w, d/dy q_1)
                (*m_exp)[i]->DGDeriv(
                    1,tmp_coeffs = m_coeffs + m_coeff_offset[i],
                    elmtToTrace[i], edgeCoeffs, out_tmp);

                (*m_exp)[i]->BwdTrans(out_tmp,qrhs);
                //(*m_exp)[i]->IProductWRTDerivBase(1,qrhs,out_tmp);
                ppExp->IProductWRTDerivBase(1,qrhs,out_tmp);

                Vmath::Vadd(nm_elmt,force,1,out_tmp,1,force,1);

                // determine force[0] = (1,u)
                (*m_exp)[i]->BwdTrans(
                    tmp_coeffs = m_coeffs + m_coeff_offset[i],qrhs);
                force[0] = (*m_exp)[i]->Integral(qrhs);

                // multiply by inverse Laplacian matrix
                // get matrix inverse
                LocalRegions::MatrixKey  lapkey(StdRegions::eInvLaplacianWithUnityMean, ppExp->DetShapeType(), *ppExp);
                DNekScalMatSharedPtr lapsys = ppExp->GetLocMatrix(lapkey); 
                
                NekVector<NekDouble> in (nm_elmt,force,eWrapper);
                NekVector<NekDouble> out(nm_elmt);

                out = (*lapsys)*in;

                // Transforming back to modified basis
                Array<OneD, NekDouble> work(nq_elmt);
                ppExp->BwdTrans(out.GetPtr(), work);
                (*m_exp)[i]->FwdTrans(work, tmp_coeffs = outarray + m_coeff_offset[i]);
            }
        }
    } // end of namespace
} //end of namespace
