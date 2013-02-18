///////////////////////////////////////////////////////////////////////////////
//
// File: DiffusionLDGNS.cpp
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
// Description: LDGNS diffusion class.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Diffusion/DiffusionLDGNS.h>

namespace Nektar
{
    namespace SolverUtils
    {
        std::string DiffusionLDGNS::type = GetDiffusionFactory().
            RegisterCreatorFunction("LDGNS", DiffusionLDGNS::create);
        
        DiffusionLDGNS::DiffusionLDGNS()
        {
        }
        
        void DiffusionLDGNS::v_InitObject(
            LibUtilities::SessionReaderSharedPtr pSession)
        {
            m_session = pSession;
        }
        
        void DiffusionLDGNS::v_Diffuse(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray)
        {
            int i, j;
            int nDim      = fields[0]->GetCoordim(0);
            int nScalars  = inarray.num_elements();
            int nPts      = fields[0]->GetTotPoints();
            int nCoeffs   = fields[0]->GetNcoeffs();
            int nTracePts = fields[0]->GetTrace()->GetTotPoints();
            
            Array<OneD, NekDouble>                             qcoeffs(nCoeffs);
            Array<OneD, Array<OneD, NekDouble> >               fluxvector(nDim);
            Array<OneD, Array<OneD, NekDouble> >               tmp_out(nConvectiveFields);
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > flux(nDim);
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > qfield(nDim);
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > viscousFlux(nDim);
            
            for (j = 0; j < nDim; ++j)
            {
                fluxvector[j]  = Array<OneD, NekDouble>(nPts, 0.0);
                qfield[j]      = Array<OneD, Array<OneD, NekDouble> >(nConvectiveFields);
                viscousFlux[j] = Array<OneD, Array<OneD, NekDouble> >(nConvectiveFields);
                flux[j]        = Array<OneD, Array<OneD, NekDouble> >(nConvectiveFields);
                
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    qfield[j][i]      = Array<OneD, NekDouble>(nPts, 0.0);
                    viscousFlux[j][i] = Array<OneD, NekDouble>(nPts, 0.0);
                    flux[j][i]        = Array<OneD, NekDouble>(nTracePts, 0.0);
                }
            }
            
            // Compute q_{\eta} and q_{\xi}
            // Obtain numerical fluxes
            v_NumFluxforScalar(fields, inarray, flux);
            
            for (j = 0; j < nDim; ++j)
            {
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    fields[i]->IProductWRTDerivBase(j, inarray[i], qcoeffs);
                    Vmath::Neg(nCoeffs, qcoeffs, 1);
                    fields[i]->AddTraceIntegral(flux[j][i], qcoeffs);
                    fields[i]->SetPhysState(false);
                    fields[i]->MultiplyByElmtInvMass(qcoeffs, qcoeffs);
                    fields[i]->BwdTrans(qcoeffs, qfield[j][i]);
                }
            }
            
            // Get the ith component of the  flux vector in the 
            // physical space
            for (i = 0; i < nConvectiveFields; ++i)
            {
                m_fluxVectorNS(i, inarray, qfield, viscousFlux); 
            }
            
            // Compute u from q_{\eta} and q_{\xi}
            // Obtain numerical fluxes
            v_NumFluxforVector(fields, inarray, viscousFlux, flux[0]);
            
            for (i = 0; i < nConvectiveFields; ++i)
            {
                tmp_out[i] = Array<OneD, NekDouble>(nCoeffs, 0.0);
                
                for (j = 0; j < nDim; ++j)
                {
                    Vmath::Vcopy(nPts, qfield[j][i], 1, fluxvector[j], 1);
                    
                    fields[i]->IProductWRTDerivBase(j, fluxvector[j], qcoeffs);
                    
                    Vmath::Vadd(nCoeffs, qcoeffs, 1, 
                                tmp_out[i], 1, 
                                tmp_out[i], 1);
                }
                
                // Evaulate  <\phi, \hat{F}\cdot n> - outarray[i]
                Vmath::Neg(nCoeffs, tmp_out[i], 1);
                fields[i]->AddTraceIntegral(flux[0][i], tmp_out[i]);
                fields[i]->SetPhysState(false);
                fields[i]->MultiplyByElmtInvMass(tmp_out[i], tmp_out[i]);
                fields[i]->BwdTrans(tmp_out[i], outarray[i]);
            }
        }
        
        void DiffusionLDGNS::v_NumFluxforScalar(
            const Array<OneD, MultiRegions::ExpListSharedPtr>        &fields,
            const Array<OneD, Array<OneD, NekDouble> >               &ufield,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &uflux)
        {
            // 1. ==============================================================
            // Initialise useful variables -------------------------------------
            int i, j;
            int nTracePts = fields[0]->GetTrace()->GetTotPoints();
            int nvariables      = fields.num_elements();
            int nDim            = fields[0]->GetCoordim(0);
            NekDouble time      = 0.0;
            
            //Array<OneD, NekDouble > Fwd     (nTracePts);
            //Array<OneD, NekDouble > Bwd     (nTracePts);
            //Array<OneD, NekDouble > Vn      (nTracePts, 0.0);
            Array<OneD, NekDouble > fluxtemp(nTracePts, 0.0);
            //------------------------------------------------------------------
            // =================================================================
            
            
            
            // 2. ==============================================================
            // Setting up the trace normals ------------------------------------
            m_traceNormals = Array<OneD, Array<OneD, NekDouble> >(nDim);
            for(i = 0; i < nDim; ++i)
            {
                m_traceNormals[i] = Array<OneD, NekDouble>(nTracePts);
            }
            fields[0]->GetTrace()->GetNormals(m_traceNormals);
            // -----------------------------------------------------------------
            // =================================================================
            
            // Store forwards/backwards space along trace space.
            Array<OneD, Array<OneD, NekDouble> > Fwd    (nvariables);
            Array<OneD, Array<OneD, NekDouble> > Bwd    (nvariables);
            Array<OneD, Array<OneD, NekDouble> > numflux(nvariables);
            
            for(i = 0; i < nvariables; ++i)
            {
                Fwd[i]     = Array<OneD, NekDouble>(nTracePts);
                Bwd[i]     = Array<OneD, NekDouble>(nTracePts);
                numflux[i] = Array<OneD, NekDouble>(nTracePts);
                fields[i]->GetFwdBwdTracePhys(ufield[i], Fwd[i], Bwd[i]);
            }
            
            m_riemann->Solve(Fwd, Bwd, numflux);

            // 3. ==============================================================
            // Evaulate Riemann flux: uflux = \hat{u} \phi \cdot u = u^{(+,-)} n
            for (j = 0; j < nDim; ++j)
            {
                for (i = 0; i < nvariables ; ++i)
                {
                                        
                    // Modify the values in case of boundary interfaces
                    if(fields[0]->GetBndCondExpansions().num_elements())
                    {
                        v_WeakPenaltyforScalar(fields, i, ufield[i], 
                                               numflux[i], time);
                    }
                    
                    // Multiply the Riemann flux by the trace normals and 
                    // store it into uflux
                    Vmath::Vmul(nTracePts, 
                                m_traceNormals[j], 1, 
                                fluxtemp, 1, uflux[j][i], 1);
                }
            }
            //==================================================================
        }        
        
        void DiffusionLDGNS::v_WeakPenaltyforScalar(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const int                                          var,
            const Array<OneD, const NekDouble>                &ufield,
                  Array<OneD,       NekDouble>                &penaltyflux,
            NekDouble                                          time)
        {
            // 1. ==============================================================
            // Initialise useful variables -------------------------------------
            int cnt  = 0;
            int nBndEdgePts, nBndEdges;
            int i, e, npoints, id1, id2;
            int nDim = fields[0]->GetCoordim(0);
            int nTracePoints = fields[0]->GetTrace()->GetTotPoints();
            int nBndRegions  = fields[var]->GetBndCondExpansions().num_elements();
            Array<OneD, NekDouble > uplus(nTracePoints);
            //------------------------------------------------------------------
            // =================================================================
            
            
            
            // 2. ==============================================================
            // Extract the physical values of the solution at the boundaries ---
            fields[var]->ExtractTracePhys(ufield, uplus);
            //------------------------------------------------------------------
            // =================================================================
            
            
            
            // 3. ============================================================== 
            // Loop on the boundary regions to apply appropriate bcs -----------
            for (i = 0; i < nBndRegions; ++i)
            {
                // 3.1) Number of boundary edges related to region 'i'
                nBndEdges = fields[var]->
                    GetBndCondExpansions()[i]->GetExpSize();
                
                // 3.2) Weakly impose bcs by modifying flux values
                for (e = 0; e < nBndEdges; ++e)
                {
                    nBndEdgePts = fields[var]->
                    GetBndCondExpansions()[i]->GetExp(e)->GetNumPoints(0);
                    
                    id1 = fields[var]->
                    GetBndCondExpansions()[i]->GetPhys_Offset(e);
                    
                    id2 = fields[0]->GetTrace()->
                    GetPhys_Offset(fields[0]->GetTraceMap()->
                                   GetBndCondTraceToGlobalTraceMap(cnt++));

                    // 3.3) In case of Dirichlet bcs: uflux = g_D
                    if (fields[var]->GetBndConditions()[i]->
                    GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                    {   
                        Vmath::Vcopy(nBndEdgePts, 
                                     &(fields[var]->
                                       GetBndCondExpansions()[i]->
                                       UpdatePhys())[id1], 1, 
                                     &penaltyflux[id2], 1);
                        
                    }
                    // 3.4) In case of Neumann bcs: uflux = u+
                    else if ((fields[var]->GetBndConditions()[i])->
                    GetBoundaryConditionType() == SpatialDomains::eNeumann)
                    {
                        Vmath::Vcopy(nBndEdgePts, 
                                     &uplus[id2], 1, 
                                     &penaltyflux[id2], 1);
                    }
                }
            }
            //------------------------------------------------------------------
            // =================================================================
        }
        
        
        void DiffusionLDGNS::v_NumFluxforVector(
            const Array<OneD, MultiRegions::ExpListSharedPtr>        &fields,
            const Array<OneD, Array<OneD, NekDouble> >               &ufield,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &qfield,
                  Array<OneD, Array<OneD, NekDouble> >               &qflux)
        {
            // 1. ==============================================================
            // Initialise useful variables -------------------------------------
            int i, j;
            int nTracePoints = fields[0]->GetTrace()->GetTotPoints();
            int nVariables   = fields.num_elements();
            int nDim         = fields[0]->GetCoordim(0);
            NekDouble time   = 0.0;
            
            NekDouble C11 = 1.0;
            Array<OneD, NekDouble > Fwd(nTracePoints);
            Array<OneD, NekDouble > Bwd(nTracePoints);
            Array<OneD, NekDouble > Vn (nTracePoints, 0.0);
            
            Array<OneD, NekDouble > qFwd     (nTracePoints);
            Array<OneD, NekDouble > qBwd     (nTracePoints);
            Array<OneD, NekDouble > qfluxtemp(nTracePoints, 0.0);
            
            Array<OneD, NekDouble > uterm(nTracePoints);
            //------------------------------------------------------------------
            // =================================================================
            
            
            
            // 2. ==============================================================
            // Setting up the normals ------------------------------------------
            m_traceNormals = Array<OneD, Array<OneD, NekDouble> >(nDim);
            for(i = 0; i < nDim; ++i)
            {
                m_traceNormals[i] = Array<OneD, NekDouble> (nTracePoints);
            }
            fields[0]->GetTrace()->GetNormals(m_traceNormals);
            //------------------------------------------------------------------
            // =================================================================
            
            
            
            // 3. ==============================================================
            // Evaulate Riemann flux 
            // qflux = \hat{q} \cdot u = q \cdot n - C_(11)*(u^+ - u^-)
            for (i = 0; i < nVariables; ++i)
            {
                qflux[i] = Array<OneD, NekDouble> (nTracePoints, 0.0);
                for (j = 0; j < nDim; ++j)
                {
                    // 3.1) Compute qFwd and qBwd value of qfield in 
                    // position 'ji'
                    fields[i]->GetFwdBwdTracePhys(qfield[j][i], qFwd, qBwd);
                    
                    // 3.2) Get Riemann flux of qflux --> LDG implies upwind
                    fields[i]->GetTrace()->Upwind(m_traceNormals[j], 
                                                  qBwd, qFwd, 
                                                  qfluxtemp);
                    
                    // 3.3) Multiply the Riemann flux by the trace normals
                    Vmath::Vmul(nTracePoints, 
                                m_traceNormals[j], 1, 
                                qfluxtemp, 1, 
                                qfluxtemp, 1);
                    
                    // 3.4) Compute Fwd and Bwd value of ufield for the 
                    // stability term
                    fields[i]->GetFwdBwdTracePhys(ufield[i], Fwd, Bwd);
                    
                    // 3.5) Compute the stability term = - C11( u- - u+ )
                    Vmath::Vsub(nTracePoints, 
                                Fwd, 1, 
                                Bwd, 1, 
                                uterm, 1);
                    
                    Vmath::Smul(nTracePoints, 
                                -1.0 * C11, 
                                uterm, 1, 
                                uterm, 1);
                    
                    // 3.6) Compute the flux
                    // Flux = {Fwd, Bwd} * (nx, ny, nz) + uterm * (nx, ny)
                    Vmath::Vadd(nTracePoints, 
                                uterm, 1, 
                                qfluxtemp, 1, 
                                qfluxtemp, 1);
                    
                    // 3.7) Impose weak boundary condition with flux
                    if (fields[0]->GetBndCondExpansions().num_elements())
                    {
                        v_WeakPenaltyforVector(fields, i, j,
                                             qfield[j][i], 
                                             qfluxtemp, 
                                             C11, time);
                    }
                    
                    // 3.8) Store the final flux into qflux
                    Vmath::Vadd(nTracePoints, 
                                qfluxtemp, 1, 
                                qflux[i], 1, 
                                qflux[i], 1);
                }
            }
            //------------------------------------------------------------------
            // =================================================================
        }
        
        
        
        /**
         * Diffusion: Imposing weak boundary condition for q with flux
         *  uflux = g_D  on Dirichlet boundary condition
         *  uflux = u_Fwd  on Neumann boundary condition
         */
        void DiffusionLDGNS::v_WeakPenaltyforVector(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const int                                          var,
            const int                                          dir,
            const Array<OneD, const NekDouble>                &qfield,
                  Array<OneD,       NekDouble>                &penaltyflux,
            NekDouble                                          C11,
            NekDouble                                          time)
        {
            // 1. ==============================================================
            // Initialise useful variables -------------------------------------
            int cnt = 0;
            int nBndEdges, nBndEdgePts;
            int i, j, e, npoints, id1, id2;
            int nDim = fields[0]->GetCoordim(0);
            int nTracePoints = fields[0]->GetTrace()->GetTotPoints();
            int nBndRegions = fields[var]->GetBndCondExpansions().num_elements();
            Array<OneD, NekDouble > uterm(nTracePoints);
            Array<OneD, NekDouble > qtemp(nTracePoints);
            //------------------------------------------------------------------
            // =================================================================
            
            
            
            // 2. ==============================================================
            // Setting up the normals ------------------------------------------
            m_traceNormals = Array<OneD, Array<OneD, NekDouble> >(nDim);
            for(i = 0; i < nDim; ++i)
            {
                m_traceNormals[i] = Array<OneD, NekDouble> (nTracePoints);
            }
            fields[0]->GetTrace()->GetNormals(m_traceNormals);
            //------------------------------------------------------------------
            // =================================================================


            
            // 2. ==============================================================
            // Extract the physical values of the solution at the boundaries ---
            fields[var]->ExtractTracePhys(qfield, qtemp);
            //------------------------------------------------------------------
            // =================================================================
            
            
            
            // 3. ============================================================== 
            // Loop on the boundary regions to apply appropriate bcs -----------
            for (i = 0; i < nBndRegions; ++i)
            {
                // 3.1) Number of boundary regions related to region 'i'
                nBndEdges = fields[var]->
                    GetBndCondExpansions()[i]->GetExpSize();
                
                // 3.2) Weakly impose bcs by modifying flux values
                for (e = 0; e < nBndEdges; ++e)
                {
                    nBndEdgePts = fields[var]->
                    GetBndCondExpansions()[i]->GetExp(e)->GetNumPoints(0);
                    
                    id1 = fields[var]->
                    GetBndCondExpansions()[i]->GetPhys_Offset(e);
                    
                    id2 = fields[0]->GetTrace()->
                        GetPhys_Offset(fields[0]->GetTraceMap()->
                                   GetBndCondTraceToGlobalTraceMap(cnt++));
                                        
                    // 3.3) In case of Dirichlet bcs: 
                    // uflux = g_D
                    // qflux = q+ - C_11 (u+  - g_D) (nx, ny)
                    if(fields[var]->GetBndConditions()[i]->
                    GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                    {
                        Vmath::Vmul(nBndEdgePts, 
                                    &m_traceNormals[dir][id2], 1, 
                                    &qtemp[id2], 1, 
                                    &penaltyflux[id2], 1);
                    }
                    // 3.4) In case of Neumann bcs: 
                    // uflux = u+
                    else if((fields[var]->GetBndConditions()[i])->
                    GetBoundaryConditionType() == SpatialDomains::eNeumann)
                    {
                        Vmath::Vmul(nBndEdgePts, 
                                    &m_traceNormals[dir][id2], 1, 
                                    &(fields[var]->
                                      GetBndCondExpansions()[i]->
                                      UpdatePhys())[id1], 1, 
                                    &penaltyflux[id2], 1);
                    }
                }
            }
            //------------------------------------------------------------------
            // =================================================================
        }
        
    }//end of namespace SolverUtils
}//end of namespace Nektar
