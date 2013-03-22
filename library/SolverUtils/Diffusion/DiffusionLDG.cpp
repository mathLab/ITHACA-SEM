///////////////////////////////////////////////////////////////////////////////
//
// File: DiffusionLDG.cpp
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
// Description: LDG diffusion class.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Diffusion/DiffusionLDG.h>
#include <iostream>
#include <iomanip>

namespace Nektar
{
    namespace SolverUtils
    {
        std::string DiffusionLDG::type = GetDiffusionFactory().
            RegisterCreatorFunction("LDG", DiffusionLDG::create);
        
        DiffusionLDG::DiffusionLDG()
        {
        }
        
        void DiffusionLDG::v_InitObject(
            LibUtilities::SessionReaderSharedPtr        pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
        {
            m_session = pSession;
        }
        
        void DiffusionLDG::v_Diffuse(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray)
        {
            //cout<<setprecision(16);
            int i, j, k;
            int nDim      = fields[0]->GetCoordim(0);
            int nPts      = fields[0]->GetTotPoints();
            int nCoeffs   = fields[0]->GetNcoeffs();
            int nTracePts = fields[0]->GetTrace()->GetTotPoints();
            
            Array<OneD, NekDouble>  qcoeffs(nCoeffs);
            Array<OneD, NekDouble>  temp   (nCoeffs);
            
            Array<OneD, Array<OneD, NekDouble> > fluxvector(nDim);
            Array<OneD, Array<OneD, NekDouble> > tmp(nConvectiveFields);
            
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > flux  (nDim);
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > qfield(nDim);
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > qfieldStd(nDim);

            

            for (j = 0; j < nDim; ++j)
            {
                qfield[j] = 
                    Array<OneD, Array<OneD, NekDouble> >(nConvectiveFields);
                qfieldStd[j] = 
                Array<OneD, Array<OneD, NekDouble> >(nConvectiveFields);
                flux[j]   = 
                    Array<OneD, Array<OneD, NekDouble> >(nConvectiveFields);
                
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    qfield[j][i] = Array<OneD, NekDouble>(nPts, 0.0);
                    qfieldStd[j][i] = Array<OneD, NekDouble>(nPts, 0.0);
                    flux[j][i]   = Array<OneD, NekDouble>(nTracePts, 0.0);
                }
            }
            
            for (k = 0; k < nDim; ++k)
            {
                fluxvector[k] = Array<OneD, NekDouble>(nPts, 0.0);
            }
                        
            // Compute q_{\eta} and q_{\xi}
            // Obtain numerical fluxes
            v_NumFluxforScalar(fields, inarray, flux);
            
            for (j = 0; j < nDim; ++j)
            {
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    fields[i]->IProductWRTDerivBase(j, inarray[i], qcoeffs);
                    Vmath::Neg                      (nCoeffs, qcoeffs, 1);
                    fields[i]->AddTraceIntegral     (flux[j][i], qcoeffs);
                    fields[i]->SetPhysState         (false);
                    fields[i]->MultiplyByElmtInvMass(qcoeffs, qcoeffs);
                    fields[i]->BwdTrans             (qcoeffs, qfield[j][i]);
                }
            }
            
            // Compute u from q_{\eta} and q_{\xi}
            // Obtain numerical fluxes
            v_NumFluxforVector(fields, inarray, qfield, flux[0]);
            
            for (i = 0; i < nConvectiveFields; ++i)
            {
                tmp[i] = Array<OneD, NekDouble>(nCoeffs, 0.0);
                
                for (j = 0; j < nDim; ++j)
                {
                    Vmath::Vcopy(nPts, qfield[j][i], 1, fluxvector[j], 1);
                    fields[i]->IProductWRTDerivBase(j, fluxvector[j], qcoeffs);
                    Vmath::Vadd(nCoeffs, qcoeffs, 1, tmp[i], 1, tmp[i], 1);
                }
                
                // Evaulate  <\phi, \hat{F}\cdot n> - outarray[i]
                Vmath::Neg                      (nCoeffs, tmp[i], 1);
                fields[i]->AddTraceIntegral     (flux[0][i], tmp[i]);
                fields[i]->SetPhysState         (false);
                fields[i]->MultiplyByElmtInvMass(tmp[i], tmp[i]);
                fields[i]->BwdTrans             (tmp[i], outarray[i]);
            }
        }
        
        
        
        void DiffusionLDG::v_NumFluxforScalar(
            const Array<OneD, MultiRegions::ExpListSharedPtr>        &fields,
            const Array<OneD, Array<OneD, NekDouble> >               &ufield,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &uflux)
        {
            int i, j;
            int nTracePts  = fields[0]->GetTrace()->GetTotPoints();
            int nvariables = fields.num_elements();
            int nDim       = uflux.num_elements();
            NekDouble time = 0.0;
            
            Array<OneD, NekDouble > Fwd     (nTracePts);
            Array<OneD, NekDouble > Bwd     (nTracePts);
            Array<OneD, NekDouble > Vn      (nTracePts, 0.0);
            Array<OneD, NekDouble > fluxtemp(nTracePts, 0.0);
            
            // Setting up the normals
            m_traceNormals = Array<OneD, Array<OneD, NekDouble> >(nDim);
            for(i = 0; i < nDim; ++i)
            {
                m_traceNormals[i] = Array<OneD, NekDouble> (nTracePts);
            }
            fields[0]->GetTrace()->GetNormals(m_traceNormals);
            
            // Get the normal velocity Vn
            for(i = 0; i < nDim; ++i)
            {
                Vmath::Svtvp(nTracePts, 1.0, m_traceNormals[i], 1, 
                             Vn, 1, Vn, 1);
            }
            
            // Get the sign of (v \cdot n), v = an arbitrary vector
            // Evaluate upwind flux:
            // uflux = \hat{u} \phi \cdot u = u^{(+,-)} n
            for (j = 0; j < nDim; ++j)
            {
                for (i = 0; i < nvariables ; ++i)
                {
                    // Compute Fwd and Bwd value of ufield of i direction
                    fields[i]->GetFwdBwdTracePhys(ufield[i], Fwd, Bwd);
                    
                    // if Vn >= 0, flux = uFwd, i.e.,
                    // edge::eForward, if V*n>=0 <=> V*n_F>=0, pick uflux = uFwd
                    // edge::eBackward, if V*n>=0 <=> V*n_B<0, pick uflux = uFwd
                    
                    // else if Vn < 0, flux = uBwd, i.e.,
                    // edge::eForward, if V*n<0 <=> V*n_F<0, pick uflux = uBwd
                    // edge::eBackward, if V*n<0 <=> V*n_B>=0, pick uflux = uBwd
                    
                    fields[i]->GetTrace()->Upwind(/*m_traceNormals[j]*/Vn, 
                                                    Fwd, Bwd, fluxtemp);
                    
                    // Imposing weak boundary condition with flux
                    // if Vn >= 0, uflux = uBwd at Neumann, i.e.,
                    // edge::eForward, if V*n>=0 <=> V*n_F>=0, pick uflux = uBwd
                    // edge::eBackward, if V*n>=0 <=> V*n_B<0, pick uflux = uBwd
                    
                    // if Vn >= 0, uflux = uFwd at Neumann, i.e.,
                    // edge::eForward, if V*n<0 <=> V*n_F<0, pick uflux = uFwd
                    // edge::eBackward, if V*n<0 <=> V*n_B>=0, pick uflux = uFwd
                    
                    if(fields[0]->GetBndCondExpansions().num_elements())
                    {
                        v_WeakPenaltyforScalar(fields, i, ufield[i], fluxtemp);
                    }
                    
                    // if Vn >= 0, flux = uFwd*(tan_{\xi}^- \cdot \vec{n}), 
                    // i.e,
                    // edge::eForward, uFwd \(\tan_{\xi}^Fwd \cdot \vec{n})
                    // edge::eBackward, uFwd \(\tan_{\xi}^Bwd \cdot \vec{n})
                    
                    // else if Vn < 0, flux = uBwd*(tan_{\xi}^- \cdot \vec{n}), 
                    // i.e,
                    // edge::eForward, uBwd \(\tan_{\xi}^Fwd \cdot \vec{n})
                    // edge::eBackward, uBwd \(\tan_{\xi}^Bwd \cdot \vec{n})
                    
                    Vmath::Vmul(nTracePts, 
                                m_traceNormals[j], 1, 
                                fluxtemp, 1, 
                                uflux[j][i], 1);
                }
            }
        }
        
        
        
        void DiffusionLDG::v_WeakPenaltyforScalar(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const int                                          var,
            const Array<OneD, const NekDouble>                &ufield,
                  Array<OneD,       NekDouble>                &penaltyflux)
        {
            int i, j, e, id1, id2;
            
            // Number of boundary regions
            int nBndEdgePts, nBndEdges;
            int cnt         = 0;
            int nBndRegions = fields[var]->GetBndCondExpansions().num_elements();
            int nDim        = fields[0]->GetCoordim(0);
            int nTracePts   = fields[0]->GetTrace()->GetTotPoints();
            
            Array<OneD, NekDouble > uplus(nTracePts);
            
            fields[var]->ExtractTracePhys(ufield, uplus);
            for (i = 0; i < nBndRegions; ++i)
            {
                // Number of boundary expansion related to that region
                nBndEdges = fields[var]->
                GetBndCondExpansions()[i]->GetExpSize();
                                                                                
                // Weakly impose boundary conditions by modifying flux values
                for (e = 0; e < nBndEdges ; ++e)
                {
                    // Number of points on the expansion
                    nBndEdgePts = fields[var]->
                    GetBndCondExpansions()[i]->GetExp(e)->GetNumPoints(0);
                    
                    id1 = fields[var]->
                    GetBndCondExpansions()[i]->GetPhys_Offset(e);
                    
                    id2 = fields[0]->GetTrace()->
                    GetPhys_Offset(fields[0]->GetTraceMap()->
                                   GetBndCondTraceToGlobalTraceMap(cnt++));
                    
                    // For Dirichlet boundary condition: uflux = g_D
                    if (fields[var]->GetBndConditions()[i]->
                        GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                    {
                        Vmath::Vcopy(nBndEdgePts, 
                                     &(fields[var]->
                                       GetBndCondExpansions()[i]->
                                       GetPhys())[id1], 1, 
                                     &penaltyflux[id2], 1);
                    }
                    // For Neumann boundary condition: uflux = u+
                    else if ((fields[var]->GetBndConditions()[i])->
                        GetBoundaryConditionType() == SpatialDomains::eNeumann)
                    {
                        Vmath::Vcopy(nBndEdgePts, 
                                     &uplus[id2], 1, 
                                     &penaltyflux[id2], 1);
                    }
                }
            }
        }
        
        
        
        void DiffusionLDG::v_NumFluxforVector(
            const Array<OneD, MultiRegions::ExpListSharedPtr>        &fields,
            const Array<OneD, Array<OneD, NekDouble> >               &ufield,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &qfield,
                  Array<OneD, Array<OneD, NekDouble> >               &qflux)
        {
            int i, j;
            int nTracePts  = fields[0]->GetTrace()->GetTotPoints();
            int nvariables = fields.num_elements();
            int nDim       = qfield.num_elements();
            
            NekDouble C11 = 0.0;
            Array<OneD, NekDouble > Fwd(nTracePts);
            Array<OneD, NekDouble > Bwd(nTracePts);
            Array<OneD, NekDouble > Vn (nTracePts, 0.0);
            
            Array<OneD, NekDouble > qFwd     (nTracePts);
            Array<OneD, NekDouble > qBwd     (nTracePts);
            Array<OneD, NekDouble > qfluxtemp(nTracePts, 0.0);
            
            Array<OneD, NekDouble > uterm(nTracePts);
            
            // Setting up the normals
            m_traceNormals = Array<OneD, Array<OneD, NekDouble> >(nDim);
            for(i = 0; i < nDim; ++i)
            {
                m_traceNormals[i] = Array<OneD, NekDouble> (nTracePts);
            }
            fields[0]->GetTrace()->GetNormals(m_traceNormals);
            
            // Get the normal velocity Vn
            for(i = 0; i < nDim; ++i)
            {
                Vmath::Svtvp(nTracePts, 1.0, m_traceNormals[i], 1, 
                             Vn, 1, Vn, 1);
            }
            
            // Evaulate upwind flux:
            // qflux = \hat{q} \cdot u = q \cdot n - C_(11)*(u^+ - u^-)
            for (i = 0; i < nvariables; ++i)
            {
                qflux[i] = Array<OneD, NekDouble> (nTracePts, 0.0);
                for (j = 0; j < nDim; ++j)
                {
                    //  Compute Fwd and Bwd value of ufield of jth direction
                    fields[i]->GetFwdBwdTracePhys(qfield[j][i], qFwd, qBwd);
                    
                    // if Vn >= 0, flux = uFwd, i.e.,
                    // edge::eForward, if V*n>=0 <=> V*n_F>=0, pick 
                    // qflux = qBwd = q+
                    // edge::eBackward, if V*n>=0 <=> V*n_B<0, pick 
                    // qflux = qBwd = q-
                    
                    // else if Vn < 0, flux = uBwd, i.e.,
                    // edge::eForward, if V*n<0 <=> V*n_F<0, pick 
                    // qflux = qFwd = q-
                    // edge::eBackward, if V*n<0 <=> V*n_B>=0, pick 
                    // qflux = qFwd = q+
                    
                    fields[i]->GetTrace()->Upwind(/*m_traceNormals[j]*/Vn, 
                                                  qBwd, qFwd, 
                                                  qfluxtemp);
                    
                    Vmath::Vmul(nTracePts, 
                                m_traceNormals[j], 1, 
                                qfluxtemp, 1, 
                                qfluxtemp, 1);
                    
                    // Generate Stability term = - C11 ( u- - u+ )
                    fields[i]->GetFwdBwdTracePhys(ufield[i], Fwd, Bwd);
                    
                    Vmath::Vsub(nTracePts, 
                                Fwd, 1, Bwd, 1, 
                                uterm, 1);
                    
                    Vmath::Smul(nTracePts, 
                                -1.0 * C11, uterm, 1, 
                                uterm, 1);
                    
                    // Flux = {Fwd, Bwd} * (nx, ny, nz) + uterm * (nx, ny)
                    Vmath::Vadd(nTracePts, 
                                uterm, 1, 
                                qfluxtemp, 1, 
                                qfluxtemp, 1);
                    
                    // Imposing weak boundary condition with flux
                    if (fields[0]->GetBndCondExpansions().num_elements())
                    {
                        v_WeakPenaltyforVector(fields, i, j, 
                                               qfield[j][i], 
                                               qfluxtemp, C11);
                    }
                    
                    // q_hat \cdot n = (q_xi \cdot n_xi) or (q_eta \cdot n_eta)
                    // n_xi = n_x * tan_xi_x + n_y * tan_xi_y + n_z * tan_xi_z
                    // n_xi = n_x * tan_eta_x + n_y * tan_eta_y + n_z*tan_eta_z
                    Vmath::Vadd(nTracePts, 
                                qfluxtemp, 1, 
                                qflux[i], 1, 
                                qflux[i], 1);
                }
            }
        }
         
        
        
        /**
         * Diffusion: Imposing weak boundary condition for q with flux
         *  uflux = g_D  on Dirichlet boundary condition
         *  uflux = u_Fwd  on Neumann boundary condition
         */
        void DiffusionLDG::v_WeakPenaltyforVector(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const int                                          var,
            const int                                          dir,
            const Array<OneD, const NekDouble>                &qfield,
                  Array<OneD,       NekDouble>                &penaltyflux,
            NekDouble                                          C11)
        {
            int i, j, e, id1, id2;
            int nBndEdges, nBndEdgePts;
            int nBndRegions = fields[var]->GetBndCondExpansions().num_elements();
            int nDim        = fields[0]->GetCoordim(0);
            int nTracePts   = fields[0]->GetTrace()->GetTotPoints();
            
            Array<OneD, NekDouble > uterm(nTracePts);
            Array<OneD, NekDouble > qtemp(nTracePts);
            int cnt = 0;
            
            // Setting up the normals
            m_traceNormals = Array<OneD, Array<OneD, NekDouble> >(nDim);
            for(i = 0; i < nDim; ++i)
            {
                m_traceNormals[i] = Array<OneD, NekDouble> (nTracePts);
            }
            fields[0]->GetTrace()->GetNormals(m_traceNormals);
            
            fields[var]->ExtractTracePhys(qfield, qtemp);
            
            for (i = 0; i < nBndRegions; ++i)
            {
                nBndEdges = fields[var]->
                    GetBndCondExpansions()[i]->GetExpSize();
                                
                // Weakly impose boundary conditions by modifying flux values
                for (e = 0; e < nBndEdges ; ++e)
                {
                    nBndEdgePts = fields[var]->
                    GetBndCondExpansions()[i]->GetExp(e)->GetNumPoints(0);
                    
                    id1 = fields[var]->
                    GetBndCondExpansions()[i]->GetPhys_Offset(e);
                    
                    id2 = fields[0]->GetTrace()->
                    GetPhys_Offset(fields[0]->GetTraceMap()->
                                   GetBndCondTraceToGlobalTraceMap(cnt++));
                    
                    // For Dirichlet boundary condition: 
                    //qflux = q+ - C_11 (u+ -    g_D) (nx, ny)
                    if(fields[var]->GetBndConditions()[i]->
                    GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                    {
                        Vmath::Vmul(nBndEdgePts, 
                                    &m_traceNormals[dir][id2], 1, 
                                    &qtemp[id2], 1, 
                                    &penaltyflux[id2], 1);
                    }
                    // For Neumann boundary condition: qflux = g_N
                    else if((fields[var]->GetBndConditions()[i])->
                    GetBoundaryConditionType() == SpatialDomains::eNeumann)
                    {
                        Vmath::Vmul(nBndEdgePts,
                                    &m_traceNormals[dir][id2], 1, 
                                    &(fields[var]->
                                      GetBndCondExpansions()[i]->
                                      GetPhys())[id1], 1, 
                                    &penaltyflux[id2], 1);
                    }
                }
            }
        }
        
    }
}
