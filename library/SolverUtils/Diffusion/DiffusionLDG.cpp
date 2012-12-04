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
            LibUtilities::SessionReaderSharedPtr        pSession)
        {
            m_session = pSession;
        }
        
        void DiffusionLDG::v_Diffuse(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray)
        {
            int i, j, k;
            int nVelDim         = fields[0]->GetCoordim(0);
            int nPointsTot      = fields[0]->GetTotPoints();
            int nCoeffs         = fields[0]->GetNcoeffs();
            int nTracePointsTot = fields[0]->GetTrace()->GetTotPoints();
            int nqvar           = 2;
            
            Array<OneD, NekDouble>  qcoeffs(nCoeffs);
            Array<OneD, NekDouble>  temp   (nCoeffs);
            
            Array<OneD, Array<OneD, NekDouble> > fluxvector(nVelDim);
            Array<OneD, Array<OneD, NekDouble> > 
                                        temp_outarray(nConvectiveFields);
            
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > flux  (nqvar);
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > qfield(nqvar);
            
            for (j = 0; j < nqvar; ++j)
            {
                qfield[j] = Array<OneD, Array<OneD, NekDouble> >(nqvar);
                flux[j]   = Array<OneD, Array<OneD, NekDouble> >(nqvar);
                
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    qfield[j][i] = Array<OneD, NekDouble>(nPointsTot, 0.0);
                    flux[j][i]   = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
                }
            }
            
            for (k = 0; k < nVelDim; ++k)
            {
                fluxvector[k] = Array<OneD, NekDouble>(nPointsTot, 0.0);
            }
                        
            // Compute q_{\eta} and q_{\xi}
            // Obtain numerical fluxes
            v_NumFluxforScalar(fields, inarray, flux);
            
            for (j = 0; j < nqvar; ++j)
            {
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    // Get the ith component of the  flux vector
                    // fluxvector = m_tanbasis * u 
                    // where m_tanbasis = 2 * nVelDim * nPointsTot
                    if (m_tanbasis.num_elements())
                    {
                        for (k = 0; k < nVelDim; ++k)
                        {
                            Vmath::Vmul(nPointsTot, 
                                        m_tanbasis[j][k], 1, 
                                        inarray[i], 1, 
                                        fluxvector[k], 1);
                        }
                    }
                    else
                    {
                        // Get the ith component of the  flux vector in the 
                        // physical space
                        m_fluxVector(i, j, inarray, fluxvector);                    
                    }
                    
                    // Calculate the i^th value of (\grad_i \phi, F)
                    v_WeakAdvectionGreensDivergenceForm(fields, 
                                                        fluxvector, 
                                                        qcoeffs);
                    
                    Vmath::Neg(nCoeffs, qcoeffs, 1);
                    fields[i]->AddTraceIntegral(flux[j][i], qcoeffs);
                    fields[i]->SetPhysState(false);
                    fields[i]->MultiplyByElmtInvMass(qcoeffs, qcoeffs);
                    fields[i]->BwdTrans(qcoeffs, qfield[j][i]);
                }
            }
            
            
            // Compute u from q_{\eta} and q_{\xi}
            // Obtain numerical fluxes
            v_NumFluxforVector(fields, inarray, qfield, flux[0]);
            
            for (i = 0; i < nConvectiveFields; ++i)
            {
                // L = L(tan_eta) q_eta + L(tan_xi) q_xi
                temp_outarray[i] = Array<OneD, NekDouble>(nCoeffs, 0.0);
                temp             = Array<OneD, NekDouble>(nCoeffs, 0.0);
                
                if (m_tanbasis.num_elements())
                {
                    for (j = 0; j < nqvar; ++j)
                    {
                        for (k = 0; k < nVelDim; ++k)
                        {
                            Vmath::Vmul(nPointsTot, m_tanbasis[j][k], 1,
                                        qfield[j][i], 1, fluxvector[k], 1);
                        }
                        
                        v_WeakAdvectionGreensDivergenceForm(fields, 
                                                            fluxvector, 
                                                            temp);
                        
                        Vmath::Vadd(nCoeffs, temp, 1, temp_outarray[i], 1,
                                    temp_outarray[i], 1);
                    }
                }
                else
                {
                    for (k = 0; k < nVelDim; ++k)
                    {
                        Vmath::Vcopy(nPointsTot, 
                                     qfield[k][i], 1, 
                                     fluxvector[k], 1);
                    }
                    
                    v_WeakAdvectionGreensDivergenceForm(
                                                    fields, 
                                                    fluxvector, 
                                                    temp_outarray[i]);
                }
                
                // Evaulate  <\phi, \hat{F}\cdot n> - outarray[i]
                Vmath::Neg                      (nCoeffs, temp_outarray[i], 1);
                fields[i]->AddTraceIntegral     (flux[0][i], temp_outarray[i]);
                fields[i]->SetPhysState(false);
                fields[i]->MultiplyByElmtInvMass(temp_outarray[i], 
                                                 temp_outarray[i]);
                fields[i]->BwdTrans             (temp_outarray[i], outarray[i]);
            }
        }
        
        void DiffusionLDG::v_NumFluxforScalar(
            const Array<OneD, MultiRegions::ExpListSharedPtr>        &fields,
            const Array<OneD, Array<OneD, NekDouble> >               &ufield,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &uflux)
        {
            int i, j;
            int nTracePointsTot = fields[0]->GetTrace()->GetTotPoints();
            int nvariables      = fields.num_elements();
            int nqvar           = uflux.num_elements();
            int nDimensions     = fields[0]->GetCoordim(0);
            NekDouble time = 0.0;
            
            Array<OneD, NekDouble > Fwd     (nTracePointsTot);
            Array<OneD, NekDouble > Bwd     (nTracePointsTot);
            Array<OneD, NekDouble > Vn      (nTracePointsTot, 0.0);
            Array<OneD, NekDouble > fluxtemp(nTracePointsTot, 0.0);
            
            // Setting up the normals
            m_traceNormals = Array<OneD, Array<OneD, NekDouble> >(nDimensions);
            for(i = 0; i < nDimensions; ++i)
            {
                m_traceNormals[i] = Array<OneD, NekDouble> (nTracePointsTot);
            }
            fields[0]->GetTrace()->GetNormals(m_traceNormals);
            
            // Get the sign of (v \cdot n), v = an arbitrary vector
            
            // Evaulate upwind flux:
            // uflux = \hat{u} \phi \cdot u = u^{(+,-)} n
            for (j = 0; j < nqvar; ++j)
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
                    
                    fields[i]->GetTrace()->Upwind(m_traceNormals[j], 
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
                        v_WeakPenaltyforScalar(fields, i, ufield[i], fluxtemp, 
                                               time);
                    }
                    
                    // if Vn >= 0, flux = uFwd*(tan_{\xi}^- \cdot \vec{n}), 
                    // i.e,
                    // edge::eForward, uFwd \(\tan_{\xi}^Fwd \cdot \vec{n})
                    // edge::eBackward, uFwd \(\tan_{\xi}^Bwd \cdot \vec{n})
                    
                    // else if Vn < 0, flux = uBwd*(tan_{\xi}^- \cdot \vec{n}), 
                    // i.e,
                    // edge::eForward, uBwd \(\tan_{\xi}^Fwd \cdot \vec{n})
                    // edge::eBackward, uBwd \(\tan_{\xi}^Bwd \cdot \vec{n})
                    
                    Vmath::Vmul(nTracePointsTot, 
                                m_traceNormals[j], 1, 
                                fluxtemp, 1, 
                                uflux[j][i], 1);
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
            int nTracePointsTot = fields[0]->GetTrace()->GetTotPoints();
            int nvariables      = fields.num_elements();
            int nqvar           = qfield.num_elements();
            int nDimensions     = fields[0]->GetCoordim(0);
            NekDouble time      = 0.0;
            
            NekDouble C11 = 1.0;
            Array<OneD, NekDouble > Fwd(nTracePointsTot);
            Array<OneD, NekDouble > Bwd(nTracePointsTot);
            Array<OneD, NekDouble > Vn (nTracePointsTot, 0.0);
            
            Array<OneD, NekDouble > qFwd     (nTracePointsTot);
            Array<OneD, NekDouble > qBwd     (nTracePointsTot);
            Array<OneD, NekDouble > qfluxtemp(nTracePointsTot, 0.0);
            
            Array<OneD, NekDouble > uterm(nTracePointsTot);
            
            // Setting up the normals
            m_traceNormals = Array<OneD, Array<OneD, NekDouble> >(nDimensions);
            for(i = 0; i < nDimensions; ++i)
            {
                m_traceNormals[i] = Array<OneD, NekDouble> (nTracePointsTot);
            }
            fields[0]->GetTrace()->GetNormals(m_traceNormals);
            
            // Evaulate upwind flux:
            // qflux = \hat{q} \cdot u = q \cdot n - C_(11)*(u^+ - u^-)
            for (i = 0; i < nvariables; ++i)
            {
                qflux[i] = Array<OneD, NekDouble> (nTracePointsTot, 0.0);
                for (j = 0; j < nqvar; ++j)
                {
                    //  Compute Fwd and Bwd value of ufield of jth direction
                    fields[i]->GetFwdBwdTracePhys(qfield[j][i],qFwd,qBwd);
                    
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
                    
                    fields[i]->GetTrace()->Upwind(m_traceNormals[j], 
                                                    qBwd, qFwd, 
                                                    qfluxtemp);
                    
                    Vmath::Vmul(nTracePointsTot, 
                                m_traceNormals[j], 1, 
                                qfluxtemp, 1, 
                                qfluxtemp, 1);
                    
                    // Generate Stability term = - C11 ( u- - u+ )
                    fields[i]->GetFwdBwdTracePhys(ufield[i], Fwd, Bwd);
                    
                    Vmath::Vsub(nTracePointsTot, 
                                Fwd, 1, Bwd, 1, 
                                uterm, 1);
                    
                    Vmath::Smul(nTracePointsTot, 
                                -1.0 * C11, uterm, 1, 
                                uterm, 1);
                    
                    // Flux = {Fwd, Bwd} * (nx, ny, nz) + uterm * (nx, ny)
                    Vmath::Vadd(nTracePointsTot, 
                                uterm, 1, 
                                qfluxtemp, 1, 
                                qfluxtemp, 1);
                    
                    // Imposing weak boundary condition with flux
                    if (fields[0]->GetBndCondExpansions().num_elements())
                    {
                        WeakPenaltyforVector(fields, 
                                             i, j, 
                                             qfield[j][i], 
                                             qfluxtemp, 
                                             C11,
                                             time);
                    }
                    
                    // q_hat \cdot n = (q_xi \cdot n_xi) or (q_eta \cdot n_eta)
                    // n_xi = n_x * tan_xi_x + n_y * tan_xi_y + n_z * tan_xi_z
                    // n_xi = n_x * tan_eta_x + n_y * tan_eta_y + n_z*tan_eta_z
                    Vmath::Vadd(nTracePointsTot, 
                                qfluxtemp, 1, 
                                qflux[i], 1, 
                                qflux[i], 1);
                }
            }
        }
        
        void DiffusionLDG::v_WeakPenaltyforScalar(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const int                                          var,
            const Array<OneD, const NekDouble>                &physfield,
                  Array<OneD,       NekDouble>                &penaltyflux,
            NekDouble                                          time)
        {
            int i, j, e, npoints, id1, id2;
            
            // Number of boundary regions
            int cnt = 0;
            int nbnd = fields[var]->GetBndCondExpansions().num_elements();
            int Nfps, numBDEdge;
            int nDimensions     = fields[0]->GetCoordim(0);
            int nTracePointsTot = fields[0]->GetTrace()->GetTotPoints();
            
            Array<OneD, NekDouble > uplus(nTracePointsTot);
            
            fields[var]->ExtractTracePhys(physfield, uplus);
            for (i = 0; i < nbnd; ++i)
            {
                // Number of boundary expansion related to that region
                numBDEdge = fields[var]->
                    GetBndCondExpansions()[i]->GetExpSize();
                
                // Evaluate boundary values g_D or g_N from input files                
                LibUtilities::EquationSharedPtr ifunc = 
                    m_session->GetFunction("InitialConditions", 0);
                
                npoints = fields[var]->
                    GetBndCondExpansions()[i]->GetNpoints();
                
                Array<OneD,NekDouble> BDphysics(npoints);
                Array<OneD,NekDouble> x0(npoints, 0.0);
                Array<OneD,NekDouble> x1(npoints, 0.0);
                Array<OneD,NekDouble> x2(npoints, 0.0);
                
                fields[var]->GetBndCondExpansions()[i]->GetCoords(x0, x1, x2);
                ifunc->Evaluate(x0, x1, x2, time, BDphysics);
                
                // Weakly impose boundary conditions by modifying flux values
                for (e = 0; e < numBDEdge ; ++e)
                {
                    // Number of points on the expansion
                    Nfps = fields[var]->
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
                        Vmath::Vcopy(Nfps, 
                                     &BDphysics[id1], 1, 
                                     &penaltyflux[id2], 1);
                    }
                    
                    // For Neumann boundary condition: uflux = u+
                    else if ((fields[var]->GetBndConditions()[i])->
                    GetBoundaryConditionType() == SpatialDomains::eNeumann)
                    {
                        Vmath::Vcopy(Nfps, 
                                     &uplus[id2], 1, 
                                     &penaltyflux[id2], 1);
                    }
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
            const Array<OneD, const NekDouble>                &physfield,
                  Array<OneD,       NekDouble>                &penaltyflux,
            NekDouble                                          C11,
            NekDouble                                          time)
        {
            int i, j, e, npoints, id1, id2;
            int nbnd = fields[var]->GetBndCondExpansions().num_elements();
            int numBDEdge, Nfps;
            int nDimensions     = fields[0]->GetCoordim(0);
            int nTracePointsTot = fields[0]->GetTrace()->GetTotPoints();
            
            Array<OneD, NekDouble > uterm(nTracePointsTot);
            Array<OneD, NekDouble > qtemp(nTracePointsTot);
            int cnt = 0;
            
            // Setting up the normals
            m_traceNormals = Array<OneD, Array<OneD, NekDouble> >(nDimensions);
            for(i = 0; i < nDimensions; ++i)
            {
                m_traceNormals[i] = Array<OneD, NekDouble> (nTracePointsTot);
            }
            fields[0]->GetTrace()->GetNormals(m_traceNormals);
            
            fields[var]->ExtractTracePhys(physfield,qtemp);
            
            for (i = 0; i < nbnd; ++i)
            {
                numBDEdge = fields[var]->
                    GetBndCondExpansions()[i]->GetExpSize();
                
                // Evaluate boundary values g_D or g_N from input files
                LibUtilities::EquationSharedPtr ifunc = 
                m_session->GetFunction("InitialConditions", 0);
                
                npoints = fields[var]->
                    GetBndCondExpansions()[i]->GetNpoints();
                
                Array<OneD,NekDouble> BDphysics(npoints);
                Array<OneD,NekDouble> x0(npoints, 0.0);
                Array<OneD,NekDouble> x1(npoints, 0.0);
                Array<OneD,NekDouble> x2(npoints, 0.0);
                
                fields[var]->GetBndCondExpansions()[i]->GetCoords(x0, x1, x2);
                ifunc->Evaluate(x0, x1, x2, time, BDphysics);
                
                // Weakly impose boundary conditions by modifying flux values
                for (e = 0; e < numBDEdge ; ++e)
                {
                    Nfps = fields[var]->
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
                        Vmath::Vmul(Nfps, 
                                    &m_traceNormals[dir][id2], 1, 
                                    &qtemp[id2], 1, 
                                    &penaltyflux[id2], 1);
                    }
                    // For Neumann boundary condition: qflux = g_N
                    else if((fields[var]->GetBndConditions()[i])->
                    GetBoundaryConditionType() == SpatialDomains::eNeumann)
                    {
                        Vmath::Vmul(Nfps,
                                    &m_traceNormals[dir][id2], 1, 
                                    &BDphysics[id1], 1, 
                                    &penaltyflux[id2], 1);
                    }
                }
            }
        }
        
        /**
         * Computes the weak Green form of advection terms (without boundary
         * integral), i.e. \f$ (\nabla \phi \cdot F) \f$ where for example
         * \f$ F=uV \f$.
         * @param   fields      Fields.
         * @param   F           Given data.
         * @param   outarray    Storage for result.
         *
         * \note Assuming all fields are of the same expansion and order so  
         * that we can use the parameters of m_fields[0].
         */
        void DiffusionLDG::v_WeakAdvectionGreensDivergenceForm(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &F,
                  Array<OneD, NekDouble>                      &Fout)
        {
            // Use dimension of velocity vector to dictate operation dimensions
            int ndim    = F.num_elements();
            int nCoeffs = fields[0]->GetNcoeffs();
            
            Array<OneD, NekDouble> iprod(nCoeffs);
            Vmath::Zero(nCoeffs, Fout, 1);
            
            for (int i = 0; i < ndim; ++i)
            {
                fields[0]->IProductWRTDerivBase(i, F[i], iprod);
                Vmath::Vadd(nCoeffs, iprod, 1, Fout, 1, Fout, 1);
            }
        }

    }
}
