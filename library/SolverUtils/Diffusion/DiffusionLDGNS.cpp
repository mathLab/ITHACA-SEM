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
            int i, j, k;
            int nDim            = fields[0]->GetCoordim(0);
            int nPointsTot      = fields[0]->GetTotPoints();
            int nCoeffs         = fields[0]->GetNcoeffs();
            int nTracePointsTot = fields[0]->GetTrace()->GetTotPoints();
            
            Array<OneD, NekDouble> qcoeffs(nCoeffs);
            Array<OneD, Array<OneD, NekDouble> > fluxvector(nDim);
            Array<OneD, Array<OneD, NekDouble> > tmp_out(nConvectiveFields);
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >flux(nDim);
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >qfield(nDim);
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >viscousFlux(nDim);
            
            for (j = 0; j < nDim; ++j)
            {
                qfield[j] = 
                    Array<OneD, Array<OneD, NekDouble> >(nConvectiveFields);
                viscousFlux[j] = 
                    Array<OneD, Array<OneD, NekDouble> >(nConvectiveFields);
                flux[j] = 
                    Array<OneD, Array<OneD, NekDouble> >(nConvectiveFields);
                
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    qfield[j][i] = Array<OneD, NekDouble>(nPointsTot, 0.0);
                    viscousFlux[j][i] = Array<OneD, NekDouble>(nPointsTot, 0.0);
                    flux[j][i] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
                }
            }
            
            for (k = 0; k < nDim; ++k)
            {
                fluxvector[k] = Array<OneD, NekDouble>(nPointsTot, 0.0);
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
                    Vmath::Vcopy(nPointsTot, qfield[j][i], 1, fluxvector[j], 1);
                    
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
            int nTracePointsTot = fields[0]->GetTrace()->GetTotPoints();
            int nvariables      = fields.num_elements();
            int nDim            = fields[0]->GetCoordim(0);
            NekDouble time      = 0.0;
            
            Array<OneD, NekDouble > Fwd     (nTracePointsTot);
            Array<OneD, NekDouble > Bwd     (nTracePointsTot);
            Array<OneD, NekDouble > Vn      (nTracePointsTot, 0.0);
            Array<OneD, NekDouble > fluxtemp(nTracePointsTot, 0.0);
            //------------------------------------------------------------------
            // =================================================================
            
            
            
            // 2. ==============================================================
            // Setting up the trace normals ------------------------------------
            m_traceNormals = Array<OneD, Array<OneD, NekDouble> >(nDim);
            for(i = 0; i < nDim; ++i)
            {
                m_traceNormals[i] = Array<OneD, NekDouble>(nTracePointsTot);
            }
            fields[0]->GetTrace()->GetNormals(m_traceNormals);
            // -----------------------------------------------------------------
            // =================================================================
            
            
            
            // 3. ==============================================================
            // Evaulate Riemann flux: uflux = \hat{u} \phi \cdot u = u^{(+,-)} n
            for (j = 0; j < nDim; ++j)
            {
                for (i = 0; i < nvariables ; ++i)
                {
                    // 3.1) Compute Fwd and Bwd value of ufield of 'i' equation
                    fields[i]->GetFwdBwdTracePhys(ufield[i], Fwd, Bwd);
                    
                    // 3.2) Take the Riemann flux
                    fields[i]->GetTrace()->Upwind(m_traceNormals[j], 
                                                  Fwd, Bwd, fluxtemp);
                                        
                    // 3.3) Modify the values in case of boundary interfaces
                    if(fields[0]->GetBndCondExpansions().num_elements())
                    {
                        v_WeakPenaltyforScalar(fields, i, ufield[i], 
                                               fluxtemp, time);
                    }
                    
                    // 3.4) Multiply the Riemann flux by the trace normals 
                    // and store it into uflux
                    Vmath::Vmul(nTracePointsTot, 
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
            int Nfps, numBDEdge;
            int i, j, e, npoints, id1, id2;
            int cnt  = 0;
            int nbnd = fields[var]->GetBndCondExpansions().num_elements();
            int nDim = fields[0]->GetCoordim(0);
            int nTracePoints = fields[0]->GetTrace()->GetTotPoints();
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
            for (i = 0; i < nbnd; ++i)
            {
                // 3.1) Number of boundary edges related to region 'i'
                numBDEdge = fields[var]->
                    GetBndCondExpansions()[i]->GetExpSize();
                
                // 3.2) Evaluate boundary values of variable 'var' using the 
                // initial condition tag (to be consistent with BC)
                LibUtilities::EquationSharedPtr ifunc = 
                    m_session->GetFunction("InitialConditions", /*0*/ var);
                
                // 3.3) Total number of boundary points related to region 'i'
                npoints = fields[var]->
                    GetBndCondExpansions()[i]->GetNpoints();
                
                Array<OneD, NekDouble> BDphysics(npoints);
                Array<OneD, NekDouble> ZDphysics(npoints);

                Array<OneD, NekDouble> x0(npoints, 0.0);
                Array<OneD, NekDouble> x1(npoints, 0.0);
                Array<OneD, NekDouble> x2(npoints, 0.0);
                
                // 3.4) Get coordinates of the boundary points related to 
                // region 'i'
                fields[var]->GetBndCondExpansions()[i]->GetCoords(x0, x1, x2);
                
                // 3.5) Evaluate Dirichlet boundary conditions using consistent 
                // initial conditions
                ifunc->Evaluate(x0, x1, x2, time, BDphysics);

                // 3.6) Weakly impose bcs by modifying flux values
                for (e = 0; e < numBDEdge; ++e)
                {
                    // 3.7) Number of points on the expansion
                    Nfps = fields[var]->
                    GetBndCondExpansions()[i]->GetExp(e)->GetNumPoints(0);
                    
                    // 3.8) Id of physical points related to the boundary 
                    // expansion
                    id1 = fields[var]->
                    GetBndCondExpansions()[i]->GetPhys_Offset(e);
                    
                    // 3.9) Id of trace points related to the boundary expansion
                    id2 = fields[0]->GetTrace()->
                    GetPhys_Offset(fields[0]->GetTraceMap()->
                                   GetBndCondTraceToGlobalTraceMap(cnt++));

                    // 3.10) In case of Dirichlet bcs: uflux = g_D
                    if (fields[var]->GetBndConditions()[i]->
                    GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                    {
                        if (var == 0 || var == 3)
                        {
                            Vmath::Vcopy(Nfps, 
                                         &BDphysics[id1], 1, 
                                         &penaltyflux[id2], 1);
                        }
                        else
                        {
                            Vmath::Zero(npoints, ZDphysics, 1);
                            Vmath::Vcopy(Nfps, &ZDphysics[id1], 1, 
                                         &penaltyflux[id2], 1);
                        }
                    }
                    // 3.11) In case of Neumann bcs: uflux = u+
                    else if ((fields[var]->GetBndConditions()[i])->
                    GetBoundaryConditionType() == SpatialDomains::eNeumann)
                    {
                        Vmath::Vcopy(Nfps, 
                                     &uplus[id2], 1, 
                                     &penaltyflux[id2], 1);
                    }
                    /*else if ((fields[var]->GetBndConditions()[i])->
                    GetUserDefined() == SpatialDomains::eWallViscous)
                    {
                        if (var == 0 || var == 3)
                        {
                            Vmath::Vcopy(Nfps, 
                                         &BDphysics[id1], 1, 
                                         &penaltyflux[id2], 1);
                        }
                        else if (var == 1 || var == 2)
                        {
                            // Zero at the wall
                            Vmath::Zero(npoints, ZDphysics, 1);
                            Vmath::Vcopy(Nfps, &ZDphysics[id1], 1, 
                                         &penaltyflux[id2], 1);
                        }
                    }*/
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
                    
                    // 3.2) Get Riemann flux of qflux
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
                        WeakPenaltyforVector(fields, i, j,
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
            int numBDEdge, Nfps;
            int i, j, e, npoints, id1, id2;
            int nbnd = fields[var]->GetBndCondExpansions().num_elements();
            int nDim = fields[0]->GetCoordim(0);
            int nTracePoints = fields[0]->GetTrace()->GetTotPoints();
            
            Array<OneD, NekDouble > uterm(nTracePoints);
            Array<OneD, NekDouble > qtemp(nTracePoints);
            int cnt = 0;
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
            for (i = 0; i < nbnd; ++i)
            {
                // 3.1) Number of boundary regions related to region 'i'
                numBDEdge = fields[var]->
                    GetBndCondExpansions()[i]->GetExpSize();
                
                // 3.2) Evaluate boundary values of variable 'var' using the 
                // initial condition tag (to be consistent with BC)
                LibUtilities::EquationSharedPtr ifunc = 
                    m_session->GetFunction("InitialConditions", /*0*/ var);
                
                // 3.3) Total number of boundary points related to region 'i'
                npoints = fields[var]->
                    GetBndCondExpansions()[i]->GetNpoints();
                
                Array<OneD,NekDouble> BDphysics(npoints); 
                Array<OneD,NekDouble> ZDphysics(npoints);                
                Array<OneD,NekDouble> x0(npoints, 0.0);
                Array<OneD,NekDouble> x1(npoints, 0.0);
                Array<OneD,NekDouble> x2(npoints, 0.0);
                
                // 3.4) Get coordinates of the boundary points related to 
                // region 'i'
                fields[var]->GetBndCondExpansions()[i]->GetCoords(x0, x1, x2);
                
                // 3.5) Evaluate Dirichlet boundary conditions using consistent 
                // initial conditions
                ifunc->Evaluate(x0, x1, x2, time, BDphysics);
                
                // 3.6) Weakly impose bcs by modifying flux values
                for (e = 0; e < numBDEdge; ++e)
                {
                    // 3.7) Number of points on the expansion
                    Nfps = fields[var]->
                    GetBndCondExpansions()[i]->GetExp(e)->GetNumPoints(0);
                    
                    // 3.8) Id of physical points related to the boundary expansion
                    id1 = fields[var]->
                    GetBndCondExpansions()[i]->GetPhys_Offset(e);
                    
                    // 3.9) Id of trace points related to the boundary expansion
                    id2 = fields[0]->GetTrace()->
                        GetPhys_Offset(fields[0]->GetTraceMap()->
                                   GetBndCondTraceToGlobalTraceMap(cnt++));
                    
                    // 3.10) In case of Dirichlet bcs: uflux = g_D
                    //qflux = q+ - C_11 (u+  - g_D) (nx, ny)
                    if(fields[var]->GetBndConditions()[i]->
                    GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                    {
                        if (var == 0 || var == 3)
                        {
                            Vmath::Vmul(Nfps, 
                                        &m_traceNormals[dir][id2], 1, 
                                        &qtemp[id2], 1, 
                                        &penaltyflux[id2], 1);
                        }
                        else
                        {
                            // Zero at the wall
                            Vmath::Zero(npoints, ZDphysics, 1);
                            Vmath::Vcopy(Nfps, &ZDphysics[id1], 1, 
                                         &penaltyflux[id2], 1);
                        }
                    }
                    // 3.11) In case of Neumann bcs: uflux = u+
                    else if((fields[var]->GetBndConditions()[i])->
                    GetBoundaryConditionType() == SpatialDomains::eNeumann)
                    {
                        if (var == 0 || var == 3)
                        {
                            Vmath::Vmul(Nfps,
                                        &m_traceNormals[dir][id2], 1, 
                                        &BDphysics[id1], 1, 
                                        &penaltyflux[id2], 1);
                        }
                        else
                        {
                            // Zero at the wall
                            Vmath::Zero(npoints, ZDphysics, 1);
                            Vmath::Vcopy(Nfps, &ZDphysics[id1], 1, 
                                         &penaltyflux[id2], 1);
                        }
                    }
                    /*else if ((fields[var]->GetBndConditions()[i])->
                             GetUserDefined() == SpatialDomains::eWallViscous)
                    {
                        if (var == 0 || var == 3)
                        {
                            Vmath::Vmul(Nfps, 
                                        &m_traceNormals[dir][id2], 1, 
                                        &qtemp[id2], 1, 
                                        &penaltyflux[id2], 1);
                        }
                        else if (var == 1 || var == 2)
                        {
                            // Zero at the wall
                            Vmath::Zero(npoints, ZDphysics, 1);
                            Vmath::Vcopy(Nfps, &ZDphysics[id1], 1, 
                                         &penaltyflux[id2], 1);
                        }
                    }
                     */
                }
            }
            //------------------------------------------------------------------
            // =================================================================
        }
        
    }//end of namespace SolverUtils
}//end of namespace Nektar
