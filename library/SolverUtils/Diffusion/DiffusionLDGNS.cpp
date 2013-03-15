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
#include <iostream>
#include <iomanip>

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
            LibUtilities::SessionReaderSharedPtr        pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
        {
            m_session = pSession;
            m_session->LoadParameter ("Gamma",         m_gamma, 1.4);
            m_session->LoadParameter ("GasConstant",   m_gasConstant, 287.058);
            m_session->LoadParameter ("Twall",         m_Twall, 300.15);
            m_session->LoadSolverInfo("ViscosityType", m_ViscosityType, 
                                      "Constant");
            m_session->LoadParameter ("mu",            m_mu, 1.78e-05);
            m_session->LoadParameter ("thermalConductivity",
                                      m_thermalConductivity, 0.0257);
            m_session->LoadParameter ("rhoInf",        m_rhoInf, 1.225);
            m_session->LoadParameter ("pInf",          m_pInf, 101325);
        }
        
        /**
         * @brief Calculate weak DG Diffusion in the LDG form for the 
         * Navier-Stokes (NS) equations:
         *
         * \f$ \langle\psi, \hat{u}\cdot n\rangle 
         *   - \langle\nabla\psi \cdot u\rangle
         *     \langle\phi, \hat{q}\cdot n\rangle - 
         *     (\nabla \phi \cdot q) \rangle \f$
         *
         * The equations that need a diffusion operator are those related 
         * with the velocities and with the energy.
         *
         */
        void DiffusionLDGNS::v_Diffuse(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray)
        {   
            //cout<<setprecision(16);
            int i, j;
            int nDim      = fields[0]->GetCoordim(0);
            int nScalars  = inarray.num_elements();
            int nPts      = fields[0]->GetTotPoints();
            int nCoeffs   = fields[0]->GetNcoeffs();
            int nTracePts = fields[0]->GetTrace()->GetTotPoints();
            
            Array<OneD, NekDouble>               tmp1(nCoeffs);
            Array<OneD, Array<OneD, NekDouble> > tmp2(nConvectiveFields);
            
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > 
                                                    numericalFluxO1(nDim);
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > 
                                                    derivativesO1(nDim);
            
            Array<OneD, Array<OneD, NekDouble> > fluxvector(nDim);            
            
            for (j = 0; j < nDim; ++j)
            {
                numericalFluxO1[j] = Array<OneD, Array<OneD, NekDouble> >(
                                                                    nScalars);
                derivativesO1[j]   = Array<OneD, Array<OneD, NekDouble> >(
                                                                    nScalars);
                
                for (i = 0; i < nScalars; ++i)
                {
                    numericalFluxO1[j][i] = Array<OneD, NekDouble>(
                                                                nTracePts, 0.0);
                    derivativesO1[j][i]   = Array<OneD, NekDouble>(nPts, 0.0);
                    
                }
            }
            
            // Compute the numerical fluxes for the first order derivatives
            v_NumericalFluxO1(fields, inarray, numericalFluxO1);
            
            for (j = 0; j < nDim; ++j)
            {
                for (i = 0; i < nScalars; ++i)
                {
                    fields[i]->IProductWRTDerivBase (j, inarray[i], tmp1);
                    Vmath::Neg                      (nCoeffs, tmp1, 1);
                    fields[i]->AddTraceIntegral     (numericalFluxO1[j][i], 
                                                     tmp1);
                    
                    fields[i]->SetPhysState         (false);
                    fields[i]->MultiplyByElmtInvMass(tmp1, tmp1);
                    fields[i]->BwdTrans             (tmp1, derivativesO1[j][i]);
                }
            }            
            
            // Initialisation viscous tensor
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > viscousTensor(
                                                                        nDim);
            Array<OneD, Array<OneD, NekDouble> > viscousFlux(nConvectiveFields);
            
            for (j = 0; j < nDim; ++j)
            {
                viscousTensor[j] = Array<OneD, Array<OneD, NekDouble> >(
                                                                    nScalars+1);
                for (i = 0; i < nScalars+1; ++i)
                {
                    viscousTensor[j][i] = Array<OneD, NekDouble>(nPts, 0.0);
                }
            }
            
            for (i = 0; i < nConvectiveFields; ++i)
            {
                viscousFlux[i] = Array<OneD, NekDouble>(nPts, 0.0);
            }
            
            m_fluxVectorNS(inarray, derivativesO1, viscousTensor); 
            
             // Compute u from q_{\eta} and q_{\xi}
             // Obtain numerical fluxes
             v_NumericalFluxO2(fields, inarray, viscousTensor, viscousFlux);
             
             for (i = 0; i < nConvectiveFields; ++i)
             {
                 tmp2[i] = Array<OneD, NekDouble>(nCoeffs, 0.0);
             
                 for (j = 0; j < nDim; ++j)
                 {
                     //Vmath::Vcopy(nPts, qfield[j][i], 1, fluxvector[j], 1);
             
                     fields[i]->IProductWRTDerivBase(j, viscousTensor[j][i], 
                                                     tmp1);
             
                     Vmath::Vadd(nCoeffs, tmp1, 1, tmp2[i], 1, tmp2[i], 1);
                 }
             
             // Evaulate  <\phi, \hat{F}\cdot n> - outarray[i]
             Vmath::Neg(nCoeffs, tmp2[i], 1);
             fields[i]->AddTraceIntegral(viscousFlux[i], tmp2[i]);
             fields[i]->SetPhysState(false);
             fields[i]->MultiplyByElmtInvMass(tmp2[i], tmp2[i]);
             fields[i]->BwdTrans(tmp2[i], outarray[i]);
             }
        }
        
        /**
         * @brief Builds the numerical flux for the 1st order derivatives
         *
         */
        void DiffusionLDGNS::v_NumericalFluxO1(
            const Array<OneD, MultiRegions::ExpListSharedPtr>        &fields,
            const Array<OneD, Array<OneD, NekDouble> >               &inarray,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > > 
                                                            &numericalFluxO1)
        {
            int i, j;
            int nTracePts  = fields[0]->GetTrace()->GetTotPoints();
            int nScalars   = inarray.num_elements();
            int nDim       = fields[0]->GetCoordim(0);
            
            Array<OneD, NekDouble > Vn      (nTracePts, 0.0);
            Array<OneD, NekDouble > fluxtemp(nTracePts, 0.0);
            Array<OneD, Array<OneD, NekDouble> > traceVel(nDim);
            
            m_traceNormals = Array<OneD, Array<OneD, NekDouble> >(nDim);
            for(i = 0; i < nDim; ++i)
            {
                m_traceNormals[i] = Array<OneD, NekDouble>(nTracePts);
                traceVel[i]       = Array<OneD, NekDouble>(nTracePts, 0.0);
            }
            fields[0]->GetTrace()->GetNormals(m_traceNormals);
            
            // Get the normal velocity Vn
            for(i = 0; i < nDim; ++i)
            {
                fields[0]->ExtractTracePhys(inarray[i], traceVel[i]);
                Vmath::Vvtvp(nTracePts, m_traceNormals[i], 1, 
                             traceVel[i], 1, Vn, 1, Vn, 1);
            }
            
            // Store forwards/backwards space along trace space
            Array<OneD, Array<OneD, NekDouble> > Fwd    (nScalars);
            Array<OneD, Array<OneD, NekDouble> > Bwd    (nScalars);
            Array<OneD, Array<OneD, NekDouble> > numflux(nScalars);
            
            for (i = 0; i < nScalars; ++i)
            {
                Fwd[i]     = Array<OneD, NekDouble>(nTracePts);
                Bwd[i]     = Array<OneD, NekDouble>(nTracePts);
                numflux[i] = Array<OneD, NekDouble>(nTracePts);
                fields[i]->GetFwdBwdTracePhys(inarray[i], Fwd[i], Bwd[i]);
                fields[0]->GetTrace()->Upwind(Vn, Fwd[i], Bwd[i], numflux[i]);
            }
             
            // Modify the values in case of boundary interfaces
            if (fields[0]->GetBndCondExpansions().num_elements())
            {
                v_WeakPenaltyO1(fields, inarray, numflux);
            }
            
            // Splitting the numerical flux into the dimensions
            for (j = 0; j < nDim; ++j)
            {
                for (i = 0; i < nScalars; ++i)
                {
                    Vmath::Vmul(nTracePts, m_traceNormals[j], 1, 
                                numflux[i], 1, numericalFluxO1[j][i], 1);
                }
            }
        }        
        
        /**
         * @brief Imposes appropriate bcs for the 1st order derivatives
         *
         */
        void DiffusionLDGNS::v_WeakPenaltyO1(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &penaltyfluxO1)
        {            
            int cnt;
            int i, j, e;            
            int id1, id2;
            
            int nBndEdgePts, nBndEdges, nBndRegions;
            
            int nDim              = fields[0]->GetCoordim(0);
            int nTracePts         = fields[0]->GetTrace()->GetTotPoints();
            int nConvectiveFields = fields.num_elements();
            int nScalars          = inarray.num_elements();
            
            Array<OneD, NekDouble> tmp1(nTracePts, 0.0);
            Array<OneD, NekDouble> tmp2(nTracePts, 0.0);
            Array<OneD, NekDouble> Tw(nTracePts, m_Twall);
            
            Array< OneD, Array<OneD, NekDouble > > scalarVariables(nScalars);
            Array< OneD, Array<OneD, NekDouble > > uplus(nScalars);
            
            // Extract internal values of the scalar variables for Neumann bcs
            for (i = 0; i < nScalars; ++i)
            {
                scalarVariables[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
                
                uplus[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
                fields[i]->ExtractTracePhys(inarray[i], uplus[i]);
            }
            
            // Compute boundary conditions for velocities
            for (i = 0; i < nScalars-1; ++i)
            {
                // Note that cnt has to loop on nBndRegions and nBndEdges
                // and has to be reset to zero per each equation
                cnt = 0;
                nBndRegions = fields[i+1]->
                GetBndCondExpansions().num_elements();
                for (j = 0; j < nBndRegions; ++j)
                {
                    nBndEdges = fields[i+1]->
                    GetBndCondExpansions()[j]->GetExpSize();
                    for (e = 0; e < nBndEdges; ++e)
                    {
                        nBndEdgePts = fields[i+1]->
                        GetBndCondExpansions()[j]->GetExp(e)->GetNumPoints(0);
                        
                        id1 = fields[i+1]->
                        GetBndCondExpansions()[j]->GetPhys_Offset(e);
                        
                        id2 = fields[0]->GetTrace()->
                        GetPhys_Offset(fields[0]->GetTraceMap()->
                                       GetBndCondTraceToGlobalTraceMap(cnt++));

                        // Reinforcing bcs for velocity in case of Wall bcs
                        if (fields[i]->GetBndConditions()[j]->
                            GetUserDefined() == 
                            SpatialDomains::eWallViscous)
                        {
                            Vmath::Zero(nBndEdgePts, 
                                        &scalarVariables[i][id2], 1);
                        }

                        // Imposing velocity bcs if not Wall
                        else if (fields[i]->GetBndConditions()[j]->
                                 GetBoundaryConditionType() == 
                                 SpatialDomains::eDirichlet)
                        {
                        Vmath::Vdiv(nBndEdgePts,
                                    &(fields[i+1]->GetBndCondExpansions()[j]->
                                      UpdatePhys())[id1], 1,
                                    &(fields[0]->GetBndCondExpansions()[j]->
                                      UpdatePhys())[id1], 1,
                                    &scalarVariables[i][id2], 1);
                        }
                        
                        // For Dirichlet boundary condition: uflux = u_bcs
                        if (fields[i]->GetBndConditions()[j]->
                            GetBoundaryConditionType() == 
                            SpatialDomains::eDirichlet)
                        {
                            Vmath::Vcopy(nBndEdgePts, 
                                         &scalarVariables[i][id2], 1, 
                                         &penaltyfluxO1[i][id2], 1);
                        }
                        
                        // For Neumann boundary condition: uflux = u_+
                        else if ((fields[i]->GetBndConditions()[j])->
                                 GetBoundaryConditionType() == 
                                 SpatialDomains::eNeumann)
                        {
                            Vmath::Vcopy(nBndEdgePts, 
                                         &uplus[i][id2], 1, 
                                         &penaltyfluxO1[i][id2], 1);
                        }
                        
                        // Building kinetic energy to be used for T bcs
                        Vmath::Vmul(nBndEdgePts, 
                                    &scalarVariables[i][id2], 1,
                                    &scalarVariables[i][id2], 1,
                                    &tmp1[id2], 1);
                        
                        Vmath::Smul(nBndEdgePts, 0.5,
                                    &tmp1[id2], 1,
                                    &tmp1[id2], 1);
                        
                        Vmath::Vadd(nBndEdgePts,
                                    &tmp2[id2], 1,
                                    &tmp1[id2], 1,
                                    &tmp2[id2], 1);
                    }
                }
            }
            
            // Compute boundary conditions  for temperature
            cnt = 0;
            nBndRegions = fields[nScalars]->
            GetBndCondExpansions().num_elements();
            for (j = 0; j < nBndRegions; ++j)
            {
                //cout<<"bcRegion = "<< j << endl;
                nBndEdges = fields[nScalars]->
                GetBndCondExpansions()[j]->GetExpSize();
                for (e = 0; e < nBndEdges; ++e)
                {
                    //cout<<"bcEdge = "<< e << endl;
                    nBndEdgePts = fields[nScalars]->
                    GetBndCondExpansions()[j]->GetExp(e)->GetNumPoints(0);
                    
                    id1 = fields[nScalars]->
                    GetBndCondExpansions()[j]->GetPhys_Offset(e);
                    
                    id2 = fields[0]->GetTrace()->
                    GetPhys_Offset(fields[0]->GetTraceMap()->
                                   GetBndCondTraceToGlobalTraceMap(cnt++));
                    
                    // Imposing Temperature Twall at the wall 
                    if (fields[i]->GetBndConditions()[j]->
                        GetUserDefined() == 
                        SpatialDomains::eWallViscous)
                    {                        
                        Vmath::Vcopy(nBndEdgePts, 
                                     &Tw[0], 1, 
                                     &scalarVariables[nScalars-1][id2], 1);
                    }                    
                    // Imposing Temperature through condition on the Energy
                    // for no wall boundaries (e.g. farfield)
                    else if (fields[i]->GetBndConditions()[j]->
                             GetBoundaryConditionType() == 
                             SpatialDomains::eDirichlet)
                    {
                        // Divide E by rho
                        Vmath::Vdiv(nBndEdgePts,
                                    &(fields[nScalars]->
                                      GetBndCondExpansions()[j]->
                                      GetPhys())[id1], 1, 
                                    &(fields[0]->
                                      GetBndCondExpansions()[j]->
                                      GetPhys())[id1], 1,
                                    &scalarVariables[nScalars-1][id2], 1);
                        
                        /*
                        // Subtract kinetic energy to E/rho
                        Vmath::Vsub(nBndEdgePts, 
                                    &scalarVariables[nScalars-1][id2], 1,
                                    &tmp2[id2], 1,
                                    &scalarVariables[nScalars-1][id2], 1);
                        */
                         
                        // Multiply by constant factor (gamma-1)/R 
                        Vmath::Smul(nBndEdgePts, (m_gamma - 1)/m_gasConstant,
                                    &scalarVariables[nScalars-1][id2], 1,
                                    &scalarVariables[nScalars-1][id2], 1);
                    }
                    /*
                    for (int bb = 0; bb < nBndEdgePts; ++bb)
                    {
                        cout<<"T-bcs = "<<scalarVariables[nScalars-1][id2+bb]<<endl;
                    }
                    int num;
                    cin>>num;
                    */
                    // For Dirichlet boundary condition: uflux = u_bcs
                    if (fields[nScalars]->GetBndConditions()[j]->
                        GetBoundaryConditionType() == 
                        SpatialDomains::eDirichlet)
                    {
                        Vmath::Vcopy(nBndEdgePts, 
                                     &scalarVariables[nScalars-1][id2], 1, 
                                     &penaltyfluxO1[nScalars-1][id2], 1);
                    }
                    
                    // For Neumann boundary condition: uflux = u_+
                    else if ((fields[nScalars]->GetBndConditions()[j])->
                             GetBoundaryConditionType() == 
                             SpatialDomains::eNeumann)
                    {
                        Vmath::Vcopy(nBndEdgePts, 
                                     &uplus[nScalars-1][id2], 1, 
                                     &penaltyfluxO1[nScalars-1][id2], 1);
                    }
                    /*
                    for (int bb = 0; bb < nBndEdgePts; ++bb)
                    {
                        cout<<"PenaltyFlux = "<<penaltyfluxO1[nScalars-1][id2+bb]<<endl;
                    }
                    cin>>num;
                     */
                }
            }
        }
        
        /**
         * @brief Build the numerical flux for the 2nd order derivatives
         *
         */
        void DiffusionLDGNS::v_NumericalFluxO2(
            const Array<OneD, MultiRegions::ExpListSharedPtr>        &fields,
            const Array<OneD, Array<OneD, NekDouble> >               &ufield,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &qfield,
                  Array<OneD, Array<OneD, NekDouble> >               &qflux)
        {
            int i, j;
            int nTracePts = fields[0]->GetTrace()->GetTotPoints();
            int nVariables   = fields.num_elements();
            int nDim         = fields[0]->GetCoordim(0);
            
            Array<OneD, NekDouble > Fwd(nTracePts);
            Array<OneD, NekDouble > Bwd(nTracePts);
            Array<OneD, NekDouble > Vn (nTracePts, 0.0);
            
            Array<OneD, NekDouble > qFwd     (nTracePts);
            Array<OneD, NekDouble > qBwd     (nTracePts);
            Array<OneD, NekDouble > qfluxtemp(nTracePts, 0.0);
                        
            Array<OneD, Array<OneD, NekDouble> > traceVel(nDim);
            m_traceNormals = Array<OneD, Array<OneD, NekDouble> >(nDim);
            for(i = 0; i < nDim; ++i)
            {
                m_traceNormals[i] = Array<OneD, NekDouble> (nTracePts);
                traceVel[i]       = Array<OneD, NekDouble>(nTracePts, 0.0);
            }
            fields[0]->GetTrace()->GetNormals(m_traceNormals);
            
            // Get the normal velocity Vn
            for(i = 0; i < nDim; ++i)
            {
                fields[0]->ExtractTracePhys(ufield[i], traceVel[i]);
                Vmath::Vvtvp(nTracePts, m_traceNormals[i], 1, 
                             traceVel[i], 1, Vn, 1, Vn, 1);
            }
                        
            // Evaulate Riemann flux 
            // qflux = \hat{q} \cdot u = q \cdot n 
            // Notice: i = 1 (first row of the viscous tensor is zero)
            for (i = 1; i < nVariables; ++i)
            {
                qflux[i] = Array<OneD, NekDouble> (nTracePts, 0.0);
                for (j = 0; j < nDim; ++j)
                {
                    // Compute qFwd and qBwd value of qfield in position 'ji'
                    fields[i]->GetFwdBwdTracePhys(qfield[j][i], qFwd, qBwd);
                    
                    // Get Riemann flux of qflux --> LDG implies upwind
                    fields[i]->GetTrace()->Upwind(/*m_traceNormals[j]*/Vn, 
                                                  qBwd, qFwd, 
                                                  qfluxtemp);
                    
                    // Multiply the Riemann flux by the trace normals
                    Vmath::Vmul(nTracePts, 
                                m_traceNormals[j], 1, 
                                qfluxtemp, 1, 
                                qfluxtemp, 1);
                                                            
                    // Impose weak boundary condition with flux
                    if (fields[0]->GetBndCondExpansions().num_elements())
                    {
                        v_WeakPenaltyO2(fields, i, j,
                                        qfield[j][i], 
                                        qfluxtemp);
                    }
                    
                    // Store the final flux into qflux
                    Vmath::Vadd(nTracePts, 
                                qfluxtemp, 1, 
                                qflux[i], 1, 
                                qflux[i], 1);
                }
            }
        }
        
        
        /**
         * @brief Imposes appropriate bcs for the 2nd order derivatives
         *
         */
        void DiffusionLDGNS::v_WeakPenaltyO2(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const int                                          var,
            const int                                          dir,
            const Array<OneD, const NekDouble>                &qfield,
                  Array<OneD,       NekDouble>                &penaltyflux)
        {
            int cnt = 0;
            int nBndEdges, nBndEdgePts;
            int i, e; 
            int id1, id2;
            
            int nDim        = fields[0]->GetCoordim(0);
            int nTracePts   = fields[0]->GetTrace()->GetTotPoints();
            int nBndRegions = fields[var]->GetBndCondExpansions().num_elements();
            
            Array<OneD, NekDouble > uterm(nTracePts);
            Array<OneD, NekDouble > qtemp(nTracePts);

            m_traceNormals = Array<OneD, Array<OneD, NekDouble> >(nDim);
            for(i = 0; i < nDim; ++i)
            {
                m_traceNormals[i] = Array<OneD, NekDouble> (nTracePts);
            }
            fields[0]->GetTrace()->GetNormals(m_traceNormals);

            // Extract the physical values of the solution at the boundaries
            fields[var]->ExtractTracePhys(qfield, qtemp);
            
            //cout<<"var = "<< var << endl;
            // Loop on the boundary regions to apply appropriate bcs
            for (i = 0; i < nBndRegions; ++i)
            {
                //cout<<"bcRegion = "<< i << endl;
                // Number of boundary regions related to region 'i'
                nBndEdges = fields[var]->
                GetBndCondExpansions()[i]->GetExpSize();
                
                // Weakly impose bcs by modifying flux values
                for (e = 0; e < nBndEdges; ++e)
                {
                    //cout<<"bcEdge = "<< e << endl;
                    nBndEdgePts = fields[var]->
                    GetBndCondExpansions()[i]->GetExp(e)->GetNumPoints(0);
                    
                    id1 = fields[var]->
                    GetBndCondExpansions()[i]->GetPhys_Offset(e);
                    
                    id2 = fields[0]->GetTrace()->
                    GetPhys_Offset(fields[0]->GetTraceMap()->
                                   GetBndCondTraceToGlobalTraceMap(cnt++));
                    
                    // In case of Dirichlet bcs: 
                    // uflux = gD
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
                        ASSERTL0(false, 
                                 "Neumann bcs not implemented for LDGNS");
                        
                        /*
                        Vmath::Vmul(nBndEdgePts, 
                                    &m_traceNormals[dir][id2], 1, 
                                    &(fields[var]->
                                      GetBndCondExpansions()[i]->
                                      UpdatePhys())[id1], 1, 
                                    &penaltyflux[id2], 1);
                         */
                    }
                    /*
                    for (int bb = 0; bb < nBndEdgePts; ++bb)
                    {
                        cout<<"2nd - bcs = "<<penaltyflux[id2+bb]<<endl;
                    }
                    int num;
                    cin>>num;
                     */
                }
            }
        }
        
    }//end of namespace SolverUtils
}//end of namespace Nektar
