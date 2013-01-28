///////////////////////////////////////////////////////////////////////////////
//
// File CompressibleFlowSystem.cpp
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
// Description: Auxiliary functions for the compressible flow system
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/EquationSystems/CompressibleFlowSystem.h>
#include <LocalRegions/TriExp.h>
#include <LocalRegions/QuadExp.h>
#include <MultiRegions/ExpList.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>

namespace Nektar
{
    string CompressibleFlowSystem::className = 
        SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
            "CompressibleFlowSystem", 
            CompressibleFlowSystem::create, 
            "Auxiliary functions for the compressible flow system.");
    
    CompressibleFlowSystem::CompressibleFlowSystem(
        const LibUtilities::SessionReaderSharedPtr& pSession)
        : UnsteadySystem(pSession)
    {
    }
    
    /**
     * @brief Initialization object for CompressibleFlowSystem class.
     */
    void CompressibleFlowSystem::v_InitObject()
    {
        UnsteadySystem::v_InitObject();
        
        ASSERTL0(m_session->DefinesSolverInfo("UPWINDTYPE"),
                 "No UPWINDTYPE defined in session.");
        
        // Set up locations of velocity vector.
        m_velLoc = Array<OneD, NekDouble>(m_expdim);
        for (int i = 0; i < m_expdim; ++i)
        {
            m_velLoc[i] = i+1;
        }
        
        // Get gamma parameter from session file.
        ASSERTL0(m_session->DefinesParameter("Gamma"),
                 "Compressible flow sessions must define a Gamma parameter.");
        m_session->LoadParameter("Gamma", m_gamma, 1.4);
        
        // Get rho0 parameter from session file.
        ASSERTL0(m_session->DefinesParameter("rho0"),
                 "Compressible flow sessions must define a rho0 parameter.");
        m_session->LoadParameter("rho0", m_rho0, 1.225);
        
        // Get rhou0 parameter from session file.
        ASSERTL0(m_session->DefinesParameter("rhou0"),
                 "Compressible flow sessions must define a rhou0 parameter.");
        m_session->LoadParameter("rhou0", m_rhou0, 0.1225);
        
        // Get rhov0 parameter from session file.
        if (m_expdim == 2 || m_expdim == 3)
        {
            ASSERTL0(m_session->DefinesParameter("rhov0"),
                 "Compressible flow sessions must define a rhov0 parameter.");        
            m_session->LoadParameter("rhov0", m_rhov0, 0.0);
        }
        
        // Get rhow0 parameter from session file.
        if (m_expdim == 3)
        {
            ASSERTL0(m_session->DefinesParameter("rhow0"),
                 "Compressible flow sessions must define a rhow0 parameter.");
            m_session->LoadParameter("rhow0", m_rhow0, 0.0);
        }
        
        // Get E0 parameter from session file.
        ASSERTL0(m_session->DefinesParameter("E0"),
                 "Compressible flow sessions must define a E0 parameter.");
        m_session->LoadParameter("E0", m_E0, 0.149875);
        
        m_session->LoadParameter("GasConstant", m_gasConstant, 287.058);
        
        // Type of advection class to be used
        switch(m_projectionType)
        {
            // Continuous field 
            case MultiRegions::eGalerkin:
            {
                ASSERTL0(false, "Continuous field not supported.");
                break;
            }
            // Discontinuous field 
            case MultiRegions::eDiscontinuous:
            {
                string advName;
                string riemName;
                
                m_session->LoadSolverInfo("AdvectionType", advName, "WeakDG");
                
                m_advection = SolverUtils::GetAdvectionFactory()
                                            .CreateInstance(advName, advName);
                
                m_advection->SetFluxVector(&CompressibleFlowSystem::
                                                         GetFluxVector, this);
                
                m_session->LoadSolverInfo("UpwindType", riemName, "Exact");
                
                m_riemannSolver = SolverUtils::GetRiemannSolverFactory()
                                                    .CreateInstance(riemName);

                m_riemannSolver->AddParam (
                                    "gamma",  
                                    &CompressibleFlowSystem::GetGamma,   this);
                m_riemannSolver->AddScalar(
                                    "velLoc", 
                                    &CompressibleFlowSystem::GetVelLoc,  this);
                m_riemannSolver->AddVector(
                                    "N",
                                    &CompressibleFlowSystem::GetNormals, this);
                
                m_advection->SetRiemannSolver(m_riemannSolver);
                m_advection->InitObject      (m_session, m_fields);
                break;
            }
            default:
            {
                ASSERTL0(false, "Unsupported projection type.");
                break;
            }
        }
    }
    
    /**
     * @brief Destructor for CompressibleFlowSystem class.
     */
    CompressibleFlowSystem::~CompressibleFlowSystem()
    {
        
    }
    
    /**
     * @brief Print out a summary with some relevant information.
     */
    void CompressibleFlowSystem::v_PrintSummary(std::ostream &out)
    {
        UnsteadySystem::v_PrintSummary(out);
    }
    
    /**
     * @brief Wall boundary conditions for compressible flow problems.
     */
    void CompressibleFlowSystem::WallBoundary(
        int                                   bcRegion,
        int                                   cnt, 
        Array<OneD, Array<OneD, NekDouble> > &physarray)
    { 
        int i;
        int nTraceNumPoints = GetTraceTotPoints();
        int nvariables      = physarray.num_elements();
        
        // get physical values of the forward trace
        Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
        for (i = 0; i < nvariables; ++i)
        {
            Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
            m_fields[i]->ExtractTracePhys(physarray[i],Fwd[i]);
        }
        
        // Adjust the physical values of the trace to take 
        // user defined boundaries into account
        int e, id1, id2, npts;
        
        for(e = 0; e < m_fields[0]->
            GetBndCondExpansions()[bcRegion]->GetExpSize(); ++e)
        {
            npts = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                GetExp(e)->GetNumPoints(0);
            id1  = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                GetPhys_Offset(e);
            id2  = m_fields[0]->GetTrace()->GetPhys_Offset(
                        m_fields[0]->GetTraceMap()->
                                    GetBndCondCoeffsToGlobalCoeffsMap(cnt+e));
            
            switch(m_expdim)
            {
                // Special case for 2D
                case 2:
                {
                    Array<OneD, NekDouble> tmp_n(npts);
                    Array<OneD, NekDouble> tmp_t(npts);
                    
                    Vmath::Vmul (npts, 
                                 &Fwd[1][id2], 1, 
                                 &m_traceNormals[0][id2], 1,
                                 &tmp_n[0], 1);
                    
                    Vmath::Vvtvp(npts, 
                                 &Fwd[2][id2], 1, 
                                 &m_traceNormals[1][id2], 1,
                                 &tmp_n[0], 1, 
                                 &tmp_n[0], 1);
                    
                    Vmath::Vmul (npts, 
                                 &Fwd[1][id2], 1, 
                                 &m_traceNormals[1][id2], 1,
                                 &tmp_t[0], 1);
                    
                    Vmath::Vvtvm(npts, 
                                 &Fwd[2][id2], 1, 
                                 &m_traceNormals[0][id2], 1,
                                 &tmp_t[0], 1, &tmp_t[0], 1);
                    
                    // negate the normal flux
                    Vmath::Neg  (npts, tmp_n, 1);		      
                    
                    // rotate back to Cartesian
                    Vmath::Vmul (npts,
                                 &tmp_t[0], 1, 
                                 &m_traceNormals[1][id2], 1,
                                 &Fwd[1][id2], 1);
                    
                    Vmath::Vvtvm(npts, 
                                 &tmp_n[0], 1, 
                                 &m_traceNormals[0][id2], 1,
                                 &Fwd[1][id2], 1, 
                                 &Fwd[1][id2], 1);
                    
                    Vmath::Vmul (npts, 
                                 &tmp_t[0], 1,
                                 &m_traceNormals[0][id2], 1,
                                 &Fwd[2][id2], 1);
                    
                    Vmath::Vvtvp(npts, 
                                 &tmp_n[0], 1,
                                 &m_traceNormals[1][id2], 1,
                                 &Fwd[2][id2], 1,
                                 &Fwd[2][id2], 1);
                    break;
                }
                    
                // For 1D/3D, define: v* = v - (v.n)n so that v*.n = 0
                case 1:
                case 3:
                {
                    Array<OneD,NekDouble> tmp(npts);
                    
                    // Calculate (v.n)
                    for (i = 0; i < m_expdim; ++i)
                    {
                        Vmath::Vvtvp(npts,
                                     &Fwd[1+i][id2], 1, 
                                     &m_traceNormals[i][id2], 1,
                                     &tmp[0], 1,
                                     &tmp[0], 1);
                    }
                    
                    for (i = 0; i < m_expdim; ++i)
                    {
                        Vmath::Vvtvm(npts, 
                                     &tmp[0], 1,
                                     &m_traceNormals[i][id2], 1,
                                     &Fwd[1+i][id2], 1,
                                     &Fwd[1+i][id2], 1);
                        
                        Vmath::Neg  (npts, &Fwd[1+i][id2], 1);
                    }
                    
                    break;
                }
                default:
                    ASSERTL0(false,"Illegal expansion dimension");
                    break;
            }
            
            // copy boundary adjusted values into the boundary expansion
            for (i = 0; i < nvariables; ++i)
            {
                Vmath::Vcopy(npts, &Fwd[i][id2], 1, 
                             &(m_fields[i]->GetBndCondExpansions()[bcRegion]->
                             UpdatePhys())[id1], 1);
            }
        }
    }
    
    /**
     * @brief Simmetry boundary conditions for compressible flow problems.
     */
    void CompressibleFlowSystem::SymmetryBoundary(
        int                                      bcRegion, 
        int                                      cnt, 
        Array<OneD, Array<OneD, NekDouble> >    &physarray)
    {  
        int i;
        int nTraceNumPoints = GetTraceTotPoints();
        int nvariables      = physarray.num_elements();
        
        // get physical values of the forward trace (from exp to phys)
        Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
        for (i = 0; i < nvariables; ++i)
        {
            Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
            m_fields[i]->ExtractTracePhys(physarray[i], Fwd[i]);
        }
        
        int e, id1, id2, npts;
        
        for(e = 0; e < m_fields[0]->
            GetBndCondExpansions()[bcRegion]->GetExpSize(); ++e)
        {
            npts = m_fields[0]->GetBndCondExpansions()[bcRegion]->
            GetExp(e)->GetNumPoints(0);
            id1  = m_fields[0]->GetBndCondExpansions()[bcRegion]->
            GetPhys_Offset(e);
            id2  = m_fields[0]->GetTrace()->GetPhys_Offset(m_fields[0]->
                    GetTraceMap()->GetBndCondCoeffsToGlobalCoeffsMap(cnt+e));
            
            switch(m_expdim)
            {
                case 1:
                {
                    ASSERTL0(false,
                             "1D not yet implemented for Compressible " 
                             "Flow Equations");
                    break;
                }
                case 2:
                {
                    Array<OneD, NekDouble> tmp_t(npts);
                    
                    Vmath::Vmul(npts, 
                                &Fwd[1][id2], 1, 
                                &m_traceNormals[1][id2], 1, 
                                &tmp_t[0], 1);
                    
                    Vmath::Vvtvm(npts, 
                                 &Fwd[2][id2], 1, 
                                 &m_traceNormals[0][id2], 1, 
                                 &tmp_t[0], 1, 
                                 &tmp_t[0], 1);
                    
                    Array<OneD, NekDouble> tmp_n(npts, 0.0);
                    
                    // rotate back to Cartesian
                    Vmath::Vmul(npts, 
                                &tmp_t[0], 1, 
                                &m_traceNormals[1][id2], 1, 
                                &Fwd[1][id2], 1);
                    
                    Vmath::Vvtvm(npts, 
                                 &tmp_n[0], 1,
                                 &m_traceNormals[0][id2], 1, 
                                 &Fwd[1][id2], 1, 
                                 &Fwd[1][id2], 1);
                    
                    Vmath::Vmul(npts, 
                                &tmp_t[0], 1, 
                                &m_traceNormals[0][id2], 1,
                                &Fwd[2][id2], 1);
                    Vmath::Vvtvp(npts, 
                                 &tmp_n[0], 1, 
                                 &m_traceNormals[1][id2], 1, 
                                 &Fwd[2][id2], 1, 
                                 &Fwd[2][id2], 1);
                    break;
                }
                case 3:
                {
                    ASSERTL0(false,
                             "3D not yet implemented for Compressible "
                             "Flow Equations");
                    break;
                }
                default:
                {
                    ASSERTL0(false, "Illegal expansion dimension");
                }
            }
            
            // copy boundary adjusted values into the boundary expansion
            for (i = 0; i < nvariables; ++i)
            {
                Vmath::Vcopy(npts, 
                             &Fwd[i][id2], 1, 
                             &(m_fields[i]->GetBndCondExpansions()[bcRegion]->
                               UpdatePhys())[id1], 1);	
            }
        }
    }
    
    /**
     * @brief Inflow boundary conditions for compressible flow problems based
     * on the eigenvalues of inviscid term of the compressible Euler equation.
     */
    void CompressibleFlowSystem::InflowCFE(
        int                                   bcRegion,
        int                                   cnt, 
        Array<OneD, Array<OneD, NekDouble> > &physarray)
    { 
        int i, e;
        int id1, id2;
        int npts, pnt;
        
        int nTraceNumPoints = GetTraceTotPoints();
        int nvariables      = physarray.num_elements();
        
        NekDouble cPlus, cMinus;
        NekDouble rPlus, rMinus;
        NekDouble Vnorm, Vdelta;
        NekDouble soundSpeedCorrection;
        NekDouble uCorrection;
        NekDouble vCorrection;
        NekDouble pressureCorrection;
        NekDouble entropyCorrection;
        
        NekDouble rhoFinalCorrection;
        NekDouble rhouFinalCorrection;
        NekDouble rhovFinalCorrection;
        NekDouble energyFinalCorrection;
        
        // Forward trace arrays
        Array<OneD, NekDouble > pressurefwd  (nTraceNumPoints, 0.0);
        Array<OneD, NekDouble > soundSpeedfwd(nTraceNumPoints, 0.0);
        Array<OneD, NekDouble > Machfwd      (nTraceNumPoints, 0.0);
        Array<OneD, NekDouble > Velfwd       (nTraceNumPoints, 0.0);
        Array<OneD, NekDouble > Vnfwd        (nTraceNumPoints, 0.0);
        
        // Boundary values
        NekDouble pressureB;
        NekDouble soundSpeedB;
        NekDouble MachB;
        NekDouble VelB;
        NekDouble VnB;

        // Get physical values of the forward trace
        Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
        for (i = 0; i < nvariables; ++i)
        {
            Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
            m_fields[i]->ExtractTracePhys(physarray[i], Fwd[i]);
        }
                
        // Switch on the dimensions of the problem
        switch(m_expdim)
        {
            // 1D case 
            case 1:
            {
                ASSERTL0(false, "1D InflowCFE boundary condition" 
                                "not implemented yet")
                break;
            }
            // 2D case
            case 2:
            {
                // Get forward trace variables
                for(i = 0; i < nTraceNumPoints; i++)
                {
                    // Forward trace normal velocity
                    Vnfwd[i] = Fwd[1][i] / Fwd[0][i] * m_traceNormals[0][i] + 
                               Fwd[2][i] / Fwd[0][i] * m_traceNormals[1][i];
                    
                    // Forward trace pressure
                    pressurefwd[i] = (m_gamma - 1) * (Fwd[3][i] - 
                                    0.5 * (Fwd[1][i] * Fwd[1][i] / Fwd[0][i] + 
                                           Fwd[2][i] * Fwd[2][i] / Fwd[0][i]));
                        
                    // Forward trace speed of sound
                    soundSpeedfwd[i] = sqrt(m_gamma * pressurefwd[i] / 
                                            Fwd[0][i]);
                }
                
                // Boundary normal velocity
                VnB = sqrt((m_rhou0 * m_rhou0) / (m_rho0 * m_rho0) + 
                           (m_rhov0 * m_rhov0) / (m_rho0 * m_rho0));
                
                // Boundary pressure
                pressureB = (m_gamma - 1) * (m_E0 - 0.5 * 
                                             ((m_rhou0 * m_rhou0) / (m_rho0) +  
                                              (m_rhov0 * m_rhov0) / (m_rho0)));
                
                // Boundary speed of sound
                soundSpeedB   = sqrt(m_gamma * pressureB / m_rho0);
                 
                // Forward trace Mach number
                Vmath::Vdiv(nTraceNumPoints, Velfwd, 1, 
                            soundSpeedfwd, 1, Machfwd, 1);
                
                // Boundary Mach number
                MachB = VnB / soundSpeedB;
                
                // Loop on bcRegions
                for(e = 0; e < m_fields[0]->
                    GetBndCondExpansions()[bcRegion]->GetExpSize(); ++e)
                {
                    npts = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                        GetExp(e)->GetNumPoints(0);
                    
                    id1  = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                        GetPhys_Offset(e) ;
                    
                    id2  = m_fields[0]->GetTrace()->
                        GetPhys_Offset(m_fields[0]->GetTraceMap()->
                                       GetBndCondTraceToGlobalTraceMap(cnt++));
                                            
                    // Loop on the points of the bcRegion
                    for (i = 0; i < npts; i++)
                    {
                        pnt = id2 + i;
                        if (Machfwd[pnt] < 0.99)
                        {
                            // + Characteristic from boundary
                            cPlus = sqrt(m_gamma * pressureB / m_rho0);
                            rPlus = VnB + 2.0 * cPlus / (m_gamma - 1);
                            
                            // - Characteristic from inside
                            cMinus = sqrt(m_gamma * pressurefwd[pnt] / 
                                     Fwd[0][pnt]);
                            rMinus = Vnfwd[pnt] - 2.0 * cMinus / (m_gamma - 1);
                            
                            // Corrections for the boundary coditions
                            Vnorm  = 0.5 * (rPlus + rMinus);
                            Vdelta = Vnorm - VnB;
                            
                            // Speed of sound correction
                            soundSpeedCorrection = 0.25 * (m_gamma - 1) * 
                                                          (rPlus - rMinus);
                            
                            // Velocity corrections
                            uCorrection = (VnB + Vdelta)*m_traceNormals[0][pnt];
                            vCorrection = (VnB + Vdelta)*m_traceNormals[1][pnt];
                            
                            // Entropy correction
                            entropyCorrection = pressureB / 
                                                (pow(m_rho0, m_gamma));
                            
                            // rho final correction
                            rhoFinalCorrection = pow((soundSpeedCorrection * 
                                                      soundSpeedCorrection) / 
                                                     (m_gamma * 
                                                      entropyCorrection), 
                                                      1.0 / (m_gamma - 1));
                            
                            // Pressure correction
                            pressureCorrection = rhoFinalCorrection * 
                                                 soundSpeedCorrection * 
                                                 soundSpeedCorrection / m_gamma;
                            
                            // rhou final correction
                            rhouFinalCorrection = rhoFinalCorrection * 
                                                    uCorrection;
                            
                            // rhov final correction
                            rhovFinalCorrection = rhoFinalCorrection * 
                                                    vCorrection;
                            
                            energyFinalCorrection = pressureCorrection / 
                                                    (m_gamma - 1) + 0.5 * 
                                                    (rhouFinalCorrection * 
                                                     uCorrection + 
                                                     rhovFinalCorrection * 
                                                     vCorrection);
                            
                            // Imposing characteristic farfield bcs
                            (m_fields[0]->GetBndCondExpansions()[bcRegion]->
                             UpdatePhys())[id1+i] = 
                                2.0 * rhoFinalCorrection - Fwd[0][pnt];
                            (m_fields[1]->GetBndCondExpansions()[bcRegion]->
                             UpdatePhys())[id1+i] = 
                                2.0 * rhouFinalCorrection - Fwd[1][pnt];
                            (m_fields[2]->GetBndCondExpansions()[bcRegion]->
                             UpdatePhys())[id1+i] = 
                                2.0 * rhovFinalCorrection - Fwd[2][pnt];
                            (m_fields[3]->GetBndCondExpansions()[bcRegion]->
                             UpdatePhys())[id1+i] = 
                                2.0 * energyFinalCorrection - Fwd[3][pnt];
                        }
                        else
                        {
                            (m_fields[0]->GetBndCondExpansions()[bcRegion]->
                             UpdatePhys())[id1+i] = m_rho0;
                            (m_fields[1]->GetBndCondExpansions()[bcRegion]->
                             UpdatePhys())[id1+i] = m_rho0 * VnB;
                            (m_fields[2]->GetBndCondExpansions()[bcRegion]->
                             UpdatePhys())[id1+i] = m_rho0 * VnB;
                            (m_fields[3]->GetBndCondExpansions()[bcRegion]->
                             UpdatePhys())[id1+i] = 
                                pressureB / (m_gamma - 1.0) + 
                                0.5 * m_rhou0 * m_rhou0 / m_rho0;
                        }
                    }
                }
                break;
            }
            // 3D case
            case 3:
            {
                ASSERTL0(false, "3D InflowCFE boundary condition" 
                                "not implemented yet")
                break;
            }
            default:
            {
                ASSERTL0(false, "Illegal expansion dimension");
                break;
            }
        }
    }

    /**
     * @brief Outflow boundary conditions for compressible flow problems based
     * on the eigenvalues of inviscid term of the compressible Euler equation.
     */
    void CompressibleFlowSystem::OutflowCFE(
        int                                   bcRegion,
        int                                   cnt, 
        Array<OneD, Array<OneD, NekDouble> > &physarray)
    { 
        int i, e;
        int id1, id2;
        int npts, pnt;
        
        int nTraceNumPoints = GetTraceTotPoints();
        int nvariables      = physarray.num_elements();
        
        NekDouble cPlus, cMinus;
        NekDouble rPlus, rMinus;
        NekDouble Vnorm, Vdelta;
        NekDouble soundSpeedCorrection;
        NekDouble uCorrection;
        NekDouble vCorrection;
        NekDouble pressureCorrection;
        NekDouble entropyCorrection;
        
        NekDouble rhoFinalCorrection;
        NekDouble rhouFinalCorrection;
        NekDouble rhovFinalCorrection;
        NekDouble energyFinalCorrection;
        
        // Forward trace arrays
        Array<OneD, NekDouble > pressurefwd  (nTraceNumPoints, 0.0);
        Array<OneD, NekDouble > soundSpeedfwd(nTraceNumPoints, 0.0);
        Array<OneD, NekDouble > Machfwd      (nTraceNumPoints, 0.0);
        Array<OneD, NekDouble > Velfwd       (nTraceNumPoints, 0.0);
        Array<OneD, NekDouble > Vnfwd        (nTraceNumPoints, 0.0);
        
        // Boundary values
        NekDouble pressureB;
        NekDouble soundSpeedB;
        NekDouble MachB;
        NekDouble VelB;
        NekDouble VnB;
        
        // Get physical values of the forward trace
        Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
        for (i = 0; i < nvariables; ++i)
        {
            Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
            m_fields[i]->ExtractTracePhys(physarray[i], Fwd[i]);
        }
        
        // Switch on the dimensions of the problem
        switch(m_expdim)
        {
            // 1D case 
            case 1:
            {
                ASSERTL0(false, "1D InflowCFE boundary condition" 
                         "not implemented yet")
                break;
            }
            // 2D case
            case 2:
            {
                // Get forward trace variables
                for(i = 0; i < nTraceNumPoints; i++)
                {
                    // Forward trace normal velocity
                    Vnfwd[i] = Fwd[1][i] / Fwd[0][i] * m_traceNormals[0][i] + 
                    Fwd[2][i] / Fwd[0][i] * m_traceNormals[1][i];
                    
                    // Forward trace pressure
                    pressurefwd[i] = 
                        (m_gamma - 1) * (Fwd[3][i] - 
                                    0.5 * (Fwd[1][i] * Fwd[1][i] / Fwd[0][i] + 
                                           Fwd[2][i] * Fwd[2][i] / Fwd[0][i]));
                    
                    // Forward trace speed of sound
                    soundSpeedfwd[i] = sqrt(m_gamma * pressurefwd[i] / 
                                            Fwd[0][i]);
                }
                
                // Boundary normal velocity
                VnB = sqrt((m_rhou0 * m_rhou0) / (m_rho0 * m_rho0) + 
                           (m_rhov0 * m_rhov0) / (m_rho0 * m_rho0));
                
                // Boundary pressure
                pressureB = (m_gamma - 1) * (m_E0 - 0.5 * 
                                             ((m_rhou0 * m_rhou0) / (m_rho0) +  
                                              (m_rhov0 * m_rhov0) / (m_rho0)));
                
                // Boundary speed of sound
                soundSpeedB   = sqrt(m_gamma * pressureB / m_rho0);
                
                // Forward trace Mach number
                Vmath::Vdiv(nTraceNumPoints, Velfwd, 1, 
                            soundSpeedfwd, 1, Machfwd, 1);
                
                // Boundary Mach number
                MachB = VnB / soundSpeedB;
                
                // Loop on bcRegions
                for(e = 0; e < m_fields[0]->
                    GetBndCondExpansions()[bcRegion]->GetExpSize(); ++e)
                {
                    npts = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                    GetExp(e)->GetNumPoints(0);
                    
                    id1  = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                    GetPhys_Offset(e) ;
                    
                    id2  = m_fields[0]->GetTrace()->
                    GetPhys_Offset(m_fields[0]->GetTraceMap()->
                                   GetBndCondTraceToGlobalTraceMap(cnt++));
                    
                    // Loop on the points of the bcRegion
                    for (i = 0; i < npts; i++)
                    {
                        pnt = id2 + i;
                        
                        // Subsonic flow charcteristic bcs
                        if (Machfwd[pnt] < 0.99)
                        {
                            // + Characteristic from boundary
                            cPlus = sqrt(m_gamma * pressurefwd[pnt] / 
                                    Fwd[0][pnt]);
                            rPlus = Vnfwd[pnt] + 2.0 * cPlus / (m_gamma - 1);
                            
                            // - Characteristic from inside
                            cMinus = sqrt(m_gamma * pressureB / m_rho0);
                            rMinus = VnB - 2.0 * cMinus / (m_gamma - 1);
                            
                            // Corrections for the boundary coditions
                            Vnorm   = 0.5 * (rPlus + rMinus);
                            Vdelta  = Vnorm - Vnfwd[pnt];

                            // Speed of sound correction
                            soundSpeedCorrection = 0.25 * (m_gamma - 1) * 
                            (rPlus - rMinus);
                            
                            // Velocity corrections
                            uCorrection = Fwd[1][pnt]/Fwd[0][pnt] + 
                            Vdelta * m_traceNormals[0][pnt];
                            vCorrection = Fwd[2][pnt]/Fwd[0][pnt] + 
                            Vdelta * m_traceNormals[1][pnt];
                            
                            // Entropy correction
                            entropyCorrection = pressurefwd[pnt] / 
                            (pow(Fwd[0][pnt], m_gamma));
                            
                            // Eq.1) Rho final correction
                            rhoFinalCorrection = pow((soundSpeedCorrection * 
                                                      soundSpeedCorrection) / 
                                                     (m_gamma * 
                                                      entropyCorrection), 
                                                     1.0 / (m_gamma - 1));
                            
                            // Pressure correction
                            pressureCorrection = rhoFinalCorrection * 
                            soundSpeedCorrection * 
                            soundSpeedCorrection / m_gamma;
                            
                            // Eq.2) Rhou final correction
                            rhouFinalCorrection = rhoFinalCorrection * 
                            uCorrection;
                            
                            // Eq.3) Rhov final correction
                            rhovFinalCorrection = rhoFinalCorrection * 
                            vCorrection;
                            
                            // Eq.4) E final correction
                            energyFinalCorrection = pressureCorrection / 
                            (m_gamma - 1) + 0.5 * 
                            (rhouFinalCorrection * 
                             uCorrection + 
                             rhovFinalCorrection * 
                             vCorrection);
                            
                            // Imposing characteristic farfield bcs
                            (m_fields[0]->GetBndCondExpansions()[bcRegion]->
                             UpdatePhys())[id1+i] = 
                            2.0 * rhoFinalCorrection - Fwd[0][pnt];
                            (m_fields[1]->GetBndCondExpansions()[bcRegion]->
                             UpdatePhys())[id1+i] = 
                            2.0 * rhouFinalCorrection - Fwd[1][pnt];
                            (m_fields[2]->GetBndCondExpansions()[bcRegion]->
                             UpdatePhys())[id1+i] = 
                            2.0 * rhovFinalCorrection - Fwd[2][pnt];
                            (m_fields[3]->GetBndCondExpansions()[bcRegion]->
                             UpdatePhys())[id1+i] = 
                            2.0 * energyFinalCorrection - Fwd[3][pnt];
                        }
                        // Supersonic flow characteristic bcs
                        else
                        {
                            (m_fields[0]->GetBndCondExpansions()[bcRegion]->
                             UpdatePhys())[id1+i] = Fwd[0][pnt];
                            (m_fields[1]->GetBndCondExpansions()[bcRegion]->
                             UpdatePhys())[id1+i] = Fwd[1][pnt];
                            (m_fields[2]->GetBndCondExpansions()[bcRegion]->
                             UpdatePhys())[id1+i] = Fwd[2][pnt];
                            (m_fields[3]->GetBndCondExpansions()[bcRegion]->
                             UpdatePhys())[id1+i] = Fwd[3][pnt];
                        }
                    }
                }
                break;
            }
                // 3D case
            case 3:
            {
                ASSERTL0(false, "3D InflowCFE boundary condition" 
                         "not implemented yet")
                break;
            }
            default:
            {
                ASSERTL0(false, "Illegal expansion dimension");
                break;
            }
        }
    }

    /**
     * @brief Return the flux vector for the compressible Euler equations.
     * 
     * @param i           Component of the flux vector to calculate.
     * @param physfield   Fields.
     * @param flux        Resulting flux.
     */
    void CompressibleFlowSystem::GetFluxVector(
        const int                                   i,
        const Array<OneD, Array<OneD, NekDouble> > &physfield,
              Array<OneD, Array<OneD, NekDouble> > &flux)
    {
        int j, nq = m_fields[0]->GetTotPoints();
        
        if (i == 0)
        {
            // Flux vector for the rho equation.
            for (j = 0; j < m_expdim; ++j)
            {
                Vmath::Vcopy(nq, physfield[j+1], 1, flux[j], 1);
            }
        } 
        else if (i >= 1 && i <= m_expdim)
        {
            // Flux vector for the velocity fields.
            Array<OneD, NekDouble> pressure(nq);
            GetPressure(physfield, pressure);
            
            for (j = 0; j < m_expdim; ++j)
            {
                Vmath::Vmul(nq, physfield[j+1], 1, physfield[i], 1, flux[j], 1);
                Vmath::Vdiv(nq, flux[j], 1, physfield[0], 1, flux[j], 1);
            }
            
            // Add pressure to appropriate field
            Vmath::Vadd(nq, flux[i-1], 1, pressure, 1, flux[i-1], 1);
        }
        else if (i == m_expdim+1) 
        {
            // Flux vector for the total energy field.
            Array<OneD, NekDouble> pressure(nq);
            GetPressure(physfield, pressure);
            Vmath::Vadd(nq, physfield[m_expdim+1], 1, pressure, 1, pressure, 1);
            
            for (j = 0; j < m_expdim; ++j)
            {
                Vmath::Vdiv(nq, physfield[j+1], 1, physfield[0], 1, flux[j], 1);
                Vmath::Vmul(nq, flux[j], 1, pressure, 1, flux[j], 1);
            }
        }
        else
        {
            ASSERTL0(false, "Invalid vector index.");
        }
    }

    /**
     * @brief Calculate the pressure field \f$ p =
     * (\gamma-1)(E-\frac{1}{2}\rho\| \mathbf{v} \|^2) \f$ assuming an ideal 
     * gas law.
     * 
     * @param physfield  Input momentum.
     * @param pressure   Computed pressure field.
     */
    void CompressibleFlowSystem::GetPressure(
        const Array<OneD, const Array<OneD, NekDouble> > &physfield,
              Array<OneD,                   NekDouble>   &pressure)
    {
        int       npts  = m_fields[0]->GetTotPoints();
        double    alpha = -0.5;
        
        // Calculate ||rho v||^2.
        Vmath::Zero(npts, pressure, 1);
        for (int i = 0; i < m_expdim; ++i)
        {
            Vmath::Vvtvp(npts, physfield[1+i], 1, physfield[1+i], 1, 
                               pressure,       1, pressure,       1);
        }
        // Divide by rho to get rho*||v||^2.
        Vmath::Vdiv (npts, pressure, 1, physfield[0], 1, pressure, 1);
        // pressure <- E - 0.5*pressure
        Vmath::Svtvp(npts,     alpha, 
                     pressure, 1, physfield[m_expdim+1], 1, pressure, 1);
        // Multiply by (gamma-1).
        Vmath::Smul (npts, m_gamma-1, pressure, 1, pressure, 1);
    }
    
    /**
     * @brief Compute the velocity field \f$ \mathbf{v} \f$ given the momentum
     * \f$ \rho\mathbf{v} \f$.
     * 
     * @param physfield  Momentum field.
     * @param velocity   Velocity field.
     */
    void CompressibleFlowSystem::GetVelocityVector(
        const Array<OneD, Array<OneD, NekDouble> > &physfield,
              Array<OneD, Array<OneD, NekDouble> > &velocity)
    {
        const int npts = m_fields[0]->GetTotPoints();
        
        for (int i = 0; i < m_expdim; ++i)
        {
            Vmath::Vdiv(npts, physfield[1+i], 1, physfield[0], 1, 
                              velocity[i],    1);
        }
    }
  
    /**
     * @brief Compute the temperature \f$ T = p/\rho R \f$.
     * 
     * @param physfield    Input physical field.
     * @param pressure     The pressure field corresponding to physfield.
     * @param temperature  The resulting temperature \f$ T \f$.
     */
    void CompressibleFlowSystem::GetTemperature(
        Array<OneD, Array<OneD, NekDouble> > &physfield,
        Array<OneD,             NekDouble  > &pressure,
        Array<OneD,             NekDouble  > &temperature)
    {
        for (int i = 0; i < m_fields[0]->GetTotPoints(); ++i)
        {
            temperature[i] = pressure[i]/(physfield[0][i]*m_gasConstant);
        }
    }
    
    /**
     * @brief Compute the sound speed \f$ c = sqrt(\gamma p/\rho) \f$.
     * 
     * @param physfield    Input physical field.
     * @param pressure     The pressure field corresponding to physfield.
     * @param soundspeed   The resulting sound speed \f$ c \f$.
     */
    void CompressibleFlowSystem::GetSoundSpeed(
        const Array<OneD, Array<OneD, NekDouble> > &physfield,
              Array<OneD,             NekDouble  > &pressure,
              Array<OneD,             NekDouble  > &soundspeed)
    {
        for (int i = 0; i < m_fields[0]->GetTotPoints(); ++i)
        {
            soundspeed[i] = sqrt(m_gamma*pressure[i]/physfield[0][i]);
        }
    }
    
    /**
     * @brief Compute the mach number \f$ M = \| \mathbf{v} \|^2 / c \f$.
     * 
     * @param physfield    Input physical field.
     * @param soundfield   The speed of sound corresponding to physfield.
     * @param mach         The resulting mach number \f$ M \f$.
     */
    void CompressibleFlowSystem::GetMach(
        Array<OneD, Array<OneD, NekDouble> > &physfield,
        Array<OneD,             NekDouble  > &soundspeed,
        Array<OneD,             NekDouble  > &mach)
    {
        const int npts = m_fields[0]->GetTotPoints();

        Vmath::Zero(npts, mach, 1);
        
        for (int i = 0; i < m_expdim; ++i)
        {
            Vmath::Vvtvp(npts, physfield[1+i], 1, physfield[1+i], 1, 
                               mach,           1, mach,           1);
        }
        Vmath::Vdiv(npts, mach, 1, physfield[0], 1, mach, 1);
        Vmath::Vdiv(npts, mach, 1, physfield[0], 1, mach, 1);
        Vmath::Vdiv(npts, mach, 1, soundspeed,   1, mach, 1);
    }
    
    /**
     * @brief Calculate the maximum timestep subject to CFL restrictions.
     */
    NekDouble CompressibleFlowSystem::v_GetTimeStep(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray)
    { 
        int i, n;
        int nvariables     = m_fields.num_elements();
        int nTotQuadPoints = GetTotPoints();
        int nElements      = m_fields[0]->GetExpSize(); 
        const Array<OneD, int> ExpOrder = GetNumExpModesPerExp();
        
        Array<OneD, NekDouble> tstep      (nElements, 0.0);
        Array<OneD, NekDouble> stdVelocity(nElements);
        
        // Get standard velocity to compute the time-step limit
        GetStdVelocity(inarray, stdVelocity);
        
        // Factors to compute the time-step limit
        NekDouble minLength;        
        NekDouble alpha   = MaxTimeStepEstimator();
        NekDouble cLambda = 0.2; // Spencer book-317
        
        // Loop over elements to compute the time-step limit for each element
        for(n = 0; n < nElements; ++n)
        {
            int npoints = m_fields[0]->GetExp(n)->GetTotPoints();
            Array<OneD, NekDouble> one2D(npoints, 1.0);
            NekDouble Area = m_fields[0]->GetExp(n)->Integral(one2D);
            
            if (boost::dynamic_pointer_cast<LocalRegions::TriExp>(
                    m_fields[0]->GetExp(n)))
            {
                minLength = 2.0 * sqrt(Area);
            }
            
            else if (boost::dynamic_pointer_cast<LocalRegions::QuadExp>(
                         m_fields[0]->GetExp(n)))
            {
                minLength = sqrt(Area);
            }
            
            tstep[n] = m_cflSafetyFactor * alpha * minLength 
                     / (stdVelocity[n] * cLambda 
                        * (ExpOrder[n] - 1) * (ExpOrder[n] - 1));
        }
        
        // Get the minimum time-step limit and return the time-step
        NekDouble TimeStep = Vmath::Vmin(nElements, tstep, 1);
        m_comm->AllReduce(TimeStep, LibUtilities::ReduceMin);
        return TimeStep;
    }
    
    /**
     * @brief Compute the advection velocity in the standard space 
     * for each element of the expansion.
     * 
     * @param inarray    Momentum field.
     * @param stdV       Standard velocity field.
     */
    void CompressibleFlowSystem::GetStdVelocity(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,                   NekDouble>   &stdV)
    {
        int nTotQuadPoints = GetTotPoints();
        int n_element      = m_fields[0]->GetExpSize();
        int npts           = 0;

        // Getting the velocity vector on the 2D normal space
        Array<OneD, Array<OneD, NekDouble> > velocity   (m_expdim);
        Array<OneD, Array<OneD, NekDouble> > stdVelocity(m_expdim);
        Array<OneD, NekDouble>               pressure   (nTotQuadPoints);
        Array<OneD, NekDouble>               soundspeed (nTotQuadPoints);
        
        // Zero output array
        Vmath::Zero(stdV.num_elements(), stdV, 1);
        
        for (int i = 0; i < m_expdim; ++i)
        {
            velocity   [i] = Array<OneD, NekDouble>(nTotQuadPoints);
            stdVelocity[i] = Array<OneD, NekDouble>(nTotQuadPoints, 0.0);
        }
        GetVelocityVector(inarray, velocity);
        GetPressure      (inarray, pressure);
        GetSoundSpeed    (inarray, pressure, soundspeed);
        
        for(int el = 0; el < n_element; ++el)
        { 
            const Array<OneD, const NekDouble> &jac  = 
                m_fields[0]->GetExp(el)->GetGeom()->GetJac();
            const Array<TwoD, const NekDouble> &gmat = 
                m_fields[0]->GetExp(el)->GetGeom()->GetGmat();
            
            int nqtot = m_fields[0]->GetExp(el)->GetTotPoints();
            
            if(m_fields[0]->GetExp(el)->GetGeom()->GetGtype() == 
                   SpatialDomains::eDeformed)
            {
                // d xi/ dx = gmat = 1/J * d x/d xi
                for (int i = 0; i < m_expdim; ++i)
                {
                    Vmath::Zero(nqtot, stdVelocity[i], 1);
                    for (int j = 0; j < m_expdim; ++j)
                    {
                        Vmath::Vvtvp(nqtot, gmat[m_expdim*j+i], 1, velocity[j], 
                                     1, stdVelocity[i], 1, stdVelocity[i], 1);
                    }
                }
            }
            else
            {
                for (int i = 0; i < m_expdim; ++i)
                {
                    Vmath::Zero(nqtot, stdVelocity[i], 1);
                    for (int j = 0; j < m_expdim; ++j)
                    {
                        Vmath::Svtvp(nqtot, gmat[m_expdim*j+i][0], velocity[j], 
                                     1, stdVelocity[i], 1, stdVelocity[i], 1);
                    }
                }
            }
            
            for(int i = 0; i < nqtot; ++i)
            {
                NekDouble pntVelocity = 0.0;
                for (int j = 0; j < m_expdim; ++j)
                {
                    pntVelocity += stdVelocity[j][i]*stdVelocity[j][i];
                }
                pntVelocity = sqrt(pntVelocity) + soundspeed[npts];
                if (pntVelocity > stdV[el])
                {
                    stdV[el] = pntVelocity;
                }
                npts++;
            }
        }
    }
    
    /**
     * @brief Set the denominator to compute the time step when a cfl  
     * control is employed. This function is no longer used but is still 
     * here for being utilised in the future.
     *
     * @param n   Order of expansion element by element.
     */
    NekDouble CompressibleFlowSystem::GetStabilityLimit(int n)
    {
        if (n > 20)
        {
            ASSERTL0(false,
                     "Illegal modes dimension for CFL calculation "
                     "(P has to be less then 20)");
        }
		
        NekDouble CFLDG[21] = {  2.0000,   6.0000,  11.8424,  19.1569, 
                                27.8419,  37.8247,  49.0518,  61.4815, 
                                75.0797,  89.8181, 105.6700, 122.6200,
                               140.6400, 159.7300, 179.8500, 201.0100,
                               223.1800, 246.3600, 270.5300, 295.6900,
                               321.8300}; //CFLDG 1D [0-20]
        
        NekDouble CFLCG[2]  = {1.0, 1.0};
        NekDouble CFL;
		
        if (m_projectionType == MultiRegions::eDiscontinuous)
        {
            CFL = CFLDG[n];
        }
        else 
        {
            ASSERTL0(false,
                     "Continuos Galerkin stability coefficients "
                     "not introduced yet.");
        }
		
        return CFL;
    }
	
    /**
     * @brief Compute the vector of denominators to compute the time step  
     * when a cfl control is employed. This function is no longer used but 
     * is still here for being utilised in the future.
     *
     * @param ExpOrder   Order of expansion element by element.
     */
    Array<OneD, NekDouble> CompressibleFlowSystem::GetStabilityLimitVector(
        const Array<OneD,int> &ExpOrder)
    {
        int i;
        Array<OneD,NekDouble> returnval(m_fields[0]->GetExpSize(), 0.0);
        for (i =0; i<m_fields[0]->GetExpSize(); i++)
        {
            returnval[i] = GetStabilityLimit(ExpOrder[i]);
        }
        return returnval;
    }
}
