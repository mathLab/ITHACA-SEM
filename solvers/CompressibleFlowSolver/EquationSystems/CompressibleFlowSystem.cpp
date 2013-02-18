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
#include <LocalRegions/HexExp.h>
#include <MultiRegions/ExpList.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <solvers/ADRSolver/EquationSystems/UnsteadyDiffusion.cpp>
#include <solvers/ADRSolver/EquationSystems/UnsteadyDiffusion.h>

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
        
        // Get E0 parameter from session file.
        ASSERTL0(m_session->DefinesParameter("pInf"),
                 "Compressible flow sessions must define a pInf parameter.");
        m_session->LoadParameter("pInf", m_pInf, 101325);
        
        // Get rhoInf parameter from session file.
        ASSERTL0(m_session->DefinesParameter("rhoInf"),
                 "Compressible flow sessions must define a rhoInf parameter.");
        m_session->LoadParameter("rhoInf", m_rhoInf, 1.225);
        
        // Get uInf parameter from session file.
        ASSERTL0(m_session->DefinesParameter("uInf"),
                 "Compressible flow sessions must define a uInf parameter.");
        m_session->LoadParameter("uInf", m_uInf, 0.1);
        
        // Get vInf parameter from session file.
        if (m_expdim == 2 || m_expdim == 3)
        {
            ASSERTL0(m_session->DefinesParameter("vInf"),
                     "Compressible flow sessions must define a vInf parameter"
                     "for 2D/3D problems.");
            m_session->LoadParameter("vInf", m_vInf, 0.0);
        }
        
        // Get wInf parameter from session file.
        if (m_expdim == 3)
        {
            ASSERTL0(m_session->DefinesParameter("wInf"),
                     "Compressible flow sessions must define a wInf parameter"
                     "for 3D problems.");
            m_session->LoadParameter("wInf", m_wInf, 0.0);
        }
                
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
                string diffName;
                string riemName;

                m_session->LoadSolverInfo("AdvectionType", advName, "WeakDG");
                m_session->LoadSolverInfo("DiffusionType", diffName, "LDGNS");

                m_advection = SolverUtils::GetAdvectionFactory()
                                            .CreateInstance(advName, advName);
                
                m_diffusion = SolverUtils::GetDiffusionFactory()
                                            .CreateInstance(diffName, diffName);

                m_advection->SetFluxVector(&CompressibleFlowSystem::
                                            GetFluxVector, this);
                
                m_diffusion->SetFluxVectorNS(&CompressibleFlowSystem::
                                             GetViscousFluxVector, this);

                m_session->LoadSolverInfo("UpwindType", riemName, "Exact");
                
                m_riemannSolver = SolverUtils::GetRiemannSolverFactory()
                                            .CreateInstance(riemName);
                
                m_riemannSolverLDG = SolverUtils::GetRiemannSolverFactory()
                                                .CreateInstance("Upwind");

                m_riemannSolver->AddParam (
                                    "gamma",  
                                    &CompressibleFlowSystem::GetGamma,   this);
                m_riemannSolver->AddScalar(
                                    "velLoc", 
                                    &CompressibleFlowSystem::GetVelLoc,  this);
                m_riemannSolver->AddVector(
                                    "N",
                                    &CompressibleFlowSystem::GetNormals, this);
                
                m_riemannSolverLDG->AddParam (
                                    "gamma",  
                                    &CompressibleFlowSystem::GetGamma,   this);
                m_riemannSolverLDG->AddScalar(
                                    "velLoc", 
                                    &CompressibleFlowSystem::GetVelLoc,  this);
                m_riemannSolverLDG->AddVector(
                                    "N",
                                    &CompressibleFlowSystem::GetNormals, this);
                
                
                m_advection->SetRiemannSolver   (m_riemannSolver);
                m_diffusion->SetRiemannSolver   (m_riemannSolverLDG);
                m_advection->InitObject         (m_session, m_fields);
                m_diffusion->InitObject         (m_session);
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
            m_fields[i]->ExtractTracePhys(physarray[i], Fwd[i]);
        }
        
        // Adjust the physical values of the trace to take 
        // user defined boundaries into account
        int e, id1, id2, npts;
        
        for (e = 0; e < m_fields[0]->GetBndCondExpansions()[bcRegion]
                 ->GetExpSize(); ++e)
        {
            npts = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                GetExp(e)->GetTotPoints();
            id1  = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                GetPhys_Offset(e);
            id2  = m_fields[0]->GetTrace()->GetPhys_Offset(
                        m_fields[0]->GetTraceMap()->
                                    GetBndCondCoeffsToGlobalCoeffsMap(cnt+e));
            
            // For 2D/3D, define: v* = v - 2(v.n)n
            Array<OneD,NekDouble> tmp(npts, 0.0);

            // Calculate (v.n)
            for (i = 0; i < m_expdim; ++i)
            {
                Vmath::Vvtvp(npts,
                             &Fwd[1+i][id2], 1,
                             &m_traceNormals[i][id2], 1,
                             &tmp[0], 1,
                             &tmp[0], 1);
            }

            // Calculate 2.0(v.n)
            Vmath::Smul(npts, -2.0, &tmp[0], 1, &tmp[0], 1);
            
            // Calculate v* = v - 2.0(v.n)n
            for (i = 0; i < m_expdim; ++i)
            {
                Vmath::Vvtvp(npts,
                             &tmp[0], 1,
                             &m_traceNormals[i][id2], 1,
                             &Fwd[1+i][id2], 1,
                             &Fwd[1+i][id2], 1);
                
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
     * @brief Wall boundary conditions for viscous compressible flow problems.
     */
    void CompressibleFlowSystem::WallBoundaryViscous(
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
            m_fields[i]->ExtractTracePhys(physarray[i], Fwd[i]);
        }
        
        // Adjust the physical values of the trace to take 
        // user defined boundaries into account
        int e, id1, id2, npts;
        
        for (e = 0; e < m_fields[0]->
            GetBndCondExpansions()[bcRegion]->GetExpSize(); ++e)
        {
            npts = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                GetExp(e)->GetTotPoints();
            id1  = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                GetPhys_Offset(e);
            id2  = m_fields[0]->GetTrace()->GetPhys_Offset(
                m_fields[0]->GetTraceMap()->
                    GetBndCondCoeffsToGlobalCoeffsMap(cnt+e));
            
            for (i = 0; i < m_expdim ; i++)
            {
                Vmath::Neg(npts,&Fwd[i+1][id2], 1);
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
                GetExp(e)->GetTotPoints();
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
     * @brief Return the flux vector for the LDG diffusion problem.
     * \todo Complete the viscous flux vector
     */
    void CompressibleFlowSystem::GetViscousFluxVector(
        const int i, 
        const Array<OneD, Array<OneD, NekDouble> >               &physfield,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &derivatives,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &viscousTensor)
    {
        int k, j, d;
        
        int nq                       = m_fields[0]->GetTotPoints();
        NekDouble lambda             = -0.66666;
        NekDouble thermalDiffusivity = 0.000019;;

        
        Array<OneD, NekDouble> mu(nq);
        Array<OneD, NekDouble> tmp_mu(nq);
        Array<OneD, NekDouble> tmp(nq);
        
        Array<OneD, Array<OneD, NekDouble> > velocities(m_expdim);
        Array<OneD, Array<OneD, NekDouble> > du(m_expdim);
        Array<OneD, Array<OneD, NekDouble> > dv(m_expdim);
        Array<OneD, Array<OneD, NekDouble> > dw(m_expdim);
        Array<OneD, Array<OneD, NekDouble> > dT(m_expdim);
        
        Array<OneD, NekDouble > pressure   (nq, 0.0);
        Array<OneD, NekDouble > temperature(nq, 0.0);
        
        GetPressure(physfield, pressure);
        GetTemperature(physfield, pressure, temperature);
        GetDynamicViscosity(physfield, mu);
        
        for (k = 0; k < m_expdim; ++k)
        {
            velocities[k] = Array<OneD, NekDouble >(nq, 0.0);
            du[k] = Array<OneD, NekDouble >(nq, 0.0);
            dv[k] = Array<OneD, NekDouble >(nq, 0.0);
            dw[k] = Array<OneD, NekDouble >(nq, 0.0);
            dT[k] = Array<OneD, NekDouble >(nq, 0.0);
        }
        
        // Building the velocities
        for (k = 0; k < m_expdim; ++k)
        {
            Vmath::Vdiv(nq, physfield[k+1], 1, 
                        physfield[0], 1, 
                        velocities[k], 1);
        }
         
        // Building the proper derivatives
        if (m_expdim == 1)
        {
            for (d = 0; d < m_expdim; ++d)
            {
                Vmath::Zero(nq, tmp, 1);
                Vmath::Vmul(nq, velocities[0], 1, derivatives[d][0], 1, tmp, 1);
                Vmath::Vsub(nq, derivatives[d][1], 1, tmp, 1, tmp, 1);
                for (j = 0; j < nq; ++j)
                {
                    du[d][j] = (1.0 / physfield[0][j]) * tmp[j];
                }
            }
        }
        else if (m_expdim == 2)
        {
            for (d = 0; d < m_expdim; ++d)
            {
                Vmath::Zero(nq, tmp, 1);
                Vmath::Vmul(nq, velocities[0], 1, derivatives[d][0], 1, tmp, 1);
                Vmath::Vsub(nq, derivatives[d][1], 1, tmp, 1, tmp, 1);
                for (j = 0; j < nq; ++j)
                {
                    du[d][j] = (1.0 / physfield[0][j]) * tmp[j];
                }
            }
            for (d = 0; d < m_expdim; ++d)
            {
                Vmath::Zero(nq, tmp, 1);
                Vmath::Vmul(nq, velocities[1], 1, derivatives[d][0], 1, tmp, 1);
                Vmath::Vsub(nq, derivatives[d][2], 1, tmp, 1, tmp, 1);
                for (j = 0; j < nq; ++j)
                {
                    dv[d][j] = (1.0 / physfield[0][j]) * tmp[j];
                }
            }
            
            // At the moment for 2D only
            for (d = 0; d < m_expdim; ++d)
            {
                for (j = 0; j < nq; ++j)
                {
                    dT[d][j] = 1.0 / physfield[0][j] * derivatives[d][3][j] -
                    (temperature[j]/physfield[0][j]) * derivatives[d][0][j] - 
                    (0.5 * m_gasConstant / (m_gamma - 1)) * 
                    (2.0 * velocities[0][j] * derivatives[d][1][j] + 
                     ((velocities[0][j] * velocities[0][j] +
                       velocities[1][j] * velocities[1][j]) / 
                      physfield[0][j]) * derivatives[d][0][j]);
                }
            }
        }
         
        else if (m_expdim == 3)
        {
            for (d = 0; d < m_expdim; ++d)
            {
                Vmath::Zero(nq, tmp, 1);
                Vmath::Vmul(nq, velocities[0], 1, derivatives[d][0], 1, tmp, 1);
                Vmath::Vsub(nq, derivatives[d][1], 1, tmp, 1, tmp, 1);
                for (j = 0; j < nq; ++j)
                {
                    du[d][j] = (1.0 / physfield[0][j]) * tmp[j];
                }
            }
            for (d = 0; d < m_expdim; ++d)
            {
                Vmath::Zero(nq, tmp, 1);
                Vmath::Vmul(nq, velocities[1], 1, derivatives[d][0], 1, tmp, 1);
                Vmath::Vsub(nq, derivatives[d][2], 1, tmp, 1, tmp, 1);
                for (j = 0; j < nq; ++j)
                {
                    dv[d][j] = (1.0 / physfield[0][j]) * tmp[j];
                }
            }
            for (d = 0; d < m_expdim; ++d)
            {
                Vmath::Zero(nq, tmp, 1);
                Vmath::Vmul(nq, velocities[2], 1, derivatives[d][0], 1, tmp, 1);
                Vmath::Vsub(nq, derivatives[d][3], 1, tmp, 1, tmp, 1);
                for (j = 0; j < nq; ++j)
                {
                    dw[d][j] = (1.0 / physfield[0][j]) * tmp[j];
                }
            }
        }
        
        // Building the viscous flux vector
        if (i == 0)
        {
            // Viscous flux vector for the rho equation
            for (k = 0; k < m_expdim; ++k)
            {
                Vmath::Zero(nq, viscousTensor[k][i], 1);
            }
        }
        
        if (m_expdim == 1)
        {
            // to be completed
        }
        else if (m_expdim == 2)
        {
            if (i == 1)
            {
                for (j = 0; j < nq; ++j)
                {
                    viscousTensor[0][i][j] = 2.0 * mu[j] * du[0][j] 
                                            + lambda * (du[0][j] + 
                                                        dv[0][j]);
                    
                    viscousTensor[1][i][j] = dv[0][j] + du[1][j];
                }
            }
            else if (i == 2)
            {
                for (j = 0; j < nq; ++j)
                {
                    viscousTensor[0][i][j] = dv[0][j] + du[1][j];
                
                    viscousTensor[1][i][j] = 2.0 * mu[j] * dv[1][j] 
                                            + lambda * (du[0][j] + 
                                                        dv[0][j]);
                }
            }
            else if (i == 3)
            {
                for (j = 0; j < nq; ++j)
                {
                    viscousTensor[0][i][j] = 
                        velocities[0][j] * viscousTensor[0][0][j] + 
                        velocities[1][j] * viscousTensor[1][1][j] + 
                        (thermalDiffusivity / mu[j]) * (dT[0][j]);
                
                    viscousTensor[1][i][j] = 
                        velocities[1][j] * viscousTensor[1][2][j] + 
                        velocities[0][j] * viscousTensor[1][1][j] + 
                        (thermalDiffusivity / mu[j]) * (dT[1][j]);
                }
            }
        }
        else if (m_expdim == 3)
        {
            if (i == 1)
            {
                // to be completed
            }
            else if (i == 2)
            {
                // to be completed
            }
            else if (i == 3)
            {
                // to be completed
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
        NekDouble    alpha = -0.5;
        
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
        const Array<OneD, const Array<OneD, NekDouble> > &physfield,
        Array<OneD,                         NekDouble  > &pressure,
        Array<OneD,                         NekDouble  > &temperature)
    {
        for (int i = 0; i < m_fields[0]->GetTotPoints(); ++i)
        {
            temperature[i] = pressure[i] / (physfield[0][i] * m_gasConstant);
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
            soundspeed[i] = sqrt(m_gamma * pressure[i] / physfield[0][i]);
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
     * @brief Compute the dynamic viscosity using the Sutherland's law 
     * \f$ \mu = \mu_star * (T / T_star)^3/2 * (T_star + 110) / (T + 110) \f$,
     * where:   \mu_star = 1.7894 * 10^-5 Kg / (m * s)
     *          T_star   = 288.15 K
     * 
     * @param physfield    Input physical field.
     * @param mu           The resulting dynamic viscosity.
     */
    void CompressibleFlowSystem::GetDynamicViscosity(
        const Array<OneD, const Array<OneD, NekDouble> > &physfield,
             Array<OneD,                    NekDouble  > &mu)
    {
        const int npts       = m_fields[0]->GetTotPoints();
        const double mu_star = 0.00001794;
        const double T_star  = 288.15;
        Vmath::Zero(npts, mu, 1);
        
        Array<OneD, NekDouble > pressure         (npts, 0.0);
        Array<OneD, NekDouble > temperature      (npts, 0.0);
        Array<OneD, NekDouble > temperature_ratio(npts, 0.0);
        
        GetPressure(physfield, pressure);
        GetTemperature(physfield, pressure, temperature);
        
        for (int i = 0; i < npts; ++i)
        {
            temperature_ratio[i] = temperature[i] / T_star;
            mu[i] = mu_star * pow(temperature_ratio[i], 1.50) * 
                    (T_star + 110.0) / (temperature[i] + 110.0);
        }
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
            else if (boost::dynamic_pointer_cast<LocalRegions::HexExp>(
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
