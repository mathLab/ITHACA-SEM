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
        m_velLoc = Array<OneD, NekDouble>(m_spacedim);
        for (int i = 0; i < m_spacedim; ++i)
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
        if (m_spacedim == 2 || m_spacedim == 3)
        {
            ASSERTL0(m_session->DefinesParameter("vInf"),
                     "Compressible flow sessions must define a vInf parameter"
                     "for 2D/3D problems.");
            m_session->LoadParameter("vInf", m_vInf, 0.0);
        }
        
        // Get wInf parameter from session file.
        if (m_spacedim == 3)
        {
            ASSERTL0(m_session->DefinesParameter("wInf"),
                     "Compressible flow sessions must define a wInf parameter"
                     "for 3D problems.");
            m_session->LoadParameter("wInf", m_wInf, 0.0);
        }
                
        m_session->LoadParameter ("GasConstant",   m_gasConstant,   287.058);
        m_session->LoadParameter ("Twall",         m_Twall,         300.15);
        m_session->LoadSolverInfo("ViscosityType", m_ViscosityType, "Constant");
        m_session->LoadParameter ("mu",            m_mu,            1.78e-05);
        m_session->LoadParameter ("thermalConductivity",
                                  m_thermalConductivity, 0.0257);
        
        m_Cp      = m_gamma / (m_gamma - 1.0) * m_gasConstant;
        m_Prandtl = m_Cp * m_mu / m_thermalConductivity;


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
                
                // Setting up advection and diffusion operators
                m_session->LoadSolverInfo("AdvectionType", advName, "WeakDG");
                m_session->LoadSolverInfo("DiffusionType", diffName, "LDGNS");
                m_advection = SolverUtils::GetAdvectionFactory()
                                            .CreateInstance(advName, advName);
                m_diffusion = SolverUtils::GetDiffusionFactory()
                                            .CreateInstance(diffName, diffName);

                // Setting up flux vector for advection operator
                m_advection->SetFluxVector(&CompressibleFlowSystem::
                                            GetFluxVector, this);
                
                // Setting up flux vector for diffusion operator
                m_diffusion->SetFluxVectorNS(&CompressibleFlowSystem::
                                             GetViscousFluxVector, this);

                // Setting up Riemann solver for advection operator
                m_session->LoadSolverInfo("UpwindType", riemName, "Average");
                m_riemannSolver = SolverUtils::GetRiemannSolverFactory()
                                            .CreateInstance(riemName);
                
                // Setting up upwind solver for diffusion operator
                m_riemannSolverLDG = SolverUtils::GetRiemannSolverFactory()
                                                .CreateInstance("UpwindLDG");

                // Setting up parameters for advection operator Riemann solver 
                m_riemannSolver->AddParam (
                                    "gamma",  
                                    &CompressibleFlowSystem::GetGamma,   this);
                m_riemannSolver->AddScalar(
                                    "velLoc", 
                                    &CompressibleFlowSystem::GetVelLoc,  this);
                m_riemannSolver->AddVector(
                                    "N",
                                    &CompressibleFlowSystem::GetNormals, this);
                
                // Setting up parameters for diffusion operator Riemann solver
                m_riemannSolverLDG->AddParam (
                                    "gamma",  
                                    &CompressibleFlowSystem::GetGamma,   this);
                m_riemannSolverLDG->AddScalar(
                                    "velLoc", 
                                    &CompressibleFlowSystem::GetVelLoc,  this);
                m_riemannSolverLDG->AddVector(
                                    "N",
                                    &CompressibleFlowSystem::GetNormals, this);
                
                // Concluding initialisation of advection / diffusion operators
                m_advection->SetRiemannSolver   (m_riemannSolver);
                m_diffusion->SetRiemannSolver   (m_riemannSolverLDG);
                m_advection->InitObject         (m_session, m_fields);
                m_diffusion->InitObject         (m_session, m_fields);
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
        int nTracePts = GetTraceTotPoints();
        int nvariables      = physarray.num_elements();
        
        // get physical values of the forward trace
        Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
        for (i = 0; i < nvariables; ++i)
        {
            Fwd[i] = Array<OneD, NekDouble>(nTracePts);
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
            Array<OneD, NekDouble> tmp(npts, 0.0);

            // Calculate (v.n)
            for (i = 0; i < m_spacedim; ++i)
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
            for (i = 0; i < m_spacedim; ++i)
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
        int nTracePts = GetTraceTotPoints();
        int nvariables      = physarray.num_elements();
        
        // get physical values of the forward trace
        Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
        for (i = 0; i < nvariables; ++i)
        {
            Fwd[i] = Array<OneD, NekDouble>(nTracePts);
            m_fields[i]->ExtractTracePhys(physarray[i], Fwd[i]);
        }
        
        // Adjust the physical values of the trace to 
        // take user defined boundaries into account
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
            
            for (i = 0; i < m_spacedim; i++)
            {
                Vmath::Neg(npts, &Fwd[i+1][id2], 1);
            }
                        
            // Imposition of the temperature
            //Fwd[nvariables-1][id2] = m_gasConstant * Fwd[0][id2] * m_Twall / 
            //                            (m_gamma - 1);
            
            
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
        int nTracePts = GetTraceTotPoints();
        int nvariables      = physarray.num_elements();
        
        // Get physical values of the forward trace (from exp to phys)
        Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
        for (i = 0; i < nvariables; ++i)
        {
            Fwd[i] = Array<OneD, NekDouble>(nTracePts);
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
            
            switch(m_spacedim)
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
     * @brief Inflow characteristic boundary conditions for compressible 
     * flow problems.
     */
    void CompressibleFlowSystem::InflowCFSBoundary(
        int                                   bcRegion, 
        int                                   cnt, 
        Array<OneD, Array<OneD, NekDouble> > &physarray)
    {
        int i;
        int nTracePts              = GetTraceTotPoints();
        int nvariables             = physarray.num_elements();
        int nvel                   = m_spacedim;
        
        NekDouble gamma            = m_gamma;
        NekDouble gammaInv         = 1.0 / gamma;
        NekDouble gammaMinusOne    = gamma - 1.0;
        NekDouble gammaMinusOneInv = 1.0 / gammaMinusOne;
        
        Array<OneD, NekDouble> tmp1 (nTracePts, 0.0);
        Array<OneD, NekDouble> tmp2 (nTracePts, 0.0);
        Array<OneD, NekDouble> VnInf(nTracePts, 0.0);

        // Computing the normal velocity for characteristics coming 
        // from outside the computational domain
        Vmath::Smul(nTracePts, m_uInf, m_traceNormals[0], 1, VnInf, 1);
        
        if (nvel == 2 || nvel == 3)
        {
            Vmath::Smul(nTracePts, m_vInf, m_traceNormals[0], 1, tmp1, 1);
            Vmath::Vadd(nTracePts, VnInf, 1, tmp1, 1, VnInf, 1);
        }
        if (nvel == 3)
        {
            Vmath::Smul(nTracePts, m_wInf, m_traceNormals[0], 1, tmp2, 1);
            Vmath::Vadd(nTracePts, VnInf, 1, tmp2, 1, VnInf, 1);

        }
        
        // Get physical values of the forward trace
        Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
        for (i = 0; i < nvariables; ++i)
        {
            Fwd[i] = Array<OneD, NekDouble>(nTracePts);
            m_fields[i]->ExtractTracePhys(physarray[i], Fwd[i]);
        }
        
        // Computing the normal velocity for characteristics coming 
        // from inside the computational domain 
        Array<OneD, NekDouble > Vn (nTracePts, 0.0);
        Array<OneD, NekDouble > Vel(nTracePts, 0.0);
        for (i = 0; i < nvel; ++i)
        {
            Vmath::Vdiv(nTracePts, Fwd[i+1], 1, Fwd[0], 1, Vel, 1);
            Vmath::Vvtvp(nTracePts, m_traceNormals[i], 1, Vel, 1, Vn, 1, Vn, 1);
        }
        
        // Computing the absolute value of the velocity
        Array<OneD, NekDouble > absVel(nTracePts, 0.0);
        for (i = 0; i < nvel; ++i)
        {
            Vmath::Vdiv(nTracePts, Fwd[i+1], 1, Fwd[0], 1, tmp1, 1);
            Vmath::Vmul(nTracePts, tmp1, 1, tmp1, 1, tmp1, 1);
            Vmath::Vadd(nTracePts, tmp1, 1, absVel, 1, absVel, 1);
        }        
        Vmath::Vsqrt(nTracePts, absVel, 1, absVel, 1);

        // Get speed of sound
        Array<OneD, NekDouble > SoundSpeed(nTracePts);
        Array<OneD, NekDouble > pressure  (nTracePts);
        
        for (i = 0; i < nTracePts; i++)
        {
            pressure[i] = (gammaMinusOne) * (Fwd[3][i] - 
                                0.5 * (Fwd[1][i] * Fwd[1][i] / Fwd[0][i] + 
                                       Fwd[2][i] * Fwd[2][i] / Fwd[0][i]));
            
            SoundSpeed[i] = sqrt(gamma * pressure[i] / Fwd[0][i]);
        }
        
        // Get Mach
        Array<OneD, NekDouble > Mach(nTracePts);
        
        // This looks wrong
        //Vmath::Vdiv(nTracePts, Vel, 1, SoundSpeed, 1, Mach, 1);
        
        // The Mach number you should get here is the total Mach number
        Vmath::Vdiv(nTracePts, absVel, 1, SoundSpeed, 1, Mach, 1);
                
        // Auxiliary variables
        int e, id1, id2, npts, pnt;
        NekDouble cPlus, rPlus, cMinus, rMinus;
        NekDouble VelNorm, VelNormRef, VDelta;
        NekDouble rhob, rhoub, rhovb, rhoeb;
        NekDouble ub, vb, cb, sb, pb;
        
        // Loop on bcRegions
        for (e = 0; e < m_fields[0]->GetBndCondExpansions()[bcRegion]->
             GetExpSize(); ++e)
        {
            npts = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                GetExp(e)->GetNumPoints(0);
            
            id1 = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                GetPhys_Offset(e);
            id2 = m_fields[0]->GetTrace()->
                GetPhys_Offset(m_fields[0]->GetTraceMap()->
                    GetBndCondTraceToGlobalTraceMap(cnt++));
            
            // Loop on the points of the bcRegion
            for (i = 0; i < npts; i++)
            {
                pnt = id2+i;
                if(Mach[pnt] < 0.99)
                {
                    // + Characteristic from boundary
                    cPlus = sqrt(gamma * /*m_pressure*/m_pInf / /*m_density*/m_rhoInf);
                    rPlus = /*m_velocity*/VnInf[pnt] + 2.0 * cPlus * gammaMinusOneInv;
                    
                    // - Characteristic from inside
                    cMinus = sqrt(gamma * pressure[pnt] / Fwd[0][pnt]);
                    rMinus = Vn[pnt] - 2.0 * cMinus * gammaMinusOneInv;
                    
                    // Boundary Variables
                    VelNorm = 0.5 * (rPlus + rMinus);
                    cb      = 0.25 * gammaMinusOne * (rPlus - rMinus);
                    
                    // Boundary Velocity
                    VDelta  = VelNorm - /*m_velocity*/VnInf[pnt];
                    
                    ub      = /*m_velocity*/VnInf[pnt] * m_traceNormals[0][pnt] + 
                                VDelta * m_traceNormals[0][pnt];
                    vb      = /*m_velocity*/VnInf[pnt] * m_traceNormals[1][pnt] + 
                                VDelta * m_traceNormals[1][pnt];
                    
                    sb      = /*m_pressure*/m_pInf / (pow(/*m_density*/m_rhoInf, gamma));
                    rhob    = pow((cb * cb) / (gamma * sb), gammaMinusOneInv);
                    pb      = rhob * cb * cb * gammaInv;
                    
                    rhoub   = rhob * ub;
                    rhovb   = rhob * vb;
                    
                    rhoeb   = pb * gammaMinusOneInv + 
                                0.5 * (rhoub * ub + rhovb * vb);
                    
                    // Extrapolation
                    (m_fields[0]->GetBndCondExpansions()[bcRegion]->
                            UpdatePhys())[id1+i] = 2.0 * rhob - Fwd[0][pnt];
                    
                    (m_fields[1]->GetBndCondExpansions()[bcRegion]->
                            UpdatePhys())[id1+i] = 2.0 * rhoub - Fwd[1][pnt];
                    
                    (m_fields[2]->GetBndCondExpansions()[bcRegion]->
                            UpdatePhys())[id1+i] = 2.0 * rhovb - Fwd[2][pnt];
                    
                    (m_fields[3]->GetBndCondExpansions()[bcRegion]->
                            UpdatePhys())[id1+i] = 2.0 * rhoeb - Fwd[3][pnt];
                }
                else
                {
                    (m_fields[0]->GetBndCondExpansions()[bcRegion]->
                            UpdatePhys())[id1+i] = /*m_density*/m_rhoInf;
                    (m_fields[1]->GetBndCondExpansions()[bcRegion]->
                            UpdatePhys())[id1+i] = /*m_density * m_velocity*/m_rhoInf * m_uInf;
                    (m_fields[2]->GetBndCondExpansions()[bcRegion]->
                            UpdatePhys())[id1+i] = /*0.0*/ m_rhoInf * m_vInf;
                    (m_fields[3]->GetBndCondExpansions()[bcRegion]->
                            UpdatePhys())[id1+i] = /*m_pressure / (m_gamma-1.0) + 
                                0.5 * m_density * m_velocity * m_velocity*/m_pInf / (gamma - 1) + 
                                0.5 * m_rhoInf * (m_uInf * m_uInf + m_vInf * m_vInf);
                }
            }
        }
    }
    
    
    /**
     * @brief Outflow characteristic boundary conditions for compressible 
     * flow problems.
     */
    void CompressibleFlowSystem::OutflowCFSBoundary(
        int                                   bcRegion, 
        int                                   cnt, 
        Array<OneD, Array<OneD, NekDouble> > &physarray)
    {
        int i;
        int nTracePts               = GetTraceTotPoints();
        int nvariables              = physarray.num_elements();
        int nvel                    = m_spacedim;
        NekDouble gamma             = m_gamma;
        NekDouble gammaInv          = 1.0 / gamma;
        NekDouble gammaMinusOne     = gamma - 1.0;
        NekDouble gammaMinusOneInv  = 1.0 / gammaMinusOne;
        
        Array<OneD, NekDouble> tmp1 (nTracePts, 0.0);
        Array<OneD, NekDouble> tmp2 (nTracePts, 0.0);
        Array<OneD, NekDouble> VnInf(nTracePts, 0.0);
        
        // Computing the normal velocity for characteristics coming 
        // from outside the computational domain
        Vmath::Smul(nTracePts, m_uInf, m_traceNormals[0], 1, VnInf, 1);
        
        if (nvel == 2 || nvel == 3)
        {
            Vmath::Smul(nTracePts, m_vInf, m_traceNormals[0], 1, tmp1, 1);
            Vmath::Vadd(nTracePts, VnInf, 1, tmp1, 1, VnInf, 1);
        }
        if (nvel == 3)
        {
            Vmath::Smul(nTracePts, m_wInf, m_traceNormals[0], 1, tmp2, 1);
            Vmath::Vadd(nTracePts, VnInf, 1, tmp2, 1, VnInf, 1);
            
        }
        
        // Get physical values of the forward trace
        Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
        for (i = 0; i < nvariables; ++i)
        {
            Fwd[i] = Array<OneD, NekDouble>(nTracePts);
            m_fields[i]->ExtractTracePhys(physarray[i], Fwd[i]);
        }
        
        // Get normal velocity
        Array<OneD, NekDouble > Vn(nTracePts, 0.0);
        Array<OneD, NekDouble > Vel(nTracePts, 0.0);
        for (i = 0; i < nvel; ++i)
        {
            Vmath::Vdiv(nTracePts, Fwd[i+1], 1, Fwd[0], 1, Vel, 1);
            Vmath::Vvtvp(nTracePts, m_traceNormals[i], 1, Vel, 1, Vn, 1, Vn, 1);
        }
        
        // Get speed of sound
        Array<OneD, NekDouble > SoundSpeed(nTracePts);
        Array<OneD, NekDouble > pressure  (nTracePts);
        for (i = 0; i < nTracePts; i++)
        {
            pressure[i] = (gammaMinusOne) 
            * (Fwd[3][i] - 0.5 * (Fwd[1][i] * Fwd[1][i] / Fwd[0][i] + 
                                  Fwd[2][i] * Fwd[2][i] / Fwd[0][i]));
            
            SoundSpeed[i] = sqrt(gamma * pressure[i] / Fwd[0][i]);
        }
        
        // Get Mach
        Array<OneD, NekDouble > Mach(nTracePts);
        Vmath::Vdiv(nTracePts, Vn, 1, SoundSpeed, 1, Mach, 1);
        
        // Auxiliary variables 
        int e, id1, id2, npts, pnt;
        NekDouble cPlus, rPlus, cMinus, rMinus;
        NekDouble VelNorm, VelNormRef, VDelta;
        NekDouble rhob, rhoub, rhovb, rhoeb;
        NekDouble ub, vb, cb, sb, pb;
        
        // Loop on the bcRegions
        for (e = 0; e < m_fields[0]->GetBndCondExpansions()[bcRegion]->
             GetExpSize(); ++e)
        {
            npts = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                GetExp(e)->GetNumPoints(0);
            id1  = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                GetPhys_Offset(e) ;
            id2  = m_fields[0]->GetTrace()->
                GetPhys_Offset(m_fields[0]->GetTraceMap()->
                    GetBndCondTraceToGlobalTraceMap(cnt++));
            
            // Loop on points of bcRegion 'e'
            for (i = 0; i < npts; i++)
            {
                pnt = id2+i;
                /*
                // Specific case for Couette flow
                if(m_fields[0]->GetBndConditions()[bcRegion]->
                   GetUserDefined().GetEquation() == "OutFlowCouette")
                {
                    NekDouble rho0      = 1.0;
                    NekDouble Re        = 15;
                    NekDouble mu        = m_viscLam;
                    NekDouble u_average = 15.0*mu/4.0;
                    NekDouble T0        = 300;
                    NekDouble P         = m_GasConstant*rho0*T0;
                    P                   = P - 1.5 / Re * rho0 * u_average * 
                    u_average / 0.4 * 4.0;
                    
                    NekDouble sigma = 0.25;
                    NekDouble K     = sigma*(1.0-0.1*0.1)*0.0006675/4.0;
                    NekDouble L1    = K*(pressure[pnt]-P);
                    
                    (m_fields[0]->GetBndCondExpansions()[bcRegion]->
                     UpdatePhys())[id1+i] = Fwd[0][pnt];
                    (m_fields[1]->GetBndCondExpansions()[bcRegion]->
                     UpdatePhys())[id1+i] = Fwd[1][pnt];
                    (m_fields[2]->GetBndCondExpansions()[bcRegion]->
                     UpdatePhys())[id1+i] = Fwd[2][pnt];
                    P = P - L1;
                    NekDouble rhoe = P/(m_gamma-1.0) + 0.5*(Fwd[1][pnt]*Fwd[1][pnt] 
                                                            + Fwd[2][pnt]*Fwd[2][pnt])/Fwd[0][pnt];
                    (m_fields[3]->GetBndCondExpansions()[bcRegion]->
                     UpdatePhys())[id1+i] = rhoe;
                }
                 */
                // Subsonic flows
                if (Mach[pnt] < 0.99 /*&& m_fields[0]->
                        GetBndConditions()[bcRegion]->
                        GetUserDefined().GetEquation() != "OutFlowSuperSonic"*/)
                {
                    // + Characteristic from inside
                    cPlus = sqrt(gamma * pressure[pnt] / Fwd[0][pnt]);
                    rPlus = Vn[pnt] + 2.0 * cPlus * gammaMinusOneInv;
                    
                    // - Characteristic from boundary
                    cMinus = sqrt(gamma * /*m_pressure*/m_pInf / /*m_density*/m_rhoInf);
                    rMinus = /*m_velocity*/ VnInf[pnt] - 2.0 * cMinus * gammaMinusOneInv;
                    
                    // Boundary Variables
                    VelNorm = 0.5 * (rPlus + rMinus);
                    cb      = 0.25 * gammaMinusOne * (rPlus - rMinus);
                    sb      = pressure[pnt] / (pow(Fwd[0][pnt], gamma));
                    rhob    = pow((cb * cb) / (gamma * sb), gammaMinusOneInv);
                    pb      = rhob * cb * cb * gammaInv;
                    
                    // Boundary Velocity
                    VDelta  = VelNorm - Vn[pnt];
                    ub      = Fwd[1][pnt] / Fwd[0][pnt] + VDelta * m_traceNormals[0][pnt];
                    vb      = Fwd[2][pnt] / Fwd[0][pnt] + VDelta * m_traceNormals[1][pnt];
                    
                    rhoub   = rhob * ub;
                    rhovb   = rhob * vb;
                    rhoeb   = pb * gammaMinusOneInv + 0.5 * (rhoub * ub + rhovb * vb);
                    
                    // Extrapolation
                    (m_fields[0]->GetBndCondExpansions()[bcRegion]->
                            UpdatePhys())[id1+i] = Fwd[0][pnt];//2.0*rhob - Fwd[0][pnt];
                    (m_fields[1]->GetBndCondExpansions()[bcRegion]->
                            UpdatePhys())[id1+i] = Fwd[1][pnt];//2.0*rhoub - Fwd[1][pnt];
                    (m_fields[2]->GetBndCondExpansions()[bcRegion]->
                            UpdatePhys())[id1+i] = Fwd[2][pnt];//2.0*rhovb - Fwd[2][pnt];
                    
                    rhoeb = /*m_pressure*/m_pInf * gammaMinusOneInv + 0.5 * 
                            (Fwd[1][pnt] * Fwd[1][pnt] + 
                             Fwd[2][pnt] * Fwd[2][pnt])/Fwd[0][pnt];
                    
                    (m_fields[3]->GetBndCondExpansions()[bcRegion]->
                            UpdatePhys())[id1+i] = 2.0 * rhoeb - Fwd[3][pnt];
                }
                // Supersonic flows
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
    }
    
    /** 
     * @brief Extrapolation of order 0 for all the variables such that, 
     * at the boundaries, a trivial Riemann problem is solved.
     */
    void CompressibleFlowSystem::ExtrapOrder0Boundary(
        int                                   bcRegion, 
        int                                   cnt, 
        Array<OneD, Array<OneD, NekDouble> > &physarray)
    {
        int i, j;
        int e, pnt;
        int id1, id2, npts;
        
        int nTracePts  = GetTraceTotPoints();
        int nvariables = physarray.num_elements();
        int nvel       = m_spacedim;
                
        // Get physical values of the forward trace
        Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
        for (i = 0; i < nvariables; ++i)
        {
            Fwd[i] = Array<OneD, NekDouble>(nTracePts);
            m_fields[i]->ExtractTracePhys(physarray[i], Fwd[i]);
        }

        // Loop on bcRegions
        for (e = 0; e < m_fields[0]->GetBndCondExpansions()[bcRegion]->
             GetExpSize(); ++e)
        {
            npts = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                GetExp(e)->GetNumPoints(0);
            id1  = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                GetPhys_Offset(e) ;
            id2  = m_fields[0]->GetTrace()->
                GetPhys_Offset(m_fields[0]->GetTraceMap()->
                    GetBndCondTraceToGlobalTraceMap(cnt++));
            
            // Loop on points of bcRegion 'e'
            for (i = 0; i < npts; i++)
            {
                pnt = id2+i;
                
                // Setting up bcs for density
                (m_fields[0]->GetBndCondExpansions()[bcRegion]->
                    UpdatePhys())[id1+i] = Fwd[0][pnt];
                
                // Setting up bcs for velocities
                for (j = 1; j <=nvel; ++j)
                {
                    (m_fields[j]->GetBndCondExpansions()[bcRegion]->
                     UpdatePhys())[id1+i] = Fwd[j][pnt];
                }
                
                // Setting up bcs for energy
                (m_fields[nvariables-1]->GetBndCondExpansions()[bcRegion]->
                    UpdatePhys())[id1+i] = Fwd[nvariables-1][pnt];
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
        const Array<OneD, Array<OneD, NekDouble> >               &physfield,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux)
    {
        int i, j, nq = m_fields[0]->GetTotPoints();

        Array<OneD, NekDouble> pressure(nq);
        Array<OneD, Array<OneD, NekDouble> > velocity(m_spacedim);
        
        // Flux vector for the rho equation.
        for (i = 0; i < m_spacedim; ++i)
        {
            velocity[i] = Array<OneD, NekDouble>(nq);
            Vmath::Vcopy(nq, physfield[i+1], 1, flux[0][i], 1);
        }

        GetVelocityVector(physfield, velocity);
        GetPressure      (physfield, velocity, pressure);

        // Flux vector for the velocity fields.
        for (i = 0; i < m_spacedim; ++i)
        {
            for (j = 0; j < m_spacedim; ++j)
            {
                Vmath::Vmul(nq, velocity[j], 1, physfield[i+1], 1, flux[i+1][j], 1);
            }
            
            // Add pressure to appropriate field
            Vmath::Vadd(nq, flux[i+1][i], 1, pressure, 1, flux[i+1][i], 1);
        }

        // Flux vector for energy.
        Vmath::Vadd(nq, physfield[m_spacedim+1], 1, pressure, 1, pressure, 1);
        for (j = 0; j < m_spacedim; ++j)
        {
            Vmath::Vmul(nq, velocity[j], 1, pressure, 1, flux[m_spacedim+1][j], 1);
        }
    }
    
    /** 
     * @brief Return the flux vector for the LDG diffusion problem.
     * \todo Complete the viscous flux vector
     */
    void CompressibleFlowSystem::GetViscousFluxVector(
        const Array<OneD, Array<OneD, NekDouble> >               &physfield,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &derivativesO1,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &viscousTensor)
    {
        int i, j, k;
        int nvariables = m_fields.num_elements();
        int nPts = m_fields[0]->GetTotPoints();
        
        // Note: the value below is referred to 20˚C
        //NekDouble thermalConductivity = /*0.0257*/0.0000257;
        
        // Stokes hypotesis
        NekDouble lambda = -0.66666;
        
        // Auxiliary variables
        Array<OneD, NekDouble > mu          (nPts, 0.0);
        Array<OneD, NekDouble > mu2         (nPts, 0.0);
        Array<OneD, NekDouble > divVel      (nPts, 0.0);
        Array<OneD, NekDouble > pressure    (nPts, 0.0);
        Array<OneD, NekDouble > temperature (nPts, 0.0);
                
        // Set up wrapper to fields data storage
        Array<OneD, Array<OneD, NekDouble> > fields(nvariables);
        
        // Reorder storage to list time-integrated fields first
        for(i = 0; i < nvariables; ++i)
        {
            fields[i] = m_fields[i]->UpdatePhys();
        }
        
        // Thermodynamic related quantities
        GetPressure(fields, pressure);
        GetTemperature(fields, pressure, temperature);
        
        // Variable viscosity through the Sutherland's law
        if (m_ViscosityType == "Variable")
        {
            GetDynamicViscosity(fields, mu);
        }
        else
        {
            Vmath::Sadd(nPts, m_mu, &mu[0], 1, &mu[0], 1);
        }
        
        // Computing diagonal terms of viscous stress tensor
        Array<OneD, Array<OneD, NekDouble> > tmp(m_spacedim);
        Array<OneD, Array<OneD, NekDouble> > Sgg(m_spacedim);
        
        // mu2 = 2 * mu
        Vmath::Smul(nPts, 2.0, &mu[0], 1, &mu2[0], 1);
        
        // Velocity divergence
        for (j = 0; j < m_spacedim; ++j)
        {
            Vmath::Vadd(nPts, &divVel[0], 1, 
                        &derivativesO1[j][j][0], 1, 
                        &divVel[0], 1);
        }
        
        // Velocity divergence scaled by lambda * mu
        Vmath::Smul(nPts, lambda, &divVel[0], 1, &divVel[0], 1);
        Vmath::Vmul(nPts, &mu[0], 1, &divVel[0], 1, &divVel[0], 1);
        
        // Digonal terms of viscous stress tensor (Sxx, Syy, Szz)
        // Sjj = 2 * mu * du_j/dx_j - (2 / 3) * mu * sum_j(du_j/dx_j)
        for (j = 0; j < m_spacedim; ++j)
        {
            tmp[j] = Array<OneD, NekDouble>(nPts, 0.0);
            Sgg[j] = Array<OneD, NekDouble>(nPts, 0.0);
            
            Vmath::Vmul(nPts, &mu2[0], 1, 
                        &derivativesO1[j][j][0], 1, 
                        &tmp[j][0], 1);
            
            Vmath::Vadd(nPts, &tmp[j][0], 1, &divVel[0], 1, &Sgg[j][0], 1);
        }
        
        // Extra diagonal terms of viscous stress tensor (Sxy, Sxz, Syz)
        // Note: they exist for 2D and 3D problems only
        Array<OneD, NekDouble > Sxy(nPts, 0.0);
        Array<OneD, NekDouble > Sxz(nPts, 0.0);
        Array<OneD, NekDouble > Syz(nPts, 0.0);
        
        if (m_spacedim == 2)
        {            
            // Sxy = (du/dy + dv/dx)
            Vmath::Vadd(nPts, &derivativesO1[0][1][0], 1,
                        &derivativesO1[1][0][0], 1,
                        &Sxy[0], 1);
            
            // Sxy = mu * (du/dy + dv/dx)
            Vmath::Vmul(nPts, &mu[0], 1, &Sxy[0], 1, &Sxy[0], 1);
        }
        else if (m_spacedim == 3)
        {
            // Sxy = (du/dy + dv/dx)
            Vmath::Vadd(nPts, &derivativesO1[0][1][0], 1,
                        &derivativesO1[1][0][0], 1,
                        &Sxy[0], 1);
            
            // Sxz = (du/dz + dw/dx)
            Vmath::Vadd(nPts, &derivativesO1[0][2][0], 1,
                        &derivativesO1[2][0][0], 1,
                        &Sxz[0], 1);
            
            // Syz = (dv/dz + dw/dy)
            Vmath::Vadd(nPts, &derivativesO1[1][2][0], 1,
                        &derivativesO1[2][1][0], 1,
                        &Syz[0], 1);
            
            // Sxy = mu * (du/dy + dv/dx)
            Vmath::Vmul(nPts, &mu[0], 1, &Sxy[0], 1, &Sxy[0], 1);
            
            // Sxz = mu * (du/dy + dv/dx)
            Vmath::Vmul(nPts, &mu[0], 1, &Sxz[0], 1, &Sxz[0], 1);
            
            // Syz = mu * (du/dy + dv/dx)
            Vmath::Vmul(nPts, &mu[0], 1, &Syz[0], 1, &Syz[0], 1);
        }
        
        // Energy-related terms
        Array<OneD, NekDouble > STx(nPts, 0.0);
        Array<OneD, NekDouble > STy(nPts, 0.0);
        Array<OneD, NekDouble > STz(nPts, 0.0);
        // Building the viscous flux vector
        if (i == 0)
        {
            // Viscous flux vector for the rho equation
            for (k = 0; k < m_spacedim; ++k)
            {
                Vmath::Zero(nPts, viscousTensor[k][i], 1);
            }
        }
        
        if (m_spacedim == 1)
        {
            Array<OneD, NekDouble > tmp1(nPts, 0.0);
            
            // u * Sxx
            Vmath::Vmul(nPts, &physfield[0][0], 1,
                        &Sgg[0][0], 1, &STx[0], 1);
            
            // k * dT/dx
            Vmath::Smul(nPts, m_thermalConductivity, 
                        &derivativesO1[0][2][0], 1, 
                        &tmp1[0], 1);
            
            // STx = u * Sxx + (K / mu) * dT/dx
            Vmath::Vadd(nPts, &STx[0], 1, &tmp1[0], 1, &STx[0], 1);
        }
        else if (m_spacedim == 2)
        {
            Array<OneD, NekDouble > tmp1(nPts, 0.0);
            Array<OneD, NekDouble > tmp2(nPts, 0.0);
            
            // Computation of STx
            
            // u * Sxx
            Vmath::Vmul(nPts, &physfield[0][0], 1,
                        &Sgg[0][0], 1, &STx[0], 1);
            
            // v * Sxy
            Vmath::Vmul(nPts, &physfield[1][0], 1,
                        &Sxy[0], 1, &tmp1[0], 1);
            
            // k * dT/dx
            Vmath::Smul(nPts, m_thermalConductivity, 
                        &derivativesO1[0][2][0], 1, 
                        &tmp2[0], 1);
            
            // STx = u * Sxx + v * Sxy + K * dT/dx
            Vmath::Vadd(nPts, &STx[0], 1, &tmp1[0], 1, &STx[0], 1);
            Vmath::Vadd(nPts, &STx[0], 1, &tmp2[0], 1, &STx[0], 1);
            
            // Computation of STy
            
            // Re-initialise temporary arrays
            Vmath::Zero(nPts, &tmp1[0], 1);
            Vmath::Zero(nPts, &tmp2[0], 1);
            
            // v * Syy
            Vmath::Vmul(nPts, &physfield[1][0], 1,
                        &Sgg[1][0], 1, &STy[0], 1);
            
            // u * Sxy
            Vmath::Vmul(nPts, &physfield[0][0], 1,
                        &Sxy[0], 1, &tmp1[0], 1);
                        
            // k * dT/dy
            Vmath::Smul(nPts, m_thermalConductivity, 
                        &derivativesO1[1][2][0], 1, 
                        &tmp2[0], 1);
            
            // STy = v * Syy + u * Sxy + K * dT/dy
            Vmath::Vadd(nPts, &STy[0], 1, &tmp1[0], 1, &STy[0], 1);
            Vmath::Vadd(nPts, &STy[0], 1, &tmp2[0], 1, &STy[0], 1);
        }
        else if (m_spacedim == 3)
        {
            Array<OneD, NekDouble > tmp1(nPts, 0.0);
            Array<OneD, NekDouble > tmp2(nPts, 0.0);
            Array<OneD, NekDouble > tmp3(nPts, 0.0);
            
            // Computation of STx
            
            // u * Sxx
            Vmath::Vmul(nPts, &physfield[0][0], 1,
                        &Sgg[0][0], 1, &STx[0], 1);
            
            // v * Sxy
            Vmath::Vmul(nPts, &physfield[1][0], 1,
                        &Sxy[0], 1, &tmp1[0], 1);
            
            // v * Sxy
            Vmath::Vmul(nPts, &physfield[2][0], 1,
                        &Sxz[0], 1, &tmp2[0], 1);
            
            // k * dT/dx
            Vmath::Smul(nPts, m_thermalConductivity, 
                        &derivativesO1[0][2][0], 1, 
                        &tmp3[0], 1);
            
            // STx = u * Sxx + v * Sxy + w * Sxz + (K / mu) * dT/dx
            Vmath::Vadd(nPts, &STx[0], 1, &tmp1[0], 1, &STx[0], 1);
            Vmath::Vadd(nPts, &STx[0], 1, &tmp2[0], 1, &STx[0], 1);
            Vmath::Vadd(nPts, &STx[0], 1, &tmp3[0], 1, &STx[0], 1);
            
            // Computation of STy
            
            // Re-initialise temporary arrays
            Vmath::Zero(nPts, &tmp1[0], 1);
            Vmath::Zero(nPts, &tmp2[0], 1);
            Vmath::Zero(nPts, &tmp3[0], 1);
            
            // v * Syy
            Vmath::Vmul(nPts, &physfield[1][0], 1,
                        &Sgg[1][0], 1, &STy[0], 1);
            
            // u * Sxy
            Vmath::Vmul(nPts, &physfield[0][0], 1,
                        &Sxy[0], 1, &tmp1[0], 1);
            
            // w * Syz
            Vmath::Vmul(nPts, &physfield[2][0], 1,
                        &Syz[0], 1, &tmp2[0], 1);
                        
            // k * dT/dy
            Vmath::Smul(nPts, m_thermalConductivity, 
                        &derivativesO1[1][2][0], 1, 
                        &tmp3[0], 1);
            
            // STy = v * Syy + u * Sxy + w * Syz + K * dT/dy
            Vmath::Vadd(nPts, &STy[0], 1, &tmp1[0], 1, &STy[0], 1);
            Vmath::Vadd(nPts, &STy[0], 1, &tmp2[0], 1, &STy[0], 1);
            Vmath::Vadd(nPts, &STy[0], 1, &tmp3[0], 1, &STy[0], 1);
            
            // Computation of STz
            
            // Re-initialise temporary arrays
            Vmath::Zero(nPts, &tmp1[0], 1);
            Vmath::Zero(nPts, &tmp2[0], 1);
            Vmath::Zero(nPts, &tmp3[0], 1);
            
            // w * Szz
            Vmath::Vmul(nPts, &physfield[2][0], 1,
                        &Sgg[2][0], 1, &STz[0], 1);
            
            // u * Sxz
            Vmath::Vmul(nPts, &physfield[0][0], 1,
                        &Sxz[0], 1, &tmp1[0], 1);
            
            // v * Syz
            Vmath::Vmul(nPts, &physfield[1][0], 1,
                        &Syz[0], 1, &tmp2[0], 1);
                        
            // k * dT/dz
            Vmath::Smul(nPts, m_thermalConductivity, 
                        &derivativesO1[2][2][0], 1, 
                        &tmp3[0], 1);
            
            // STz = w * Szz + u * Sxz + v * Syz + K * dT/dz
            Vmath::Vadd(nPts, &STz[0], 1, &tmp1[0], 1, &STz[0], 1);
            Vmath::Vadd(nPts, &STz[0], 1, &tmp2[0], 1, &STz[0], 1);
            Vmath::Vadd(nPts, &STz[0], 1, &tmp3[0], 1, &STz[0], 1);
        }
        
        switch(m_spacedim)
        {
            case 1:
            {
                // f_11v = f_rho = 0
                Vmath::Zero(nPts, &viscousTensor[0][0][0], 1);
                
                // f_21v = f_rhou 
                Vmath::Vcopy(nPts, &Sgg[0][0], 1, &viscousTensor[0][1][0], 1);
                
                // f_31v = f_E
                Vmath::Vcopy(nPts, &STx[0], 1, &viscousTensor[0][2][0], 1);
                break;
            }
            case 2:
            {
                // f_11v = f_rho1 = 0
                Vmath::Zero(nPts, &viscousTensor[0][0][0], 1);
                // f_12v = f_rho2 = 0
                Vmath::Zero(nPts, &viscousTensor[1][0][0], 1);
                
                // f_21v = f_rhou1
                Vmath::Vcopy(nPts, &Sgg[0][0], 1, &viscousTensor[0][1][0], 1);
                // f_22v = f_rhou2
                Vmath::Vcopy(nPts, &Sxy[0],    1, &viscousTensor[1][1][0], 1);
                
                // f_31v = f_rhov1
                Vmath::Vcopy(nPts, &Sxy[0],    1, &viscousTensor[0][2][0], 1);
                // f_32v = f_rhov2
                Vmath::Vcopy(nPts, &Sgg[1][0], 1, &viscousTensor[1][2][0], 1);
                
                // f_41v = f_E1
                Vmath::Vcopy(nPts, &STx[0], 1, &viscousTensor[0][3][0], 1);
                // f_42v = f_E2
                Vmath::Vcopy(nPts, &STy[0], 1, &viscousTensor[1][3][0], 1);
                break;
            }
            case 3:
            {
                // f_11v = f_rho1 = 0
                Vmath::Zero(nPts, &viscousTensor[0][0][0], 1);
                // f_12v = f_rho2 = 0
                Vmath::Zero(nPts, &viscousTensor[1][0][0], 1);
                // f_13v = f_rho3 = 0
                Vmath::Zero(nPts, &viscousTensor[2][0][0], 1);
                
                // f_21v = f_rhou1
                Vmath::Vcopy(nPts, &Sgg[0][0], 1, &viscousTensor[0][1][0], 1);
                // f_22v = f_rhou2
                Vmath::Vcopy(nPts, &Sxy[0],    1, &viscousTensor[1][1][0], 1);
                // f_23v = f_rhou3
                Vmath::Vcopy(nPts, &Sxz[0],    1, &viscousTensor[2][1][0], 1);
                
                // f_31v = f_rhov1
                Vmath::Vcopy(nPts, &Sxy[0],    1, &viscousTensor[0][2][0], 1);
                // f_32v = f_rhov2
                Vmath::Vcopy(nPts, &Sgg[1][0], 1, &viscousTensor[1][2][0], 1);
                // f_33v = f_rhov3
                Vmath::Vcopy(nPts, &Syz[0],    1, &viscousTensor[2][2][0], 1);
                
                // f_31v = f_rhow1
                Vmath::Vcopy(nPts, &Sxz[0],    1, &viscousTensor[0][3][0], 1);
                // f_32v = f_rhow2
                Vmath::Vcopy(nPts, &Syz[0],    1, &viscousTensor[1][3][0], 1);
                // f_33v = f_rhow3
                Vmath::Vcopy(nPts, &Sgg[2][0], 1, &viscousTensor[2][3][0], 1);
                
                // f_41v = f_E1
                Vmath::Vcopy(nPts, &STx[0], 1, &viscousTensor[0][4][0], 1);
                // f_42v = f_E2
                Vmath::Vcopy(nPts, &STy[0], 1, &viscousTensor[1][4][0], 1);
                // f_43v = f_E3
                Vmath::Vcopy(nPts, &STz[0], 1, &viscousTensor[2][4][0], 1);
                break;
            }
            default:
            {
                ASSERTL0(false, "Illegal expansion dimension");
            }
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
        NekDouble alpha = -0.5;
        
        // Calculate ||rho v||^2.
        Vmath::Vmul(npts, physfield[1], 1, physfield[1], 1, pressure, 1);
        for (int i = 1; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(npts, physfield[1+i], 1, physfield[1+i], 1, 
                               pressure,       1, pressure,       1);
        }
        // Divide by rho to get rho*||v||^2.
        Vmath::Vdiv (npts, pressure, 1, physfield[0], 1, pressure, 1);
        // pressure <- E - 0.5*pressure
        Vmath::Svtvp(npts,     alpha, 
                     pressure, 1, physfield[m_spacedim+1], 1, pressure, 1);
        // Multiply by (gamma-1).
        Vmath::Smul (npts, m_gamma-1, pressure, 1, pressure, 1);
    }
    
    /**
     * @brief Calculate the pressure field \f$ p =
     * (\gamma-1)(E-\frac{1}{2}\rho\| \mathbf{v} \|^2) \f$ assuming an ideal 
     * gas law.
     *
     * This is a slightly optimised way to calculate the pressure field which
     * avoids division by the density field if the velocity field has already
     * been calculated.
     * 
     * @param physfield  Input momentum.
     * @param velocity   Velocity vector.
     * @param pressure   Computed pressure field.
     */
    void CompressibleFlowSystem::GetPressure(
        const Array<OneD, const Array<OneD, NekDouble> > &physfield,
        const Array<OneD, const Array<OneD, NekDouble> > &velocity,
              Array<OneD,                   NekDouble>   &pressure)
    {
        int       npts  = m_fields[0]->GetTotPoints();
        NekDouble alpha = -0.5;
        
        // Calculate ||\rho v||^2.
        Vmath::Vmul (npts, velocity[0], 1, physfield[1], 1, pressure, 1);
        for (int i = 1; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(npts, velocity[i], 1, physfield[1+i], 1, 
                               pressure,    1, pressure,       1);
        }
        // pressure <- E - 0.5*pressure
        Vmath::Svtvp(npts,     alpha, 
                     pressure, 1, physfield[m_spacedim+1], 1, pressure, 1);
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
        
        for (int i = 0; i < m_spacedim; ++i)
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
        const int nq = m_fields[0]->GetTotPoints();
        Vmath::Vdiv(nq, pressure, 1, physfield[0], 1, temperature, 1);
        Vmath::Smul(nq, 1.0/m_gasConstant, temperature, 1, temperature, 1);
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
        const int nq = m_fields[0]->GetTotPoints();
        Vmath::Vdiv (nq, pressure, 1, physfield[0], 1, soundspeed, 1);
        Vmath::Smul (nq, m_gamma, soundspeed, 1, soundspeed, 1);
        Vmath::Vsqrt(nq, soundspeed, 1, soundspeed, 1);
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

        Vmath::Vmul(npts, physfield[1], 1, physfield[1], 1, mach, 1);

        for (int i = 1; i < m_spacedim; ++i)
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
        
        GetPressure   (physfield, pressure);
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
        int n;
        int nElements = m_fields[0]->GetExpSize(); 
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
        Array<OneD, Array<OneD, NekDouble> > velocity   (m_spacedim);
        Array<OneD, Array<OneD, NekDouble> > stdVelocity(m_spacedim);
        Array<OneD, NekDouble>               pressure   (nTotQuadPoints);
        Array<OneD, NekDouble>               soundspeed (nTotQuadPoints);
        
        // Zero output array
        Vmath::Zero(stdV.num_elements(), stdV, 1);
        
        for (int i = 0; i < m_spacedim; ++i)
        {
            velocity   [i] = Array<OneD, NekDouble>(nTotQuadPoints);
            stdVelocity[i] = Array<OneD, NekDouble>(nTotQuadPoints, 0.0);
        }
        GetVelocityVector(inarray, velocity);
        GetPressure      (inarray, velocity, pressure);
        GetSoundSpeed    (inarray, pressure, soundspeed);
        
        for(int el = 0; el < n_element; ++el)
        { 
            // Possible bug: not multiply by jacobian??
            const Array<TwoD, const NekDouble> &gmat = 
                m_fields[0]->GetExp(el)->GetGeom()->GetGmat();
            
            int nq = m_fields[0]->GetExp(el)->GetTotPoints();
            
            if(m_fields[0]->GetExp(el)->GetGeom()->GetGtype() == 
                   SpatialDomains::eDeformed)
            {
                // d xi/ dx = gmat = 1/J * d x/d xi
                for (int i = 0; i < m_spacedim; ++i)
                {
                    Vmath::Vmul(nq, gmat[i], 1, velocity[0], 1,
                                stdVelocity[i], 1);
                    for (int j = 1; j < m_spacedim; ++j)
                    {
                        Vmath::Vvtvp(nq, gmat[m_spacedim*j+i], 1, velocity[j],
                                     1, stdVelocity[i], 1, stdVelocity[i], 1);
                    }
                }
            }
            else
            {
                for (int i = 0; i < m_spacedim; ++i)
                {
                    Vmath::Smul(nq, gmat[i][0], velocity[0], 1,
                                stdVelocity[i], 1);
                    for (int j = 1; j < m_spacedim; ++j)
                    {
                        Vmath::Svtvp(nq, gmat[m_spacedim*j+i][0], velocity[j], 
                                     1, stdVelocity[i], 1, stdVelocity[i], 1);
                    }
                }
            }
            
            for (int i = 0; i < nq; ++i)
            {
                NekDouble pntVelocity = 0.0;
                for (int j = 0; j < m_spacedim; ++j)
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
        ASSERTL0(n <= 20, "Illegal modes dimension for CFL calculation "
                          "(P has to be less then 20)");

        NekDouble CFLDG[21] = {  2.0000,   6.0000,  11.8424,  19.1569, 
                                27.8419,  37.8247,  49.0518,  61.4815, 
                                75.0797,  89.8181, 105.6700, 122.6200,
                               140.6400, 159.7300, 179.8500, 201.0100,
                               223.1800, 246.3600, 270.5300, 295.6900,
                               321.8300}; //CFLDG 1D [0-20]
        NekDouble CFL;
		
        if (m_projectionType == MultiRegions::eDiscontinuous)
        {
            CFL = CFLDG[n];
        }
        else 
        {
            ASSERTL0(false, "Continuous Galerkin stability coefficients "
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
