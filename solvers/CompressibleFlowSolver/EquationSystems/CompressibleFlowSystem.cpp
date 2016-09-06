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
// Description: Compressible flow system base class with auxiliary functions
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/EquationSystems/CompressibleFlowSystem.h>
#include <LocalRegions/TriExp.h>
#include <MultiRegions/ExpList.h>


using namespace std;

namespace Nektar
{
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

        m_varConv = MemoryManager<VariableConverter>::AllocateSharedPtr(
                    m_session, m_spacedim);

        ASSERTL0(m_session->DefinesSolverInfo("UPWINDTYPE"),
                 "No UPWINDTYPE defined in session.");

        // Do not forwards transform initial condition
        m_homoInitialFwd = false;

        // Set up locations of velocity vector.
        m_vecLocs = Array<OneD, Array<OneD, NekDouble> >(1);
        m_vecLocs[0] = Array<OneD, NekDouble>(m_spacedim);
        for (int i = 0; i < m_spacedim; ++i)
        {
            m_vecLocs[0][i] = 1 + i;
        }

        // Loading parameters from session file
        InitialiseParameters();

        // Setting up advection and diffusion operators
        InitAdvectionDiffusion();

        // Forcing terms for the sponge region
        m_forcing = SolverUtils::Forcing::Load(m_session, m_fields,
                                               m_fields.num_elements());

        // User-defined boundary conditions
        int cnt = 0;
        for (int n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
        {
            std::string type = 
                    m_fields[0]->GetBndConditions()[n]->GetUserDefined();
            if(!type.empty())
            {
                m_bndConds.push_back(GetCFSBndCondFactory().CreateInstance(
                        type,
                        m_session,
                        m_fields,
                        m_traceNormals,
                        m_spacedim,
                        n,
                        cnt));
            }
            cnt += m_fields[0]->GetBndCondExpansions()[n]->GetExpSize();
        }

        if (m_explicitAdvection)
        {
            m_ode.DefineOdeRhs    (&CompressibleFlowSystem::DoOdeRhs, this);
            m_ode.DefineProjection(&CompressibleFlowSystem::DoOdeProjection, this);
        }
        else
        {
            ASSERTL0(false, "Implicit CFS not set up.");
        }
    }

    /**
     * @brief Destructor for CompressibleFlowSystem class.
     */
    CompressibleFlowSystem::~CompressibleFlowSystem()
    {

    }

    /**
     * @brief Load CFS parameters from the session file
     */
    void CompressibleFlowSystem::InitialiseParameters()
    {
        NekDouble velInf, gasConstant;

        // Get gamma parameter from session file.
        m_session->LoadParameter("Gamma", m_gamma, 1.4);

        // Get gas constant from session file and compute Cp
        m_session->LoadParameter ("GasConstant",   gasConstant,   287.058);
        m_Cp      = m_gamma / (m_gamma - 1.0) * gasConstant;

        // Get pInf parameter from session file.
        m_session->LoadParameter("pInf", m_pInf, 101325);

        // Get rhoInf parameter from session file.
        m_session->LoadParameter("rhoInf", m_rhoInf, 1.225);

        // Get uInf parameter from session file.
        m_session->LoadParameter("uInf", velInf, 0.1);

        m_UInf = velInf*velInf;

        // Get vInf parameter from session file.
        if (m_spacedim == 2 || m_spacedim == 3)
        {
            m_session->LoadParameter("vInf", velInf, 0.0);
            m_UInf += velInf*velInf;
        }

        // Get wInf parameter from session file.
        if (m_spacedim == 3)
        {
            m_session->LoadParameter("wInf", velInf, 0.0);
            m_UInf += velInf*velInf;
        }
        m_UInf = sqrt(m_UInf);

        // Viscosity
        m_session->LoadSolverInfo("ViscosityType", m_ViscosityType, "Constant");
        m_session->LoadParameter ("mu",            m_mu,            1.78e-05);

        // Thermal conductivity or Prandtl
        if( m_session->DefinesParameter("thermalConductivity"))
        {
            ASSERTL0( !m_session->DefinesParameter("Pr"),
                 "Cannot define both Pr and thermalConductivity.");

            m_session->LoadParameter ("thermalConductivity",
                                        m_thermalConductivity);
            m_Prandtl = m_Cp * m_mu / m_thermalConductivity;
        }
        else
        {
            m_session->LoadParameter ("Pr",
                                        m_Prandtl, 0.72);
            m_thermalConductivity = m_Cp * m_mu / m_Prandtl;
        }

        // Parameters for sensor
        m_session->LoadParameter ("Skappa",        m_Skappa,        -2.048);
        m_session->LoadParameter ("Kappa",         m_Kappa,         0.0);
        m_session->LoadParameter ("mu0",           m_mu0,           1.0);

        // Steady state tolerance
        m_session->LoadParameter("SteadyStateTol", m_steadyStateTol, 0.0);
    }

    /**
     * @brief Create advection and diffusion objects for CFS
     */
    void CompressibleFlowSystem::InitAdvectionDiffusion()
    {
        // Check if projection type is correct
        ASSERTL0(m_projectionType == MultiRegions::eDiscontinuous,
                "Unsupported projection type.");

        string advName, diffName, riemName;
        m_session->LoadSolverInfo("AdvectionType", advName, "WeakDG");
        m_session->LoadSolverInfo("DiffusionType", diffName, "LDGNS");

        m_advection = SolverUtils::GetAdvectionFactory()
                                    .CreateInstance(advName, advName);
        m_diffusion = SolverUtils::GetDiffusionFactory()
                                    .CreateInstance(diffName, diffName);

        if (m_specHP_dealiasing)
        {
            m_advection->SetFluxVector(&CompressibleFlowSystem::
                                       GetFluxVectorDeAlias, this);
            m_diffusion->SetFluxVectorNS(
                &CompressibleFlowSystem::GetViscousFluxVectorDeAlias,
                this);
        }
        else
        {
            m_advection->SetFluxVector  (&CompressibleFlowSystem::
                                          GetFluxVector, this);
            m_diffusion->SetFluxVectorNS(&CompressibleFlowSystem::
                                          GetViscousFluxVector, this);
        }

        // Setting up Riemann solver for advection operator
        m_session->LoadSolverInfo("UpwindType", riemName, "Average");

        SolverUtils::RiemannSolverSharedPtr riemannSolver;
        riemannSolver = SolverUtils::GetRiemannSolverFactory()
                                    .CreateInstance(riemName);

        // Setting up upwind solver for diffusion operator
        SolverUtils::RiemannSolverSharedPtr riemannSolverLDG;
        riemannSolverLDG = SolverUtils::GetRiemannSolverFactory()
                                        .CreateInstance("UpwindLDG");

        // Setting up parameters for advection operator Riemann solver
        riemannSolver->SetParam (
            "gamma",   &CompressibleFlowSystem::GetGamma,   this);
        riemannSolver->SetAuxVec(
            "vecLocs", &CompressibleFlowSystem::GetVecLocs, this);
        riemannSolver->SetVector(
            "N",       &CompressibleFlowSystem::GetNormals, this);

        // Setting up parameters for diffusion operator Riemann solver
        riemannSolverLDG->SetParam (
            "gamma",   &CompressibleFlowSystem::GetGamma,   this);
        riemannSolverLDG->SetVector(
            "vecLocs", &CompressibleFlowSystem::GetVecLocs, this);
        riemannSolverLDG->SetVector(
            "N",       &CompressibleFlowSystem::GetNormals, this);

        // Concluding initialisation of advection / diffusion operators
        m_advection->SetRiemannSolver   (riemannSolver);
        m_diffusion->SetRiemannSolver   (riemannSolverLDG);
        m_advection->InitObject         (m_session, m_fields);
        m_diffusion->InitObject         (m_session, m_fields);
    }

    /**
     * @brief Compute the right-hand side.
     */
    void CompressibleFlowSystem::DoOdeRhs(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
        const NekDouble                                   time)
    {
        int i;
        int nvariables = inarray.num_elements();
        int npoints    = GetNpoints();
        int nTracePts  = GetTraceTotPoints();

        // Store forwards/backwards space along trace space
        Array<OneD, Array<OneD, NekDouble> > Fwd    (nvariables);
        Array<OneD, Array<OneD, NekDouble> > Bwd    (nvariables);

        if (m_HomogeneousType == eHomogeneous1D)
        {
            Fwd = NullNekDoubleArrayofArray;
            Bwd = NullNekDoubleArrayofArray;
        }
        else
        {
            for(i = 0; i < nvariables; ++i)
            {
                Fwd[i]     = Array<OneD, NekDouble>(nTracePts, 0.0);
                Bwd[i]     = Array<OneD, NekDouble>(nTracePts, 0.0);
                m_fields[i]->GetFwdBwdTracePhys(inarray[i], Fwd[i], Bwd[i]);
            }
        }

        // Calculate advection
        DoAdvection(inarray, outarray, time, Fwd, Bwd);

        // Negate results
        for (i = 0; i < nvariables; ++i)
        {
            Vmath::Neg(npoints, outarray[i], 1);
        }

        // Add diffusion terms
        DoDiffusion(inarray, outarray, Fwd, Bwd);

        // Add forcing terms
        std::vector<SolverUtils::ForcingSharedPtr>::const_iterator x;
        for (x = m_forcing.begin(); x != m_forcing.end(); ++x)
        {
            (*x)->Apply(m_fields, inarray, outarray, time);
        }
    }

    /**
     * @brief Compute the projection and call the method for imposing the 
     * boundary conditions in case of discontinuous projection.
     */
    void CompressibleFlowSystem::DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
        const NekDouble                                   time)
    {
        int i;
        int nvariables = inarray.num_elements();

        switch(m_projectionType)
        {
            case MultiRegions::eDiscontinuous:
            {
                // Just copy over array
                int npoints = GetNpoints();

                for(i = 0; i < nvariables; ++i)
                {
                    Vmath::Vcopy(npoints, inarray[i], 1, outarray[i], 1);
                }
                SetBoundaryConditions(outarray, time);
                break;
            }
            case MultiRegions::eGalerkin:
            case MultiRegions::eMixed_CG_Discontinuous:
            {
                ASSERTL0(false, "No Continuous Galerkin for full compressible "
                                "Navier-Stokes equations");
                break;
            }
            default:
                ASSERTL0(false, "Unknown projection scheme");
                break;
        }
    }

    /**
     * @brief Compute the advection terms for the right-hand side
     */
    void CompressibleFlowSystem::DoAdvection(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
        const NekDouble                                   time,
        const Array<OneD, Array<OneD, NekDouble> >       &pFwd,
        const Array<OneD, Array<OneD, NekDouble> >       &pBwd)
    {
        int nvariables = inarray.num_elements();
        Array<OneD, Array<OneD, NekDouble> > advVel(m_spacedim);

        m_advection->Advect(nvariables, m_fields, advVel, inarray,
                            outarray, time, pFwd, pBwd);
    }

    /**
     * @brief Add the diffusions terms to the right-hand side
     */
    void CompressibleFlowSystem::DoDiffusion(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const Array<OneD, Array<OneD, NekDouble> >   &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >   &pBwd)
    {
        v_DoDiffusion(inarray, outarray, pFwd, pBwd);
    }

    void CompressibleFlowSystem::SetBoundaryConditions(
            Array<OneD, Array<OneD, NekDouble> >             &physarray,
            NekDouble                                         time)
    {
        if (m_bndConds.size())
        {
            int nTracePts  = GetTraceTotPoints();
            int nvariables = physarray.num_elements();

            Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
            for (int i = 0; i < nvariables; ++i)
            {
                Fwd[i] = Array<OneD, NekDouble>(nTracePts);
                m_fields[i]->ExtractTracePhys(physarray[i], Fwd[i]);
            }

            // Loop over user-defined boundary conditions
            std::vector<CFSBndCondSharedPtr>::iterator x;
            for (x = m_bndConds.begin(); x != m_bndConds.end(); ++x)
            {
                (*x)->Apply(Fwd, physarray, time);
            }
        }
    }

    /**
     * @brief Return the flux vector for the compressible Euler equations.
     *
     * @param physfield   Fields.
     * @param flux        Resulting flux.
     */
    void CompressibleFlowSystem::GetFluxVector(
        const Array<OneD, Array<OneD, NekDouble> >               &physfield,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux)
    {
        int i, j;
        int nq = physfield[0].num_elements();
        int nVariables = m_fields.num_elements();

        Array<OneD, NekDouble> pressure(nq);
        Array<OneD, Array<OneD, NekDouble> > velocity(m_spacedim);

        // Flux vector for the rho equation
        for (i = 0; i < m_spacedim; ++i)
        {
            velocity[i] = Array<OneD, NekDouble>(nq);
            Vmath::Vcopy(nq, physfield[i+1], 1, flux[0][i], 1);
        }

        m_varConv->GetVelocityVector(physfield, velocity);
        m_varConv->GetPressure(physfield, velocity, pressure);

        // Flux vector for the velocity fields
        for (i = 0; i < m_spacedim; ++i)
        {
            for (j = 0; j < m_spacedim; ++j)
            {
                Vmath::Vmul(nq, velocity[j], 1, physfield[i+1], 1,
                            flux[i+1][j], 1);
            }

            // Add pressure to appropriate field
            Vmath::Vadd(nq, flux[i+1][i], 1, pressure, 1, flux[i+1][i], 1);
        }

        // Flux vector for energy.
        Vmath::Vadd(nq, physfield[m_spacedim+1], 1, pressure, 1,
                    pressure, 1);

        for (j = 0; j < m_spacedim; ++j)
        {
            Vmath::Vmul(nq, velocity[j], 1, pressure, 1,
                        flux[m_spacedim+1][j], 1);
        }

        // For the smooth viscosity model
        if (nVariables == m_spacedim+3)
        {
            // Add a zero row for the advective fluxes
            for (j = 0; j < m_spacedim; ++j)
            {
                Vmath::Zero(nq, flux[m_spacedim+2][j], 1);
            }
        }
    }

    /**
     * @brief Return the flux vector for the compressible Euler equations
     * by using the de-aliasing technique.
     *
     * @param physfield   Fields.
     * @param flux        Resulting flux.
     */
    void CompressibleFlowSystem::GetFluxVectorDeAlias(
        const Array<OneD, Array<OneD, NekDouble> >               &physfield,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux)
    {
        int i, j;
        int nq = physfield[0].num_elements();
        int nVariables = m_fields.num_elements();

        // Factor to rescale 1d points in dealiasing
        NekDouble OneDptscale = 2;
        nq = m_fields[0]->Get1DScaledTotPoints(OneDptscale);

        Array<OneD, NekDouble> pressure(nq);
        Array<OneD, Array<OneD, NekDouble> > velocity(m_spacedim);

        Array<OneD, Array<OneD, NekDouble> > physfield_interp(nVariables);
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > flux_interp(
                                                            nVariables);

        for (i = 0; i < nVariables; ++ i)
        {
            physfield_interp[i] = Array<OneD, NekDouble>(nq);
            flux_interp[i] = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
            m_fields[0]->PhysInterp1DScaled(
                OneDptscale, physfield[i], physfield_interp[i]);

            for (j = 0; j < m_spacedim; ++j)
            {
                flux_interp[i][j] = Array<OneD, NekDouble>(nq);
            }
        }

        // Flux vector for the rho equation
        for (i = 0; i < m_spacedim; ++i)
        {
            velocity[i] = Array<OneD, NekDouble>(nq);

            // Galerkin project solution back to original space
            m_fields[0]->PhysGalerkinProjection1DScaled(
                OneDptscale, physfield_interp[i+1], flux[0][i]);
        }

        m_varConv->GetVelocityVector(physfield_interp, velocity);
        m_varConv->GetPressure      (physfield_interp, velocity, pressure);

        // Evaluation of flux vector for the velocity fields
        for (i = 0; i < m_spacedim; ++i)
        {
            for (j = 0; j < m_spacedim; ++j)
            {
                Vmath::Vmul(nq, velocity[j], 1, physfield_interp[i+1], 1,
                            flux_interp[i+1][j], 1);
            }

            // Add pressure to appropriate field
            Vmath::Vadd(nq, flux_interp[i+1][i], 1, pressure,1,
                        flux_interp[i+1][i], 1);
        }

        // Galerkin project solution back to origianl space
        for (i = 0; i < m_spacedim; ++i)
        {
            for (j = 0; j < m_spacedim; ++j)
            {
                m_fields[0]->PhysGalerkinProjection1DScaled(
                    OneDptscale, flux_interp[i+1][j], flux[i+1][j]);
            }
        }

        // Evaluation of flux vector for energy
        Vmath::Vadd(nq, physfield_interp[m_spacedim+1], 1, pressure, 1,
                    pressure, 1);

        for (j = 0; j < m_spacedim; ++j)
        {
            Vmath::Vmul(nq, velocity[j], 1, pressure, 1,
                        flux_interp[m_spacedim+1][j], 1);

            // Galerkin project solution back to origianl space
            m_fields[0]->PhysGalerkinProjection1DScaled(
                OneDptscale,
                flux_interp[m_spacedim+1][j],
                flux[m_spacedim+1][j]);
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
        int i, j;
        int nVariables = m_fields.num_elements();
        int nPts       = physfield[0].num_elements();

        // Stokes hypothesis
        const NekDouble lambda = -2.0/3.0;

        // Auxiliary variables
        Array<OneD, NekDouble > mu                 (nPts, 0.0);
        Array<OneD, NekDouble > thermalConductivity(nPts, 0.0);
        Array<OneD, NekDouble > divVel             (nPts, 0.0);

        // Variable viscosity through the Sutherland's law
        if (m_ViscosityType == "Variable")
        {
            m_varConv->GetDynamicViscosity(physfield[nVariables-2], mu);
            NekDouble tRa = m_Cp / m_Prandtl;
            Vmath::Smul(nPts, tRa, mu, 1, thermalConductivity, 1);
        }
        else
        {
            Vmath::Fill(nPts, m_mu, mu, 1);
            Vmath::Fill(nPts, m_thermalConductivity,
                        thermalConductivity, 1);
        }

        // Velocity divergence
        for (j = 0; j < m_spacedim; ++j)
        {
            Vmath::Vadd(nPts, divVel, 1, derivativesO1[j][j], 1,
                        divVel, 1);
        }

        // Velocity divergence scaled by lambda * mu
        Vmath::Smul(nPts, lambda, divVel, 1, divVel, 1);
        Vmath::Vmul(nPts, mu,  1, divVel, 1, divVel, 1);

        // Viscous flux vector for the rho equation = 0
        for (i = 0; i < m_spacedim; ++i)
        {
            Vmath::Zero(nPts, viscousTensor[i][0], 1);
        }

        // Viscous stress tensor (for the momentum equations)
        for (i = 0; i < m_spacedim; ++i)
        {
            for (j = i; j < m_spacedim; ++j)
            {
                Vmath::Vadd(nPts, derivativesO1[i][j], 1,
                                  derivativesO1[j][i], 1,
                                  viscousTensor[i][j+1], 1);

                Vmath::Vmul(nPts, mu, 1,
                                  viscousTensor[i][j+1], 1,
                                  viscousTensor[i][j+1], 1);

                if (i == j)
                {
                    // Add divergence term to diagonal
                    Vmath::Vadd(nPts, viscousTensor[i][j+1], 1,
                                  divVel, 1,
                                  viscousTensor[i][j+1], 1);
                }
                else
                {
                    // Copy to make symmetric
                    Vmath::Vcopy(nPts, viscousTensor[i][j+1], 1,
                                       viscousTensor[j][i+1], 1);
                }
            }
        }

        // Terms for the energy equation
        for (i = 0; i < m_spacedim; ++i)
        {
            Vmath::Zero(nPts, viscousTensor[i][m_spacedim+1], 1);
            // u_j * tau_ij
            for (j = 0; j < m_spacedim; ++j)
            {
                Vmath::Vvtvp(nPts, physfield[j], 1,
                               viscousTensor[i][j+1], 1,
                               viscousTensor[i][m_spacedim+1], 1,
                               viscousTensor[i][m_spacedim+1], 1);
            }
            // Add k*T_i
            Vmath::Vvtvp(nPts, thermalConductivity, 1,
                               derivativesO1[i][m_spacedim], 1,
                               viscousTensor[i][m_spacedim+1], 1,
                               viscousTensor[i][m_spacedim+1], 1);
        }
    }

    /**
     * @brief Return the flux vector for the LDG diffusion problem.
     * \todo Complete the viscous flux vector
     */
    void CompressibleFlowSystem::GetViscousFluxVectorDeAlias(
        const Array<OneD, Array<OneD, NekDouble> >               &physfield,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &derivativesO1,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &viscousTensor)
    {
        int i, j;
        int nVariables = m_fields.num_elements();
        // Factor to rescale 1d points in dealiasing.
        NekDouble OneDptscale = 2;
        // Get number of points to dealias a cubic non-linearity
        int nPts      = m_fields[0]->Get1DScaledTotPoints(OneDptscale);
        int nPts_orig = physfield[0].num_elements();

        // Stokes hypothesis
        const NekDouble lambda = -2.0/3.0;

        // Auxiliary variables
        Array<OneD, NekDouble > mu                 (nPts, 0.0);
        Array<OneD, NekDouble > thermalConductivity(nPts, 0.0);
        Array<OneD, NekDouble > divVel             (nPts, 0.0);

        // Variable viscosity through the Sutherland's law
        if (m_ViscosityType == "Variable")
        {
            m_varConv->GetDynamicViscosity(physfield[nVariables-2], mu);
            NekDouble tRa = m_Cp / m_Prandtl;
            Vmath::Smul(nPts, tRa, mu, 1, thermalConductivity, 1);
        }
        else
        {
            Vmath::Fill(nPts, m_mu, mu, 1);
            Vmath::Fill(nPts, m_thermalConductivity,
                        thermalConductivity, 1);
        }

        // Interpolate inputs and initialise interpolated output
        Array<OneD, Array<OneD, NekDouble> > vel_interp(m_spacedim);
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                                             deriv_interp(m_spacedim);
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                                             out_interp(m_spacedim);
        for (i = 0; i < m_spacedim; ++i)
        {
            // Interpolate velocity
            vel_interp[i]   = Array<OneD, NekDouble> (nPts);
            m_fields[0]->PhysInterp1DScaled(
                OneDptscale, physfield[i], vel_interp[i]);

            // Interpolate derivatives
            deriv_interp[i] = Array<OneD,Array<OneD,NekDouble> > (m_spacedim+1);
            for (j = 0; j < m_spacedim+1; ++j)
            {
                deriv_interp[i][j] = Array<OneD, NekDouble> (nPts);
                m_fields[0]->PhysInterp1DScaled(
                    OneDptscale, derivativesO1[i][j], deriv_interp[i][j]);
            }

            // Output (start from j=1 since flux is zero for rho)
            out_interp[i] = Array<OneD,Array<OneD,NekDouble> > (m_spacedim+2);
            for (j = 1; j < m_spacedim+2; ++j)
            {
                out_interp[i][j] = Array<OneD, NekDouble> (nPts);
            }
        }

        // Velocity divergence
        for (j = 0; j < m_spacedim; ++j)
        {
            Vmath::Vadd(nPts, divVel, 1, deriv_interp[j][j], 1,
                        divVel, 1);
        }

        // Velocity divergence scaled by lambda * mu
        Vmath::Smul(nPts, lambda, divVel, 1, divVel, 1);
        Vmath::Vmul(nPts, mu,  1, divVel, 1, divVel, 1);

        // Viscous flux vector for the rho equation = 0 (no need to dealias)
        for (i = 0; i < m_spacedim; ++i)
        {
            Vmath::Zero(nPts_orig, viscousTensor[i][0], 1);
        }

        // Viscous stress tensor (for the momentum equations)
        for (i = 0; i < m_spacedim; ++i)
        {
            for (j = i; j < m_spacedim; ++j)
            {
                Vmath::Vadd(nPts, deriv_interp[i][j], 1,
                                  deriv_interp[j][i], 1,
                                  out_interp[i][j+1], 1);

                Vmath::Vmul(nPts, mu, 1,
                                  out_interp[i][j+1], 1,
                                  out_interp[i][j+1], 1);

                if (i == j)
                {
                    // Add divergence term to diagonal
                    Vmath::Vadd(nPts, out_interp[i][j+1], 1,
                                  divVel, 1,
                                  out_interp[i][j+1], 1);
                }
                else
                {
                    // Make symmetric
                    out_interp[j][i+1] = out_interp[i][j+1];
                }
            }
        }

        // Terms for the energy equation
        for (i = 0; i < m_spacedim; ++i)
        {
            Vmath::Zero(nPts, out_interp[i][m_spacedim+1], 1);
            // u_j * tau_ij
            for (j = 0; j < m_spacedim; ++j)
            {
                Vmath::Vvtvp(nPts, vel_interp[j], 1,
                               out_interp[i][j+1], 1,
                               out_interp[i][m_spacedim+1], 1,
                               out_interp[i][m_spacedim+1], 1);
            }
            // Add k*T_i
            Vmath::Vvtvp(nPts, thermalConductivity, 1,
                               deriv_interp[i][m_spacedim], 1,
                               out_interp[i][m_spacedim+1], 1,
                               out_interp[i][m_spacedim+1], 1);
        }

        // Project to original space
        for (i = 0; i < m_spacedim; ++i)
        {
            for (j = 1; j < m_spacedim+2; ++j)
            {
                m_fields[0]->PhysGalerkinProjection1DScaled(
                    OneDptscale,
                    out_interp[i][j],
                    viscousTensor[i][j]);
            }
        }
    }

    /**
     * @brief Perform post-integration checks, presently just to check steady
     * state behaviour.
     */
    bool CompressibleFlowSystem::v_PostIntegrate(int step)
    {
        if (m_steadyStateTol > 0.0)
        {
            bool doOutput = step % m_infosteps == 0;
            if (CalcSteadyState(doOutput))
            {
                if (m_comm->GetRank() == 0)
                {
                    cout << "Reached Steady State to tolerance "
                         << m_steadyStateTol << endl;
                }
                return true;
            }
        }
        return false;
    }

    /**
     * @brief Calculate whether the system has reached a steady state by
     * observing residuals to a user-defined tolerance.
     */
    bool CompressibleFlowSystem::CalcSteadyState(bool output)
    {
        const int nPoints = GetTotPoints();
        const int nFields = m_fields.num_elements();

        // Holds L2 errors.
        Array<OneD, NekDouble> L2      (nFields);
        Array<OneD, NekDouble> residual(nFields);

        for (int i = 0; i < nFields; ++i)
        {
            Array<OneD, NekDouble> diff(nPoints);

            Vmath::Vsub(nPoints, m_fields[i]->GetPhys(), 1, m_un[i], 1, diff, 1);
            Vmath::Vmul(nPoints, diff, 1, diff, 1, diff, 1);
            residual[i] = Vmath::Vsum(nPoints, diff, 1);
        }

        m_comm->AllReduce(residual, LibUtilities::ReduceSum);

        // L2 error
        L2[0] = sqrt(residual[0]) / m_rhoInf;

        for (int i = 1; i < nFields-1; ++i)
        {
            L2[i] = sqrt(residual[i]) / m_UInf / m_rhoInf;
        }

        NekDouble Einf = m_pInf / (m_gamma-1.0) + 0.5 * m_rhoInf * m_UInf;
        L2[nFields-1] = sqrt(residual[nFields-1]) / Einf;

        if (m_comm->GetRank() == 0 && output)
        {
            // Output time
            m_errFile << setprecision(8) << setw(17) << scientific << m_time;

            // Output residuals
            for (int i = 0; i < nFields; ++i)
            {
                m_errFile << setprecision(11) << setw(22) << scientific
                          << L2[i];
            }

            m_errFile << endl;
        }

        // Calculate maximum L2 error
        NekDouble maxL2 = Vmath::Vmax(nFields, L2, 1);

        if (m_session->DefinesCmdLineArgument("verbose") &&
            m_comm->GetRank() == 0 && output)
        {
            cout << "-- Maximum L^2 residual: " << maxL2 << endl;
        }

        if (maxL2 <= m_steadyStateTol)
        {
            return true;
        }

        for (int i = 0; i < m_fields.num_elements(); ++i)
        {
            Vmath::Vcopy(nPoints, m_fields[i]->GetPhys(), 1, m_un[i], 1);
        }

        return false;
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
        NekDouble minLength = 0.0;
        NekDouble alpha     = MaxTimeStepEstimator();
        NekDouble cLambda   = 0.2; // Spencer book-317

        // Loop over elements to compute the time-step limit for each element
        for(n = 0; n < nElements; ++n)
        {
            int npoints = m_fields[0]->GetExp(n)->GetTotPoints();
            Array<OneD, NekDouble> one2D(npoints, 1.0);
            NekDouble Area = m_fields[0]->GetExp(n)->Integral(one2D);

            minLength = sqrt(Area);
            if (m_fields[0]->GetExp(n)->as<LocalRegions::TriExp>())
            {
                minLength *= 2.0;
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
     * @brief Set up logic for residual calculation.
     */
    void CompressibleFlowSystem::v_SetInitialConditions(
        NekDouble initialtime,
        bool      dumpInitialConditions,
        const int domain)
    {
        EquationSystem::v_SetInitialConditions(initialtime, false);

        // insert white noise in initial condition
        NekDouble Noise;
        int phystot = m_fields[0]->GetTotPoints();
        Array<OneD, NekDouble> noise(phystot);

        m_session->LoadParameter("Noise", Noise,0.0);
        int m_nConvectiveFields =  m_fields.num_elements();

        if (Noise > 0.0)
        {
            for (int i = 0; i < m_nConvectiveFields; i++)
            {
                Vmath::FillWhiteNoise(phystot, Noise, noise, 1,
                                      m_comm->GetColumnComm()->GetRank()+1);
                Vmath::Vadd(phystot, m_fields[i]->GetPhys(), 1,
                            noise, 1, m_fields[i]->UpdatePhys(), 1);
                m_fields[i]->FwdTrans_IterPerExp(m_fields[i]->GetPhys(),
                                                 m_fields[i]->UpdateCoeffs());
            }
        }

        InitializeSteadyState();

        if (dumpInitialConditions)
        {
            // Dump initial conditions to file
            Checkpoint_Output(0);
        }
    }

    void CompressibleFlowSystem::InitializeSteadyState()
    {
        if (m_session->DefinesParameter("SteadyStateTol"))
        {
            const int nPoints = m_fields[0]->GetTotPoints();
            m_un = Array<OneD, Array<OneD, NekDouble> > (
                m_fields.num_elements());

            for (int i = 0; i < m_fields.num_elements(); ++i)
            {
                m_un[i] = Array<OneD, NekDouble>(nPoints);
                Vmath::Vcopy(nPoints, m_fields[i]->GetPhys(), 1, m_un[i], 1);
            }

            if (m_comm->GetRank() == 0)
            {
                std::string fName = m_session->GetSessionName() +
                    std::string(".res");
                m_errFile.open(fName.c_str());
                m_errFile << "# "
                          << setw(15) << left << "Time"
                          << setw(22) << left << "rho";

                std::string velFields[3] = {"u", "v", "w"};

                for (int i = 0; i < m_fields.num_elements()-2; ++i)
                {
                    m_errFile << setw(22) << "rho"+velFields[i];
                }

                m_errFile << setw(22) << left << "E" << endl;
            }
        }
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
        int nBCEdgePts           = 0;

        // Getting the velocity vector on the 2D normal space
        Array<OneD, Array<OneD, NekDouble> > velocity   (m_spacedim);
        Array<OneD, Array<OneD, NekDouble> > stdVelocity(m_spacedim);
        Array<OneD, NekDouble>               pressure   (nTotQuadPoints);
        Array<OneD, NekDouble>               soundspeed (nTotQuadPoints);
        LibUtilities::PointsKeyVector        ptsKeys;

        // Zero output array
        Vmath::Zero(stdV.num_elements(), stdV, 1);

        for (int i = 0; i < m_spacedim; ++i)
        {
            velocity   [i] = Array<OneD, NekDouble>(nTotQuadPoints);
            stdVelocity[i] = Array<OneD, NekDouble>(nTotQuadPoints, 0.0);
        }

        m_varConv->GetVelocityVector(inarray, velocity);
        m_varConv->GetPressure      (inarray, velocity, pressure);
        m_varConv->GetSoundSpeed    (inarray, pressure, soundspeed);

        for(int el = 0; el < n_element; ++el)
        {
            ptsKeys = m_fields[0]->GetExp(el)->GetPointsKeys();

            // Possible bug: not multiply by jacobian??
            const SpatialDomains::GeomFactorsSharedPtr metricInfo =
                m_fields[0]->GetExp(el)->GetGeom()->GetMetricInfo();
            const Array<TwoD, const NekDouble> &gmat =
                m_fields[0]->GetExp(el)->GetGeom()->GetMetricInfo()
                                                  ->GetDerivFactors(ptsKeys);

            int nq = m_fields[0]->GetExp(el)->GetTotPoints();

            if(metricInfo->GetGtype() == SpatialDomains::eDeformed)
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
                pntVelocity = sqrt(pntVelocity) + soundspeed[nBCEdgePts];
                if (pntVelocity > stdV[el])
                {
                    stdV[el] = pntVelocity;
                }
                nBCEdgePts++;
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
        NekDouble CFL = 0.0;

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

    void CompressibleFlowSystem::GetSensor(
        const Array<OneD, const Array<OneD, NekDouble> > &physarray,
              Array<OneD,                   NekDouble>   &Sensor,
              Array<OneD,                   NekDouble>   &SensorKappa)
    {
        int e, NumModesElement, nQuadPointsElement;
        int nTotQuadPoints  = GetTotPoints();
        int nElements       = m_fields[0]->GetExpSize();

        // Find solution (SolP) at p = P;
        // The input array (physarray) is the solution at p = P;

        Array<OneD,int> ExpOrderElement = GetNumExpModesPerExp();

        Array<OneD, NekDouble> SolP    (nTotQuadPoints, 0.0);
        Array<OneD, NekDouble> SolPmOne(nTotQuadPoints, 0.0);
        Array<OneD, NekDouble> SolNorm (nTotQuadPoints, 0.0);

        Vmath::Vcopy(nTotQuadPoints, physarray[0], 1, SolP, 1);

        int CoeffsCount = 0;

        for (e = 0; e < nElements; e++)
        {
            NumModesElement        = ExpOrderElement[e];
            int nQuadPointsElement = m_fields[0]->GetExp(e)->GetTotPoints();
            int nCoeffsElement     = m_fields[0]->GetExp(e)->GetNcoeffs();
            int numCutOff          = NumModesElement - 1;

            // Set-up of the Orthogonal basis for a Quadrilateral element which
            // is needed to obtain thesolution at P =  p - 1;

            Array<OneD, NekDouble> SolPElementPhys  (nQuadPointsElement, 0.0);
            Array<OneD, NekDouble> SolPElementCoeffs(nCoeffsElement,     0.0);

            Array<OneD, NekDouble> SolPmOneElementPhys(nQuadPointsElement, 0.0);
            Array<OneD, NekDouble> SolPmOneElementCoeffs(nCoeffsElement, 0.0);

            // create vector the save the solution points per element at P = p;

            for (int i = 0; i < nQuadPointsElement; i++)
            {
                SolPElementPhys[i] = SolP[CoeffsCount+i];
            }

            m_fields[0]->GetExp(e)->FwdTrans(SolPElementPhys,
                                             SolPElementCoeffs);

            // ReduceOrderCoeffs reduces the polynomial order of the solution
            // that is represented by the coeffs given as an inarray. This is
            // done by projecting the higher order solution onto the orthogonal
            // basis and padding the higher order coefficients with zeros.

            m_fields[0]->GetExp(e)->ReduceOrderCoeffs(numCutOff,
                                                      SolPElementCoeffs,
                                                      SolPmOneElementCoeffs);

            m_fields[0]->GetExp(e)->BwdTrans(SolPmOneElementCoeffs,
                                             SolPmOneElementPhys);

            for (int i = 0; i < nQuadPointsElement; i++)
            {
                SolPmOne[CoeffsCount+i] = SolPmOneElementPhys[i];
            }

            NekDouble SolPmeanNumerator   = 0.0;
            NekDouble SolPmeanDenumerator = 0.0;

            // Determining the norm of the numerator of the Sensor

            Vmath::Vsub(nQuadPointsElement,
                        SolPElementPhys, 1,
                        SolPmOneElementPhys, 1,
                        SolNorm, 1);

            Vmath::Vmul(nQuadPointsElement,
                        SolNorm, 1,
                        SolNorm, 1,
                        SolNorm, 1);

            for (int i = 0; i < nQuadPointsElement; i++)
            {
                SolPmeanNumerator   += SolNorm[i];
                SolPmeanDenumerator += SolPElementPhys[i];
            }

            for (int i = 0; i < nQuadPointsElement; ++i)
            {
                Sensor[CoeffsCount+i] =
                    sqrt(SolPmeanNumerator / nQuadPointsElement) /
                    sqrt(SolPmeanDenumerator / nQuadPointsElement);

                Sensor[CoeffsCount+i] = log10(Sensor[CoeffsCount+i]);
            }
            CoeffsCount += nQuadPointsElement;
        }

        CoeffsCount = 0.0;

        for (e = 0; e < nElements; e++)
        {
            NumModesElement    = ExpOrderElement[e];
            NekDouble ThetaS   = m_mu0;
            NekDouble Phi0     = m_Skappa;
            NekDouble DeltaPhi = m_Kappa;
            nQuadPointsElement = m_fields[0]->GetExp(e)->GetTotPoints();

            for (int i = 0; i < nQuadPointsElement; i++)
            {
                if (Sensor[CoeffsCount+i] <= (Phi0 - DeltaPhi))
                {
                    SensorKappa[CoeffsCount+i] = 0;
                }
                else if(Sensor[CoeffsCount+i] >= (Phi0 + DeltaPhi))
                {
                    SensorKappa[CoeffsCount+i] = ThetaS;
                }
                else if(abs(Sensor[CoeffsCount+i]-Phi0) < DeltaPhi)
                {
                    SensorKappa[CoeffsCount+i] =
                        ThetaS / 2 * (1 + sin(M_PI * (Sensor[CoeffsCount+i] -
                                                      Phi0) / (2 * DeltaPhi)));
                }
            }

            CoeffsCount += nQuadPointsElement;
        }

    }

    void CompressibleFlowSystem::v_ExtraFldOutput(
        std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
        std::vector<std::string>             &variables)
    {
        bool extraFields;
        m_session->MatchSolverInfo("OutputExtraFields","True",
                                   extraFields, true);
        if (extraFields)
        {
            const int nPhys   = m_fields[0]->GetNpoints();
            const int nCoeffs = m_fields[0]->GetNcoeffs();
            Array<OneD, Array<OneD, NekDouble> > tmp(m_fields.num_elements());

            for (int i = 0; i < m_fields.num_elements(); ++i)
            {
                tmp[i] = m_fields[i]->GetPhys();
            }

            Array<OneD, NekDouble> pressure(nPhys), soundspeed(nPhys), mach(nPhys);
            Array<OneD, NekDouble> sensor(nPhys), SensorKappa(nPhys);

            m_varConv->GetPressure  (tmp, pressure);
            m_varConv->GetSoundSpeed(tmp, pressure, soundspeed);
            m_varConv->GetMach      (tmp, soundspeed, mach);
            GetSensor    (tmp, sensor, SensorKappa);

            Array<OneD, NekDouble> pFwd(nCoeffs), sFwd(nCoeffs), mFwd(nCoeffs);
            Array<OneD, NekDouble> sensFwd(nCoeffs);

            m_fields[0]->FwdTrans(pressure,   pFwd);
            m_fields[0]->FwdTrans(soundspeed, sFwd);
            m_fields[0]->FwdTrans(mach,       mFwd);
            m_fields[0]->FwdTrans(sensor,     sensFwd);

            variables.push_back  ("p");
            variables.push_back  ("a");
            variables.push_back  ("Mach");
            variables.push_back  ("Sensor");
            fieldcoeffs.push_back(pFwd);
            fieldcoeffs.push_back(sFwd);
            fieldcoeffs.push_back(mFwd);
            fieldcoeffs.push_back(sensFwd);
        }
    }
}
