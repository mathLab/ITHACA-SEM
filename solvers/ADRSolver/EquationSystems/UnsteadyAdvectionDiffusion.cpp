///////////////////////////////////////////////////////////////////////////////
//
// File UnsteadyAdvectionDiffusion.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: Unsteady advection-diffusion solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <boost/core/ignore_unused.hpp>

#include <ADRSolver/EquationSystems/UnsteadyAdvectionDiffusion.h>

#include <LibUtilities/TimeIntegration/TimeIntegrationSolution.h>

using namespace std;

namespace Nektar
{
    string UnsteadyAdvectionDiffusion::className
        = SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
                "UnsteadyAdvectionDiffusion",
                UnsteadyAdvectionDiffusion::create);

    UnsteadyAdvectionDiffusion::UnsteadyAdvectionDiffusion(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const SpatialDomains::MeshGraphSharedPtr& pGraph)
        : UnsteadySystem(pSession, pGraph),
          AdvectionSystem(pSession, pGraph)
    {
        m_planeNumber = 0;
    }

    /**
     * @brief Initialisation object for the unsteady linear advection
     * diffusion equation.
     */
    void UnsteadyAdvectionDiffusion::v_InitObject()
    {
        AdvectionSystem::v_InitObject();

        m_session->LoadParameter("wavefreq",   m_waveFreq, 0.0);
        m_session->LoadParameter("epsilon",    m_epsilon,  0.0);

        // turn on substepping
        m_session->MatchSolverInfo("Extrapolation", "SubStepping",
                                   m_subSteppingScheme, false);

        // Define Velocity fields
        m_velocity = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
        std::vector<std::string> vel;
        vel.push_back("Vx");
        vel.push_back("Vy");
        vel.push_back("Vz");
        vel.resize(m_spacedim);

        GetFunction( "AdvectionVelocity")->Evaluate(vel,  m_velocity);

        m_session->MatchSolverInfo(
            "SpectralVanishingViscosity", "True", m_useSpecVanVisc, false);

        if(m_useSpecVanVisc)
        {
            m_session->LoadParameter("SVVCutoffRatio",m_sVVCutoffRatio,0.75);
            m_session->LoadParameter("SVVDiffCoeff",m_sVVDiffCoeff,0.1);
        }

        // Type of advection and diffusion classes to be used
        switch(m_projectionType)
        {
            // Discontinuous field
            case MultiRegions::eDiscontinuous:
            {
                // Do not forwards transform initial condition
                m_homoInitialFwd = false;

                // Advection term
                string advName;
                string riemName;
                m_session->LoadSolverInfo("AdvectionType", advName, "WeakDG");
                m_advObject = SolverUtils::GetAdvectionFactory().
                    CreateInstance(advName, advName);
                m_advObject->SetFluxVector(&UnsteadyAdvectionDiffusion::
                                           GetFluxVectorAdv, this);
                m_session->LoadSolverInfo("UpwindType", riemName, "Upwind");
                m_riemannSolver = SolverUtils::GetRiemannSolverFactory().
                    CreateInstance(riemName, m_session);
                m_riemannSolver->SetScalar("Vn", &UnsteadyAdvectionDiffusion::
                                           GetNormalVelocity, this);
                m_advObject->SetRiemannSolver(m_riemannSolver);
                m_advObject->InitObject      (m_session, m_fields);

                // Diffusion term
                if (m_explicitDiffusion)
                {
                    std::string diffName;
                    m_session->LoadSolverInfo("DiffusionType", diffName, "LDG");
                    m_diffusion = SolverUtils::GetDiffusionFactory().
                        CreateInstance(diffName, diffName);
                    m_diffusion->SetFluxVector(&UnsteadyAdvectionDiffusion::
                                               GetFluxVectorDiff, this);
                    m_diffusion->InitObject(m_session, m_fields);
                }

                ASSERTL0(m_subSteppingScheme == false,"SubSteppingScheme is not set up for DG projection");
                break;
            }
            // Continuous field
            case MultiRegions::eGalerkin:
            case MultiRegions::eMixed_CG_Discontinuous:
            {
                // Advection term
                std::string advName;
                m_session->LoadSolverInfo("AdvectionType", advName,
                                          "NonConservative");
                m_advObject = SolverUtils::GetAdvectionFactory().
                    CreateInstance(advName, advName);
                m_advObject->SetFluxVector(&UnsteadyAdvectionDiffusion::
                                           GetFluxVectorAdv, this);

                if(advName.compare("WeakDG") == 0)
                {
                    string riemName;
                    m_session->LoadSolverInfo("UpwindType", riemName, "Upwind");
                    m_riemannSolver = SolverUtils::GetRiemannSolverFactory().
                        CreateInstance(riemName, m_session);
                    m_riemannSolver->SetScalar("Vn",
                                               &UnsteadyAdvectionDiffusion::
                                               GetNormalVelocity, this);
                    m_advObject->SetRiemannSolver(m_riemannSolver);
                    m_advObject->InitObject      (m_session, m_fields);
                }

                // In case of Galerkin explicit diffusion gives an error
                if (m_explicitDiffusion)
                {
                    ASSERTL0(false, "Explicit Galerkin diffusion not set up.");
                }
                // In case of Galerkin implicit diffusion: do nothing
                break;
            }
            default:
            {
                ASSERTL0(false, "Unsupported projection type.");
                break;
            }
        }

        m_ode.DefineImplicitSolve (
            &UnsteadyAdvectionDiffusion::DoImplicitSolve, this);
        m_ode.DefineProjection(
            &UnsteadyAdvectionDiffusion::DoOdeProjection, this);
        m_ode.DefineOdeRhs(&UnsteadyAdvectionDiffusion::DoOdeRhs, this);

        if(m_subSteppingScheme) // Substepping
        {
            ASSERTL0(m_projectionType == MultiRegions::eMixed_CG_Discontinuous,
                     "Projection must be set to Mixed_CG_Discontinuous for "
                     "substepping");
            SetUpSubSteppingTimeIntegration(m_intScheme);
        }
    }

    /**
     * @brief Unsteady linear advection diffusion equation destructor.
     */
    UnsteadyAdvectionDiffusion::~UnsteadyAdvectionDiffusion()
    {
    }

    /**
     * @brief Get the normal velocity for the unsteady linear advection
     * diffusion equation.
     */
    Array<OneD, NekDouble> &UnsteadyAdvectionDiffusion::GetNormalVelocity()
    {
        return GetNormalVel(m_velocity);
    }


    Array<OneD, NekDouble> &UnsteadyAdvectionDiffusion::GetNormalVel(
        const Array<OneD, const Array<OneD, NekDouble> > &velfield)
    {
        // Number of trace (interface) points
        int i;
        int nTracePts = GetTraceNpoints();

        // Auxiliary variable to compute the normal velocity
        Array<OneD, NekDouble> tmp(nTracePts);
        m_traceVn = Array<OneD, NekDouble>(nTracePts, 0.0);

        // Reset the normal velocity
        Vmath::Zero(nTracePts, m_traceVn, 1);

        for (i = 0; i < velfield.size(); ++i)
        {
            m_fields[0]->ExtractTracePhys(velfield[i], tmp);

            Vmath::Vvtvp(nTracePts,
                         m_traceNormals[i], 1,
                         tmp, 1,
                         m_traceVn, 1,
                         m_traceVn, 1);
        }

        return m_traceVn;
    }

    /**
     * @brief Compute the right-hand side for the unsteady linear advection
     * diffusion problem.
     *
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     */
    void UnsteadyAdvectionDiffusion::DoOdeRhs(
        const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
              Array<OneD,        Array<OneD, NekDouble> >&outarray,
        const NekDouble time)
    {
        // Number of fields (variables of the problem)
        int nVariables = inarray.size();

        // Number of solution points
        int nSolutionPts = GetNpoints();

        // RHS computation using the new advection base class
        m_advObject->Advect(nVariables, m_fields, m_velocity,
                            inarray, outarray, time);

        // Negate the RHS
        for (int i = 0; i < nVariables; ++i)
        {
            Vmath::Neg(nSolutionPts, outarray[i], 1);
        }

        if (m_explicitDiffusion)
        {
            Array<OneD, Array<OneD, NekDouble> > outarrayDiff(nVariables);
            for (int i = 0; i < nVariables; ++i)
            {
                outarrayDiff[i] = Array<OneD, NekDouble>(nSolutionPts, 0.0);
            }

            m_diffusion->Diffuse(nVariables, m_fields, inarray, outarrayDiff);

            for (int i = 0; i < nVariables; ++i)
            {
                Vmath::Vadd(nSolutionPts, &outarrayDiff[i][0], 1,
                    &outarray[i][0], 1, &outarray[i][0], 1);
            }
        }

    }

    /**
     * @brief Compute the projection for the unsteady advection
     * diffusion problem.
     *
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     */
    void UnsteadyAdvectionDiffusion::DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
        const NekDouble time)
    {
        int i;
        int nvariables = inarray.size();
        SetBoundaryConditions(time);
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
                break;
            }
            case MultiRegions::eGalerkin:
            case MultiRegions::eMixed_CG_Discontinuous:
            {
                Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs());

                for(i = 0; i < nvariables; ++i)
                {
                    m_fields[i]->FwdTrans(inarray[i], coeffs);
                    m_fields[i]->BwdTrans_IterPerExp(coeffs, outarray[i]);
                }
                break;
            }
            default:
            {
                ASSERTL0(false, "Unknown projection scheme");
                break;
            }
        }
    }


    /* @brief Compute the diffusion term implicitly.
     *
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     * @param lambda     Diffusion coefficient.
     */
    void UnsteadyAdvectionDiffusion::DoImplicitSolve(
        const Array<OneD, const Array<OneD, NekDouble> >&inarray,
              Array<OneD,       Array<OneD, NekDouble> >&outarray,
        const NekDouble time,
        const NekDouble lambda)
    {
        int nvariables = inarray.size();
        int nq = m_fields[0]->GetNpoints();

        StdRegions::ConstFactorMap factors;
        factors[StdRegions::eFactorLambda] = 1.0/lambda/m_epsilon;

        if(m_useSpecVanVisc)
        {
            factors[StdRegions::eFactorSVVCutoffRatio] = m_sVVCutoffRatio;
            factors[StdRegions::eFactorSVVDiffCoeff]   = m_sVVDiffCoeff/m_epsilon;
        }
        if(m_projectionType == MultiRegions::eDiscontinuous)
        {
            factors[StdRegions::eFactorTau] = 1.0;
        }

        Array<OneD, Array< OneD, NekDouble> > F(nvariables);
        for (int n = 0; n < nvariables; ++n)
        {
            F[n] = Array<OneD, NekDouble> (nq);
        }

        // We solve ( \nabla^2 - HHlambda ) Y[i] = rhs [i]
        // inarray = input: \hat{rhs} -> output: \hat{Y}
        // outarray = output: nabla^2 \hat{Y}
        // where \hat = modal coeffs
        for (int i = 0; i < nvariables; ++i)
        {
            // Multiply 1.0/timestep/lambda
            Vmath::Smul(nq, -factors[StdRegions::eFactorLambda],
                        inarray[i], 1, F[i], 1);
        }

        //Setting boundary conditions
        SetBoundaryConditions(time);

        for (int i = 0; i < nvariables; ++i)
        {
            // Solve a system of equations with Helmholtz solver
            m_fields[i]->HelmSolve(F[i], m_fields[i]->UpdateCoeffs(), factors);

            m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(), outarray[i]);
        }
    }

    /**
     * @brief Return the flux vector for the advection part.
     *
     * @param physfield   Fields.
     * @param flux        Resulting flux.
     */
    void UnsteadyAdvectionDiffusion::GetFluxVectorAdv(
        const Array<OneD, Array<OneD, NekDouble> >               &physfield,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux)
    {
        ASSERTL1(flux[0].size() == m_velocity.size(),
                 "Dimension of flux array and velocity array do not match");

        const int nq = m_fields[0]->GetNpoints();

        for (int i = 0; i < flux.size(); ++i)
        {
            for (int j = 0; j < flux[0].size(); ++j)
            {
                Vmath::Vmul(nq, physfield[i], 1, m_velocity[j], 1,
                            flux[i][j], 1);
            }
        }
    }

    /**
     * @brief Return the flux vector for the diffusion part.
     *
     * @param i           Equation number.
     * @param j           Spatial direction.
     * @param physfield   Fields.
     * @param derivatives First order derivatives.
     * @param flux        Resulting flux.
     */
    void UnsteadyAdvectionDiffusion::GetFluxVectorDiff(
        const Array<OneD, Array<OneD, NekDouble> > &inarray,
        const Array<OneD, Array<OneD, Array<OneD, NekDouble> > >&qfield,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > >&viscousTensor)
    {
        boost::ignore_unused(inarray);

        unsigned int nDim = qfield.size();
        unsigned int nConvectiveFields = qfield[0].size();
        unsigned int nPts = qfield[0][0].size();
        for (unsigned int j = 0; j < nDim; ++j)
        {
            for (unsigned int i = 0; i < nConvectiveFields; ++i)
            {
                Vmath::Smul(nPts, m_epsilon, qfield[j][i], 1,
                    viscousTensor[j][i], 1 );
            }
        }
    }

    void UnsteadyAdvectionDiffusion::v_GenerateSummary(
            SolverUtils::SummaryList& s)
    {
        AdvectionSystem::v_GenerateSummary(s);
    }


    /**
     * Perform the extrapolation.
     */
    bool UnsteadyAdvectionDiffusion::v_PreIntegrate(int step)
    {
        if(m_subSteppingScheme)
        {
            SubStepAdvance(m_intSoln,step,m_time);
        }

        return false;
    }


    /**
     *
     */
    void UnsteadyAdvectionDiffusion::SubStepAdvance(
        const LibUtilities::TimeIntegrationScheme::TimeIntegrationSolutionSharedPtr &integrationSoln,
        int nstep,
        NekDouble time)
    {
        int n;
        int nsubsteps;

        NekDouble dt;

        Array<OneD, Array<OneD, NekDouble> > fields, velfields;

        static int ncalls = 1;
        int  nint         = min(ncalls++, m_intSteps);

        Array<OneD, NekDouble> CFL(m_fields[0]->GetExpSize(),
                                   m_cflSafetyFactor);

        LibUtilities::CommSharedPtr comm = m_session->GetComm();

        // Get the proper time step with CFL control
        dt = GetSubstepTimeStep();

        nsubsteps = (m_timestep > dt)? ((int)(m_timestep/dt)+1):1;
        nsubsteps = max(m_minsubsteps, nsubsteps);

        dt = m_timestep/nsubsteps;

        if (m_infosteps && !((nstep+1)%m_infosteps) && comm->GetRank() == 0)
        {
            cout << "Sub-integrating using "<< nsubsteps
                 << " steps over Dt = "     << m_timestep
                 << " (SubStep CFL="        << m_cflSafetyFactor << ")"<< endl;
        }

        for (int m = 0; m < nint; ++m)
        {
            // We need to update the fields held by the m_integrationSoln
            fields = integrationSoln->UpdateSolutionVector()[m];

            // Initialise NS solver which is set up to use a GLM method
            // with calls to EvaluateAdvection_SetPressureBCs and
            // SolveUnsteadyStokesSystem
            LibUtilities::TimeIntegrationScheme::TimeIntegrationSolutionSharedPtr
                SubIntegrationSoln = m_subStepIntegrationScheme->
                InitializeScheme(dt, fields, time, m_subStepIntegrationOps);

            for(n = 0; n < nsubsteps; ++n)
            {
                fields = m_subStepIntegrationScheme->TimeIntegrate(
                                                n, dt, SubIntegrationSoln,
                                                m_subStepIntegrationOps );
            }

            // Reset time integrated solution in m_integrationSoln
            integrationSoln->SetSolVector(m,fields);
        }
    }


    /**
     *
     */
    NekDouble UnsteadyAdvectionDiffusion::GetSubstepTimeStep()
    {
        int n_element      = m_fields[0]->GetExpSize();

        const Array<OneD, int> ExpOrder=m_fields[0]->EvalBasisNumModesMaxPerExp();
        Array<OneD, int> ExpOrderList (n_element, ExpOrder);

        const NekDouble cLambda = 0.2; // Spencer book pag. 317

        Array<OneD, NekDouble> tstep      (n_element, 0.0);
        Array<OneD, NekDouble> stdVelocity(n_element, 0.0);

        stdVelocity = GetMaxStdVelocity(m_velocity);

        for(int el = 0; el < n_element; ++el)
        {
            tstep[el] = m_cflSafetyFactor /
                (stdVelocity[el] * cLambda *
                 (ExpOrder[el]-1) * (ExpOrder[el]-1));
        }

        NekDouble TimeStep = Vmath::Vmin(n_element, tstep, 1);
        m_session->GetComm()->AllReduce(TimeStep,LibUtilities::ReduceMin);

        return TimeStep;
    }

    void UnsteadyAdvectionDiffusion::SetUpSubSteppingTimeIntegration(
        const LibUtilities::TimeIntegrationSchemeSharedPtr &IntegrationScheme)
    {
        // Set to 1 for first step and it will then be increased in
        // time advance routines
        unsigned int order = IntegrationScheme->GetOrder();

        // Set to 1 for first step and it will then be increased in
        // time advance routines
        if( IntegrationScheme->GetName() == "BackwardEuler" ||
            (IntegrationScheme->GetName() == "BDFImplicit" &&
             (order == 1 || order == 2)) )
        {
            // Note RK first order SSP is just Forward Euler.
            m_subStepIntegrationScheme =
                LibUtilities::GetTimeIntegrationSchemeFactory()
	            .CreateInstance( "RungeKutta", "SSP", order,
				     std::vector<NekDouble>());
        }
        else
        {
            NEKERROR(ErrorUtil::efatal, "Integration method not suitable: "
                     "Options include BackwardEuler or BDFImplicitOrder1");
        }

        m_intSteps = IntegrationScheme->GetNumIntegrationPhases();

        // set explicit time-integration class operators
        m_subStepIntegrationOps.DefineOdeRhs(&UnsteadyAdvectionDiffusion::SubStepAdvection, this);
        m_subStepIntegrationOps.DefineProjection(&UnsteadyAdvectionDiffusion::SubStepProjection, this);
    }

    /**
     * Explicit Advection terms used by SubStepAdvance time integration
     */
    void UnsteadyAdvectionDiffusion::SubStepAdvection(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
        Array<OneD, Array<OneD,       NekDouble> > &outarray,
        const NekDouble time)
    {
        int i;
        int nVariables     = inarray.size();

        /// Get the number of coefficients
        int ncoeffs = m_fields[0]->GetNcoeffs();

        /// Define an auxiliary variable to compute the RHS
        Array<OneD, Array<OneD, NekDouble> > WeakAdv(nVariables);
        WeakAdv[0] = Array<OneD, NekDouble> (ncoeffs*nVariables);
        for(i = 1; i < nVariables; ++i)
        {
            WeakAdv[i] = WeakAdv[i-1] + ncoeffs;
        }

        // Currently assume velocity field is time independent and does not therefore
        // need extrapolating.
        // RHS computation using the advection base class
        m_advObject->Advect(nVariables, m_fields, m_velocity,
                            inarray, outarray, time);

        for(i = 0; i < nVariables; ++i)
        {
            m_fields[i]->IProductWRTBase(outarray[i],WeakAdv[i]);
            // negation requried due to sign of DoAdvection term to be consistent
            Vmath::Neg(ncoeffs, WeakAdv[i], 1);
        }

        AddAdvectionPenaltyFlux(m_velocity, inarray, WeakAdv);


        /// Operations to compute the RHS
        for(i = 0; i < nVariables; ++i)
        {
            // Negate the RHS
            Vmath::Neg(ncoeffs, WeakAdv[i], 1);

            /// Multiply the flux by the inverse of the mass matrix
            m_fields[i]->MultiplyByElmtInvMass(WeakAdv[i], WeakAdv[i]);

            /// Store in outarray the physical values of the RHS
            m_fields[i]->BwdTrans(WeakAdv[i], outarray[i]);
        }
    }

    /**
     * Projection used by SubStepAdvance time integration
     */
    void UnsteadyAdvectionDiffusion::SubStepProjection(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
        Array<OneD, Array<OneD, NekDouble> > &outarray,
        const NekDouble time)
    {
        boost::ignore_unused(time);

        ASSERTL1(inarray.size() == outarray.size(),"Inarray and outarray of different sizes ");

        for(int i = 0; i < inarray.size(); ++i)
        {
            Vmath::Vcopy(inarray[i].size(),inarray[i],1,outarray[i],1);
        }
    }

    void UnsteadyAdvectionDiffusion::AddAdvectionPenaltyFlux(
                                const Array<OneD, const Array<OneD, NekDouble> > &velfield,
                                const Array<OneD, const Array<OneD, NekDouble> > &physfield,
                                Array<OneD, Array<OneD, NekDouble> > &Outarray)
    {
        ASSERTL1(physfield.size() == Outarray.size(),
                 "Physfield and outarray are of different dimensions");

        int i;

        /// Number of trace points
        int nTracePts   = m_fields[0]->GetTrace()->GetNpoints();

        /// Forward state array
        Array<OneD, NekDouble> Fwd(3*nTracePts);

        /// Backward state array
        Array<OneD, NekDouble> Bwd = Fwd + nTracePts;

        /// upwind numerical flux state array
        Array<OneD, NekDouble> numflux = Bwd + nTracePts;

        /// Normal velocity array
        Array<OneD, NekDouble> Vn  = GetNormalVel(velfield);

        for(i = 0; i < physfield.size(); ++i)
        {
            /// Extract forwards/backwards trace spaces
            /// Note: Needs to have correct i value to get boundary conditions
            m_fields[i]->GetFwdBwdTracePhys(physfield[i], Fwd, Bwd);

            /// Upwind between elements
            m_fields[0]->GetTrace()->Upwind(Vn, Fwd, Bwd, numflux);

            /// Construct difference between numflux and Fwd,Bwd
            Vmath::Vsub(nTracePts, numflux, 1, Fwd, 1, Fwd, 1);
            Vmath::Vsub(nTracePts, numflux, 1, Bwd, 1, Bwd, 1);

            /// Calculate the numerical fluxes multipling Fwd, Bwd and
            /// numflux by the normal advection velocity
            Vmath::Vmul(nTracePts, Fwd, 1, Vn, 1, Fwd, 1);
            Vmath::Vmul(nTracePts, Bwd, 1, Vn, 1, Bwd, 1);

            m_fields[0]->AddFwdBwdTraceIntegral(Fwd,Bwd,Outarray[i]);
        }
    }


    Array<OneD, NekDouble> UnsteadyAdvectionDiffusion::GetMaxStdVelocity(
        const Array<OneD, Array<OneD,NekDouble> > inarray)
    {

        int n_points_0      = m_fields[0]->GetExp(0)->GetTotPoints();
        int n_element       = m_fields[0]->GetExpSize();
        int nvel            = inarray.size();
        int cnt;

        ASSERTL0(nvel >= 2, "Method not implemented for 1D");

        NekDouble pntVelocity;

        // Getting the standard velocity vector on the 2D normal space
        Array<OneD, Array<OneD, NekDouble> > stdVelocity(nvel);
        Array<OneD, NekDouble> maxV(n_element, 0.0);
        LibUtilities::PointsKeyVector ptsKeys;

        for (int i = 0; i < nvel; ++i)
        {
            stdVelocity[i] = Array<OneD, NekDouble>(n_points_0);
        }

        if (nvel == 2)
        {
            cnt = 0.0;
            for (int el = 0; el < n_element; ++el)
            {
                int n_points = m_fields[0]->GetExp(el)->GetTotPoints();
                ptsKeys = m_fields[0]->GetExp(el)->GetPointsKeys();

                // reset local space if necessary
                if(n_points != n_points_0)
                {
                    for (int j = 0; j < nvel; ++j)
                    {
                        stdVelocity[j] = Array<OneD, NekDouble>(n_points);
                    }
                    n_points_0 = n_points;
                }

                Array<TwoD, const NekDouble> gmat =
                    m_fields[0]->GetExp(el)->GetGeom()->GetMetricInfo()->GetDerivFactors(ptsKeys);

                if (m_fields[0]->GetExp(el)->GetGeom()->GetMetricInfo()->GetGtype()
                    == SpatialDomains::eDeformed)
                {
                    for (int i = 0; i < n_points; i++)
                    {
                        stdVelocity[0][i] = gmat[0][i]*inarray[0][i+cnt]
                            + gmat[2][i]*inarray[1][i+cnt];

                        stdVelocity[1][i] = gmat[1][i]*inarray[0][i+cnt]
                            + gmat[3][i]*inarray[1][i+cnt];
                    }
                }
                else
                {
                    for (int i = 0; i < n_points; i++)
                    {
                        stdVelocity[0][i] = gmat[0][0]*inarray[0][i+cnt]
                            + gmat[2][0]*inarray[1][i+cnt];

                        stdVelocity[1][i] = gmat[1][0]*inarray[0][i+cnt]
                            + gmat[3][0]*inarray[1][i+cnt];
                    }
                }

                cnt += n_points;


                for (int i = 0; i < n_points; i++)
                {
                    pntVelocity = stdVelocity[0][i]*stdVelocity[0][i]
                        + stdVelocity[1][i]*stdVelocity[1][i];

                    if (pntVelocity>maxV[el])
                    {
                        maxV[el] = pntVelocity;
                    }
                }
                maxV[el] = sqrt(maxV[el]);
            }
        }
        else
        {
            cnt = 0;
            for (int el = 0; el < n_element; ++el)
            {

                int n_points = m_fields[0]->GetExp(el)->GetTotPoints();
                ptsKeys = m_fields[0]->GetExp(el)->GetPointsKeys();

                // reset local space if necessary
                if(n_points != n_points_0)
                {
                    for (int j = 0; j < nvel; ++j)
                    {
                        stdVelocity[j] = Array<OneD, NekDouble>(n_points);
                    }
                    n_points_0 = n_points;
                }

                Array<TwoD, const NekDouble> gmat =
                    m_fields[0]->GetExp(el)->GetGeom()->GetMetricInfo()->GetDerivFactors(ptsKeys);

                if (m_fields[0]->GetExp(el)->GetGeom()->GetMetricInfo()->GetGtype()
                    == SpatialDomains::eDeformed)
                {
                    for (int i = 0; i < n_points; i++)
                    {
                        stdVelocity[0][i] = gmat[0][i]*inarray[0][i+cnt]
                            + gmat[3][i]*inarray[1][i+cnt]
                            + gmat[6][i]*inarray[2][i+cnt];

                        stdVelocity[1][i] = gmat[1][i]*inarray[0][i+cnt]
                            + gmat[4][i]*inarray[1][i+cnt]
                            + gmat[7][i]*inarray[2][i+cnt];

                        stdVelocity[2][i] = gmat[2][i]*inarray[0][i+cnt]
                            + gmat[5][i]*inarray[1][i+cnt]
                            + gmat[8][i]*inarray[2][i+cnt];
                    }
                }
                else
                {
                    for (int i = 0; i < n_points; i++)
                    {
                        stdVelocity[0][i] = gmat[0][0]*inarray[0][i+cnt]
                            + gmat[3][0]*inarray[1][i+cnt]
                            + gmat[6][0]*inarray[2][i+cnt];

                        stdVelocity[1][i] = gmat[1][0]*inarray[0][i+cnt]
                            + gmat[4][0]*inarray[1][i+cnt]
                            + gmat[7][0]*inarray[2][i+cnt];

                        stdVelocity[2][i] = gmat[2][0]*inarray[0][i+cnt]
                            + gmat[5][0]*inarray[1][i+cnt]
                            + gmat[8][0]*inarray[2][i+cnt];
                    }
                }

                cnt += n_points;

                for (int i = 0; i < n_points; i++)
                {
                    pntVelocity = stdVelocity[0][i]*stdVelocity[0][i]
                        + stdVelocity[1][i]*stdVelocity[1][i]
                        + stdVelocity[2][i]*stdVelocity[2][i];

                    if (pntVelocity > maxV[el])
                    {
                        maxV[el] = pntVelocity;
                    }
                }

                maxV[el] = sqrt(maxV[el]);
                //cout << maxV[el]*maxV[el] << endl;
            }
        }

        return maxV;
    }
}
