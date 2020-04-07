/////////////////////////////////////////////////////////////////////////////
//
// File MMFAdvection.cpp
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
// Description: MMF solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include <iomanip>
#include <iostream>

#include <boost/core/ignore_unused.hpp>
#include <boost/algorithm/string.hpp>

#include <ADRSolver/EquationSystems/MMFAdvection.h>
#include <LibUtilities/BasicUtils/Timer.h>

namespace Nektar
{
std::string MMFAdvection::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "MMFAdvection", MMFAdvection::create, "MMFAdvection equation.");

MMFAdvection::MMFAdvection(const LibUtilities::SessionReaderSharedPtr &pSession,
            const SpatialDomains::MeshGraphSharedPtr& pGraph)
    : UnsteadySystem(pSession, pGraph), MMFSystem(pSession, pGraph),
      AdvectionSystem(pSession, pGraph)
{
    m_planeNumber = 0;
}

/**
 * @brief Initialisation object for the unsteady linear advection equation.
 */
void MMFAdvection::v_InitObject()
{
    // Call to the initialisation object
    UnsteadySystem::v_InitObject();

    int nq       = m_fields[0]->GetNpoints();
    int shapedim = m_fields[0]->GetShapeDimension();
    Array<OneD, Array<OneD, NekDouble>> Anisotropy(shapedim);
    for (int j = 0; j < shapedim; ++j)
    {
        Anisotropy[j] = Array<OneD, NekDouble>(nq, 1.0);
    }

    MMFSystem::MMFInitObject(Anisotropy);

    // Define TestType
    ASSERTL0(m_session->DefinesSolverInfo("TESTTYPE"),
             "No TESTTYPE defined in session.");
    std::string TestTypeStr = m_session->GetSolverInfo("TESTTYPE");
    for (int i = 0; i < (int)SIZE_TestType; ++i)
    {
        if (boost::iequals(TestTypeMap[i], TestTypeStr))
        {
            m_TestType = (TestType)i;
            break;
        }
    }

    m_session->LoadParameter("Angular Frequency", m_waveFreq, m_pi);
    m_session->LoadParameter("Rotational Angle", m_RotAngle, 0.0);
    m_session->LoadParameter("Velocity Projection", m_VelProjection, 0);

    // Read the advection velocities from session file
    m_session->LoadParameter("advx", m_advx, 1.0);
    m_session->LoadParameter("advy", m_advy, 1.0);
    m_session->LoadParameter("advz", m_advz, 1.0);

    std::vector<std::string> vel;
    vel.push_back("Vx");
    vel.push_back("Vy");
    vel.push_back("Vz");

    // Resize the advection velocities vector to dimension of the problem
    vel.resize(m_spacedim);

    // Store in the global variable m_velocity the advection velocities

    m_velocity = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
    for (int k = 0; k < m_spacedim; ++k)
    {
        m_velocity[k] = Array<OneD, NekDouble>(nq);
    }

    switch (m_surfaceType)
    {
        case SolverUtils::eSphere:
        case SolverUtils::eTRSphere:
        case SolverUtils::eIrregular:
        case SolverUtils::eNonconvex:
        {
            // true = project velocity onto the tangent plane
            EvaluateAdvectionVelocity(m_velocity);
        }
        break;

        case SolverUtils::ePlane:
        case SolverUtils::eCube:
        {
            GetFunction("AdvectionVelocity")->Evaluate(vel, m_velocity);
        }
        break;

        default:
            break;
    }

    std::cout << "|Velocity vector| = ( " << RootMeanSquare(m_velocity[0])
             << " , " << RootMeanSquare(m_velocity[1]) << " , "
             << RootMeanSquare(m_velocity[2]) << " ) " << std::endl;

    // Define the normal velocity fields
    if (m_fields[0]->GetTrace())
    {
        m_traceVn = Array<OneD, NekDouble>(GetTraceNpoints());
    }

    std::string advName;
    std::string riemName;
    m_session->LoadSolverInfo(
        "AdvectionType", advName, "WeakDG");
    m_advObject = SolverUtils::
        GetAdvectionFactory().CreateInstance(advName, advName);
    m_advObject->SetFluxVector(
        &MMFAdvection::GetFluxVector, this);
    m_session->LoadSolverInfo(
        "UpwindType", riemName, "Upwind");
    m_riemannSolver = SolverUtils::
        GetRiemannSolverFactory().CreateInstance(riemName, m_session);
    m_riemannSolver->SetScalar(
        "Vn", &MMFAdvection::GetNormalVelocity, this);

    m_advObject->SetRiemannSolver(m_riemannSolver);
    m_advObject->InitObject(m_session, m_fields);

    // Compute m_traceVn = n \cdot v
    GetNormalVelocity();

    // Compute m_vellc = nabal a^j \cdot m_vel
    ComputeNablaCdotVelocity(m_vellc);
    std::cout << "m_vellc is generated with mag = " << AvgInt(m_vellc)
             << std::endl;

    // Compute vel \cdot MF
    ComputeveldotMF(m_veldotMF);

    // Modify e^i as v^i e^i
    for (int j = 0; j < m_shapedim; j++)
    {
        for (int k = 0; k < m_spacedim; k++)
        {
            Vmath::Vmul(nq, &m_veldotMF[j][0], 1, &m_movingframes[j][k * nq], 1,
                        &m_movingframes[j][k * nq], 1);
        }
    }

    // Reflect it into m_ncdotMFFwd and Bwd
    ComputencdotMF();

    // If explicit it computes RHS and PROJECTION for the time integration
    if (m_explicitAdvection)
    {
        m_ode.DefineOdeRhs(&MMFAdvection::DoOdeRhs, this);
        m_ode.DefineProjection(&MMFAdvection::DoOdeProjection, this);
    }
    // Otherwise it gives an error (no implicit integration)
    else
    {
        ASSERTL0(false, "Implicit unsteady Advection not set up.");
    }
}

/**
 * @brief Unsteady linear advection equation destructor.
 */
MMFAdvection::~MMFAdvection()
{
}

void MMFAdvection::v_DoSolve()
{
    ASSERTL0(m_intScheme != 0, "No time integration scheme.");

    int i, nchk = 1;
    int nvariables = 0;
    int nfields    = m_fields.size();
    int nq         = m_fields[0]->GetNpoints();

    if (m_intVariables.empty())
    {
        for (i = 0; i < nfields; ++i)
        {
            m_intVariables.push_back(i);
        }
        nvariables = nfields;
    }
    else
    {
        nvariables = m_intVariables.size();
    }

    // Set up wrapper to fields data storage.
    Array<OneD, Array<OneD, NekDouble>> fields(nvariables);
    Array<OneD, Array<OneD, NekDouble>> tmp(nvariables);

    // Order storage to list time-integrated fields first.
    for (i = 0; i < nvariables; ++i)
    {
        fields[i] = m_fields[m_intVariables[i]]->GetPhys();
        m_fields[m_intVariables[i]]->SetPhysState(false);
    }

    // Initialise time integration scheme
    m_intSoln = m_intScheme->InitializeScheme( m_timestep, fields, m_time, m_ode );

    // Check uniqueness of checkpoint output
    ASSERTL0((m_checktime == 0.0 && m_checksteps == 0) ||
                 (m_checktime > 0.0 && m_checksteps == 0) ||
                 (m_checktime == 0.0 && m_checksteps > 0),
             "Only one of IO_CheckTime and IO_CheckSteps "
             "should be set!");

    LibUtilities::Timer timer;
    bool doCheckTime  = false;
    int step          = 0;
    NekDouble intTime = 0.0;
    NekDouble cpuTime = 0.0;
    NekDouble elapsed = 0.0;

    int Ntot, indx;
    // Perform integration in time.
    Ntot = m_steps / m_checksteps + 1;

    Array<OneD, NekDouble> dMass(Ntot);

    Array<OneD, NekDouble> zeta(nq);
    Array<OneD, Array<OneD, NekDouble>> fieldsprimitive(nvariables);
    for (int i = 0; i < nvariables; ++i)
    {
        fieldsprimitive[i] = Array<OneD, NekDouble>(nq);
    }

    while (step < m_steps || m_time < m_fintime - NekConstants::kNekZeroTol)
    {
        timer.Start();
        fields = m_intScheme->TimeIntegrate(step, m_timestep, m_intSoln, m_ode);
        timer.Stop();

        m_time += m_timestep;
        elapsed = timer.TimePerTest(1);
        intTime += elapsed;
        cpuTime += elapsed;

        // Write out status information
        if (m_session->GetComm()->GetRank() == 0 && !((step + 1) % m_infosteps))
        {
            std::cout << "Steps: " << std::setw(8) << std::left << step + 1 << " "
                     << "Time: " << std::setw(12) << std::left << m_time;

            std::stringstream ss;
            ss << cpuTime << "s";
            std::cout << " CPU Time: " << std::setw(8) << std::left << ss.str()
                     << std::endl;

            // Masss = h^*
            indx = (step + 1) / m_checksteps;
            dMass[indx] =
                (m_fields[0]->PhysIntegral(fields[0]) - m_Mass0) / m_Mass0;

            std::cout << "dMass = " << std::setw(8) << std::left << dMass[indx]
                     << std::endl;

            cpuTime = 0.0;
        }

        // Transform data into coefficient space
        for (i = 0; i < nvariables; ++i)
        {
            m_fields[m_intVariables[i]]->SetPhys(fields[i]);
            m_fields[m_intVariables[i]]->FwdTrans_IterPerExp(
                fields[i], m_fields[m_intVariables[i]]->UpdateCoeffs());
            m_fields[m_intVariables[i]]->SetPhysState(false);
        }

        // Write out checkpoint files
        if ((m_checksteps && step && !((step + 1) % m_checksteps)) ||
            doCheckTime)
        {
            Checkpoint_Output(nchk++);
            doCheckTime = false;
        }

        // Step advance
        ++step;
    }

    std::cout << "dMass = ";
    for (i = 0; i < Ntot; ++i)
    {
        std::cout << dMass[i] << " , ";
    }
    std::cout << std::endl << std::endl;

    // Print out summary statistics
    if (m_session->GetComm()->GetRank() == 0)
    {
        if (m_cflSafetyFactor > 0.0)
        {
            std::cout << "CFL safety factor : " << m_cflSafetyFactor << std::endl
                     << "CFL time-step     : " << m_timestep << std::endl;
        }

        if (m_session->GetSolverInfo("Driver") != "SteadyState")
        {
            std::cout << "Time-integration  : " << intTime << "s" << std::endl;
        }
    }

    for (i = 0; i < nvariables; ++i)
    {
        m_fields[m_intVariables[i]]->SetPhys(fields[i]);
        m_fields[m_intVariables[i]]->SetPhysState(true);
    }

    for (i = 0; i < nvariables; ++i)
    {
        m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),
                              m_fields[i]->UpdateCoeffs());
    }
}

/**
 * @brief Get the normal velocity for the linear advection equation.
 */
Array<OneD, NekDouble> &MMFAdvection::GetNormalVelocity()
{
    // Number of trace (interface) points
    int i;
    int nTracePts = GetTraceNpoints();

    // Auxiliary variable to compute the normal velocity
    Array<OneD, NekDouble> tmp(nTracePts);

    // Reset the normal velocity
    Vmath::Zero(nTracePts, m_traceVn, 1);

    for (i = 0; i < m_velocity.size(); ++i)
    {
        m_fields[0]->ExtractTracePhys(m_velocity[i], tmp);

        Vmath::Vvtvp(nTracePts, m_traceNormals[i], 1, tmp, 1, m_traceVn, 1,
                     m_traceVn, 1);
    }

    return m_traceVn;
}

/**
 * @brief Compute the right-hand side for the linear advection equation.
 *
 * @param inarray    Given fields.
 * @param outarray   Calculated solution.
 * @param time       Time.
 */
void MMFAdvection::DoOdeRhs(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    boost::ignore_unused(time);

    int i;
    int nvariables = inarray.size();
    int npoints    = GetNpoints();

    switch (m_projectionType)
    {
        case MultiRegions::eDiscontinuous:
        {
            int ncoeffs = inarray[0].size();

            if (m_spacedim == 3)
            {
                Array<OneD, Array<OneD, NekDouble>> WeakAdv(nvariables);

                WeakAdv[0] = Array<OneD, NekDouble>(ncoeffs * nvariables);
                for (i = 1; i < nvariables; ++i)
                {
                    WeakAdv[i] = WeakAdv[i - 1] + ncoeffs;
                }

                // Compute \nabla \cdot \vel u according to MMF scheme
                WeakDGDirectionalAdvection(inarray, WeakAdv);

                for (i = 0; i < nvariables; ++i)
                {
                    m_fields[i]->MultiplyByElmtInvMass(WeakAdv[i], WeakAdv[i]);
                    m_fields[i]->BwdTrans(WeakAdv[i], outarray[i]);

                    // Add  m_vellc * inarray[i] = \nabla v^m \cdot e^m to
                    // outarray[i]
                    // Vmath::Vvtvp(npoints, &m_vellc[0], 1, &inarray[i][0], 1,
                    // &outarray[i][0], 1, &outarray[i][0], 1);
                    Vmath::Neg(npoints, outarray[i], 1);
                }
            }
            else
            {
                m_advObject->Advect(2, m_fields, m_velocity, inarray,
                                    outarray, 0.0);

                for (i = 0; i < nvariables; ++i)
                {
                    Vmath::Neg(npoints, outarray[i], 1);
                }
            }


        }
        break;

        default:
        {
            ASSERTL0(false, "Unknown projection scheme");
        }
        break;
    }
}
/**
 * @brief Compute the projection for the linear advection equation.
 *
 * @param inarray    Given fields.
 * @param outarray   Calculated solution.
 * @param time       Time.
 */
void MMFAdvection::DoOdeProjection(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    // Counter variable
    int i;

    // Number of fields (variables of the problem)
    int nVariables = inarray.size();

    // Set the boundary conditions
    SetBoundaryConditions(time);

    // Switch on the projection type (Discontinuous or Continuous)
    switch (m_projectionType)
    {
        // Discontinuous projection
        case MultiRegions::eDiscontinuous:
        {
            // Number of quadrature points
            int nQuadraturePts = GetNpoints();

            // Just copy over array
            for (i = 0; i < nVariables; ++i)
            {
                Vmath::Vcopy(nQuadraturePts, inarray[i], 1, outarray[i], 1);
            }
            break;
        }

        // Continuous projection
        case MultiRegions::eGalerkin:
        case MultiRegions::eMixed_CG_Discontinuous:
        {
            Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs(), 0.0);
            for (i = 0; i < nVariables; ++i)
            {
                m_fields[i]->FwdTrans(inarray[i], coeffs);
                m_fields[i]->BwdTrans_IterPerExp(coeffs, outarray[i]);
            }
            break;
        }

        default:
            ASSERTL0(false, "Unknown projection scheme");
            break;
    }
}

/**
 * @brief Return the flux vector for the linear advection equation.
 *
 * @param i           Component of the flux vector to calculate.
 * @param physfield   Fields.
 * @param flux        Resulting flux.
 */
void MMFAdvection::GetFluxVector(
    const Array<OneD, Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux)
{
    ASSERTL1(flux[0].size() == m_velocity.size(),
             "Dimension of flux array and velocity array do not match");

    int i, j;
    int nq = physfield[0].size();

    for (i = 0; i < flux.size(); ++i)
    {
        for (j = 0; j < flux[0].size(); ++j)
        {
            Vmath::Vmul(nq, physfield[i], 1, m_velocity[j], 1, flux[i][j], 1);
        }
    }
}

void MMFAdvection::WeakDGDirectionalAdvection(
    const Array<OneD, Array<OneD, NekDouble>> &InField,
    Array<OneD, Array<OneD, NekDouble>> &OutField)
{
    int i, j;
    int ncoeffs         = GetNcoeffs();
    int nTracePointsTot = GetTraceNpoints();
    int nvariables      = m_fields.size();

    Array<OneD, Array<OneD, NekDouble>> physfield(nvariables);

    // Get the variables in physical space
    // already in physical space
    for (i = 0; i < nvariables; ++i)
    {
        physfield[i] = InField[i];
    }

    Array<OneD, Array<OneD, NekDouble>> WeakDeriv(m_shapedim);
    for (i = 0; i < nvariables; ++i)
    {
        for (j = 0; j < m_shapedim; ++j)
        {
            WeakDeriv[j] = Array<OneD, NekDouble>(ncoeffs, 0.0);

            // Directional derivation with respect to the j'th moving frame
            // tmp[j] = \nabla \physfield[i] \cdot \mathbf{e}^j
            // Implemented at TriExp::v_IProductWRTDirectionalDerivBase_SumFa
            m_fields[0]->IProductWRTDirectionalDerivBase(
                m_movingframes[j], physfield[i], WeakDeriv[j]);
        }

        // Get the numerical flux and add to the modal coeffs
        // if the NumericalFluxs function already includes the
        // normal in the output

        Array<OneD, NekDouble> Fwd(nTracePointsTot);
        Array<OneD, NekDouble> Bwd(nTracePointsTot);

        Array<OneD, NekDouble> flux(nTracePointsTot, 0.0);
        Array<OneD, NekDouble> fluxFwd(nTracePointsTot);
        Array<OneD, NekDouble> fluxBwd(nTracePointsTot);

        // Evaluate numerical flux in physical space which may in
        // general couple all component of vectors
        m_fields[i]->GetFwdBwdTracePhys(physfield[i], Fwd, Bwd);

        // evaulate upwinded m_fields[i]
        m_fields[i]->GetTrace()->Upwind(m_traceVn, Fwd, Bwd, flux);

        OutField[i] = Array<OneD, NekDouble>(ncoeffs, 0.0);
        for (j = 0; j < m_shapedim; ++j)
        {
            // calculate numflux = (n \cdot MF)*flux
            Vmath::Vmul(nTracePointsTot, &flux[0], 1, &m_ncdotMFFwd[j][0], 1,
                        &fluxFwd[0], 1);
            Vmath::Vmul(nTracePointsTot, &flux[0], 1, &m_ncdotMFBwd[j][0], 1,
                        &fluxBwd[0], 1);
            Vmath::Neg(ncoeffs, WeakDeriv[j], 1);

            // FwdBwdtegral because generallize (N \cdot MF)_{FWD} \neq -(N
            // \cdot MF)_{BWD}
            m_fields[i]->AddFwdBwdTraceIntegral(fluxFwd, fluxBwd, WeakDeriv[j]);
            m_fields[i]->SetPhysState(false);

            Vmath::Vadd(ncoeffs, &WeakDeriv[j][0], 1, &OutField[i][0], 1,
                        &OutField[i][0], 1);
        }
    }
}

void MMFAdvection::EvaluateAdvectionVelocity(
    Array<OneD, Array<OneD, NekDouble>> &velocity)
{
    int nq = m_fields[0]->GetNpoints();

    NekDouble vel_phi, vel_theta, sin_varphi, cos_varphi, sin_theta, cos_theta;
    NekDouble x0j, x1j, x2j;

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    // theta = a*sin(z/r),  phi = a*tan(y/x);
    for (int j = 0; j < nq; j++)
    {
        x0j = x0[j];
        x1j = x1[j];
        x2j = x2[j];

        CartesianToSpherical(x0j, x1j, x2j, sin_varphi, cos_varphi, sin_theta,
                             cos_theta);

        vel_phi = m_waveFreq * (cos_theta * cos(m_RotAngle) +
                                sin_theta * cos_varphi * sin(m_RotAngle));
        vel_theta = -1.0 * m_waveFreq * sin_theta * sin(m_RotAngle);

        velocity[0][j] =
            -1.0 * vel_phi * sin_varphi - vel_theta * sin_theta * cos_varphi;
        velocity[1][j] =
            vel_phi * cos_varphi - vel_theta * sin_theta * sin_varphi;
        velocity[2][j] = vel_theta * cos_theta;
    }

    // Project the veloicty on the tangent plane

    if (m_VelProjection)
    {
        // Check MovingFrames \cdot SurfaceNormals = 0
        Array<OneD, Array<OneD, NekDouble>> newvelocity(m_spacedim);

        Array<OneD, Array<OneD, NekDouble>> MF1(m_spacedim);
        Array<OneD, Array<OneD, NekDouble>> MF2(m_spacedim);
        Array<OneD, Array<OneD, NekDouble>> SN(m_spacedim);

        for (int k = 0; k < m_spacedim; ++k)
        {
            newvelocity[k] = Array<OneD, NekDouble>(nq);
            MF1[k]         = Array<OneD, NekDouble>(nq);
            MF2[k]         = Array<OneD, NekDouble>(nq);
            SN[k]          = Array<OneD, NekDouble>(nq);

            Vmath::Vcopy(nq, &m_movingframes[0][k * nq], 1, &MF1[k][0], 1);
            Vmath::Vcopy(nq, &m_movingframes[1][k * nq], 1, &MF2[k][0], 1);
        }

        VectorCrossProd(MF1, MF2, SN);
        GramSchumitz(SN, m_velocity, newvelocity, true);

        Array<OneD, NekDouble> tmp(nq, 0.0);
        Array<OneD, NekDouble> Projection(nq, 0.0);

        for (int k = 0; k < m_spacedim; ++k)
        {
            Vmath::Vsub(nq, &m_velocity[0][0], 1, &newvelocity[0][0], 1,
                        &tmp[0], 1);
            Vmath::Vmul(nq, &tmp[0], 1, &tmp[0], 1, &tmp[0], 1);
            Vmath::Vadd(nq, tmp, 1, Projection, 1, Projection, 1);
        }

        std::cout
            << "Velocity vector is projected onto the tangent plane: Diff = "
            << RootMeanSquare(Projection) << std::endl;

        for (int k = 0; k < m_spacedim; ++k)
        {
            Vmath::Vcopy(nq, &newvelocity[k][0], 1, &m_velocity[k][0], 1);
        }
    }
}

/*
void MMFAdvection::EvaluateAdvectionVelocity(Array<OneD, Array<OneD, NekDouble>
> &velocity)
  {
    int nq = m_fields[0]->GetNpoints();

    NekDouble vel_phi, sin_phi, cos_phi;
    NekDouble vel_theta, sin_theta, cos_theta;
    NekDouble radius, xyrad, Tol = 0.00001;

    Array<OneD,NekDouble> x0(nq);
    Array<OneD,NekDouble> x1(nq);
    Array<OneD,NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0,x1,x2);

    // theta = a*sin(z/r),  phi = a*tan(y/x);
    for (int j = 0; j < nq; j++)
      {
        switch(m_surfaceType)
          {
          case SolverUtils::eSphere:
          case SolverUtils::eTRSphere:
            {
              radius = sqrt( x0[j]*x0[j] + x1[j]*x1[j] + x2[j]*x2[j] );
            }
            break;

          case SolverUtils::eIrregular:
            {
              radius = sqrt( 2.0*x0[j]*x0[j] + x1[j]*x1[j]*x1[j]*x1[j] +
x1[j]*x1[j]
                             + x2[j]*x2[j]*x2[j]*x2[j] + x2[j]*x2[j] );
            }
            break;

            // 2 x^2 + 2(y^4 - y ) + z^4 + z^2 = 2.0
          case SolverUtils::eNonconvex:
            {
              radius = sqrt( 2.0*x0[j]*x0[j] + 2.0*( x1[j]*x1[j]*x1[j]*x1[j] -
x1[j]*x1[j] )
                             + x2[j]*x2[j]*x2[j]*x2[j] + x2[j]*x2[j] );
            }
            break;

          default:
            break;
          }

            // At North and South poles, (ax,ay,ax) = (0, Omega_f*sin(alpha),0)
            if( fabs(fabs(x2[j]) - radius) < Tol )
              {
                sin_theta = x2[j]/radius;

                velocity[0][j] = 0.0;
                velocity[1][j] = 0.0;
                velocity[2][j] = 0.0;
              }
            else
              {
                // Compute the arc length of the trajectory
                NekDouble zp, velmag0, velmag, length0, zcutoff=0.0, zmax;

                zp = fabs(x2[j]);

                zmax = Vmath::Vmax(nq, x2, 1);
                velmag = ComputeCirculatingArclength(x2[j], radius*radius);

                switch(m_surfaceType)
                  {

                  case SolverUtils::eSphere:
                  case SolverUtils::eTRSphere:
                    {
                      zcutoff = 0.7;
                      velmag0 = velmag;
                    }
                    break;

                  case SolverUtils::eIrregular:
                    {
                      zcutoff = 0.7;
                      length0 = 0.88;
                      //   velmag0 = length0*( 1.0 - (zp -
zcutoff)/(1.0-zcutoff) );
                      velmag0 = (length0/(zcutoff*zcutoff-zmax*zmax))*( zp*zp -
zmax*zmax);
                    }
                    break;

                  case SolverUtils::eNonconvex:
                    {
                      zcutoff = 0.7;
                      length0 = 1.21;
                      //   velmag0 = length0*( 1.0 - (zp -
zcutoff)/(1.0-zcutoff) );
                      velmag0 = (length0/(zcutoff*zcutoff - 1.0))*( zp*zp -
1.0);
                    }
                    break;

                  default:
                    break;
                  }

                if( zp > zcutoff )
                  {
                    velmag = velmag0;
                  }

                vel_phi = m_waveFreq*velmag;
                vel_theta = 0.0;

                xyrad = sqrt( x0[j]*x0[j] + x1[j]*x1[j] );
                if(xyrad<Tol)
                  {
                    sin_phi = 0.0;
                    cos_phi = 0.0;
                  }

                else
                  {
                    sin_phi = x1[j]/xyrad;
                    cos_phi = x0[j]/xyrad;
                  }

                // sin_theta = x2[j]/radius;
                cos_theta = sqrt( x0[j]*x0[j] + x1[j]*x1[j] )/radius;

                if(fabs(velmag) < 0.001)
                  {
                    velocity[0][j] = 0.0 ;
                    velocity[1][j] = 0.0  ;
                    velocity[2][j] = 0.0;
                  }

                else
                  {
                    velocity[0][j] = -1.0*vel_phi*sin_phi ;
                    velocity[1][j] = vel_phi*cos_phi  ;
                    velocity[2][j] = 0.0;
                  }

              } // else-loop
      }

    // Project the veloicty on the tangent plane
      int nq = m_fields[0]->GetNpoints();

      // Check MovingFrames \cdot SurfaceNormals = 0
      Array<OneD, NekDouble> temp0(nq,0.0);
      Array<OneD, NekDouble> temp1(nq,0.0);
      Array<OneD, Array<OneD, NekDouble> > newvelocity(m_spacedim);

      for(int k=0; k<m_spacedim; ++k)
        {
          newvelocity[k] = Array<OneD, NekDouble>(nq);
        }

        std::cout <<"=====================================================" <<
std::endl;
        std::cout << "Velocity vector is projected onto the tangent plane " <<
std::endl;
        std::cout <<"=====================================================" <<
std::endl;
        GramSchumitz(m_surfaceNormal, m_velocity, newvelocity);

      for(int k=0; k<m_spacedim; ++k)
        {
          Vmath::Vcopy(nq, &newvelocity[k][0], 1, &m_velocity[k][0], 1);
        }
  }

*/

NekDouble MMFAdvection::ComputeCirculatingArclength(const NekDouble zlevel,
                                                    const NekDouble Rhs)
{

    NekDouble Tol = 0.0001, Maxiter = 1000, N = 100;
    NekDouble newy, F, dF, y0, tmp;

    Array<OneD, NekDouble> xp(N + 1);
    Array<OneD, NekDouble> yp(N + 1);

    NekDouble intval = 0.0;
    switch (m_surfaceType)
    {
        case SolverUtils::eSphere:
        case SolverUtils::eTRSphere:
        {
            intval = sqrt(Rhs - zlevel * zlevel);
        }
        break;

        case SolverUtils::eIrregular:
        {
            intval = sqrt(0.5 * (Rhs - zlevel * zlevel * zlevel * zlevel -
                                 zlevel * zlevel));
        }
        break;

        case SolverUtils::eNonconvex:
        {
            tmp = 0.5 *
                  (Rhs - zlevel * zlevel * zlevel * zlevel - zlevel * zlevel);
            intval = sqrt(0.5 * (1.0 + sqrt(1.0 + 4.0 * tmp)));
        }
        break;

        default:
            break;
    }

    switch (m_surfaceType)
    {
        // Find the half of all the xp and yp on zlevel ....
        case SolverUtils::eSphere:
        case SolverUtils::eTRSphere:
        case SolverUtils::eIrregular:
        {
            for (int j = 0; j < N + 1; ++j)
            {
                xp[j] = j * 2.0 * intval / N - intval;

                y0 = 1.0;
                for (int i = 0; i < Maxiter; ++i)
                {
                    switch (m_surfaceType)
                    {
                        // Find the half of all the xp and yp on zlevel ....
                        case SolverUtils::eSphere:
                        case SolverUtils::eTRSphere:
                        {
                            F = xp[j] * xp[j] + y0 * y0 + zlevel * zlevel - Rhs;
                            dF = 2.0 * y0;
                        }
                        break;

                        case SolverUtils::eIrregular:
                        {
                            F = 2.0 * xp[j] * xp[j] + y0 * y0 * y0 * y0 +
                                y0 * y0 + zlevel * zlevel * zlevel * zlevel +
                                zlevel * zlevel - Rhs;
                            dF = 4.0 * y0 * y0 * y0 + 2.0 * y0;
                        }
                        break;

                        default:
                            break;
                    }

                    newy = y0 - F / dF;

                    if (fabs(F / dF) < Tol)
                    {
                        yp[j] = newy;
                        break;
                    }

                    else
                    {
                        y0 = newy;
                    }

                    ASSERTL0(i < Maxiter,
                             "Advection Velocity convergence fails");

                } // i-loop
            }
        }
        break;

        case SolverUtils::eNonconvex:
        {
            for (int j = 0; j < N + 1; ++j)
            {
                xp[j] = j * 2.0 * intval / N - intval;
                tmp   = 0.5 * Rhs -
                      0.5 * (zlevel * zlevel * zlevel * zlevel +
                             zlevel * zlevel) -
                      (xp[j] * xp[j] * xp[j] * xp[j] - xp[j] * xp[j]);
                if (tmp < 0)
                {
                    tmp = -1.0 * tmp;
                }
                yp[j] = sqrt(tmp);
            } // j-loop
        }
        break;

        default:
            break;

    } // switch-loop

    NekDouble pi        = 3.14159265358979323846;
    NekDouble arclength = 0.0;
    for (int j = 0; j < N; ++j)
    {
        arclength = arclength +
                    sqrt((yp[j + 1] - yp[j]) * (yp[j + 1] - yp[j]) +
                         (xp[j + 1] - xp[j]) * (xp[j + 1] - xp[j])) /
                        pi;
    }

    return arclength;
}

void MMFAdvection::v_SetInitialConditions(const NekDouble initialtime,
                                          bool dumpInitialConditions,
                                          const int domain)
{
    boost::ignore_unused(domain);

    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> u(nq);

    switch (m_TestType)
    {
        case eAdvectionBell:
        {
            AdvectionBellSphere(u);
            m_fields[0]->SetPhys(u);

            m_Mass0 = m_fields[0]->PhysIntegral(u);

            // forward transform to fill the modal coeffs
            for (int i = 0; i < m_fields.size(); ++i)
            {
                m_fields[i]->SetPhysState(true);
                m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),
                                      m_fields[i]->UpdateCoeffs());
            }
        }
        break;

        case eTestPlane:
        {
            Test2Dproblem(initialtime, u);
            m_fields[0]->SetPhys(u);

            m_Mass0 = m_fields[0]->PhysIntegral(u);

            // forward transform to fill the modal coeffs
            for (int i = 0; i < m_fields.size(); ++i)
            {
                m_fields[i]->SetPhysState(true);
                m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),
                                      m_fields[i]->UpdateCoeffs());
            }
        }
        break;

        case eTestPlaneMassConsv:
        {
            AdvectionBellPlane(u);
            m_fields[0]->SetPhys(u);

            m_Mass0 = m_fields[0]->PhysIntegral(u);
            std::cout << "m_Mass0 = " << m_Mass0 << std::endl;

            // forward transform to fill the modal coeffs
            for (int i = 0; i < m_fields.size(); ++i)
            {
                m_fields[i]->SetPhysState(true);
                m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),
                                      m_fields[i]->UpdateCoeffs());
            }
        }
        break;

        case eTestCube:
        {
            Test3Dproblem(initialtime, u);
            m_fields[0]->SetPhys(u);

            // forward transform to fill the modal coeffs
            for (int i = 0; i < m_fields.size(); ++i)
            {
                m_fields[i]->SetPhysState(true);
                m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),
                                      m_fields[i]->UpdateCoeffs());
            }
        }
        break;

        default:
            break;
    }

    if (dumpInitialConditions)
    {
        // dump initial conditions to file
        std::string outname = m_sessionName + "_initial.chk";
        WriteFld(outname);
    }
}

void MMFAdvection::AdvectionBellPlane(Array<OneD, NekDouble> &outfield)
{
    int nq = m_fields[0]->GetNpoints();

    NekDouble dist, m_radius_limit;

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    // Sets of parameters
    m_radius_limit = 0.5;

    NekDouble x0j, x1j;
    outfield = Array<OneD, NekDouble>(nq, 0.0);
    for (int j = 0; j < nq; ++j)
    {
        x0j = x[j];
        x1j = y[j];

        dist = sqrt(x0j * x0j + x1j * x1j);

        if (dist < m_radius_limit)
        {
            outfield[j] = 0.5 * (1.0 + cos(m_pi * dist / m_radius_limit));
        }
        else
        {
            outfield[j] = 0.0;
        }
    }
}

void MMFAdvection::AdvectionBellSphere(Array<OneD, NekDouble> &outfield)
{
    int nq = m_fields[0]->GetNpoints();

    NekDouble dist, radius, cosdiff, sin_theta, cos_theta, sin_varphi,
        cos_varphi;
    NekDouble m_theta_c, m_varphi_c, m_radius_limit, m_c0;

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    // Sets of parameters
    m_theta_c      = 0.0;
    m_varphi_c     = 3.0 * m_pi / 2.0;
    m_radius_limit = 7.0 * m_pi / 64.0;
    m_c0           = 0.0;

    NekDouble x0j, x1j, x2j;
    outfield = Array<OneD, NekDouble>(nq, 0.0);
    for (int j = 0; j < nq; ++j)
    {
        x0j = x[j];
        x1j = y[j];
        x2j = z[j];

        radius = sqrt(x0j * x0j + x1j * x1j + x2j * x2j);

        sin_varphi = x1j / sqrt(x0j * x0j + x1j * x1j);
        cos_varphi = x0j / sqrt(x0j * x0j + x1j * x1j);

        sin_theta = x2j / radius;
        cos_theta = sqrt(x0j * x0j + x1j * x1j) / radius;

        cosdiff = cos_varphi * cos(m_varphi_c) + sin_varphi * sin(m_varphi_c);
        dist    = radius * acos(sin(m_theta_c) * sin_theta +
                             cos(m_theta_c) * cos_theta * cosdiff);

        if (dist < m_radius_limit)
        {
            outfield[j] =
                0.5 * (1.0 + cos(m_pi * dist / m_radius_limit)) + m_c0;
        }
        else
        {
            outfield[j] = m_c0;
        }
    }
}

void MMFAdvection::Test2Dproblem(const NekDouble time,
                                 Array<OneD, NekDouble> &outfield)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    Array<OneD, NekDouble> u(nq);

    for (int i = 0; i < nq; ++i)
    {
        u[i] = cos(m_pi * (x0[i] - m_advx * time)) *
               cos(m_pi * (x1[i] - m_advy * time));
    }

    outfield = u;
}

void MMFAdvection::Test3Dproblem(const NekDouble time,
                                 Array<OneD, NekDouble> &outfield)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    Array<OneD, NekDouble> u(nq);

    for (int i = 0; i < nq; ++i)
    {
        u[i] = cos(m_pi * (x0[i] - m_advx * time)) *
               cos(m_pi * (x1[i] - m_advy * time)) *
               cos(m_pi * (x2[i] - m_advz * time));
    }

    outfield = u;
}

void MMFAdvection::ComputeNablaCdotVelocity(Array<OneD, NekDouble> &vellc)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> velcoeff(nq, 0.0);

    Array<OneD, NekDouble> Dtmp0(nq);
    Array<OneD, NekDouble> Dtmp1(nq);
    Array<OneD, NekDouble> Dtmp2(nq);
    Array<OneD, NekDouble> Drv(nq);

    vellc = Array<OneD, NekDouble>(nq, 0.0);

    // m_vellc = \nabla m_vel \cdot tan_i
    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> vessel(nq);

    for (int j = 0; j < m_shapedim; ++j)
    {
        Vmath::Zero(nq, velcoeff, 1);
        for (int k = 0; k < m_spacedim; ++k)
        {
            // a_j = tan_j cdot m_vel
            Vmath::Vvtvp(nq, &m_movingframes[j][k * nq], 1, &m_velocity[k][0],
                         1, &velcoeff[0], 1, &velcoeff[0], 1);
        }

        // d a_j / d x^k
        m_fields[0]->PhysDeriv(velcoeff, Dtmp0, Dtmp1, Dtmp2);

        for (int k = 0; k < m_spacedim; ++k)
        {
            // tan_j^k ( d a_j / d x^k )
            switch (k)
            {
                case (0):
                {
                    Vmath::Vvtvp(nq, &Dtmp0[0], 1, &m_movingframes[j][k * nq],
                                 1, &vellc[0], 1, &vellc[0], 1);
                }
                break;

                case (1):
                {
                    Vmath::Vvtvp(nq, &Dtmp1[0], 1, &m_movingframes[j][k * nq],
                                 1, &vellc[0], 1, &vellc[0], 1);
                }
                break;

                case (2):
                {
                    Vmath::Vvtvp(nq, &Dtmp2[0], 1, &m_movingframes[j][k * nq],
                                 1, &vellc[0], 1, &vellc[0], 1);
                }
                break;
            }
        }
    }
}

void MMFAdvection::ComputeveldotMF(
    Array<OneD, Array<OneD, NekDouble>> &veldotMF)
{
    int nq = m_fields[0]->GetNpoints();

    veldotMF = Array<OneD, Array<OneD, NekDouble>>(m_shapedim);

    Array<OneD, NekDouble> magMF(nq);
    for (int j = 0; j < m_shapedim; ++j)
    {
        veldotMF[j] = Array<OneD, NekDouble>(nq, 0.0);

        for (int k = 0; k < m_spacedim; ++k)
        {
            Vmath::Vvtvp(nq, &m_movingframes[j][k * nq], 1, &m_velocity[k][0],
                         1, &veldotMF[j][0], 1, &veldotMF[j][0], 1);
        }
    }
}

void MMFAdvection::v_EvaluateExactSolution(unsigned int field,
                                           Array<OneD, NekDouble> &outfield,
                                           const NekDouble time)
{
    boost::ignore_unused(field);

    switch (m_TestType)
    {
        case eAdvectionBell:
        {
            AdvectionBellSphere(outfield);
        }
        break;

        case eTestPlane:
        {
            Test2Dproblem(time, outfield);
        }
        break;

        case eTestPlaneMassConsv:
        {
            AdvectionBellPlane(outfield);
        }
        break;

        case eTestCube:
        {
            Test3Dproblem(time, outfield);
        }
        break;

        default:
            break;
    }
}

void MMFAdvection::v_GenerateSummary(SolverUtils::SummaryList &s)
{
    MMFSystem::v_GenerateSummary(s);
    SolverUtils::AddSummaryItem(s, "TestType", TestTypeMap[m_TestType]);
    SolverUtils::AddSummaryItem(s, "Rotation Angle", m_RotAngle);
}
}
