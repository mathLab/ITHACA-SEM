///////////////////////////////////////////////////////////////////////////////
//
// File VCSMapping.cpp
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
// Description: Velocity Correction Scheme w/ coordinate transformation
// for the Incompressible Navier Stokes equations
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/EquationSystems/VCSMapping.h>
#include <LibUtilities/BasicUtils/Timer.h>
#include <SolverUtils/Core/Misc.h>

#include <boost/algorithm/string.hpp>

using namespace std;

namespace Nektar
{
    string VCSMapping::className =
        SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
            "VCSMapping",
            VCSMapping::create);

    /**
     * Constructor. Creates ...
     *
     * \param
     * \param
     */
    VCSMapping::VCSMapping(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
        : UnsteadySystem(pSession, pGraph),
          VelocityCorrectionScheme(pSession, pGraph)
    {

    }

    void VCSMapping::v_InitObject()
    {
        VelocityCorrectionScheme::v_InitObject();

        m_mapping = GlobalMapping::Mapping::Load(m_session, m_fields);
        ASSERTL0(m_mapping,
             "Could not create mapping in VCSMapping.");

        std::string vExtrapolation = "Mapping";
        m_extrapolation = GetExtrapolateFactory().CreateInstance(
            vExtrapolation,
            m_session,
            m_fields,
            m_pressure,
            m_velocity,
            m_advObject);
        m_extrapolation->SubSteppingTimeIntegration(m_intScheme);
        m_extrapolation->GenerateHOPBCMap(m_session);

       // Storage to extrapolate pressure forcing
        int physTot = m_fields[0]->GetTotPoints();
        int intSteps = 1;

        if ( m_intScheme->GetName() == "IMEX" ||
             m_intScheme->GetName() == "IMEXGear" )
        {
            m_intSteps = m_intScheme->GetOrder();
        }
        else
        {
            NEKERROR(ErrorUtil::efatal, "Integration method not suitable: "
                     "Options include IMEXGear or IMEXOrder{1,2,3,4}");
        }

        m_presForcingCorrection= Array<OneD, Array<OneD, NekDouble> >(intSteps);
        for(int i = 0; i < m_presForcingCorrection.size(); i++)
        {
            m_presForcingCorrection[i] = Array<OneD, NekDouble>(physTot,0.0);
        }
        m_verbose = (m_session->DefinesCmdLineArgument("verbose"))? true :false;

        // Load solve parameters related to the mapping
        // Flags determining if pressure/viscous terms should be treated implicitly
        m_session->MatchSolverInfo("MappingImplicitPressure","True",
                                                    m_implicitPressure,false);
        m_session->MatchSolverInfo("MappingImplicitViscous","True",
                                                    m_implicitViscous,false);
        m_session->MatchSolverInfo("MappingNeglectViscous","True",
                                                    m_neglectViscous,false);

        if (m_neglectViscous)
        {
            m_implicitViscous = false;
        }

        // Tolerances and relaxation parameters for implicit terms
        m_session->LoadParameter("MappingPressureTolerance",
                                            m_pressureTolerance,1e-12);
        m_session->LoadParameter("MappingViscousTolerance",
                                            m_viscousTolerance,1e-12);
        m_session->LoadParameter("MappingPressureRelaxation",
                                            m_pressureRelaxation,1.0);
        m_session->LoadParameter("MappingViscousRelaxation",
                                            m_viscousRelaxation,1.0);

    }

    /**
     * Destructor
     */
    VCSMapping::~VCSMapping(void)
    {
    }

    void VCSMapping::v_DoInitialise(void)
    {
        UnsteadySystem::v_DoInitialise();

        // Set up Field Meta Data for output files
        m_fieldMetaDataMap["Kinvis"]   =
                                boost::lexical_cast<std::string>(m_kinvis);
        m_fieldMetaDataMap["TimeStep"] =
                                boost::lexical_cast<std::string>(m_timestep);

        // Correct Dirichlet boundary conditions to account for mapping
        m_mapping->UpdateBCs(0.0);
        //
        m_F = Array<OneD, Array< OneD, NekDouble> > (m_nConvectiveFields);
        for(int i = 0; i < m_nConvectiveFields; ++i)
        {
            m_fields[i]->ImposeDirichletConditions(m_fields[i]->UpdateCoeffs());
            m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),
                                  m_fields[i]->UpdatePhys());
            m_F[i] = Array< OneD, NekDouble> (m_fields[0]->GetTotPoints(), 0.0);
        }

        // Initialise m_gradP
        int physTot = m_fields[0]->GetTotPoints();
        m_gradP = Array<OneD, Array<OneD, NekDouble> >(m_nConvectiveFields);
        for(int i = 0; i < m_nConvectiveFields; ++i)
        {
            m_gradP[i] = Array<OneD, NekDouble>(physTot,0.0);
            m_pressure->PhysDeriv(MultiRegions::DirCartesianMap[i],
                                    m_pressure->GetPhys(),
                                    m_gradP[i]);
            if(m_pressure->GetWaveSpace())
            {
                m_pressure->HomogeneousBwdTrans(m_gradP[i],m_gradP[i]);
            }
        }
    }

    /**
     * Explicit part of the method - Advection, Forcing + HOPBCs
     */
    void VCSMapping::v_EvaluateAdvection_SetPressureBCs(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
        Array<OneD, Array<OneD, NekDouble> > &outarray,
        const NekDouble time)
    {
        EvaluateAdvectionTerms(inarray, outarray);

        // Smooth advection
        if(m_SmoothAdvection)
        {
            for(int i = 0; i < m_nConvectiveFields; ++i)
            {
                m_pressure->SmoothField(outarray[i]);
            }
        }

        // Add forcing terms
        for (auto &x : m_forcing)
        {
            x->Apply(m_fields, inarray, outarray, time);
        }

        // Add mapping terms
        ApplyIncNSMappingForcing( inarray, outarray);

        // Calculate High-Order pressure boundary conditions
        m_extrapolation->EvaluatePressureBCs(inarray,outarray,m_kinvis);

        // Update mapping and deal with Dirichlet boundary conditions
        if (m_mapping->IsTimeDependent())
        {
            if (m_mapping->IsFromFunction())
            {
                // If the transformation is explicitly defined, update it here
                // Otherwise, it will be done somewhere else (ForcingMovingBody)
                m_mapping->UpdateMapping(time+m_timestep);
            }
            m_mapping->UpdateBCs(time+m_timestep);
        }
    }


    /**
     * Forcing term for Poisson solver solver
     */
    void   VCSMapping::v_SetUpPressureForcing(
        const Array<OneD, const Array<OneD, NekDouble> > &fields,
        Array<OneD, Array<OneD, NekDouble> > &Forcing,
        const NekDouble aii_Dt)
    {
        if (m_mapping->HasConstantJacobian())
        {
            VelocityCorrectionScheme::v_SetUpPressureForcing(fields,
                                                            Forcing, aii_Dt);
        }
        else
        {
            int physTot = m_fields[0]->GetTotPoints();
            int nvel = m_nConvectiveFields;
            Array<OneD, NekDouble> wk(physTot, 0.0);

            Array<OneD, NekDouble> Jac(physTot,0.0);
            m_mapping->GetJacobian(Jac);

            // Calculate div(J*u/Dt)
            Vmath::Zero(physTot,Forcing[0],1);
            for(int i = 0; i < nvel; ++i)
            {
                if (m_fields[i]->GetWaveSpace())
                {
                    m_fields[i]->HomogeneousBwdTrans(fields[i],wk);
                }
                else
                {
                    Vmath::Vcopy(physTot, fields[i], 1, wk, 1);
                }
                Vmath::Vmul(physTot,wk,1,Jac,1,wk,1);
                if (m_fields[i]->GetWaveSpace())
                {
                    m_fields[i]->HomogeneousFwdTrans(wk,wk);
                }
                m_fields[i]->PhysDeriv(MultiRegions::DirCartesianMap[i],wk, wk);
                Vmath::Vadd(physTot,wk,1,Forcing[0],1,Forcing[0],1);
            }
            Vmath::Smul(physTot,1.0/aii_Dt,Forcing[0],1,Forcing[0],1);

            //
            // If the mapping viscous terms are being treated explicitly
            //        we need to apply a correction to the forcing
            if (!m_implicitViscous)
            {
                bool wavespace = m_fields[0]->GetWaveSpace();
                m_fields[0]->SetWaveSpace(false);

                //
                //  Part 1: div(J*grad(U/J . grad(J)))
                Array<OneD, Array<OneD, NekDouble> > tmp (nvel);
                Array<OneD, Array<OneD, NekDouble> > velocity (nvel);
                for(int i = 0; i < tmp.size(); i++)
                {
                    tmp[i] = Array<OneD, NekDouble>(physTot,0.0);
                    velocity[i] = Array<OneD, NekDouble>(physTot,0.0);
                    if (wavespace)
                    {
                        m_fields[0]->HomogeneousBwdTrans(m_fields[i]->GetPhys(),
                                                         velocity[i]);
                    }
                    else
                    {
                        Vmath::Vcopy(physTot, m_fields[i]->GetPhys(), 1,
                                                velocity[i], 1);
                    }
                }
                // Calculate wk = U.grad(J)
                m_mapping->DotGradJacobian(velocity, wk);
                // Calculate wk = (U.grad(J))/J
                Vmath::Vdiv(physTot, wk, 1, Jac, 1, wk, 1);
                // J*grad[(U.grad(J))/J]
                for(int i = 0; i < nvel; ++i)
                {
                    m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[i],
                                wk, tmp[i]);
                    Vmath::Vmul(physTot, Jac, 1, tmp[i], 1, tmp[i], 1);
                }
                // div(J*grad[(U.grad(J))/J])
                Vmath::Zero(physTot, wk, 1);
                for(int i = 0; i < nvel; ++i)
                {
                    m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[i],
                                tmp[i], tmp[i]);
                    Vmath::Vadd(physTot, wk, 1, tmp[i], 1, wk, 1);
                }

                // Part 2: grad(J) . curl(curl(U))
                m_mapping->CurlCurlField(velocity, tmp, m_implicitViscous);
                // dont need velocity any more, so reuse it
                m_mapping->DotGradJacobian(tmp, velocity[0]);

                // Add two parts
                Vmath::Vadd(physTot, velocity[0], 1, wk, 1, wk, 1);

                // Multiply by kinvis and prepare to extrapolate
                int nlevels = m_presForcingCorrection.size();
                Vmath::Smul(physTot, m_kinvis, wk, 1,
                                     m_presForcingCorrection[nlevels-1], 1);

                // Extrapolate correction
                m_extrapolation->ExtrapolateArray(m_presForcingCorrection);

                // Put in wavespace
                if (wavespace)
                {
                    m_fields[0]->HomogeneousFwdTrans(
                                   m_presForcingCorrection[nlevels-1],wk);
                }
                else
                {
                    Vmath::Vcopy(physTot, m_presForcingCorrection[nlevels-1], 1,
                                          wk, 1);
                }
                // Apply correction: Forcing = Forcing - correction
                Vmath::Vsub(physTot, Forcing[0], 1, wk, 1, Forcing[0], 1);

                m_fields[0]->SetWaveSpace(wavespace);
            }
        }
    }

    /**
     * Forcing term for Helmholtz solver
     */
    void   VCSMapping::v_SetUpViscousForcing(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
        Array<OneD, Array<OneD, NekDouble> > &Forcing,
        const NekDouble aii_Dt)
    {
        NekDouble aii_dtinv = 1.0/aii_Dt;
        int physTot = m_fields[0]->GetTotPoints();

        // Grad p
        m_pressure->BwdTrans(m_pressure->GetCoeffs(),m_pressure->UpdatePhys());

        int nvel = m_velocity.size();
        if(nvel == 2)
        {
            m_pressure->PhysDeriv(m_pressure->GetPhys(), Forcing[0], Forcing[1]);
        }
        else
        {
            m_pressure->PhysDeriv(m_pressure->GetPhys(), Forcing[0], Forcing[1],
                                  Forcing[2]);
        }

        // Copy grad p in physical space to m_gradP to reuse later
        if (m_pressure->GetWaveSpace())
        {
            for (int i=0; i<nvel; i++)
            {
                m_pressure->HomogeneousBwdTrans(Forcing[i],m_gradP[i]);
            }
        }
        else
        {
            for (int i=0; i<nvel; i++)
            {
                Vmath::Vcopy(physTot, Forcing[i], 1, m_gradP[i], 1);
            }
        }

        if ( (!m_mapping->HasConstantJacobian()) || m_implicitPressure)
        {
            // If pressure terms are treated explicitly, we need to divide by J
            //    if they are implicit, we need to calculate G(p)
            if (m_implicitPressure)
            {
                m_mapping->RaiseIndex(m_gradP, Forcing);
            }
            else
            {
                Array<OneD, NekDouble> Jac(physTot,0.0);
                m_mapping->GetJacobian(Jac);
                for (int i=0; i<nvel; i++)
                {
                    Vmath::Vdiv(physTot, m_gradP[i], 1, Jac, 1, Forcing[i], 1);
                }
            }
            // Transform back to wavespace
            if (m_pressure->GetWaveSpace())
            {
                for (int i=0; i<nvel; i++)
                {
                    m_pressure->HomogeneousFwdTrans(Forcing[i],Forcing[i]);
                }
            }
        }

        // Subtract inarray/(aii_dt) and divide by kinvis. Kinvis will
        // need to be updated for the convected fields.
        for(int i = 0; i < nvel; ++i)
        {
            Blas::Daxpy(physTot,-aii_dtinv,inarray[i],1,Forcing[i],1);
            Blas::Dscal(physTot,1.0/m_kinvis,&(Forcing[i])[0],1);
        }
    }

    /**
     * Solve pressure system
     */
    void   VCSMapping::v_SolvePressure(
        const Array<OneD, NekDouble>  &Forcing)
    {
        if (!m_implicitPressure)
        {
            VelocityCorrectionScheme::v_SolvePressure(Forcing);
        }
        else
        {
            int physTot = m_fields[0]->GetTotPoints();
            int nvel = m_nConvectiveFields;
            bool converged = false;          // flag to mark if system converged
            int s = 0;                       // iteration counter
            NekDouble error;                 // L2 error at current iteration
            NekDouble forcing_L2 = 0.0;      // L2 norm of F

            int maxIter;
            m_session->LoadParameter("MappingMaxIter",maxIter,5000);

            // rhs of the equation at current iteration
            Array< OneD, NekDouble> F_corrected(physTot, 0.0);
            // Pressure field at previous iteration
            Array<OneD, NekDouble>  previous_iter (physTot, 0.0);
            // Temporary variables
            Array<OneD, Array<OneD, NekDouble> > wk1(nvel);
            Array<OneD, Array<OneD, NekDouble> > wk2(nvel);
            Array<OneD, Array<OneD, NekDouble> > gradP(nvel);
            for(int i = 0; i < nvel; ++i)
            {
                wk1[i]   = Array<OneD, NekDouble> (physTot, 0.0);
                wk2[i]   = Array<OneD, NekDouble> (physTot, 0.0);
                gradP[i] = Array<OneD, NekDouble> (physTot, 0.0);
            }

            // Jacobian
            Array<OneD, NekDouble> Jac(physTot, 0.0);
            m_mapping->GetJacobian(Jac);

            // Factors for Laplacian system
            StdRegions::ConstFactorMap factors;
            factors[StdRegions::eFactorLambda] = 0.0;

            m_pressure->BwdTrans(m_pressure->GetCoeffs(),
                                    m_pressure->UpdatePhys());
            forcing_L2 = m_pressure->L2(Forcing, wk1[0]);
            while (!converged)
            {
                // Update iteration counter and set previous iteration field
                // (use previous timestep solution for first iteration)
                s++;
                ASSERTL0(s < maxIter,
                         "VCSMapping exceeded maximum number of iterations.");

                Vmath::Vcopy(physTot, m_pressure->GetPhys(), 1,
                                        previous_iter, 1);

                // Correct pressure bc to account for iteration
                m_extrapolation->CorrectPressureBCs(previous_iter);

                //
                // Calculate forcing term for this iteration
                //
                for(int i = 0; i < nvel; ++i)
                {
                    m_pressure->PhysDeriv(MultiRegions::DirCartesianMap[i],
                                                    previous_iter, gradP[i]);
                    if(m_pressure->GetWaveSpace())
                    {
                        m_pressure->HomogeneousBwdTrans(gradP[i], wk1[i]);
                    }
                    else
                    {
                        Vmath::Vcopy(physTot, gradP[i], 1, wk1[i], 1);
                    }
                }
                m_mapping->RaiseIndex(wk1, wk2);   // G(p)

                m_mapping->Divergence(wk2, F_corrected);   // div(G(p))
                if (!m_mapping->HasConstantJacobian())
                {
                    Vmath::Vmul(physTot, F_corrected, 1,
                                         Jac, 1,
                                         F_corrected, 1);
                }
                // alpha*J*div(G(p))
                Vmath::Smul(physTot, m_pressureRelaxation, F_corrected, 1,
                                                           F_corrected, 1);
                if(m_pressure->GetWaveSpace())
                {
                    m_pressure->HomogeneousFwdTrans(F_corrected, F_corrected);
                }
                // alpha*J*div(G(p)) - p_ii
                for (int i = 0; i < m_nConvectiveFields; ++i)
                {
                    m_pressure->PhysDeriv(MultiRegions::DirCartesianMap[i],
                                                            gradP[i], wk1[0]);
                    Vmath::Vsub(physTot, F_corrected, 1, wk1[0], 1,
                                                            F_corrected, 1);
                }
                // p_i,i - J*div(G(p))
                Vmath::Neg(physTot, F_corrected, 1);
                // alpha*F -  alpha*J*div(G(p)) + p_i,i
                Vmath::Smul(physTot, m_pressureRelaxation, Forcing, 1,
                                                            wk1[0], 1);
                Vmath::Vadd(physTot, wk1[0], 1, F_corrected, 1, F_corrected, 1);

                //
                // Solve system
                //
                m_pressure->HelmSolve(F_corrected, m_pressure->UpdateCoeffs(),
                                      factors);
                m_pressure->BwdTrans(m_pressure->GetCoeffs(),
                                    m_pressure->UpdatePhys());

                //
                // Test convergence
                //
                error = m_pressure->L2(m_pressure->GetPhys(), previous_iter);
                if ( forcing_L2 != 0)
                {
                    if ( (error/forcing_L2 < m_pressureTolerance))
                    {
                        converged = true;
                    }
                }
                else
                {
                    if ( error < m_pressureTolerance)
                    {
                        converged = true;
                    }
                }
            }
            if (m_verbose && m_session->GetComm()->GetRank()==0)
            {
                std::cout << " Pressure system (mapping) converged in " << s <<
                            " iterations with error = " << error << std::endl;
            }
        }
    }

    /**
     * Solve velocity system
     */
    void   VCSMapping::v_SolveViscous(
        const Array<OneD, const Array<OneD, NekDouble> > &Forcing,
        Array<OneD, Array<OneD, NekDouble> > &outarray,
        const NekDouble aii_Dt)
    {
        if(!m_implicitViscous)
        {
            VelocityCorrectionScheme::v_SolveViscous(Forcing, outarray, aii_Dt);
        }
        else
        {
            int physTot = m_fields[0]->GetTotPoints();
            int nvel = m_nConvectiveFields;
            bool converged = false;          // flag to mark if system converged
            int s = 0;                       // iteration counter
            NekDouble error, max_error;      // L2 error at current iteration

            int maxIter;
            m_session->LoadParameter("MappingMaxIter",maxIter,5000);

            //L2 norm of F
            Array<OneD, NekDouble> forcing_L2(m_nConvectiveFields,0.0);

            // rhs of the equation at current iteration
            Array<OneD, Array<OneD, NekDouble> > F_corrected(nvel);
            // Solution at previous iteration
            Array<OneD, Array<OneD, NekDouble> > previous_iter(nvel);
            // Working space
            Array<OneD, Array<OneD, NekDouble> >  wk(nvel);
            for(int i = 0; i < nvel; ++i)
            {
                F_corrected[i]   = Array<OneD, NekDouble> (physTot, 0.0);
                previous_iter[i] = Array<OneD, NekDouble> (physTot, 0.0);
                wk[i]            = Array<OneD, NekDouble> (physTot, 0.0);
            }

            // Factors for Helmholtz system
            StdRegions::ConstFactorMap factors;
            factors[StdRegions::eFactorLambda] =
                                    1.0*m_viscousRelaxation/aii_Dt/m_kinvis;
            if(m_useSpecVanVisc)
            {
                factors[StdRegions::eFactorSVVCutoffRatio] = m_sVVCutoffRatio;
                factors[StdRegions::eFactorSVVDiffCoeff]   =
                                                      m_sVVDiffCoeff/m_kinvis;
            }

            // Calculate L2-norm of F and set initial solution for iteration
            for(int i = 0; i < nvel; ++i)
            {
                forcing_L2[i] = m_fields[0]->L2(Forcing[i],wk[0]);
                m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),
                                        previous_iter[i]);
            }

            while (!converged)
            {
                converged = true;
                // Iteration counter
                s++;
                ASSERTL0(s < maxIter,
                         "VCSMapping exceeded maximum number of iterations.");

                max_error = 0.0;

                //
                // Calculate forcing term for next iteration
                //

                // Calculate L(U)- in this parts all components might be coupled
                if(m_fields[0]->GetWaveSpace())
                {
                    for (int i = 0; i < nvel; ++i)
                    {
                        m_fields[0]->HomogeneousBwdTrans(previous_iter[i],
                                                                    wk[i]);
                    }
                }
                else
                {
                    for (int i = 0; i < nvel; ++i)
                    {
                        Vmath::Vcopy(physTot, previous_iter[i], 1, wk[i], 1);
                    }
                }

                // (L(U^i) - 1/alpha*U^i_jj)
                m_mapping->VelocityLaplacian(wk, F_corrected,
                                                1.0/m_viscousRelaxation);

                if(m_fields[0]->GetWaveSpace())
                {
                    for (int i = 0; i < nvel; ++i)
                    {
                        m_fields[0]->HomogeneousFwdTrans(F_corrected[i],
                                                            F_corrected[i]);
                    }
                }
                else
                {
                    for (int i = 0; i < nvel; ++i)
                    {
                        Vmath::Vcopy(physTot, F_corrected[i], 1,
                                                F_corrected[i], 1);
                    }
                }

                // Loop velocity components
                for (int i = 0; i < nvel; ++i)
                {
                    // (-alpha*L(U^i) + U^i_jj)
                    Vmath::Smul(physTot, -1.0*m_viscousRelaxation,
                                                    F_corrected[i], 1,
                                                    F_corrected[i], 1);
                    //  F_corrected = alpha*F + (-alpha*L(U^i) + U^i_jj)
                    Vmath::Smul(physTot, m_viscousRelaxation, Forcing[i], 1,
                                                                    wk[0], 1);
                    Vmath::Vadd(physTot, wk[0], 1, F_corrected[i], 1,
                                                    F_corrected[i], 1);

                    //
                    // Solve System
                    //
                    m_fields[i]->HelmSolve(F_corrected[i],
                                    m_fields[i]->UpdateCoeffs(),factors);
                    m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),outarray[i]);

                    //
                    // Test convergence
                    //
                    error = m_fields[i]->L2(outarray[i], previous_iter[i]);

                    if ( forcing_L2[i] != 0)
                    {
                        if ( (error/forcing_L2[i] >= m_viscousTolerance))
                        {
                            converged = false;
                        }
                    }
                    else
                    {
                        if ( error >= m_viscousTolerance)
                        {
                            converged = false;
                        }
                    }
                    if (error > max_error)
                    {
                        max_error = error;
                    }

                    // Copy field to previous_iter
                    Vmath::Vcopy(physTot, outarray[i], 1, previous_iter[i], 1);
                }
            }
            if (m_verbose && m_session->GetComm()->GetRank()==0)
            {
                std::cout << " Velocity system (mapping) converged in " << s <<
                        " iterations with error = " << max_error << std::endl;
            }
        }
    }

    /**
     * Explicit terms of the mapping
     */
    void   VCSMapping::ApplyIncNSMappingForcing(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
        Array<OneD, Array<OneD, NekDouble> >             &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        Array<OneD, Array<OneD, NekDouble> >       vel(m_nConvectiveFields);
        Array<OneD, Array<OneD, NekDouble> >       velPhys(m_nConvectiveFields);
        Array<OneD, Array<OneD, NekDouble> >       Forcing(m_nConvectiveFields);
        Array<OneD, Array<OneD, NekDouble> >       tmp(m_nConvectiveFields);
        for (int i = 0; i < m_nConvectiveFields; ++i)
        {
            velPhys[i] = Array<OneD, NekDouble> (physTot, 0.0);
            Forcing[i] = Array<OneD, NekDouble> (physTot, 0.0);
            tmp[i] = Array<OneD, NekDouble> (physTot, 0.0);
        }

        // Get fields and store velocity in wavespace and physical space
        if(m_fields[0]->GetWaveSpace())
        {
            for (int i = 0; i < m_nConvectiveFields; ++i)
            {
                vel[i] = inarray[i];
                m_fields[0]->HomogeneousBwdTrans(vel[i],velPhys[i]);
            }
        }
        else
        {
            for (int i = 0; i < m_nConvectiveFields; ++i)
            {
                vel[i] = inarray[i];
                Vmath::Vcopy(physTot, inarray[i], 1, velPhys[i], 1);
            }
        }

        //Advection contribution
        MappingAdvectionCorrection(velPhys, Forcing);

        // Time-derivative contribution
        if ( m_mapping->IsTimeDependent() )
        {
            MappingAccelerationCorrection(vel, velPhys, tmp);
            for (int i = 0; i < m_nConvectiveFields; ++i)
            {
                Vmath::Vadd(physTot, tmp[i], 1, Forcing[i], 1, Forcing[i], 1);
            }
        }

        // Pressure contribution
        if (!m_implicitPressure)
        {
            MappingPressureCorrection(tmp);
            for (int i = 0; i < m_nConvectiveFields; ++i)
            {
                Vmath::Vadd(physTot, tmp[i], 1, Forcing[i], 1, Forcing[i], 1);
            }
        }
        // Viscous contribution
        if ( (!m_implicitViscous) && (!m_neglectViscous))
        {
            MappingViscousCorrection(velPhys, tmp);
            for (int i = 0; i < m_nConvectiveFields; ++i)
            {
                Vmath::Smul(physTot, m_kinvis, tmp[i], 1, tmp[i], 1);
                Vmath::Vadd(physTot, tmp[i], 1, Forcing[i], 1, Forcing[i], 1);
            }
        }

        // If necessary, transform to wavespace
        if(m_fields[0]->GetWaveSpace())
        {
            for (int i = 0; i < m_nConvectiveFields; ++i)
            {
                m_fields[0]->HomogeneousFwdTrans(Forcing[i],Forcing[i]);
            }
        }

        // Add to outarray
        for (int i = 0; i < m_nConvectiveFields; ++i)
        {
            Vmath::Vadd(physTot, outarray[i], 1, Forcing[i], 1, outarray[i], 1);
        }
    }

        void VCSMapping::MappingAdvectionCorrection(
            const Array<OneD, const Array<OneD, NekDouble> >  &velPhys,
            Array<OneD, Array<OneD, NekDouble> >              &outarray)
        {
            int physTot = m_fields[0]->GetTotPoints();
            int nvel = m_nConvectiveFields;

            Array<OneD, Array<OneD, NekDouble> > wk(nvel*nvel);

            // Apply Christoffel symbols to obtain {i,kj}vel(k)
            m_mapping->ApplyChristoffelContravar(velPhys, wk);

            // Calculate correction -U^j*{i,kj}vel(k)
            for (int i = 0; i< nvel; i++)
            {
                Vmath::Zero(physTot,outarray[i],1);
                for (int j = 0; j< nvel; j++)
                {
                        Vmath::Vvtvp(physTot,wk[i*nvel+j],1,velPhys[j],1,
                                        outarray[i],1,outarray[i],1);
                }
                Vmath::Neg(physTot, outarray[i], 1);
            }
        }

        void VCSMapping::MappingAccelerationCorrection(
            const Array<OneD, const Array<OneD, NekDouble> >  &vel,
            const Array<OneD, const Array<OneD, NekDouble> >  &velPhys,
            Array<OneD, Array<OneD, NekDouble> >              &outarray)
        {
            int physTot = m_fields[0]->GetTotPoints();
            int nvel = m_nConvectiveFields;

            Array<OneD, Array<OneD, NekDouble> > wk(nvel*nvel);
            Array<OneD, Array<OneD, NekDouble> > tmp(nvel);
            Array<OneD, Array<OneD, NekDouble> > coordVel(nvel);
            for (int i = 0; i< nvel; i++)
            {
                tmp[i] = Array<OneD, NekDouble> (physTot, 0.0);
                coordVel[i] = Array<OneD, NekDouble> (physTot, 0.0);
            }
            // Get coordinates velocity in transformed system
            m_mapping->GetCoordVelocity(tmp);
            m_mapping->ContravarFromCartesian(tmp, coordVel);

            // Calculate first term: U^j u^i,j = U^j (du^i/dx^j + {i,kj}u^k)
            m_mapping->ApplyChristoffelContravar(velPhys, wk);
            for (int i=0; i< nvel; i++)
            {
                Vmath::Zero(physTot,outarray[i],1);

                m_fields[0]->PhysDeriv(velPhys[i], tmp[0], tmp[1]);
                for (int j=0; j< nvel; j++)
                {
                    if (j == 2)
                    {
                        m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[j],
                                        vel[i], tmp[2]);
                        if (m_fields[0]->GetWaveSpace())
                        {
                            m_fields[0]->HomogeneousBwdTrans(tmp[2],tmp[2]);
                        }
                    }

                    Vmath::Vadd(physTot,wk[i*nvel+j],1,tmp[j],1,
                                                        wk[i*nvel+j], 1);

                    Vmath::Vvtvp(physTot, coordVel[j], 1, wk[i*nvel+j], 1,
                                          outarray[i], 1, outarray[i], 1);
                }
            }

            // Set wavespace to false and store current value
            bool wavespace = m_fields[0]->GetWaveSpace();
            m_fields[0]->SetWaveSpace(false);

            // Add -u^j U^i,j
            m_mapping->ApplyChristoffelContravar(coordVel, wk);
            for (int i=0; i< nvel; i++)
            {
                if(nvel == 2)
                {
                    m_fields[0]->PhysDeriv(coordVel[i], tmp[0], tmp[1]);
                }
                else
                {
                    m_fields[0]->PhysDeriv(coordVel[i], tmp[0], tmp[1], tmp[2]);
                }

                for (int j=0; j< nvel; j++)
                {
                    Vmath::Vadd(physTot,wk[i*nvel+j],1,tmp[j],1,
                                                        wk[i*nvel+j], 1);
                    Vmath::Neg(physTot, wk[i*nvel+j], 1);

                    Vmath::Vvtvp(physTot, velPhys[j], 1, wk[i*nvel+j], 1,
                                          outarray[i], 1, outarray[i], 1);
                }
            }

            // Restore value of wavespace
            m_fields[0]->SetWaveSpace(wavespace);
        }

        void VCSMapping::MappingPressureCorrection(
            Array<OneD, Array<OneD, NekDouble> >              &outarray)
        {
            int physTot = m_fields[0]->GetTotPoints();
            int nvel = m_nConvectiveFields;

            // Calculate g^(ij)p_(,j)
            m_mapping->RaiseIndex(m_gradP, outarray);

            // Calculate correction = (nabla p)/J - g^(ij)p_,j
            // (Jac is not required if it is constant)
            if ( !m_mapping->HasConstantJacobian())
            {
                Array<OneD, NekDouble> Jac(physTot, 0.0);
                m_mapping->GetJacobian(Jac);
                for(int i = 0; i < nvel; ++i)
                {
                    Vmath::Vdiv(physTot, m_gradP[i], 1, Jac, 1, m_gradP[i], 1);
                }
            }
            for(int i = 0; i < nvel; ++i)
            {
                Vmath::Vsub(physTot, m_gradP[i], 1,outarray[i], 1,
                                                    outarray[i],1);
            }
        }

        void VCSMapping::MappingViscousCorrection(
            const Array<OneD, const Array<OneD, NekDouble> >  &velPhys,
            Array<OneD, Array<OneD, NekDouble> >              &outarray)
        {
            // L(U) - 1.0*d^2(u^i)/dx^jdx^j
            m_mapping->VelocityLaplacian(velPhys, outarray, 1.0);
        }

} //end of namespace
