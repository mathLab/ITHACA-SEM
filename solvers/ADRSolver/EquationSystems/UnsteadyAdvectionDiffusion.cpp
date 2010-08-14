#include <iostream>

#include <ADRSolver/EquationSystems/UnsteadyAdvectionDiffusion.h>

namespace Nektar
{
    string UnsteadyAdvectionDiffusion::className = EquationSystemFactory::RegisterCreatorFunction("UnsteadyAdvectionDiffusion", UnsteadyAdvectionDiffusion::create);

    UnsteadyAdvectionDiffusion::UnsteadyAdvectionDiffusion(SessionReaderSharedPtr& pSession)
        : UnsteadySystem(pSession)
    {
        pSession->LoadParameter("wavefreq",   m_waveFreq, 0.0);
        pSession->LoadParameter("epsilon",    m_epsilon,  0.0);

        // Define Velocity fields
        m_velocity = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
        int nq = m_fields[0]->GetNpoints();
        std::string velStr[3] = {"Vx","Vy","Vz"};

        for(int i = 0; i < m_spacedim; ++i)
        {
            m_velocity[i] = Array<OneD, NekDouble> (nq,0.0);

            SpatialDomains::ConstUserDefinedEqnShPtr ifunc
                = m_boundaryConditions->GetUserDefinedEqn(velStr[i]);

            EvaluateFunction(m_velocity[i],ifunc);
        }

        if (m_explicitAdvection)
        {
            if (!m_explicitDiffusion)
            {
                m_ode.DefineImplicitSolve (&UnsteadyAdvectionDiffusion::DoImplicitSolve, this);
                m_ode.DefineOdeRhs        (&UnsteadyAdvectionDiffusion::DoOdeRhs,        this);
                m_ode.DefineProjection    (&UnsteadyAdvectionDiffusion::DoOdeProjection, this);
            }
            else
            {
                ASSERTL0(false, "Explicit diffusion with explicit reaction option not set up.");
            }
        }
        else
        {
            ASSERTL0(false, "Implicit reaction schemes not set up.");
        }
    }

    UnsteadyAdvectionDiffusion::~UnsteadyAdvectionDiffusion()
    {

    }

    void UnsteadyAdvectionDiffusion::DoOdeRhs(
            const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                  Array<OneD,        Array<OneD, NekDouble> >&outarray,
            const NekDouble time)
    {
        int i;
        int nvariables = inarray.num_elements();
        int npoints = GetNpoints();

        switch (m_projectionType)
        {
        case eDiscontinuousGalerkin:
            {
                int ncoeffs    = GetNcoeffs();
                Array<OneD, Array<OneD, NekDouble> > WeakAdv(nvariables);

                WeakAdv[0] = Array<OneD, NekDouble>(ncoeffs*nvariables);
                for(i = 1; i < nvariables; ++i)
                {
                    WeakAdv[i] = WeakAdv[i-1] + ncoeffs;
                }

                WeakDGAdvection(inarray, WeakAdv,true,true);

                for(i = 0; i < nvariables; ++i)
                {
                    m_fields[i]->MultiplyByElmtInvMass(WeakAdv[i],
                                                       WeakAdv[i]);
                    m_fields[i]->BwdTrans(WeakAdv[i],outarray[i]);
                    Vmath::Neg(npoints,outarray[i],1);
                }

                break;
            }
            case eGalerkin:
            {
                // Calculate -V\cdot Grad(u);
                for(i = 0; i < nvariables; ++i)
                {
                    AdvectionNonConservativeForm(m_velocity,
                                                 inarray[i],
                                                 outarray[i]);
                    Vmath::Neg(npoints,outarray[i],1);
                }
                break;
            }
        }
    }

    void UnsteadyAdvectionDiffusion::DoImplicitSolve(
            const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                  Array<OneD, Array<OneD, NekDouble> >&outarray,
            const NekDouble time,
            const NekDouble lambda)
    {
        int nvariables = inarray.num_elements();
        int nq = m_fields[0]->GetNpoints();

        // We solve ( \nabla^2 - HHlambda ) Y[i] = rhs [i]
        // inarray = input: \hat{rhs} -> output: \hat{Y}
        // outarray = output: nabla^2 \hat{Y}
        // where \hat = modal coeffs
        for (int i = 0; i < nvariables; ++i)
        {
            // Multiply 1.0/timestep/lambda
            Vmath::Smul(nq, -1.0/lambda/m_epsilon, inarray[i], 1, m_fields[i]->UpdatePhys(), 1);

            NekDouble kappa = 1.0/lambda/m_epsilon;

            // Solve a system of equations with Helmholtz solver
            m_fields[i]->HelmSolve(m_fields[i]->GetPhys(),
                                            m_fields[i]->UpdateCoeffs(),kappa);
            m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(), m_fields[i]->UpdatePhys());
            m_fields[i]->SetPhysState(false);

            // The solution is Y[i]
            outarray[i] = m_fields[i]->GetPhys();
        }
    }


    /**
     *
     */
    void UnsteadyAdvectionDiffusion::DoOdeProjection(const Array<OneD,
                                            const Array<OneD, NekDouble> >&inarray,
                                            Array<OneD,       Array<OneD, NekDouble> >&outarray,
                                            const NekDouble time)
    {
        int i;
        int nvariables = inarray.num_elements();
        SetBoundaryConditions(time);

        switch(m_projectionType)
        {
        case eDiscontinuousGalerkin:
            {
                // Just copy over array
                int npoints = GetNpoints();

                for(i = 0; i < nvariables; ++i)
                {
                    Vmath::Vcopy(npoints,inarray[i],1,outarray[i],1);
                }
            }
            break;
        case eGalerkin:
            {
                Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs());

                for(i = 0; i < nvariables; ++i)
                {
                    m_fields[i]->FwdTrans(inarray[i],coeffs,false);
                    m_fields[i]->BwdTrans_IterPerExp(coeffs,outarray[i]);
                }
                break;
            }
        default:
            ASSERTL0(false,"Unknown projection scheme");
            break;
        }
    }



    void UnsteadyAdvectionDiffusion::v_GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> > &physfield,
                           Array<OneD, Array<OneD, NekDouble> > &flux)
    {
        ASSERTL1(flux.num_elements() == mVelocity.num_elements(),"Dimension of flux array and velocity array do not match");

        for(int j = 0; j < flux.num_elements(); ++j)
          {
            Vmath::Vmul(GetNpoints(),physfield[i],1,
                m_velocity[j],1,flux[j],1);
          }
    }


    void UnsteadyAdvectionDiffusion::v_NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, Array<OneD, Array<OneD, NekDouble> > &numflux)
    {
        int i;

        int nTraceNumPoints = GetTraceNpoints();
        int nvel = m_spacedim; //m_velocity.num_elements();

        Array<OneD, NekDouble > Fwd(nTraceNumPoints);
        Array<OneD, NekDouble > Bwd(nTraceNumPoints);
        Array<OneD, NekDouble > Vn (nTraceNumPoints,0.0);

        // Get Edge Velocity - Could be stored if time independent
        for(i = 0; i < nvel; ++i)
        {
            m_fields[0]->ExtractTracePhys(m_velocity[i], Fwd);
            Vmath::Vvtvp(nTraceNumPoints,m_traceNormals[i],1,Fwd,1,Vn,1,Vn,1);
        }

        for(i = 0; i < numflux.num_elements(); ++i)
        {
            m_fields[i]->GetFwdBwdTracePhys(physfield[i],Fwd,Bwd);
            //evaulate upwinded m_fields[i]
            m_fields[i]->GetTrace()->Upwind(Vn,Fwd,Bwd,numflux[i]);
            // calculate m_fields[i]*Vn
            Vmath::Vmul(nTraceNumPoints,numflux[i],1,Vn,1,numflux[i],1);
        }
    }


    void UnsteadyAdvectionDiffusion::v_PrintSummary(std::ostream &out)
    {
        UnsteadySystem::v_PrintSummary(out);
    }
}
