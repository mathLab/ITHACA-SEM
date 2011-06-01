#include <iostream>

#include <ADRSolver/EquationSystems/UnsteadyDiffusion.h>

namespace Nektar
{
    string UnsteadyDiffusion::className = GetEquationSystemFactory().RegisterCreatorFunction("UnsteadyDiffusion", UnsteadyDiffusion::create);

    UnsteadyDiffusion::UnsteadyDiffusion(
            LibUtilities::CommSharedPtr& pComm,
            LibUtilities::SessionReaderSharedPtr& pSession)
        : UnsteadySystem(pComm, pSession)
    {
        pSession->LoadParameter("wavefreq",   m_waveFreq, 0.0);
        pSession->LoadParameter("epsilon",    m_epsilon,  0.0);

        if (m_explicitDiffusion)
        {
            m_ode.DefineOdeRhs        (&UnsteadyDiffusion::DoOdeRhs,        this);
            m_ode.DefineProjection    (&UnsteadyDiffusion::DoOdeProjection, this);
        }
        else
        {
            m_ode.DefineImplicitSolve (&UnsteadyDiffusion::DoImplicitSolve, this);
        }
    }

    UnsteadyDiffusion::~UnsteadyDiffusion()
    {

    }

    void UnsteadyDiffusion::DoOdeRhs(
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
                int ncoeffs    = inarray[0].num_elements();
                Array<OneD, Array<OneD, NekDouble> > WeakAdv(nvariables);

                WeakAdv[0] = Array<OneD, NekDouble>(ncoeffs*nvariables);
                for(i = 1; i < nvariables; ++i)
                {
                    WeakAdv[i] = WeakAdv[i-1] + ncoeffs;
                }

                WeakDGDiffusion(inarray, WeakAdv,true,true);

                for(i = 0; i < nvariables; ++i)
                {
                    m_fields[i]->MultiplyByElmtInvMass(WeakAdv[i],
                                                       WeakAdv[i]);
                    m_fields[i]->BwdTrans(WeakAdv[i],outarray[i]);
                }

                break;
            }
        case eGalerkin:
            {
                ASSERTL0(false, "Explicit Galerkin diffusion not set up.");
            }
        }
    }

    void UnsteadyDiffusion::DoImplicitSolve(
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
    void UnsteadyDiffusion::DoOdeProjection(const Array<OneD,
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


}
