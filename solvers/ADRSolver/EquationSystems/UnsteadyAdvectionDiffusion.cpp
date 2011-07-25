#include <iostream>

#include <ADRSolver/EquationSystems/UnsteadyAdvectionDiffusion.h>

namespace Nektar
{
    string UnsteadyAdvectionDiffusion::className = GetEquationSystemFactory().RegisterCreatorFunction("UnsteadyAdvectionDiffusion", UnsteadyAdvectionDiffusion::create);

    UnsteadyAdvectionDiffusion::UnsteadyAdvectionDiffusion(
            LibUtilities::CommSharedPtr& pComm,
            LibUtilities::SessionReaderSharedPtr& pSession)
        : UnsteadySystem(pComm, pSession)
    {
    }

    void UnsteadyAdvectionDiffusion::v_InitObject()
    {
        UnsteadySystem::v_InitObject();

        m_session->LoadParameter("wavefreq",   m_waveFreq, 0.0);
        m_session->LoadParameter("epsilon",    m_epsilon,  0.0);

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

		m_ode.DefineImplicitSolve (&UnsteadyAdvectionDiffusion::DoImplicitSolve, this);
		m_ode.DefineOdeRhs        (&UnsteadyAdvectionDiffusion::DoOdeRhs,        this);
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
        case MultiRegions::eDiscontinuousGalerkin:
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
            case MultiRegions::eGalerkin:
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
		
		NekDouble kappa = 1.0/lambda/m_epsilon;
		
		Array<OneD, Array< OneD, NekDouble> > F(nvariables);
		F[0] = Array<OneD, NekDouble> (nq*nvariables);
        for(int n = 1; n < nvariables; ++n)
        {
            F[n] = F[n-1] + nq;
        }

        // We solve ( \nabla^2 - HHlambda ) Y[i] = rhs [i]
        // inarray = input: \hat{rhs} -> output: \hat{Y}
        // outarray = output: nabla^2 \hat{Y}
        // where \hat = modal coeffs
        for (int i = 0; i < nvariables; ++i)
        {
            // Multiply 1.0/timestep/lambda
            Vmath::Smul(nq, -kappa, inarray[i], 1, F[i], 1);
		}
            
	    //Setting boundary conditions
		SetBoundaryConditions(time);
	    
		for (int i = 0; i < nvariables; ++i)
		{
            // Solve a system of equations with Helmholtz solver
            m_fields[i]->HelmSolve(F[i],m_fields[i]->UpdateCoeffs(),kappa);
            m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),outarray[i]);
        }
    }

    void UnsteadyAdvectionDiffusion::v_GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> > &physfield,
                           Array<OneD, Array<OneD, NekDouble> > &flux)
    {
        ASSERTL1(flux.num_elements() == m_velocity.num_elements(),"Dimension of flux array and velocity array do not match");

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
