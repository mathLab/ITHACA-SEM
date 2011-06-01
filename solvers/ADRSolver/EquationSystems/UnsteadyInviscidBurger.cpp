#include <iostream>

#include <ADRSolver/EquationSystems/UnsteadyInviscidBurger.h>

namespace Nektar
{
    string UnsteadyInviscidBurger::className = GetEquationSystemFactory().RegisterCreatorFunction("UnsteadyInviscidBurger", UnsteadyInviscidBurger::create);

    UnsteadyInviscidBurger::UnsteadyInviscidBurger(
            LibUtilities::CommSharedPtr& pComm,
            LibUtilities::SessionReaderSharedPtr& pSession)
        : UnsteadySystem(pComm,pSession)
    {
        pSession->LoadParameter("wavefreq",   m_waveFreq, 0.0);

        if (m_explicitAdvection)
        {
            m_ode.DefineOdeRhs        (&UnsteadyInviscidBurger::DoOdeRhs,        this);
            m_ode.DefineProjection (&UnsteadyInviscidBurger::DoOdeProjection, this);
        }
        else
        {
            ASSERTL0(false, "Implicit unsteady Advection not set up.");
        }
    }

    UnsteadyInviscidBurger::~UnsteadyInviscidBurger()
    {

    }

    void UnsteadyInviscidBurger::DoOdeRhs(
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

                //SetBoundaryConditions(time);
                WeakDGAdvection(inarray, WeakAdv, true, true);

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

    /**
     *
     */
    void UnsteadyInviscidBurger::DoOdeProjection(const Array<OneD,
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


    void UnsteadyInviscidBurger::v_GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> > &physfield,
                           Array<OneD, Array<OneD, NekDouble> > &flux)
    {
        for(int j = 0; j < flux.num_elements(); ++j)
          {
            Vmath::Vmul(GetNpoints(),physfield[i],1,
                physfield[i],1,flux[j],1);
            Vmath::Smul(GetNpoints(),0.5,flux[j],1,flux[j],1);
          }
    }

   // Evaulate flux = m_fields*ivel for i th component of Vu for direction j
    void UnsteadyInviscidBurger::v_GetFluxVector(const int i, const int j, Array<OneD, Array<OneD, NekDouble> > &physfield,
                           Array<OneD, Array<OneD, NekDouble> > &flux)
    {
        ASSERTL0(false,"should never arrive here ...");
    }

    void UnsteadyInviscidBurger::v_NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, Array<OneD, Array<OneD, NekDouble> > &numflux)
    {
        int i;

        int nTraceNumPoints = GetTraceNpoints();
        int nvel = m_spacedim; //m_velocity.num_elements();

        Array<OneD, NekDouble > Fwd(nTraceNumPoints);
        Array<OneD, NekDouble > Bwd(nTraceNumPoints);
        Array<OneD, NekDouble > Vn (nTraceNumPoints,0.0);

        // Get Edge Velocity - Could be stored if time independent
        m_fields[0]->ExtractTracePhys(physfield[0], Fwd);
        for(i = 0; i < nvel; ++i)
      {
        // m_fields[0]->ExtractTracePhys(physfield[0], Fwd);
        Vmath::Vvtvp(nTraceNumPoints,m_traceNormals[i],1,Fwd,1,Vn,1,Vn,1);
      }

        for(i = 0; i < numflux.num_elements(); ++i)
      {
        m_fields[i]->GetFwdBwdTracePhys(physfield[i],Fwd,Bwd);
        //evaulate upwinded m_fields[i]
        m_fields[i]->GetTrace()->Upwind(Vn,Fwd,Bwd,numflux[i]);
        // calculate m_fields[i]*Vn
        Vmath::Vmul(nTraceNumPoints,numflux[i],1,Vn,1,numflux[i],1);
        Vmath::Smul(nTraceNumPoints,0.5,numflux[i],1,numflux[i],1);
      }
    }
}
