#include <iostream>

#include <ADRSolver/EquationSystems/UnsteadyInviscidBurger.h>

namespace Nektar
{
    string UnsteadyInviscidBurger::className = EquationSystemFactory::RegisterCreatorFunction("UnsteadyInviscidBurger", UnsteadyInviscidBurger::create);

    UnsteadyInviscidBurger::UnsteadyInviscidBurger(SessionReaderSharedPtr& pSession,
            LibUtilities::TimeIntegrationSchemeOperators& pOde)
        : EquationSystem(pSession, pOde)
    {
        pSession->loadParameter("wavefreq",   mWaveFreq, 0.0);

        mTimeIntMethod = LibUtilities::eClassicalRungeKutta4;

        if (mExplicitAdvection)
        {
            pOde.DefineOdeRhs        (&UnsteadyInviscidBurger::doOdeRhs,        this);
        }
        else
        {
            ASSERTL0(false, "Implicit unsteady Advection not set up.");
        }
    }

    UnsteadyInviscidBurger::~UnsteadyInviscidBurger()
    {

    }

    void UnsteadyInviscidBurger::v_doOdeRhs(
            const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                  Array<OneD,        Array<OneD, NekDouble> >&outarray,
            const NekDouble time)
    {
        int i;
        int nvariables = inarray.num_elements();
        int ncoeffs    = inarray[0].num_elements();

        switch (m_projectionType)
        {
            case eDiscontinuousGalerkin:
            {
                Array<OneD, Array<OneD, NekDouble> > Forcing(1);

                setBoundaryConditions(time);
                WeakDGAdvection(inarray, outarray);
                for(i = 0; i < nvariables; ++i)
                {
                    m_fields[i]->MultiplyByElmtInvMass(outarray[i],
                                                       outarray[i]);
                    Vmath::Neg(ncoeffs,outarray[i],1);
                }

                if(mWaveFreq>0)
                {
                    Forcing[0] = Array<OneD, NekDouble> (ncoeffs);
                    doReaction(outarray,Forcing,time);
                    Vmath::Vadd(ncoeffs, Forcing[0], 1, outarray[0], 1, outarray[0], 1);
                }
                break;
            }
            case eGalerkin:
            {
                setBoundaryConditions(time);
                Array<OneD, NekDouble> physfield(GetNpoints());

                for(i = 0; i < nvariables; ++i)
                {
                    m_fields[i]->MultiplyByInvMassMatrix(inarray[i],
                                                    outarray[i], false);
                    // Calculate -(\phi, V\cdot Grad(u))
                    m_fields[i]->BwdTrans_IterPerExp(outarray[i],
                                                        physfield);

                    WeakAdvectionNonConservativeForm(mVelocity,
                                                physfield, outarray[i]);

                    Vmath::Neg(ncoeffs,outarray[i],1);
                }
                break;
            }
        }
    }

    void UnsteadyInviscidBurger::doReaction(
            const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                  Array<OneD, Array<OneD, NekDouble> >&outarray,
            const NekDouble time)
    {
        int i,k;
        int nvariables = inarray.num_elements();
        int ncoeffs    = inarray[0].num_elements();

        // PI*PI*exp(-1.0*PI*PI*FinTime)*sin(PI*x)*sin(PI*y)
        int nq = m_fields[0]->GetNpoints();

        Array<OneD,NekDouble> x0(nq);
        Array<OneD,NekDouble> x1(nq);
        Array<OneD,NekDouble> x2(nq);

        // get the coordinates (assuming all fields have the same
        // discretisation)
        m_fields[0]->GetCoords(x0,x1,x2);

        Array<OneD, NekDouble> physfield(nq);

        NekDouble kt, kx, ky;
        for (i=0; i<nq; ++i)
        {
              kt = mWaveFreq*time;
              kx = mWaveFreq*x0[i];
              ky = mWaveFreq*x1[i];

              // F(x,y,t) = du/dt + V \cdot \nabla u - \varepsilon \nabla^2 u
              physfield[i] = (2.0*mEpsilon*mWaveFreq*mWaveFreq + mWaveFreq*cos(kt))*exp(sin(kt))*sin(kx)*sin(ky);
        }
        m_fields[0]->FwdTrans(physfield, outarray[0]);
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
