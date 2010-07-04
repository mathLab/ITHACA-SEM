#include <iostream>

#include <ADRSolver/EquationSystems/UnsteadyDiffusion.h>

namespace Nektar
{
    string UnsteadyDiffusion::className = EquationSystemFactory::RegisterCreatorFunction("UnsteadyDiffusion", UnsteadyDiffusion::create);

    UnsteadyDiffusion::UnsteadyDiffusion(SessionReaderSharedPtr& pSession,
            LibUtilities::TimeIntegrationSchemeOperators& pOde)
        : EquationSystem(pSession, pOde)
    {
        pSession->loadParameter("wavefreq",   mWaveFreq, 0.0);

        mTimeIntMethod = LibUtilities::eClassicalRungeKutta4;

        if (mExplicitDiffusion)
        {
            pOde.DefineOdeRhs        (&UnsteadyDiffusion::doOdeRhs,        this);
        }
        else
        {
            if(mWaveFreq>1000)
            {
                mTimeIntMethod = LibUtilities::eIMEXdirk_3_4_3;
            }
            else
            {
                mTimeIntMethod = LibUtilities::eDIRKOrder3;
            }
            pOde.DefineImplicitSolve (&UnsteadyDiffusion::doImplicitSolve, this);
        }
    }

    UnsteadyDiffusion::~UnsteadyDiffusion()
    {

    }

    void UnsteadyDiffusion::v_doOdeRhs(
            const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                  Array<OneD,        Array<OneD, NekDouble> >&outarray,
            const NekDouble time)
    {
        int i;
        int nvariables = inarray.num_elements();
        int ncoeffs    = inarray[0].num_elements();

        Array<OneD, Array<OneD, NekDouble> > Forcing(1);

        // BoundaryConditions are imposed weakly at the
        // Diffusion operator

        WeakDGDiffusion(inarray,outarray);

        for(int i = 0; i < nvariables; ++i)
        {
            m_fields[i]->MultiplyByElmtInvMass(outarray[i],
                                                outarray[i]);
        }

        if(mWaveFreq>0)
        {
            Forcing[0] = Array<OneD, NekDouble> (ncoeffs);
            doReaction(outarray,Forcing,time);
            Vmath::Vadd(ncoeffs, Forcing[0], 1, outarray[0], 1, outarray[0], 1);
        }
    }

    void UnsteadyDiffusion::doReaction(
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

    void UnsteadyDiffusion::v_GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> > &physfield,
                           Array<OneD, Array<OneD, NekDouble> > &flux)
    {
        ASSERTL1(flux.num_elements() == mVelocity.num_elements(),"Dimension of flux array and velocity array do not match");

        for(int j = 0; j < flux.num_elements(); ++j)
          {
            Vmath::Vmul(GetNpoints(),physfield[i],1,
                mVelocity[j],1,flux[j],1);
          }
    }

   // Evaulate flux = m_fields*ivel for i th component of Vu for direction j
    void UnsteadyDiffusion::v_GetFluxVector(const int i, const int j, Array<OneD, Array<OneD, NekDouble> > &physfield,
                           Array<OneD, Array<OneD, NekDouble> > &flux)
    {
        ASSERTL1(flux.num_elements() == mVelocity.num_elements(),"Dimension of flux array and velocity array do not match");

        for(int k = 0; k < flux.num_elements(); ++k)
          {
            Vmath::Zero(GetNpoints(),flux[k],1);
          }
        Vmath::Vcopy(GetNpoints(),physfield[i],1,flux[j],1);
    }



}
