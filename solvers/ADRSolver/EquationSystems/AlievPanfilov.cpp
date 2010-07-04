#include <iostream>

#include <ADRSolver/EquationSystems/AlievPanfilov.h>

namespace Nektar
{
    string AlievPanfilov::className = EquationSystemFactory::RegisterCreatorFunction("AlievPanfilov", AlievPanfilov::create);

    AlievPanfilov::AlievPanfilov(SessionReaderSharedPtr& pSession,
            LibUtilities::TimeIntegrationSchemeOperators& pOde)
        : EquationSystem(pSession, pOde)
    {
        pSession->loadParameter("k",          mK,         0.0);
        pSession->loadParameter("a",          mA,         0.0);
        pSession->loadParameter("mu1",        mMu1,       0.0);
        pSession->loadParameter("mu2",        mMu2,       0.0);
        pSession->loadParameter("eps",        mEps,       0.0);

        mTimeIntMethod = LibUtilities::eIMEXdirk_3_4_3;

        pOde.DefineOdeRhs        (&AlievPanfilov::doOdeRhs,        this);
        if (!mExplicitDiffusion)
        {
            pOde.DefineImplicitSolve (&AlievPanfilov::doImplicitSolve, this);
        }
    }

    AlievPanfilov::~AlievPanfilov()
    {

    }

    void AlievPanfilov::v_doOdeRhs(
            const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                  Array<OneD,        Array<OneD, NekDouble> >&outarray,
            const NekDouble time)
    {
        int nvariables = inarray.num_elements();
        int ncoeffs    = inarray[0].num_elements();
        int npoints    = m_fields[0]->GetNpoints();

        Array<OneD, NekDouble> physfieldu(npoints);
        Array<OneD, NekDouble> physfieldv(npoints);

        Array<OneD, NekDouble> Ru(npoints,0.0);
        Array<OneD, NekDouble> Rv(npoints, 0.0);
        Array<OneD, NekDouble> u2(npoints,0.0);
        Array<OneD, NekDouble> u3(npoints,0.0);

        m_fields[0]->BwdTrans(inarray[0],physfieldu);
        m_fields[0]->SetPhysState(true);

        m_fields[1]->BwdTrans(inarray[1],physfieldv);
        m_fields[1]->SetPhysState(true);

        // u2 = u*u
        Vmath::Vmul(npoints, &physfieldu[0], 1, &physfieldu[0], 1, &u2[0], 1);
        // u3 = u*u*u
        Vmath::Vmul(npoints, &physfieldu[0], 1, &u2[0], 1, &u3[0], 1);

        if (m_spatialParameters->Exists("a"))
        {
          Vmath::Vmul(npoints,  &m_spatialParameters->GetData("a")->GetPhys()[0], 1, &physfieldu[0], 1, &Ru[0], 1);
          Vmath::Vvtvm(npoints, &m_spatialParameters->GetData("a")->GetPhys()[0], 1, &u2[0], 1, &Ru[0], 1, &Ru[0], 1);
          Vmath::Svtvm(npoints, -1.0, &u2[0], 1, &Ru[0], 1, &Ru[0], 1);
        }
        else
        {
          // Ru = au
          Vmath::Smul(npoints, mA, &physfieldu[0], 1, &Ru[0], 1);
          // Ru = (-1-a)u*u + au
          Vmath::Svtvp(npoints, (-1.0-mA), &u2[0], 1, &Ru[0], 1, &Ru[0], 1);
        }
        // Ru = u*u*u - (1+a)u*u + au
        Vmath::Vadd(npoints, &u3[0], 1, &Ru[0], 1, &Ru[0], 1);
        // Ru = k(u*u*u - (1+a)u*u + au)
        if (m_spatialParameters->Exists("k"))
        {
          Vmath::Vmul(npoints, &m_spatialParameters->GetData("k")->GetPhys()[0], 1, &Ru[0], 1, &Ru[0], 1);
        }
        else
        {
          Vmath::Smul(npoints, mK, &Ru[0], 1, &Ru[0], 1);
        }
        // Ru = k(u*u*u - (1+a)u*u + au) + uv
        Vmath::Vvtvp(npoints, &physfieldu[0], 1, &physfieldv[0], 1, &Ru[0], 1, &Ru[0], 1);
        // Ru = -k(u*u*u - (1+a)u*u + au) - uv
        Vmath::Neg(npoints, &Ru[0], 1);

        m_fields[0]->FwdTrans(Ru,outarray[0]);
        m_fields[0]->SetPhysState(false);

        // For v:
        // Rv = mu2 + u
        Vmath::Sadd(npoints, mMu2, &physfieldu[0], 1, &Rv[0], 1);
        // Rv = v/(mu2 + u)
        Vmath::Vdiv(npoints, &physfieldv[0], 1, &Rv[0], 1, &Rv[0], 1);
        // Rv = mu1*v/(mu2 + u)
        Vmath::Smul(npoints, mMu1, &Rv[0], 1, &Rv[0], 1);
        // Rv = Eps + mu1*v/(mu2+u)
        Vmath::Sadd(npoints, mEps, &Rv[0], 1, &Rv[0], 1);

        // Ru = (-a-1) + u
        if (m_spatialParameters->Exists("a"))
        {
          Vmath::Vsub(npoints, &physfieldu[0], 1, &m_spatialParameters->GetData("a")->GetPhys()[0], 1, &Ru[0], 1);
          Vmath::Sadd(npoints, -1.0, &physfieldu[0], 1, &Ru[0], 1);
        }
        else
        {
          Vmath::Sadd(npoints, (-mA-1), &physfieldu[0], 1, &Ru[0], 1);
        }
        // Ru = k(u-a-1)
        if (m_spatialParameters->Exists("k"))
        {
          Vmath::Vmul(npoints, &m_spatialParameters->GetData("k")->GetPhys()[0], 1, &Ru[0], 1, &Ru[0], 1);
        }
        else
        {
          Vmath::Smul(npoints, mK, &Ru[0], 1, &Ru[0], 1);
        }
        // Ru = ku(u-a-1) + v
        Vmath::Vvtvp(npoints, &physfieldu[0], 1, &Ru[0], 1, &physfieldv[0], 1, &Ru[0], 1);
        // Ru = -ku(u-a-1)-v
        Vmath::Neg(npoints, &Ru[0], 1);

        Vmath::Vmul(npoints, &Ru[0], 1, &Rv[0], 1, &Rv[0], 1);

        m_fields[1]->FwdTrans(Rv,outarray[1]);
        m_fields[1]->SetPhysState(false);
    }

    void AlievPanfilov::v_doImplicitSolve(
            const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                  Array<OneD, Array<OneD, NekDouble> >&outarray,
            const NekDouble time,
            const NekDouble lambda)
    {
        int nvariables = inarray.num_elements();
        int ncoeffs    = inarray[0].num_elements();
        int nq = m_fields[0]->GetNpoints();

        Array<OneD, Array<OneD, NekDouble> > Forcing(1);
        Forcing[0] = Array<OneD, NekDouble>(ncoeffs,0.0);

        // We solve ( \nabla^2 - HHlambda ) Y[i] = rhs [i]
        // inarray = input: \hat{rhs} -> output: \hat{Y}
        // outarray = output: nabla^2 \hat{Y}
        // where \hat = modal coeffs

        //MultiRegions::GlobalMatrixKey key(StdRegions::eMass);

        for (int i = 0; i < nvariables; ++i)
        {
            // Only apply diffusion to first variable.
            if (i > 0) {
                Vmath::Vcopy(ncoeffs, &inarray[i][0], 1, &outarray[i][0], 1);
                continue;
            }

            // Multiply 1.0/timestep/lambda
            Vmath::Smul(ncoeffs, -1.0/lambda/mEpsilon, inarray[i], 1, outarray[i], 1);

            // Update coeffs to m_fields
            m_fields[i]->UpdateCoeffs() = outarray[i];

            // Backward Transformation to nodal coefficients
            m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(), m_fields[i]->UpdatePhys());
            // m_fields[i]->SetPhysState(true);

            NekDouble kappa = 1.0/lambda/mEpsilon;

            // Solve a system of equations with Helmholtz solver
            m_fields[i]->HelmSolve(m_fields[i]->GetPhys(),
                                            m_fields[i]->UpdateCoeffs(),kappa);
            m_fields[i]->SetPhysState(false);

            // The solution is Y[i]
            outarray[i] = m_fields[i]->GetCoeffs();
        }
    }

    void AlievPanfilov::v_printSummary(std::ostream &out)
    {
        out << "\tepsilon         : " << mEpsilon << endl;
        out << "\tk               : " << mK << endl;
        out << "\ta               : " << mA << endl;
        out << "\teps             : " << mEps << endl;
        out << "\tmu1             : " << mMu1 << endl;
        out << "\tmu2             : " << mMu2 << endl;
    }

}
