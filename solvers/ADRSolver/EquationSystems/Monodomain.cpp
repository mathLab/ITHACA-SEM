#include <iostream>

#include <ADRSolver/EquationSystems/Monodomain.h>

namespace Nektar
{
    /**
     * @class Monodomain
     *
     * Base model of cardiac electrophysiology of the form
     * \f{align*}{
     *     \frac{\partial u}{\partial t} = \nabla^2 u + J_{ion},
     * \f}
     * where the reaction term, \f$J_{ion}\f$ is defined by a specific cell
     * model.
     *
     * This implementation, at present, treats the reaction terms explicitly
     * and the diffusive element implicitly.
     */

    /**
     * Registers the class with the Factory.
     */
    string Monodomain::className
            = EquationSystemFactory::RegisterCreatorFunction(
                "Monodomain",
                Monodomain::create,
                "Phenomological model of canine cardiac electrophysiology.");


    /**
     *
     */
    Monodomain::Monodomain(SessionReaderSharedPtr& pSession)
        : UnsteadySystem(pSession)
    {
        pSession->LoadParameter("epsilon",    m_epsilon,   1.0);

        std::string vCellModel;
        pSession->LoadSolverInfo("CELLMODEL", vCellModel, "");

        ASSERTL0(vCellModel != "", "Cell Model not specified.");

        m_cell = CellModelFactory::CreateInstance(vCellModel, pSession, m_fields[0]->GetNpoints());

        if (!m_explicitDiffusion)
        {
            m_ode.DefineImplicitSolve (&Monodomain::DoImplicitSolve, this);
        }
        m_ode.DefineOdeRhs(&Monodomain::DoOdeRhs, this);
    }


    /**
     *
     */
    Monodomain::~Monodomain()
    {

    }


    /**
     * @param   inarray         Input array.
     * @param   outarray        Output array.
     * @param   time            Current simulation time.
     * @param   lambda          Timestep.
     */
    void Monodomain::DoImplicitSolve(
            const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                  Array<OneD, Array<OneD, NekDouble> >&outarray,
            const NekDouble time,
            const NekDouble lambda)
    {
        int nvariables  = inarray.num_elements();
        int ncoeffs     = inarray[0].num_elements();
        int nq          = m_fields[0]->GetNpoints();
        NekDouble kappa = 1.0/lambda/m_epsilon;

        // We solve ( \nabla^2 - HHlambda ) Y[i] = rhs [i]
        // inarray = input: \hat{rhs} -> output: \hat{Y}
        // outarray = output: nabla^2 \hat{Y}
        // where \hat = modal coeffs
        for (int i = 0; i < nvariables; ++i)
        {
            // Only apply diffusion to first variable.
            if (i > 0) {
                Vmath::Vcopy(nq, &inarray[i][0], 1, &outarray[i][0], 1);
                continue;
            }

            // Multiply 1.0/timestep/lambda
            Vmath::Smul(nq, -1.0/lambda/m_epsilon, inarray[i], 1,
                                            m_fields[i]->UpdatePhys(), 1);

            // Solve a system of equations with Helmholtz solver and transform
            // back into physical space.
            m_fields[i]->HelmSolve(m_fields[i]->GetPhys(),
                                   m_fields[i]->UpdateCoeffs(),kappa);

            m_fields[i]->BwdTrans( m_fields[i]->GetCoeffs(),
                                   m_fields[i]->UpdatePhys());
            m_fields[i]->SetPhysState(true);

            // Copy the solution vector (required as m_fields must be set).
            outarray[i] = m_fields[i]->GetPhys();
        }
    }


    void Monodomain::DoOdeRhs(
            const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                  Array<OneD,        Array<OneD, NekDouble> >&outarray,
            const NekDouble time)
    {
        m_cell->Update(inarray, outarray, time);
    }


    void Monodomain::v_SetInitialConditions(NekDouble initialtime,
                        bool dumpInitialConditions)
    {
        ADRBase::v_SetInitialConditions(initialtime, dumpInitialConditions);
        /*
        cout << "Set initial conditions." << endl;
        int nq = m_fields[0]->GetNpoints();
        Array<OneD,NekDouble> x0(nq);
        Array<OneD,NekDouble> x1(nq);
        Array<OneD,NekDouble> x2(nq);

        m_fields[0]->GetCoords(x0,x1,x2);
        NekDouble rad;
        NekDouble mDiam = 4.0;
        NekDouble m_x0c = 14.4;
        NekDouble m_x1c = 65.0;
        NekDouble m_x2c = 15.5;
        for(int j = 0; j < nq; j++)
        {
            rad = sqrt( (x0[j]-m_x0c)*(x0[j]-m_x0c) + (x1[j]-m_x1c)*(x1[j]-m_x1c) + (x2[j]-m_x2c)*(x2[j]-m_x2c) );

            if( rad <= mDiam )
            {
                (m_fields[0]->UpdatePhys())[j] = 1.0;
            }
        }

        for(int i = 0 ; i < m_fields.num_elements(); i++)
        {
            m_fields[i]->SetPhysState(true);
            m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),m_fields[i]->UpdateCoeffs());
        }
         */
    }


    /**
     *
     */
    void Monodomain::v_PrintSummary(std::ostream &out)
    {
        UnsteadySystem::v_PrintSummary(out);
        out << "\tEpsilon         : " << m_epsilon << endl;
        m_cell->v_PrintSummary(out);
    }


    /**
     * Cell model base class constructor.
     */
    CellModel::CellModel(SessionReaderSharedPtr& pSession, const int nq)
    {
        m_spatialParameters = MemoryManager<SpatialDomains::SpatialParameters>
                                          ::AllocateSharedPtr(nq);
        //m_spatialParameters->Read(m_filename);

        //Array<OneD, NekDouble> x(nq), y(nq), z(nq);
        //m_fields[0]->GetCoords(x,y,z);
        //m_spatialParameters->EvaluateParameters(x,y,z);
    }

    /**
     * Registers the class with the Factory.
     */
    string CellModelAlievPanfilov::className
            = CellModelFactory::RegisterCreatorFunction(
                "AlievPanfilov",
                CellModelAlievPanfilov::create,
                "Phenomological model of canine cardiac electrophysiology.");

    CellModelAlievPanfilov::CellModelAlievPanfilov(
                    SessionReaderSharedPtr& pSession, const int nq)
            : CellModel(pSession, nq)
    {
        pSession->LoadParameter("k",          m_k,         0.0);
        pSession->LoadParameter("a",          m_a,         0.0);
        pSession->LoadParameter("mu1",        m_mu1,       0.0);
        pSession->LoadParameter("mu2",        m_mu2,       0.0);
        pSession->LoadParameter("eps",        m_eps,       0.0);

        m_nq   = nq;
        m_uu   = Array<OneD, NekDouble>(nq, 0.0);
        m_uuu  = Array<OneD, NekDouble>(nq, 0.0);
        m_tmp1 = Array<OneD, NekDouble>(nq, 0.0);
        m_tmp2 = Array<OneD, NekDouble>(nq, 0.0);
    }


    void CellModelAlievPanfilov::Update(
                    const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                          Array<OneD,        Array<OneD, NekDouble> >&outarray,
                    const NekDouble time)
    {
        // inarray[0] holds initial physical u values throughout
        // inarray[1] holds initial physical v values throughout
        int nvariables  = inarray.num_elements();
        int nq = m_nq;

        // compute u^2: m_u = u*u
        Vmath::Vmul(nq, &inarray[0][0], 1, &inarray[0][0], 1, &m_uu[0], 1);

        // compute u^3: m_u = u*u*u
        Vmath::Vmul(nq, &inarray[0][0], 1, &m_uu[0], 1, &m_uuu[0], 1);

        // --------------------------------------
        // Compute reaction term f(u,v)
        // --------------------------------------
        if (m_spatialParameters->Exists("a"))
        {
          Vmath::Vmul(nq,  &m_spatialParameters->GetData("a")->GetPhys()[0], 1,
                           &inarray[0][0], 1, &m_tmp1[0], 1);

          Vmath::Vvtvm(nq, &m_spatialParameters->GetData("a")->GetPhys()[0], 1,
                           &m_uu[0], 1, &m_tmp1[0], 1, &m_tmp1[0], 1);

          Vmath::Svtvm(nq, -1.0, &m_uu[0], 1, &m_tmp1[0], 1, &m_tmp1[0], 1);
        }
        else
        {
          // Ru = au
          Vmath::Smul(nq, m_a, &inarray[0][0], 1, &m_tmp1[0], 1);
          // Ru = (-1-a)u*u + au
          Vmath::Svtvp(nq, (-1.0-m_a), &m_uu[0], 1, &m_tmp1[0], 1,
                                       &m_tmp1[0], 1);
        }
        // Ru = u*u*u - (1+a)u*u + au
        Vmath::Vadd(nq, &m_uuu[0], 1, &m_tmp1[0], 1, &m_tmp1[0], 1);
        // Ru = k(u*u*u - (1+a)u*u + au)
        if (m_spatialParameters->Exists("k"))
        {
          Vmath::Vmul(nq, &m_spatialParameters->GetData("k")->GetPhys()[0], 1,
                          &m_tmp1[0], 1, &m_tmp1[0], 1);
        }
        else
        {
          Vmath::Smul(nq, m_k, &m_tmp1[0], 1, &m_tmp1[0], 1);
        }
        // Ru = k(u*u*u - (1+a)u*u + au) + uv
        Vmath::Vvtvp(nq, &inarray[0][0], 1, &inarray[1][0], 1, &m_tmp1[0], 1,
                         &outarray[0][0], 1);
        // Ru = -k(u*u*u - (1+a)u*u + au) - uv
        Vmath::Neg(nq, &outarray[0][0], 1);


        // --------------------------------------
        // Compute reaction term g(u,v)
        // --------------------------------------
        // tmp2 = mu2 + u
        Vmath::Sadd(nq, m_mu2, &inarray[0][0], 1, &m_tmp2[0], 1);

        // tmp2 = v/(mu2 + u)
        Vmath::Vdiv(nq, &inarray[1][0], 1, &m_tmp2[0], 1, &m_tmp2[0], 1);

        // tmp2 = mu1*v/(mu2 + u)
        Vmath::Smul(nq, m_mu1, &m_tmp2[0], 1, &m_tmp2[0], 1);

        // tmp1 = Eps + mu1*v/(mu2+u)
        Vmath::Sadd(nq, m_eps, &m_tmp2[0], 1, &m_tmp2[0], 1);

        // tmp1 = (-a-1) + u
        if (m_spatialParameters->Exists("a"))
        {
          Vmath::Vsub(nq, &inarray[0][0], 1,
                          &m_spatialParameters->GetData("a")->GetPhys()[0], 1,
                          &m_tmp1[0], 1);

          Vmath::Sadd(nq, -1.0, &inarray[0][0], 1, &m_tmp1[0], 1);
        }
        else
        {
          Vmath::Sadd(nq, (-m_a-1), &inarray[0][0], 1, &m_tmp1[0], 1);
        }

        // tmp1 = k(u-a-1)
        if (m_spatialParameters->Exists("k"))
        {
          Vmath::Vmul(nq, &m_spatialParameters->GetData("k")->GetPhys()[0], 1,
                          &m_tmp1[0], 1, &m_tmp1[0], 1);
        }
        else
        {
          Vmath::Smul(nq, m_k, &m_tmp1[0], 1, &m_tmp1[0], 1);
        }

        // tmp1 = ku(u-a-1) + v
        Vmath::Vvtvp(nq, &inarray[0][0], 1, &m_tmp1[0], 1, &inarray[1][0], 1,
                         &m_tmp1[0], 1);

        // tmp1 = -ku(u-a-1)-v
        Vmath::Neg(nq, &m_tmp1[0], 1);

        // outarray = [Eps + mu1*v/(mu2+u)] * [-ku(u-a-1)-v]
        Vmath::Vmul(nq, &m_tmp1[0], 1, &m_tmp2[0], 1, &outarray[1][0], 1);
    }

    /**
     *
     */
    void CellModelAlievPanfilov::v_PrintSummary(std::ostream &out)
    {
        out << "\tCell model      : Aliev-Panfilov" << endl;
        out << "\tk               : " << m_k << endl;
        out << "\ta               : " << m_a << endl;
        out << "\teps             : " << m_eps << endl;
        out << "\tmu1             : " << m_mu1 << endl;
        out << "\tmu2             : " << m_mu2 << endl;
    }


    /**
     * Registers the class with the Factory.
     */
    string CellModelFitzHughNagumo::className
            = CellModelFactory::RegisterCreatorFunction(
                "FitzHughNagumo",
                CellModelFitzHughNagumo::create,
                "Phenomological model of squid nerve cell.");

    CellModelFitzHughNagumo::CellModelFitzHughNagumo(
                    SessionReaderSharedPtr& pSession, const int nq)
            : CellModel(pSession, nq)
    {
        pSession->LoadParameter("beta",          m_beta,         0.0);
        pSession->LoadParameter("epsilon",       m_epsilon,      1.0);

        m_uuu  = Array<OneD, NekDouble>(nq, 0.0);
    }


    void CellModelFitzHughNagumo::Update(
                    const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                          Array<OneD,        Array<OneD, NekDouble> >&outarray,
                    const NekDouble time)
    {
        NekDouble m_gamma = 0.5;

        int nvariables  = inarray.num_elements();
        int nq          = m_nq;

        // compute u^2: m_u = u*u
        Vmath::Vmul(nq, &inarray[0][0], 1, &inarray[0][0], 1, &m_uuu[0], 1);

        // compute u^3: m_u = u*u*u
        Vmath::Vmul(nq, &inarray[0][0], 1, &m_uuu[0], 1, &m_uuu[0], 1);

        // For u: (1/m_epsilon)*( u*-u*u*u/3 - v )
        // physfield = u - (1.0/3.0)*u*u*u
        Vmath::Svtvp(nq, (-1.0/3.0), &m_uuu[0], 1, &inarray[0][0], 1, &outarray[1][0], 1);

        Vmath::Vsub(nq, &inarray[1][0], 1, &outarray[1][0], 1, &outarray[1][0], 1);
        Vmath::Smul(nq, -1.0/m_epsilon, &outarray[1][0], 1, &outarray[1][0], 1);

        // For v: m_epsilon*( u + m_beta - m_gamma*v )
        Vmath::Svtvp(nq, -1.0*m_gamma, &inarray[1][0], 1, &inarray[0][0], 1, &outarray[1][0], 1);
        Vmath::Sadd(nq, m_beta, &outarray[1][0], 1, &outarray[1][0], 1);
        Vmath::Smul(nq, m_epsilon, &outarray[1][0], 1, &outarray[1][0], 1);
    }

    /**
     *
     */
    void CellModelFitzHughNagumo::v_PrintSummary(std::ostream &out)
    {
        out << "\tCell model      : FitzHugh-Nagumo" << endl;
        out << "\tBeta            : " << m_beta << endl;
    }

}
