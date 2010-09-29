#include <iostream>

#include <ADRSolver/EquationSystems/MonodomainAlievPanfilov.h>

namespace Nektar
{
    /**
     * @class MonodomainAlievPanfilov
     *
     * The Aliev-Panfilov model of cardiac conduction is an improvement on the
     * FitzHugh-Nagumo model of conduction in squid axons which better
     * represents the upstroke and restitution properties of cardiac cells.
     * It is a mono-domain model of the form
     * \f{align*}{
     *     \frac{\partial u}{\partial t} &= \nabla^2 u + f(u,v), \\
     *     \frac{\partial v}{\partial t} &= g(u,v),
     * \f}
     * where the reaction terms, \f$f(u,v)\f$ and \f$g(u,v)\f$ are given by
     * \f{align*}{
     *      f(u,v) &= -ku(u-a)(1-u) - uv \\
     *      g(u,v) &= \left(\epsilon_0 + \frac{\mu_1v}{\mu_2 + u} \right)
     *                \left(-v-ku(u-a-1)\right)
     * \f}
     *
     * This implementation, at present, treats the reaction terms explicitly
     * and the diffusive element implicitly.
     */

    /**
     * Registers the class with the Factory.
     */
    string MonodomainAlievPanfilov::className
            = EquationSystemFactory::RegisterCreatorFunction(
                "AlievPanfilov",
                MonodomainAlievPanfilov::create,
                "Phenomological model of canine cardiac electrophysiology.");


    /**
     *
     */
    MonodomainAlievPanfilov::MonodomainAlievPanfilov(
                                            SessionReaderSharedPtr& pSession)
        : Monodomain(pSession)
    {
        int nq = m_fields[0]->GetNpoints();
        pSession->LoadParameter("k",          m_k,         0.0);
        pSession->LoadParameter("a",          m_a,         0.0);
        pSession->LoadParameter("mu1",        m_mu1,       0.0);
        pSession->LoadParameter("mu2",        m_mu2,       0.0);
        pSession->LoadParameter("eps",        m_eps,       0.0);

        m_uu   = Array<OneD, NekDouble>(nq, 0.0);
        m_uuu  = Array<OneD, NekDouble>(nq, 0.0);
        m_tmp1 = Array<OneD, NekDouble>(nq, 0.0);
        m_tmp2 = Array<OneD, NekDouble>(nq, 0.0);

        m_ode.DefineOdeRhs(&MonodomainAlievPanfilov::DoOdeRhs, this);
    }


    /**
     *
     */
    MonodomainAlievPanfilov::~MonodomainAlievPanfilov()
    {

    }


    /**
     * @param   inarray         Input array.
     * @param   outarray        Output array after addition of reaction terms.
     * @param   time            Current simulation time.
     */
    void MonodomainAlievPanfilov::DoOdeRhs(
            const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                  Array<OneD,        Array<OneD, NekDouble> >&outarray,
            const NekDouble time)
    {
        // inarray[0] holds initial physical u values throughout
        // inarray[1] holds initial physical v values throughout
        int nvariables  = inarray.num_elements();
        int nq          = m_fields[0]->GetNpoints();

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
    void MonodomainAlievPanfilov::v_PrintSummary(std::ostream &out)
    {
        Monodomain::v_PrintSummary(out);
        out << "\tCell model      : Aliev-Panfilov" << endl;
        out << "\tk               : " << m_k << endl;
        out << "\ta               : " << m_a << endl;
        out << "\teps             : " << m_eps << endl;
        out << "\tmu1             : " << m_mu1 << endl;
        out << "\tmu2             : " << m_mu2 << endl;
    }

}
