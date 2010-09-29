#include <iostream>

#include <ADRSolver/EquationSystems/MonodomainFitzHughNagumo.h>

namespace Nektar
{
    /**
     * @class MonodomainFitzhughNagumo
     *
     * The Fitzhugh-Nagumo model of cardiac conduction is an improvement on the
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
    string MonodomainFitzHughNagumo::className
            = EquationSystemFactory::RegisterCreatorFunction(
                    "FitzHughNagumo",
                    MonodomainFitzHughNagumo::create,
                    "Phenomological model of nerve cell electrophysiology.");


    /**
     *
     */
    MonodomainFitzHughNagumo::MonodomainFitzHughNagumo(
                                            SessionReaderSharedPtr& pSession)
        : Monodomain(pSession)
    {
        int nq = m_fields[0]->GetNpoints();

        pSession->LoadParameter("beta",          m_beta,         0.0);
        m_uuu  = Array<OneD, NekDouble>(nq, 0.0);
        m_tmp1 = Array<OneD, NekDouble>(nq, 0.0);

        m_ode.DefineOdeRhs (&MonodomainFitzHughNagumo::DoOdeRhs, this);
    }


    /**
     *
     */
    MonodomainFitzHughNagumo::~MonodomainFitzHughNagumo()
    {

    }


    /**
     * @param   inarray         Input array.
     * @param   outarray        Output array after addition of reaction terms.
     * @param   time            Current simulation time.
     */
    void MonodomainFitzHughNagumo::DoOdeRhs(
            const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                  Array<OneD,        Array<OneD, NekDouble> >&outarray,
            const NekDouble time)
    {
        NekDouble m_gamma = 0.5;

        int nvariables  = inarray.num_elements();
        int nq          = m_fields[0]->GetNpoints();

        // compute u^2: m_u = u*u
        Vmath::Vmul(nq, &inarray[0][0], 1, &inarray[0][0], 1, &m_uuu[0], 1);

        // compute u^3: m_u = u*u*u
        Vmath::Vmul(nq, &inarray[0][0], 1, &m_uuu[0], 1, &m_uuu[0], 1);

        // For u: (1/m_epsilon)*( u*-u*u*u/3 - v )
        // physfield = u - (1.0/3.0)*u*u*u
        Vmath::Svtvp(nq, (-1.0/3.0), &m_uuu[0], 1, &inarray[0][0], 1, &m_tmp1[0], 1);

        Vmath::Vsub(nq, &inarray[1][0], 1, &m_tmp1[0], 1, &m_tmp1[0], 1);
        Vmath::Smul(nq, -1.0/m_epsilon, &m_tmp1[0], 1, &m_tmp1[0], 1);

        m_fields[0]->FwdTrans(m_tmp1,outarray[0]);
        m_fields[0]->SetPhysState(false);

        // For v: m_epsilon*( u + m_beta - m_gamma*v )
        Vmath::Svtvp(nq, -1.0*m_gamma, &inarray[1][0], 1, &inarray[0][0], 1, &m_tmp1[0], 1);
        Vmath::Sadd(nq, m_beta, &m_tmp1[0], 1, &m_tmp1[0], 1);
        Vmath::Smul(nq, m_epsilon, &m_tmp1[0], 1, &m_tmp1[0], 1);

        m_fields[1]->FwdTrans(m_tmp1,outarray[1]);
        m_fields[1]->SetPhysState(false);
    }


    /**
     *
     */
    void MonodomainFitzHughNagumo::v_PrintSummary(std::ostream &out)
    {
        Monodomain::v_PrintSummary(out);
        out << "\tCell model      : FitzHugh-Nagumo" << endl;
        out << "\tBeta            : " << m_beta << endl;
    }

}
