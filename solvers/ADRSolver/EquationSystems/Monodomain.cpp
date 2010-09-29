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
     *
     */
    Monodomain::Monodomain(SessionReaderSharedPtr& pSession)
        : UnsteadySystem(pSession)
    {
        pSession->LoadParameter("epsilon",    m_epsilon,   1.0);

        if (!m_explicitDiffusion)
        {
            m_ode.DefineImplicitSolve (&Monodomain::DoImplicitSolve, this);
        }
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


    /**
     *
     */
    void Monodomain::v_PrintSummary(std::ostream &out)
    {
        UnsteadySystem::v_PrintSummary(out);
        out << "\tEpsilon         : " << m_epsilon << endl;
    }

}
