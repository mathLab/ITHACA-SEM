#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_MONODOMAINALIEVPANFILOV_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_MONODOMAINALIEVPANFILOV_H

#include <ADRSolver/EquationSystems/Monodomain.h>

namespace Nektar
{
    /// A two-variable model for cardiac conduction.
    class MonodomainAlievPanfilov : public Monodomain
    {
    public:
        /// Creates an instance of this class
        static EquationSystemSharedPtr create(SessionReaderSharedPtr& pSession)
        {
            return MemoryManager<MonodomainAlievPanfilov>::AllocateSharedPtr(pSession);
        }

        /// Name of class
        static std::string className;

        /// Constructor
        MonodomainAlievPanfilov(SessionReaderSharedPtr& pSession);

        /// Desctructor
        virtual ~MonodomainAlievPanfilov();

    protected:
        /// Computes the reaction terms \f$f(u,v)\f$ and \f$g(u,v)\f$.
        void DoOdeRhs(
                const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                      Array<OneD,        Array<OneD, NekDouble> >&outarray,
                const NekDouble time);

        /// Prints a summary of the model parameters.
        virtual void v_PrintSummary(std::ostream &out);

    private:
        /// Trigger parameter a.
        NekDouble m_a;
        /// Scaling parameter k.
        NekDouble m_k;
        /// Restitution parameter \f$\mu_1\f$.
        NekDouble m_mu1;
        /// Restitution parameter \f$\mu_2\f$.
        NekDouble m_mu2;
        /// Restitution parameter \f$\epsilon\f$.
        NekDouble m_eps;

        /// Temporary space for storing \f$u^2\f$ when computing reaction term.
        Array<OneD, NekDouble> m_uu;
        /// Temporary space for storing \f$u^3\f$ when computing reaction term.
        Array<OneD, NekDouble> m_uuu;
        /// Workspace for computing reaction term.
        Array<OneD, NekDouble> m_tmp1;
        /// Workspace for computing reaction term.
        Array<OneD, NekDouble> m_tmp2;
    };
}

#endif
