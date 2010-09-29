#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_MONODOMAINFITZHUGHNAGUMO_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_MONODOMAINFITZHUGHNAGUMO_H

#include <ADRSolver/EquationSystems/Monodomain.h>

namespace Nektar
{
    /// A two-variable model for cardiac conduction.
    class MonodomainFitzHughNagumo : public Monodomain
    {
    public:
        /// Creates an instance of this class
        static EquationSystemSharedPtr create(SessionReaderSharedPtr& pSession)
        {
            return MemoryManager<MonodomainFitzHughNagumo>::AllocateSharedPtr(pSession);
        }

        /// Name of class
        static std::string className;

        /// Constructor
        MonodomainFitzHughNagumo(SessionReaderSharedPtr& pSession);

        /// Desctructor
        virtual ~MonodomainFitzHughNagumo();

    protected:
        /// Computes the reaction terms \f$f(u,v)\f$ and \f$g(u,v)\f$.
        void DoOdeRhs(
                const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                      Array<OneD,        Array<OneD, NekDouble> >&outarray,
                const NekDouble time);

        /// Prints a summary of the model parameters.
        virtual void v_PrintSummary(std::ostream &out);

    private:
        NekDouble              m_beta;
        NekDouble              m_kr;
        /// Temporary space for storing \f$u^3\f$ when computing reaction term.
        Array<OneD, NekDouble> m_uuu;
        /// Workspace for computing reaction term.
        Array<OneD, NekDouble> m_tmp1;
    };
}

#endif
