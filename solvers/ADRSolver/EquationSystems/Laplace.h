#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_LAPLACE_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_LAPLACE_H

#include <ADRSolver/EquationSystem.h>

namespace Nektar
{
    class Laplace : public EquationSystem
    {
    public:
        /// Creates an instance of this class
        static EquationSystemSharedPtr create(SessionReaderSharedPtr& pSession,
                LibUtilities::TimeIntegrationSchemeOperators& pOde) {
            return MemoryManager<Laplace>::AllocateSharedPtr(pSession, pOde);
        }
        /// Name of class
        static std::string className;

        Laplace(SessionReaderSharedPtr& pSession,
                LibUtilities::TimeIntegrationSchemeOperators& pOde);
        virtual ~Laplace();

    protected:
        NekDouble mLambda;

        virtual void v_printSummary(std::ostream &out);
        virtual void v_doSolveHelmholtz();
        virtual bool v_isSteady();
    };
}

#endif
