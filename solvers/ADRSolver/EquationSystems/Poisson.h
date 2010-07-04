#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_POISSON_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_POISSON_H

#include <ADRSolver/EquationSystems/Laplace.h>

namespace Nektar
{
    class Poisson : public Laplace
    {
    public:
        /// Creates an instance of this class
        static EquationSystemSharedPtr create(SessionReaderSharedPtr& pSession,
                LibUtilities::TimeIntegrationSchemeOperators& pOde) {
            return MemoryManager<Poisson>::AllocateSharedPtr(pSession, pOde);
        }
        /// Name of class
        static std::string className1;
        static std::string className2;

        Poisson(SessionReaderSharedPtr& pSession,
                LibUtilities::TimeIntegrationSchemeOperators& pOde);
        virtual ~Poisson();

    protected:
        virtual void v_printSummary(std::ostream &out);
    };
}

#endif
