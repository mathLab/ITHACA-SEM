#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_HELMHOLTZ_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_HELMHOLTZ_H

#include <ADRSolver/EquationSystems/Poisson.h>

namespace Nektar
{
    class Helmholtz : public Poisson
    {
    public:
        /// Creates an instance of this class
        static EquationSystemSharedPtr create(SessionReaderSharedPtr& pSession,
                LibUtilities::TimeIntegrationSchemeOperators& pOde) {
            return MemoryManager<Helmholtz>::AllocateSharedPtr(pSession, pOde);
        }
        /// Name of class
        static std::string className1;
        static std::string className2;

        Helmholtz(SessionReaderSharedPtr& pSession,
                LibUtilities::TimeIntegrationSchemeOperators& pOde);
        virtual ~Helmholtz();

    protected:
        virtual void v_printSummary(std::ostream &out);
    };
}

#endif
