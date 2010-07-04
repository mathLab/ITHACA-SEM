#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_ALIEVPANFILOV_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_ALIEVPANFILOV_H

#include <ADRSolver/EquationSystem.h>

namespace Nektar
{
    class AlievPanfilov : public EquationSystem
    {
    public:
        /// Creates an instance of this class
        static EquationSystemSharedPtr create(SessionReaderSharedPtr& pSession,
                LibUtilities::TimeIntegrationSchemeOperators& pOde) {
            return MemoryManager<AlievPanfilov>::AllocateSharedPtr(pSession, pOde);
        }
        /// Name of class
        static std::string className;

        AlievPanfilov(SessionReaderSharedPtr& pSession,
                LibUtilities::TimeIntegrationSchemeOperators& pOde);
        virtual ~AlievPanfilov();

    protected:
        virtual void v_doOdeRhs(const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                          Array<OneD,        Array<OneD, NekDouble> >&outarray,
                    const NekDouble time);

        virtual void v_doImplicitSolve(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
              Array<OneD, Array<OneD, NekDouble> >&outarray,
              NekDouble time,
              NekDouble lambda);

        virtual void v_printSummary(std::ostream &out);

    private:
        NekDouble mA;
        NekDouble mK;
        NekDouble mMu1;
        NekDouble mMu2;
        NekDouble mEps;
    };
}

#endif
