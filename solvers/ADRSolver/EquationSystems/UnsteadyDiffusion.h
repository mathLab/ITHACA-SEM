#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_UNSTEADYDIFFUSION_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_UNSTEADYDIFFUSION_H

#include <ADRSolver/EquationSystem.h>

namespace Nektar
{
    class UnsteadyDiffusion : public EquationSystem
    {
    public:
        /// Creates an instance of this class
        static EquationSystemSharedPtr create(SessionReaderSharedPtr& pSession,
                LibUtilities::TimeIntegrationSchemeOperators& pOde) {
            return MemoryManager<UnsteadyDiffusion>::AllocateSharedPtr(pSession, pOde);
        }
        /// Name of class
        static std::string className;

        UnsteadyDiffusion(SessionReaderSharedPtr& pSession,
                LibUtilities::TimeIntegrationSchemeOperators& pOde);
        virtual ~UnsteadyDiffusion();

    protected:
        virtual void v_doOdeRhs(const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                          Array<OneD,        Array<OneD, NekDouble> >&outarray,
                    const NekDouble time);

        virtual void v_GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> > &physfield, Array<OneD, Array<OneD, NekDouble> > &flux);
        virtual void v_GetFluxVector(const int i, const int j, Array<OneD, Array<OneD, NekDouble> > &physfield, Array<OneD, Array<OneD, NekDouble> > &flux);

    private:
        NekDouble mWaveFreq;
        Array<OneD, Array<OneD, NekDouble> > mVelocity;

        void doReaction(
                    const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                          Array<OneD, Array<OneD, NekDouble> >&outarray,
                    const NekDouble time);
    };
}

#endif
