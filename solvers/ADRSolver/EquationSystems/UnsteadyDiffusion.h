#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_UNSTEADYDIFFUSION_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_UNSTEADYDIFFUSION_H

#include <Auxiliary/UnsteadySystem.h>

namespace Nektar
{
    class UnsteadyDiffusion : public UnsteadySystem
    {
    public:
        friend class MemoryManager<UnsteadyDiffusion>;

        /// Creates an instance of this class
        static EquationSystemSharedPtr create(
                LibUtilities::CommSharedPtr& pComm,
                LibUtilities::SessionReaderSharedPtr& pSession) {
            EquationSystemSharedPtr p = MemoryManager<UnsteadyDiffusion>::AllocateSharedPtr(pComm, pSession);
            p->InitObject();
            return p;
        }
        /// Name of class
        static std::string className;

        virtual ~UnsteadyDiffusion();

    protected:
        UnsteadyDiffusion(
                LibUtilities::CommSharedPtr& pComm,
                LibUtilities::SessionReaderSharedPtr& pSession);

        void DoOdeRhs(const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                          Array<OneD,        Array<OneD, NekDouble> >&outarray,
                    const NekDouble time);


        void DoOdeProjection(const Array<OneD,  const  Array<OneD, NekDouble> > &inarray,
                          Array<OneD,  Array<OneD, NekDouble> > &outarray,
                          const NekDouble time);

        virtual void DoImplicitSolve(const Array<OneD, const Array<OneD,      NekDouble> >&inarray,
                      Array<OneD, Array<OneD, NekDouble> >&outarray,
                      NekDouble time,
                      NekDouble lambda);

        virtual void v_InitObject();

    private:
        NekDouble m_waveFreq;
        NekDouble m_epsilon;
    };
}

#endif
