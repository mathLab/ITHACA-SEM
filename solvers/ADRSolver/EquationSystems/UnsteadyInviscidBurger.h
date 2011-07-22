#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_UNSTEADYINVISCIDBURGER_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_UNSTEADYINVISCIDBURGER_H

#include <Auxiliary/UnsteadySystem.h>

namespace Nektar
{
    class UnsteadyInviscidBurger : public UnsteadySystem
    {
    public:
        friend class MemoryManager<UnsteadyInviscidBurger>;

        /// Creates an instance of this class
        static EquationSystemSharedPtr create(
                LibUtilities::CommSharedPtr& pComm,
                LibUtilities::SessionReaderSharedPtr& pSession) {
            EquationSystemSharedPtr p = MemoryManager<UnsteadyInviscidBurger>::AllocateSharedPtr(pComm, pSession);
            p->InitObject();
            return p;
        }
        /// Name of class
        static std::string className;

        virtual ~UnsteadyInviscidBurger();

    protected:
        UnsteadyInviscidBurger(
                LibUtilities::CommSharedPtr& pComm,
                LibUtilities::SessionReaderSharedPtr& pSession);

        void DoOdeRhs(
                const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                      Array<OneD,        Array<OneD, NekDouble> >&outarray,
                const NekDouble time);

        void DoOdeProjection(const Array<OneD,  const  Array<OneD, NekDouble> > &inarray,
                          Array<OneD,  Array<OneD, NekDouble> > &outarray,
                          const NekDouble time);

        virtual void v_InitObject();

        virtual void v_GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> > &physfield, Array<OneD, Array<OneD, NekDouble> > &flux);
        virtual void v_GetFluxVector(const int i, const int j, Array<OneD, Array<OneD, NekDouble> > &physfield, Array<OneD, Array<OneD, NekDouble> > &flux);

        virtual void v_NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, Array<OneD, Array<OneD, NekDouble> > &numflux);

    private:
        NekDouble m_waveFreq;
        Array<OneD, Array<OneD, NekDouble> > m_velocity;
    };
}

#endif
