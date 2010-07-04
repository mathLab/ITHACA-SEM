#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_UNSTEADYADVECTIONDIFFUSION_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_UNSTEADYADVECTIONDIFFUSION_H

#include <ADRSolver/EquationSystem.h>

namespace Nektar
{
    class UnsteadyAdvectionDiffusion : public EquationSystem
    {
    public:
        /// Creates an instance of this class
        static EquationSystemSharedPtr create(SessionReaderSharedPtr& pSession,
                LibUtilities::TimeIntegrationSchemeOperators& pOde) {
            return MemoryManager<UnsteadyAdvectionDiffusion>::AllocateSharedPtr(pSession, pOde);
        }
        /// Name of class
        static std::string className;

        UnsteadyAdvectionDiffusion(SessionReaderSharedPtr& pSession,
                LibUtilities::TimeIntegrationSchemeOperators& pOde);
        virtual ~UnsteadyAdvectionDiffusion();

    protected:
        virtual void v_doOdeRhs(const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                          Array<OneD,        Array<OneD, NekDouble> >&outarray,
                    const NekDouble time);

        virtual void v_GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> > &physfield, Array<OneD, Array<OneD, NekDouble> > &flux);
        virtual void v_GetFluxVector(const int i, const int j, Array<OneD, Array<OneD, NekDouble> > &physfield, Array<OneD, Array<OneD, NekDouble> > &flux);

        virtual void v_NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, Array<OneD, Array<OneD, NekDouble> > &numflux);

        virtual void v_NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield,
                                     Array<OneD, Array<OneD, NekDouble> > &numfluxX,
                                     Array<OneD, Array<OneD, NekDouble> > &numfluxY );


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
