#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_MONODOMAIN_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_MONODOMAIN_H

#include <ADRSolver/EquationSystems/UnsteadySystem.h>

namespace Nektar
{
    /// A base model for cardiac conduction (no cell model).
    class Monodomain : public UnsteadySystem
    {
    public:
        /// Desctructor
        virtual ~Monodomain();

    protected:
        /// Constructor
        Monodomain(SessionReaderSharedPtr& pSession);

        /// Solve for the diffusion term.
        void DoImplicitSolve(
                const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                      Array<OneD, Array<OneD, NekDouble> >&outarray,
                      NekDouble time,
                      NekDouble lambda);

        /// Prints a summary of the model parameters.
        virtual void v_PrintSummary(std::ostream &out);

    private:

    };
}

#endif
