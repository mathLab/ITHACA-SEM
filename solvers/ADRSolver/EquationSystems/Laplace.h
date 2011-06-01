#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_LAPLACE_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_LAPLACE_H

#include <ADRSolver/EquationSystem.h>

namespace Nektar
{
    class Laplace : public EquationSystem
    {
    public:
        /// Creates an instance of this class
        static EquationSystemSharedPtr create(
                LibUtilities::CommSharedPtr& pComm,
                LibUtilities::SessionReaderSharedPtr& pSession)
        {
            return MemoryManager<Laplace>::AllocateSharedPtr(pComm, pSession);
        }
        
        /// Name of class
        static std::string className;

        Laplace(LibUtilities::CommSharedPtr& pComm,
                LibUtilities::SessionReaderSharedPtr& pSession);
        virtual ~Laplace();

    protected:
        NekDouble m_lambda;

        virtual void v_PrintSummary(std::ostream &out);
        virtual void v_DoSolve();
    };
}

#endif
