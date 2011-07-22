#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_LAPLACE_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_LAPLACE_H

#include <Auxiliary/EquationSystem.h>

namespace Nektar
{
    class Laplace : public EquationSystem
    {
    public:
        /// Class may only be instantiated through the MemoryManager.
        friend class MemoryManager<Laplace>;

        /// Creates an instance of this class
        static EquationSystemSharedPtr create(
                LibUtilities::CommSharedPtr& pComm,
                LibUtilities::SessionReaderSharedPtr& pSession)
        {
            EquationSystemSharedPtr p = MemoryManager<Laplace>::AllocateSharedPtr(pComm, pSession);
            p->InitObject();
            return p;
        }

        /// Name of class
        static std::string className;

    protected:
        NekDouble m_lambda;

        Laplace(LibUtilities::CommSharedPtr& pComm,
                LibUtilities::SessionReaderSharedPtr& pSession);

        virtual ~Laplace();

        virtual void v_InitObject();
        virtual void v_PrintSummary(std::ostream &out);
        virtual void v_DoSolve();

    private:
        virtual Array<OneD, bool> v_GetSystemSingularChecks();

    };
}

#endif
