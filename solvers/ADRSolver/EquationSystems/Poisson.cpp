#include <ADRSolver/EquationSystems/Poisson.h>

namespace Nektar
{
    string Poisson::className1 = EquationSystemFactory::RegisterCreatorFunction("Poisson", Poisson::create);
    string Poisson::className2 = EquationSystemFactory::RegisterCreatorFunction("SteadyDiffusion", Poisson::create);

    Poisson::Poisson(SessionReaderSharedPtr& pSession,
            LibUtilities::TimeIntegrationSchemeOperators& pOde)
        : Laplace(pSession, pOde)
    {
        SetPhysForcingFunctions(m_fields);
    }

    Poisson::~Poisson()
    {

    }

    void Poisson::v_printSummary(std::ostream &out)
    {
        Laplace::v_printSummary(out);
        for (int i = 0; i < m_fields.num_elements(); ++i)
        {
            out << "\tForcing func " << i << "  : " << m_boundaryConditions->GetForcingFunction(i)->GetEquation() << endl;
        }
    }
}
