#include <ADRSolver/EquationSystems/Helmholtz.h>

namespace Nektar
{
    string Helmholtz::className1 = EquationSystemFactory::RegisterCreatorFunction("Helmholtz", Helmholtz::create);
    string Helmholtz::className2 = EquationSystemFactory::RegisterCreatorFunction("SteadyDiffusionReaction", Helmholtz::create);

    Helmholtz::Helmholtz(SessionReaderSharedPtr& pSession,
            LibUtilities::TimeIntegrationSchemeOperators& pOde)
        : Poisson(pSession, pOde)
    {

    }

    Helmholtz::~Helmholtz()
    {

    }

    void Helmholtz::v_printSummary(std::ostream &out)
    {
        Poisson::v_printSummary(out);
    }
}
