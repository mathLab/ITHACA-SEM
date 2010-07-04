#include <ADRSolver/EquationSystems/Laplace.h>

namespace Nektar
{
    string Laplace::className = EquationSystemFactory::RegisterCreatorFunction("Laplace", Laplace::create);

    Laplace::Laplace(SessionReaderSharedPtr& pSession,
            LibUtilities::TimeIntegrationSchemeOperators& pOde)
        : EquationSystem(pSession, pOde),
          mLambda(0)
    {
        if (pSession->definesParameter("Lambda"))
        {
            mLambda = pSession->getParameter("Lambda");
        }
    }

    Laplace::~Laplace()
    {

    }

    void Laplace::v_printSummary(std::ostream &out)
    {
        out << "\tLambda          : " << mLambda << endl;
    }

    void Laplace::v_doSolveHelmholtz()
    {
        for(int i = 0; i < m_fields.num_elements(); ++i)
        {
            m_fields[i]->HelmSolve(m_fields[i]->GetPhys(),
                                   m_fields[i]->UpdateCoeffs(),
                                   mLambda);
            m_fields[i]->SetPhysState(false);
        }
    }

    bool Laplace::v_isSteady()
    {
        return true;
    }
}
