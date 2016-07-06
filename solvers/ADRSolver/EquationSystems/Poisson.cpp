///////////////////////////////////////////////////////////////////////////////
//
// File Poisson.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Poisson solve routines 
//
///////////////////////////////////////////////////////////////////////////////

#include <ADRSolver/EquationSystems/Poisson.h>

using namespace std;

namespace Nektar
{
    string Poisson::className1 = GetEquationSystemFactory().RegisterCreatorFunction("Poisson", Poisson::create);
    string Poisson::className2 = GetEquationSystemFactory().RegisterCreatorFunction("SteadyDiffusion", Poisson::create);

    Poisson::Poisson(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        : Laplace(pSession)
    {
    }

    void Poisson::v_InitObject()
    {
        Laplace::v_InitObject();

        EvaluateFunction(m_session->GetVariables(), m_fields, "Forcing");
    }

    Poisson::~Poisson()
    {

    }

    void Poisson::v_GenerateSummary(SolverUtils::SummaryList& s)
    {
        Laplace::v_GenerateSummary(s);
        for (int i = 0; i < m_fields.num_elements(); ++i)
        {
            stringstream name;
            name << "Forcing func [" << i << "]";
            SolverUtils::AddSummaryItem(s, name.str(),
                    m_session->GetFunction("Forcing", i)->GetExpression());
        }
    }

    Array<OneD, bool> Poisson::v_GetSystemSingularChecks()
    {
        return Array<OneD, bool>(m_session->GetVariables().size(), true);
    }
}
