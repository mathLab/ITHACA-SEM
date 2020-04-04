///////////////////////////////////////////////////////////////////////////////
//
// File Projection.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: Projection solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include <ADRSolver/EquationSystems/Projection.h>

using namespace std;

namespace Nektar
{
string Projection::className =
    GetEquationSystemFactory().RegisterCreatorFunction("Projection",
                                                       Projection::create);

Projection::Projection(const LibUtilities::SessionReaderSharedPtr &pSession,
                       const SpatialDomains::MeshGraphSharedPtr& pGraph)
    : EquationSystem(pSession, pGraph)
{
}

void Projection::v_InitObject()
{
    EquationSystem::v_InitObject();

    GetFunction("Forcing")->Evaluate(m_session->GetVariables(), m_fields);
}

Projection::~Projection()
{
}

void Projection::v_DoSolve()
{
    for (int i = 0; i < m_fields.size(); ++i)
    {
        // Zero field so initial conditions are zero
        Vmath::Zero(m_fields[i]->GetNcoeffs(), m_fields[i]->UpdateCoeffs(), 1);
        m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),
                              m_fields[i]->UpdateCoeffs());
        m_fields[i]->SetPhysState(false);
    }
}

void Projection::v_GenerateSummary(SolverUtils::SummaryList &s)
{
    EquationSystem::SessionSummary(s);
    for (int i = 0; i < m_fields.size(); ++i)
    {
        stringstream name;
        name << "Forcing func [" << i << "]";
        SolverUtils::AddSummaryItem(
            s, name.str(),
            m_session->GetFunction("Forcing", i)->GetExpression());
    }
}
}
