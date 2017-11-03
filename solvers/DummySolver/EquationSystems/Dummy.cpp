///////////////////////////////////////////////////////////////////////////////
//
// File Dummy.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2016 Kilian Lackhove
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
// Description: Dummy Equation System that only sends/receives fields
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <DummySolver/EquationSystems/Dummy.h>

using namespace std;

namespace Nektar
{
string Dummy::className = GetEquationSystemFactory().RegisterCreatorFunction(
    "Dummy",
    Dummy::create,
    "Dummy Equation System that only sends/receives fields");

Dummy::Dummy(const LibUtilities::SessionReaderSharedPtr &pSession)
    : UnsteadySystem(pSession)
{
}

/**
 * @brief Initialization object for the Dummy class.
 */
void Dummy::v_InitObject()
{
    UnsteadySystem::v_InitObject();

    m_nanSteps = 0;

    auto sV = m_session->GetVariables();
    for (auto const &sendVar : m_coupling->GetSendFieldNames())
    {
        int i = distance(sV.begin(), find(sV.begin(), sV.end(), sendVar));
        m_intVariables.push_back(i);
    }

    m_ode.DefineOdeRhs(&Dummy::DoOdeRhs, this);
    m_ode.DefineProjection(&Dummy::DoOdeProjection, this);

    ASSERTL0(m_session->DefinesCmdLineArgument("cwipi"),
             "This EquationSystem requires the --cwipi command line switch");
}

/**
 * @brief Destructor for Dummy class.
 */
Dummy::~Dummy()
{
}

/**
 * @brief v_PostIntegrate
 */
bool Dummy::v_PostIntegrate(int step)
{
    if (m_coupling->GetSendFieldNames().size() > 0)
    {
        Timer timer1;
        timer1.Start();

        auto sV = m_session->GetVariables();
        for (auto const &sendVar : m_coupling->GetSendFieldNames())
        {
            int i = distance(sV.begin(), find(sV.begin(), sV.end(), sendVar));
            cout << "sendVar = " << sendVar << ", i = " << i << endl;
            GetFunction("SendFields", m_fields[i])
                ->Evaluate(sendVar, m_fields[i]->UpdatePhys(), m_time);
        }

        timer1.Stop();
        if (m_session->DefinesCmdLineArgument("verbose"))
        {
            cout << "Field evaluation time: " << timer1.TimePerTest(1) << endl;
        }
    }

    for (int i = 0; i < m_session->GetVariables().size(); ++i)
    {

        m_fields[i]->FwdTrans_IterPerExp(m_fields[i]->UpdatePhys(),
                                         m_fields[i]->UpdateCoeffs());
        m_fields[i]->SetPhysState(false);
    }

    return UnsteadySystem::v_PostIntegrate(step);
}

/**
 * @brief Compute the right-hand side.
 */
void Dummy::DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                     Array<OneD, Array<OneD, NekDouble> > &outarray,
                     const NekDouble time)
{
    int nVariables = inarray.num_elements();
    int nq         = GetTotPoints();

    for (int i = 0; i < nVariables; ++i)
    {
        Vmath::Zero(nq, outarray[i], 1);
    }
}

/**
 * @brief Compute the projection and call the method for imposing the
 * boundary conditions in case of discontinuous projection.
 */
void Dummy::DoOdeProjection(
    const Array<OneD, const Array<OneD, NekDouble> > &inarray,
    Array<OneD, Array<OneD, NekDouble> > &outarray,
    const NekDouble time)
{
    int nvariables = inarray.num_elements();
    int nq         = m_fields[0]->GetNpoints();

    // deep copy
    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Vcopy(nq, inarray[i], 1, outarray[i], 1);
    }
}

} // end of namespace
