///////////////////////////////////////////////////////////////////////////////
//
// File APE_coupled.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2016 Kilian Lackhove
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
// Description: APE1/APE4 (Acoustic Perturbation Equations)
//
///////////////////////////////////////////////////////////////////////////////

#include <APESolver/EquationSystems/APE_coupled.h>

namespace Nektar
{

using namespace std;

string APE_coupled::className =
    GetEquationSystemFactory().RegisterCreatorFunction(
        "APE_coupled",
        APE_coupled::create,
        "Coupled APE1/APE4 (Acoustic Perturbation Equations)");

APE_coupled::APE_coupled(const LibUtilities::SessionReaderSharedPtr &pSession)
    : APE(pSession)
{
}

/**
 * @brief Initialization object for the APE class.
 */
void APE_coupled::v_InitObject()
{
    APE::v_InitObject();

    ASSERTL0(m_session->DefinesCmdLineArgument("cwipi"),
             "This EquationSystem requires the --cwipi command line switch");

    m_coupling = MemoryManager<CwipiCoupling>::AllocateSharedPtr(
        m_bfField, "cpl1", 0, 1.0);
}

/**
 * @brief Destructor for APE class.
 */
APE_coupled::~APE_coupled()
{
}

/**
 * @brief v_PreIntegrate
 */
bool APE_coupled::v_PreIntegrate(int step)
{
    receiveFields();

    Array<OneD, NekDouble> tmpC(GetNcoeffs());
    for (int i = 0; i < m_spacedim + 2; ++i)
    {
        // ensure the field is C0-continuous
        m_bfField->IProductWRTBase(m_bf[i], tmpC);
        m_bfField->MultiplyByElmtInvMass(tmpC, tmpC);
        m_bfField->LocalToGlobal(tmpC, tmpC);
        m_bfField->GlobalToLocal(tmpC, tmpC);
        m_bfField->BwdTrans(tmpC, m_bf[i]);
    }

    for (int i = 0; i < m_st.num_elements(); ++i)
    {
        // ensure the field is C0-continuous
        m_bfField->IProductWRTBase(m_st[i], tmpC);
        m_bfField->MultiplyByElmtInvMass(tmpC, tmpC);
        m_bfField->LocalToGlobal(tmpC, tmpC);
        m_bfField->GlobalToLocal(tmpC, tmpC);
        m_bfField->BwdTrans(tmpC, m_st[i]);
    }

    return UnsteadySystem::v_PreIntegrate(step);
}

/**
 * @brief v_PostIntegrate
 */
bool APE_coupled::v_PostIntegrate(int step)
{
    if (m_cflsteps && !((step + 1) % m_cflsteps))
    {
        NekDouble cfl = GetCFLEstimate();
        if (m_comm->GetRank() == 0)
        {
            cout << "CFL: " << cfl << endl;
        }
    }

    return UnsteadySystem::v_PostIntegrate(step);
}

void APE_coupled::receiveFields()
{
    int nq = GetTotPoints();

    Array<OneD, Array<OneD, NekDouble> > recField(6);

    recField[3] = m_bf[0]; // p0
    recField[4] = m_bf[1]; // rho0
    recField[0] = m_bf[2]; // u0
    if (m_spacedim > 1)
    {
        recField[1] = m_bf[3]; // v0
    }
    else
    {
        recField[1] = Array<OneD, NekDouble>(nq);
    }
    if (m_spacedim > 2)
    {
        recField[2] = m_bf[4]; // w0
    }
    else
    {
        recField[2] = Array<OneD, NekDouble>(nq);
    }
    recField[5] = m_st[0]; // S

    m_coupling->ReceiveFields(0, m_time, recField);

    // HACK
    for (int i = 0; i < nq; ++i)
    {
        if (m_bf[0][i] < 201164.11)
        {
            m_bf[0][i] = 201164.11;
        }
    }

    ASSERTL0(Vmath::Vmin(nq, m_bf[0], 1) > 0.0, "received p0 <= 0.0");
    ASSERTL0(Vmath::Vmin(nq, m_bf[1], 1) > 0.0, "received rho0 <= 0.0");
}

void APE_coupled::v_Output(void)
{
    Nektar::SolverUtils::EquationSystem::v_Output();

    m_coupling->FinalizeCoupling();
}

} // end of namespace
