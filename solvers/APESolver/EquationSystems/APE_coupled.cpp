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

#include <SolverUtils/CwipiExchange.h>

namespace Nektar
{
string APE_coupled::className = GetEquationSystemFactory().RegisterCreatorFunction(
                                    "APE_coupled", APE_coupled::create,
                                    "Coupled APE1/APE4 (Acoustic Perturbation Equations)");


APE_coupled::APE_coupled(
    const LibUtilities::SessionReaderSharedPtr &pSession) : APE(pSession)
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

    m_stOld = Array<OneD, NekDouble>(GetTotPoints());
    m_stNew = Array<OneD, NekDouble>(GetTotPoints());
    Vmath::Vcopy(GetTotPoints(), m_st[0]->GetPhys(), 1, m_stOld, 1);
    Vmath::Vcopy(GetTotPoints(), m_st[0]->GetPhys(), 1, m_stNew, 1);

    m_bfOld = Array<OneD, Array<OneD, NekDouble> >(m_spacedim + 2);
    m_bfNew = Array<OneD, Array<OneD, NekDouble> >(m_spacedim + 2);
    for (int i = 0; i < m_spacedim + 2; ++i)
    {
        m_bfOld[i] = Array<OneD, NekDouble>(GetTotPoints());
        m_bfNew[i] = Array<OneD, NekDouble>(GetTotPoints());
        Vmath::Vcopy(GetTotPoints(), m_bf[i]->GetPhys(), 1, m_bfOld[i], 1);
        Vmath::Vcopy(GetTotPoints(), m_bf[i]->GetPhys(), 1, m_bfNew[i], 1);
    }

    m_session->LoadParameter("EX_RecvSteps", m_recvSteps, 1);
    m_session->LoadParameter("EX_SendSteps", m_sendSteps, 1);
    NekDouble filtWidth;
    m_session->LoadParameter("EX_FiltWidth", filtWidth, 0.0);

    m_nRecvVars = 6;

    m_coupling = MemoryManager<CwipiCoupling>::AllocateSharedPtr(
                        m_bf[0], "cpl1", "precise", 0, 1.0, filtWidth);
    m_sendExchange = MemoryManager<CwipiExchange>::AllocateSharedPtr(
                            m_coupling, "ex1", m_nRecvVars);

}


/**
 * @brief Destructor for APE class.
 */
APE_coupled::~APE_coupled()
{

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

    receiveFields();

    // ensure the new fields are C0-continuous
    // HACK: normally, we would perform a FwdTrans and BwdTrans here, but for
    // hexas and parrallel this is currently broken. The following has the same
    // effect and might even be faster
    m_st[0]->IProductWRTBase(m_st[0]->GetPhys(), m_st[0]->UpdateCoeffs());
    m_st[0]->MultiplyByElmtInvMass(m_st[0]->GetCoeffs(), m_st[0]->UpdateCoeffs());
    m_st[0]->LocalToGlobal();
    m_st[0]->GlobalToLocal();
    m_st[0]->BwdTrans(m_st[0]->GetCoeffs(), m_st[0]->UpdatePhys());

    for (int i = 0; i < m_spacedim + 2; ++i)
    {
        m_bf[i]->IProductWRTBase(m_bf[i]->GetPhys(), m_bf[i]->UpdateCoeffs());
        m_bf[i]->MultiplyByElmtInvMass(m_bf[i]->GetCoeffs(), m_bf[i]->UpdateCoeffs());
        m_bf[i]->LocalToGlobal();
        m_bf[i]->GlobalToLocal();
        m_bf[i]->BwdTrans(m_bf[i]->GetCoeffs(), m_bf[i]->UpdatePhys());
    }

    return UnsteadySystem::v_PostIntegrate(step);
}


void APE_coupled::receiveFields()
{
    static NekDouble last_update = -1E23;
    int nq = GetTotPoints();

    if (m_time >= last_update + m_recvSteps * m_timestep)
    {
        last_update = m_time;

        Vmath::Vcopy(nq, m_stNew, 1, m_stOld, 1);
        for (int i = 0; i < m_spacedim + 2; ++i)
        {
            Vmath::Vcopy(nq, m_bfNew[i], 1, m_bfOld[i], 1);
        }

        Array<OneD, Array<OneD, NekDouble> > recField(m_nRecvVars);

        recField[3] = m_bfNew[0]; // p0
        recField[4] = m_bfNew[1]; // rho0
        recField[0] = m_bfNew[2]; // u0
        if (m_spacedim > 1)
        {
            recField[1] = m_bfNew[3]; // v0
        }
        else
        {
            recField[1] = Array<OneD, NekDouble>(nq);
        }
        if (m_spacedim > 2)
        {
            recField[2] = m_bfNew[4]; // w0
        }
        else
        {
            recField[2] = Array<OneD, NekDouble>(nq);
        }
        recField[5] = m_stNew; // S

        m_sendExchange->ReceiveFields(0, m_time, recField);


        ASSERTL0(Vmath::Vmin(nq, m_bfNew[0], 1) > 0.0, "received p0 <= 0.0");
        ASSERTL0(Vmath::Vmin(nq, m_bfNew[1], 1) > 0.0, "received rho0 <= 0.0");
    }

    // linear interpolation in time. We cant do this iteratively because m_sourceTerms
    // and m_bf will be changed in v_PostIntegrate()
    NekDouble fact = (m_time - last_update + m_timestep) / (m_recvSteps * m_timestep);
    Vmath::Svtsvtp(GetTotPoints(), fact, m_stNew, 1, (1 - fact), m_stOld, 1, m_st[0]->UpdatePhys(), 1);
    for (int i = 0; i < m_spacedim + 2; ++i)
    {
        Vmath::Svtsvtp(GetTotPoints(), fact, m_bfNew[i], 1, (1 - fact), m_bfOld[i], 1, m_bf[i]->UpdatePhys(), 1);
    }
}


void APE_coupled::v_Output(void)
{
    Nektar::SolverUtils::EquationSystem::v_Output();

    m_coupling->FinalizeCoupling();
}

} //end of namespace

