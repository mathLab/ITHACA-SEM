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

    m_ode.DefineOdeRhs(&Dummy::DoOdeRhs, this);
    m_ode.DefineProjection(&Dummy::DoOdeProjection, this);

    ASSERTL0(m_session->DefinesCmdLineArgument("cwipi"),
             "This EquationSystem requires the --cwipi command line switch");

    m_session->LoadParameter("EX_RecvSteps", m_recvSteps, 1);
    m_session->LoadParameter("EX_SendSteps", m_sendSteps, 1);

    m_recvFieldNames.push_back("u0");
    m_recvFieldNames.push_back("v0");
    m_recvFieldNames.push_back("w0");
    m_recvFieldNames.resize(m_spacedim);

    m_recvFieldNames.push_back("p0");
    m_recvFieldNames.push_back("rho0");
    m_recvFieldNames.push_back("S");

    // HACK
    m_nRecvVars = 6; //m_recvFieldNames.size();

    m_recFields = Array<OneD, Array<OneD, NekDouble> >(m_nRecvVars);
    for (int i = 0; i < m_recFields.num_elements(); ++i)
    {
        m_recFields[i] = Array<OneD, NekDouble>(GetTotPoints());
    }

    m_coupling = MemoryManager<CwipiCoupling>::AllocateSharedPtr(
        m_fields[0], "cpl1",  0, 1.0);
    m_sendExchange = MemoryManager<CwipiExchange>::AllocateSharedPtr(
        m_coupling, "ex1", m_nRecvVars);
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
    receiveFields();

    return UnsteadySystem::v_PostIntegrate(step);
}

/**
 * @brief Compute the right-hand side.
 */
void Dummy::DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                     Array<OneD, Array<OneD, NekDouble> > &outarray,
                     const NekDouble time)
{
    // do nothing
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
    // do nothing
}

void Dummy::receiveFields()
{
    static NekDouble last_update = -1E23;

    if (m_time >= last_update + m_recvSteps * m_timestep)
    {
        last_update = m_time;

        m_sendExchange->ReceiveFields(0, m_time, m_recFields);

        DumpFields("recFields_" + boost::lexical_cast<std::string>(m_time) +
                       ".pts",
                   m_fields[0],
                   m_recFields);

        Array<OneD, MultiRegions::ExpListSharedPtr> recvFields =
            m_coupling->GetRecvFields();
        Array<OneD, Array<OneD, NekDouble> > tmp(m_nRecvVars);
        for (int i = 0; i < m_nRecvVars; ++i)
        {
            tmp[i] = recvFields[i]->GetPhys();
        }
        DumpFields("rawFields_" + boost::lexical_cast<std::string>(m_time) +
                       ".pts",
                   recvFields[0],
                   tmp);
    }
}

void Dummy::v_Output(void)
{
    Nektar::SolverUtils::EquationSystem::v_Output();

    m_coupling->FinalizeCoupling();
}

void Dummy::DumpFields(const string filename,
                       MultiRegions::ExpListSharedPtr field,
                       Array<OneD, Array<OneD, NekDouble> > data)
{
    ASSERTL0(data.num_elements() == m_nRecvVars, "dimension mismatch");
    ASSERTL0(data[0].num_elements() == field->GetTotPoints(),
             "dimension mismatch");

    int nq = field->GetTotPoints();

    Array<OneD, Array<OneD, NekDouble> > tmp(m_nRecvVars + m_spacedim);

    for (int i = 0; i < m_spacedim; ++i)
    {
        tmp[i] = Array<OneD, NekDouble>(nq, 0.0);
    }
    field->GetCoords(tmp[0], tmp[1], tmp[2]);

    for (int i = 0; i < m_nRecvVars; ++i)
    {
        tmp[m_spacedim + i] = data[i];
    }

    LibUtilities::PtsIO ptsIO(m_session->GetComm());
    LibUtilities::PtsFieldSharedPtr rvPts =
        MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(
            m_spacedim, m_recvFieldNames, tmp);
    ptsIO.Write(filename, rvPts);
}

void Dummy::v_ExtraFldOutput(std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
                             std::vector<std::string> &variables)
{
    for (int i = 0; i < m_spacedim + 2; i++)
    {
        Array<OneD, NekDouble> tmpC(GetNcoeffs());

        m_fields[0]->FwdTrans(m_recFields[i], tmpC);
        variables.push_back(m_recvFieldNames[i]);
        fieldcoeffs.push_back(tmpC);
    }
}

} // end of namespace
