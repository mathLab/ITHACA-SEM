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
            "Dummy", Dummy::create,
            "Dummy Equation System that only sends/receives fields");


Dummy::Dummy(
        const LibUtilities::SessionReaderSharedPtr& pSession)
    : UnsteadySystem(pSession)
{
}


/**
 * @brief Initialization object for the Dummy class.
 */
void Dummy::v_InitObject()
{
    UnsteadySystem::v_InitObject();

    ASSERTL0(m_session->DefinesCmdLineArgument("cwipi"),
             "This EquationSystem requires the --cwipi command line switch");

    m_session->LoadParameter("EX_RecvSteps", m_recvSteps, 1);
    m_session->LoadParameter("EX_SendSteps", m_sendSteps, 1);

    m_nRecvVars = 6;

    m_coupling = MemoryManager<CwipiCoupling>::AllocateSharedPtr(
                        m_fields[0], "cpl1", "precise", 0, 1.0, 0.0);
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
void Dummy::DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                         Array<OneD,       Array<OneD, NekDouble> >&outarray,
                   const NekDouble time)
{
    // do nothing
}


/**
 * @brief Compute the projection and call the method for imposing the
 * boundary conditions in case of discontinuous projection.
 */
void Dummy::DoOdeProjection(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                                Array<OneD,       Array<OneD, NekDouble> >&outarray,
                          const NekDouble time)
{
    // do nothing
}


void Dummy::receiveFields()
{
    static NekDouble last_update = -1E23;
    int nq = GetTotPoints();

    if (m_time >= last_update + m_recvSteps * m_timestep)
    {
        last_update = m_time;

        Array<OneD, Array<OneD, NekDouble> > recField(m_nRecvVars);
        for (int i = 0; i < recField.num_elements(); ++i)
        {
            recField[i] = Array<OneD, NekDouble>(nq);
        }

        m_sendExchange->ReceiveFields(0, m_time, recField);
    }
}


void Dummy::v_Output(void)
{
    Nektar::SolverUtils::EquationSystem::v_Output();

    m_coupling->FinalizeCoupling();
}


} //end of namespace

