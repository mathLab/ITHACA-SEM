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

    m_ode.DefineOdeRhs(&Dummy::DoOdeRhs, this);
    m_ode.DefineProjection(&Dummy::DoOdeProjection, this);

    ASSERTL0(m_session->DefinesCmdLineArgument("cwipi"),
             "This EquationSystem requires the --cwipi command line switch");

    m_coupling = MemoryManager<CwipiCoupling>::AllocateSharedPtr(
        m_fields[0], 0, 1.0);

    m_recFields = Array<OneD, Array<OneD, NekDouble> >(
        m_coupling->GetRecvFieldNames().size());
    for (int i = 0; i < m_recFields.num_elements(); ++i)
    {
        m_recFields[i] = Array<OneD, NekDouble>(GetTotPoints(), 0.0);
    }

    m_sendFields = Array<OneD, Array<OneD, NekDouble> >(
        m_coupling->GetSendFieldNames().size());
    for (int i = 0; i < m_sendFields.num_elements(); ++i)
    {
        m_sendFields[i] = Array<OneD, NekDouble>(GetTotPoints(), 0.0);
    }
}

/**
 * @brief Destructor for Dummy class.
 */
Dummy::~Dummy()
{
}

/**
 * @brief v_PreIntegrate
 */
bool Dummy::v_PreIntegrate(int step)
{
    if (m_sendFields.num_elements() > 0)
    {
        Timer timer1;
        timer1.Start();
        GetFunction("SendFields", m_fields[0])->Evaluate(m_coupling->GetSendFieldNames(), m_sendFields, m_time);
        timer1.Stop();
        if (m_session->DefinesCmdLineArgument("verbose"))
        {
            cout << "Field evaluation time: " << timer1.TimePerTest(1) << endl;
        }
    }

    m_coupling->Send(step, m_time, m_sendFields);

    m_coupling->ReceiveInterp(step, m_time, m_recFields);

    for (int i = 0; i < m_recFields.num_elements(); ++i)
    {
        NekDouble intVal = m_fields[0]->PhysIntegral(m_recFields[i]);
        m_comm->AllReduce(intVal, LibUtilities::ReduceSum);
        cout << "Integral of received field " << i;
        cout << " = " << intVal << endl;
    }

    return UnsteadySystem::v_PreIntegrate(step);
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

void Dummy::v_Output(void)
{
    Nektar::SolverUtils::EquationSystem::v_Output();

    m_coupling->FinalizeCoupling();
}

void Dummy::v_ExtraFldOutput(std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
                             std::vector<std::string> &variables)
{
    for (int i = 0; i < m_recFields.num_elements(); i++)
    {
        Array<OneD, NekDouble> tmpC(GetNcoeffs());

        m_fields[0]->FwdTrans(m_recFields[i], tmpC);
        variables.push_back(m_coupling->GetRecvFieldNames()[i]);
        fieldcoeffs.push_back(tmpC);
    }

    for (int i = 0; i < m_sendFields.num_elements(); i++)
    {
        Array<OneD, NekDouble> tmpC(GetNcoeffs());

        m_fields[0]->FwdTrans(m_sendFields[i], tmpC);
        variables.push_back(m_coupling->GetSendFieldNames()[i]);
        fieldcoeffs.push_back(tmpC);
    }
}

} // end of namespace
