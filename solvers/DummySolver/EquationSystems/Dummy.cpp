///////////////////////////////////////////////////////////////////////////////
//
// File Dummy.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2017 Kilian Lackhove
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
// Description: Dummy Equation System that only sends/receives fields
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/BasicUtils/Timer.h>
#include <DummySolver/EquationSystems/Dummy.h>

using namespace std;

namespace Nektar
{
string Dummy::className = GetEquationSystemFactory().RegisterCreatorFunction(
    "Dummy",
    Dummy::create,
    "Dummy Equation System that only sends/receives fields");

Dummy::Dummy(const LibUtilities::SessionReaderSharedPtr &pSession,
             const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadySystem(pSession, pGraph)
{
}

/**
 * @brief Initialization object for the Dummy class.
 */
void Dummy::v_InitObject()
{
    UnsteadySystem::v_InitObject();

    m_nanSteps = 0;

    m_ode.DefineOdeRhs(&Dummy::DoOdeRhs, this);
    m_ode.DefineProjection(&Dummy::DoOdeProjection, this);

    m_forcing = SolverUtils::Forcing::Load(
        m_session, shared_from_this(), m_fields, m_fields.size());

    if (m_session->DefinesElement("Nektar/Coupling"))
    {
        TiXmlElement *vCoupling = m_session->GetElement("Nektar/Coupling");

        ASSERTL0(vCoupling->Attribute("TYPE"),
                 "Missing TYPE attribute in Coupling");
        string vType = vCoupling->Attribute("TYPE");
        ASSERTL0(!vType.empty(),
                 "TYPE attribute must be non-empty in Coupling");

        m_coupling = GetCouplingFactory().CreateInstance(vType, m_fields[0]);

        auto sV = m_session->GetVariables();
        for (auto const &sendVar : m_coupling->GetSendFieldNames())
        {
            auto match = find(sV.begin(), sV.end(), sendVar);
            if (match != sV.end())
            {
                int id = distance(sV.begin(), match);
                m_intVariables.push_back(id);
            }
        }
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
    if (m_coupling)
    {
        int numForceFields = 0;
        for (auto &x : m_forcing)
        {
            numForceFields += x->GetForces().size();
        }
        vector<string> varNames;
        Array<OneD, Array<OneD, NekDouble> > phys(m_fields.size() +
                                                  numForceFields);
        for (int i = 0; i < m_fields.size(); ++i)
        {
            varNames.push_back(m_session->GetVariable(i));
            phys[i] = m_fields[i]->UpdatePhys();
        }

        int f = 0;
        for (auto &x : m_forcing)
        {
            for (int i = 0; i < x->GetForces().size(); ++i)
            {
                phys[m_fields.size() + f + i] = x->GetForces()[i];
                varNames.push_back("F_" + boost::lexical_cast<string>(f) + "_" +
                                   m_session->GetVariable(i));
            }
            f++;
        }

        m_coupling->Send(step, m_time, phys, varNames);
        m_coupling->Receive(step, m_time, phys, varNames);
    }

    return UnsteadySystem::v_PreIntegrate(step);
}

/**
 * @brief v_PostIntegrate
 */
bool Dummy::v_PostIntegrate(int step)
{
    if (m_coupling && m_coupling->GetSendFieldNames().size() > 0)
    {
        LibUtilities::Timer timer1;
        timer1.Start();

        auto sV = m_session->GetVariables();
        for (auto const &sendVar : m_coupling->GetSendFieldNames())
        {
            auto match = find(sV.begin(), sV.end(), sendVar);
            if (match != sV.end())
            {
                int id = distance(sV.begin(), match);
                GetFunction("SendFields", m_fields[id])
                    ->Evaluate(sendVar, m_fields[id]->UpdatePhys(), m_time);
            }
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

void Dummy::v_Output()
{
    if (m_coupling)
    {
        m_coupling->Finalize();
    }

    UnsteadySystem::v_Output();

    int f = 0;
    for (auto &x : m_forcing)
    {
        for (int i = 0; i < x->GetForces().size(); ++i)
        {
            int npts = GetTotPoints();

            NekDouble l2err   = 0.0;
            NekDouble linferr = 0.0;
            for (int j = 0; j < npts; ++j)
            {
                l2err += x->GetForces()[i][j] * x->GetForces()[i][j];
                linferr = max(linferr, fabs(x->GetForces()[i][j]));
            }

            m_comm->AllReduce(l2err, LibUtilities::ReduceSum);
            m_comm->AllReduce(npts, LibUtilities::ReduceSum);
            m_comm->AllReduce(linferr, LibUtilities::ReduceMax);

            l2err /= npts;
            l2err = sqrt(l2err);

            if (m_comm->TreatAsRankZero())
            {
                cout << "L 2 error (variable "
                     << "F_" + boost::lexical_cast<string>(f) + "_" +
                            m_session->GetVariable(i)
                     << ") : " << l2err << endl;

                cout << "L inf error (variable "
                     << "F_" + boost::lexical_cast<string>(f) + "_" +
                            m_session->GetVariable(i)
                     << ") : " << linferr << endl;
            }
        }
        f++;
    }
}

/**
 * @brief Compute the right-hand side.
 */
void Dummy::DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                     Array<OneD, Array<OneD, NekDouble> > &outarray,
                     const NekDouble time)
{
    boost::ignore_unused(time);

    int nVariables = inarray.size();
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
    boost::ignore_unused(time);

    int nvariables = inarray.size();
    int nq         = m_fields[0]->GetNpoints();

    // deep copy
    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Vcopy(nq, inarray[i], 1, outarray[i], 1);
    }
}

} // namespace Nektar
