////////////////////////////////////////////////////////////////////////////////
//
// File: CouplingFile.cpp
//
// For more information, please see: http://www.nektar.info/
//
// The MIT License
//
// Copyright (c) 2017 Kilian Lackhove
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
// Description: CWIPI Exchange class
//
////////////////////////////////////////////////////////////////////////////////

#include "CouplingFile.h"

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/BasicUtils/PtsField.h>
#include <LibUtilities/BasicUtils/PtsIO.h>

#include <boost/format.hpp>

namespace Nektar
{
namespace SolverUtils
{

using namespace std;

std::string CouplingFile::className =
    GetCouplingFactory().RegisterCreatorFunction(
        "File", CouplingFile::create, "File Coupling");

CouplingFile::CouplingFile(MultiRegions::ExpListSharedPtr field)
    : Coupling(field), m_lastSend(-1E6), m_lastReceive(-1E6)
{
    m_config["RECEIVEFUNCTION"] = "CouplingIn";
    m_config["SENDFILENAME"]    = "CouplingOut_%14.8E.pts";
}

CouplingFile::~CouplingFile()
{
}

void CouplingFile::v_Init()
{
    Coupling::v_Init();

    if (m_nRecvVars > 0 && m_recvSteps > 0)
    {
        m_inputFunction = MemoryManager<SessionFunction>::AllocateSharedPtr(
            m_evalField->GetSession(),
            m_evalField,
            m_config["RECEIVEFUNCTION"],
            true);
    }
}

void CouplingFile::v_Send(
    const int step,
    const NekDouble time,
    const Array<OneD, const Array<OneD, NekDouble> > &field,
    vector<string> &varNames)
{
    if (m_nSendVars < 1 || m_sendSteps < 1)
    {
        return;
    }

    if (step < m_lastSend + m_sendSteps)
    {
        return;
    }
    m_lastSend = step;

    if (m_evalField->GetComm()->GetRank() == 0 &&
        m_evalField->GetSession()->DefinesCmdLineArgument("verbose"))
    {
        cout << "sending fields at i = " << step << ", t = " << time << endl;
    }

    vector<int> sendVarsToVars =
        GenerateVariableMapping(varNames, m_sendFieldNames);

#if (defined _WIN32 && _MSC_VER < 1900)
    // We need this to make sure boost::format has always
    // two digits in the exponents of Scientific notation.
    unsigned int old_exponent_format;
    old_exponent_format = _set_output_format(_TWO_DIGIT_EXPONENT);
    std::string filename =
        boost::str(boost::format(m_config["SENDFILENAME"]) % time);
    _set_output_format(old_exponent_format);
#else
    std::string filename =
        boost::str(boost::format(m_config["SENDFILENAME"]) % time);
#endif

    Array<OneD, Array<OneD, NekDouble> > pts(m_nSendVars + 3);
    for (int i = 0; i < 3; ++i)
    {
        pts[i] = Array<OneD, NekDouble>(m_evalField->GetTotPoints(), 0.0);
    }
    m_evalField->GetCoords(pts[0], pts[1], pts[2]);

    for (int i = 0; i < m_nSendVars; ++i)
    {
        pts[3 + i] = field[sendVarsToVars[i]];
    }

    LibUtilities::PtsIO ptsIO(m_evalField->GetSession()->GetComm());
    LibUtilities::PtsFieldSharedPtr sPts =
        MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(
            3, m_sendFieldNames, pts);
    // we first write to a temp file and rename this to make sure the
    // receiver doesnt try to read it before we finished writing
    string tmpFn = filename + ".tmp";
    ptsIO.Write(tmpFn, sPts);
    fs::rename(tmpFn, filename);
}

void CouplingFile::v_Receive(const int step,
                             const NekDouble time,
                             Array<OneD, Array<OneD, NekDouble> > &field,
                             vector<string> &varNames)
{
    if (m_nRecvVars < 1 || m_recvSteps < 1)
    {
        return;
    }

    if (step < m_lastReceive + m_recvSteps)
    {
        return;
    }
    m_lastReceive = step;

    if (m_evalField->GetComm()->GetRank() == 0 &&
        m_evalField->GetSession()->DefinesCmdLineArgument("verbose"))
    {
        cout << "receiving fields at i = " << step << ", t = " << time << endl;
    }

    string filename = m_evalField->GetSession()->GetFunctionFilename(
        m_config["RECEIVEFUNCTION"], m_recvFieldNames[0]);

#if (defined _WIN32 && _MSC_VER < 1900)
    // We need this to make sure boost::format has always
    // two digits in the exponents of Scientific notation.
    unsigned int old_exponent_format;
    old_exponent_format = _set_output_format(_TWO_DIGIT_EXPONENT);
    filename            = boost::str(boost::format(filename) % time);
    _set_output_format(old_exponent_format);
#else
    filename = boost::str(boost::format(filename) % time);
#endif

    int exists = 0;
    while (!exists)
    {
        exists = fs::exists(filename);
        m_evalField->GetComm()->AllReduce(exists, LibUtilities::ReduceMin);
    }

    Array<OneD, Array<OneD, NekDouble> > recvFields(m_nRecvVars);
    m_inputFunction->Evaluate(m_recvFieldNames, recvFields, time);

    vector<int> recvVarsToVars =
        GenerateVariableMapping(varNames, m_recvFieldNames);
    ASSERTL1(m_nRecvVars == recvVarsToVars.size(), "field size mismatch");
    for (int i = 0; i < recvVarsToVars.size(); ++i)
    {
        Vmath::Vcopy(recvFields[i].size(),
                     recvFields[i],
                     1,
                     field[recvVarsToVars[i]],
                     1);
    }
}
}
}
