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
// Description: CWIPI Exchange class
//
////////////////////////////////////////////////////////////////////////////////

#include "CouplingFile.h"

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/PtsField.h>
#include <LibUtilities/BasicUtils/PtsIO.h>

namespace Nektar
{
namespace SolverUtils
{

using namespace std;

std::string CouplingFile::className =
    GetCouplingFactory().RegisterCreatorFunction(
        "File", CouplingFile::create, "File Coupling");

CouplingFile::CouplingFile(MultiRegions::ExpListSharedPtr field)
    : Coupling(field)
{
    m_inputFunction = MemoryManager<SessionFunction>::AllocateSharedPtr(m_evalField->GetSession(), field, "CouplingIn", true);
}

CouplingFile::~CouplingFile()
{
}

void CouplingFile::v_Send(
    const int step,
    const NekDouble time,
    const Array<OneD, const Array<OneD, NekDouble> > &field,
    LibUtilities::FieldMetaDataMap &fieldMetaDataMap)
{
    int nSendVars = m_sendFieldNames.size();

    ASSERTL0(nSendVars == field.num_elements(), "size mismatch");

#ifdef _WIN32
    // We need this to make sure boost::format has always
    // two digits in the exponents of Scientific notation.
    unsigned int old_exponent_format;
    old_exponent_format = _set_output_format(_TWO_DIGIT_EXPONENT);
    filename = boost::str(boost::format("CouplingOut_%14.8E.pts") % m_time);
    _set_output_format(old_exponent_format);
#else
    std::string filename =
        boost::str(boost::format("CouplingOut_%14.8E.pts") % time);
#endif

    Array<OneD, Array<OneD, NekDouble> > tmp(nSendVars + 3);
    for (int i = 0; i < 3; ++i)
    {
        tmp[i] = Array<OneD, NekDouble>(m_evalField->GetTotPoints(), 0.0);
    }
    m_evalField->GetCoords(tmp[0], tmp[1], tmp[2]);

    for (int i = 0; i < nSendVars; ++i)
    {
        tmp[3 + i] = field[i];
    }

    LibUtilities::PtsIO ptsIO(m_evalField->GetSession()->GetComm());
    LibUtilities::PtsFieldSharedPtr sPts =
        MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(
            3, m_sendFieldNames, tmp);
    ptsIO.Write(filename, sPts);
}

void CouplingFile::v_Receive(const int step,
                             const NekDouble time,
                             Array<OneD, Array<OneD, NekDouble> > &field,
                             LibUtilities::FieldMetaDataMap &fieldMetaDataMap)
{
    vector<string> fieldNames;

    m_inputFunction->Evaluate(fieldNames, field, time);
}

}
}
