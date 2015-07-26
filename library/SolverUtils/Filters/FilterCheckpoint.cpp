///////////////////////////////////////////////////////////////////////////////
//
// File FilterCheckpoint.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
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
// Description: Outputs solution fields during time-stepping.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Filters/FilterCheckpoint.h>

namespace Nektar
{
namespace SolverUtils
{
std::string FilterCheckpoint::className =
        GetFilterFactory().RegisterCreatorFunction(
                "Checkpoint", FilterCheckpoint::create);

FilterCheckpoint::FilterCheckpoint(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const ParamMap &pParams) :
    Filter(pSession)
{
    ParamMap::const_iterator it;

    // OutputFile
    it = pParams.find("OutputFile");
    if (it == pParams.end())
    {
        m_outputFile = m_session->GetSessionName();
    }
    else
    {
        ASSERTL0(it->second.length() > 0, "Empty parameter 'OutputFile'.");
        m_outputFile = it->second;
    }

    // OutputFrequency
    ASSERTL0(it != pParams.end(), "Missing parameter 'OutputFrequency'.");
    LibUtilities::Equation equ(m_session, it->second);
    m_outputFrequency = floor(equ.Evaluate());

    m_outputIndex = 0;
    m_index       = 0;
    m_fld         = MemoryManager<LibUtilities::FieldIO>
                        ::AllocateSharedPtr(pSession->GetComm());

}

FilterCheckpoint::~FilterCheckpoint()
{

}

void FilterCheckpoint::v_Initialise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time)
{
    m_index = 0;
    m_outputIndex = 0;
}

void FilterCheckpoint::v_Update(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time)
{
    m_index++;
    if (m_index % m_outputFrequency > 0)
    {
        return;
    }

    std::stringstream vOutputFilename;
    vOutputFilename << m_outputFile << "_" << m_outputIndex << ".chk";

    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
        = pFields[0]->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

    // copy Data into FieldData and set variable
    for(int j = 0; j < pFields.num_elements(); ++j)
    {
        for(int i = 0; i < FieldDef.size(); ++i)
        {
            // Could do a search here to find correct variable
            FieldDef[i]->m_fields.push_back(m_session->GetVariable(j));
            pFields[0]->AppendFieldData(FieldDef[i],
                                        FieldData[i],
                                        pFields[j]->UpdateCoeffs());
        }
    }
    m_fld->Write(vOutputFilename.str(),FieldDef,FieldData);
    m_outputIndex++;
}

void FilterCheckpoint::v_Finalise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time)
{

}

bool FilterCheckpoint::v_IsTimeDependent()
{
    return true;
}
}
}
