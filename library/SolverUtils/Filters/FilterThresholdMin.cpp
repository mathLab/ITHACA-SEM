///////////////////////////////////////////////////////////////////////////////
//
// File FilterThresholdMin.cpp
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
// Description: Outputs time when solution first drops below a threshold value.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Filters/FilterThresholdMin.h>

namespace Nektar
{
namespace SolverUtils
{

std::string FilterThresholdMin::className = 
        GetFilterFactory().RegisterCreatorFunction(
                "ThresholdMin", FilterThresholdMin::create);


/**
 *
 */
FilterThresholdMin::FilterThresholdMin(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::map<std::string, std::string> &pParams)
    : Filter(pSession)
{
    ASSERTL0(pParams.find("ThresholdValue") != pParams.end(),
             "Missing parameter 'ThresholdValue'.");
    m_thresholdValue = atof(pParams.find("ThresholdValue")->second.c_str());
    ASSERTL0(pParams.find("InitialValue") != pParams.end(),
             "Missing parameter 'InitialValue'.");
    m_initialValue = atof(pParams.find("InitialValue")->second.c_str());
    ASSERTL0(!(pParams.find("OutputFile")->second.empty()),
             "Missing parameter 'OutputFile'.");
    m_outputFile = pParams.find("OutputFile")->second;

    if (pParams.find("StartTime") != pParams.end())
    {
        m_startTime = atof(pParams.find("StartTime")->second.c_str());
    }

    m_fld = MemoryManager<LibUtilities::FieldIO>::AllocateSharedPtr(
                                                        pSession->GetComm());

}


/**
 *
 */
FilterThresholdMin::~FilterThresholdMin()
{

}


/**
 *
 */
void FilterThresholdMin::v_Initialise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, 
        const NekDouble &time)
{
    m_threshold = Array<OneD, NekDouble> (pFields[0]->GetNpoints(),
                                          m_initialValue);
}


/**
 *
 */
void FilterThresholdMin::v_Update(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, 
        const NekDouble &time)
{
    if (time < m_startTime)
    {
        return;
    }

    int i;
    NekDouble timestep = pFields[0]->GetSession()->GetParameter("TimeStep");

    for (i = 0; i < pFields[0]->GetNpoints(); ++i)
    {
        if (m_threshold[i]           < timestep && 
            pFields[0]->GetPhys()[i] < m_thresholdValue)
        {
            m_threshold[i] = time;
        }
    }
}


/**
 *
 */
void FilterThresholdMin::v_Finalise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, 
        const NekDouble &time)
{
    std::stringstream vOutputFilename;
    vOutputFilename << m_outputFile << ".fld";

    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
        = pFields[0]->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

    Array<OneD, NekDouble> vCoeffs(pFields[0]->GetNcoeffs());
    pFields[0]->FwdTrans_IterPerExp(m_threshold, vCoeffs);

    // copy Data into FieldData and set variable
    for(int i = 0; i < FieldDef.size(); ++i)
    {
        // Could do a search here to find correct variable
        FieldDef[i]->m_fields.push_back("m");
        pFields[0]->AppendFieldData(FieldDef[i], FieldData[i], vCoeffs);
    }

    m_fld->Write(vOutputFilename.str(),FieldDef,FieldData);

}


/**
 *
 */
bool FilterThresholdMin::v_IsTimeDependent()
{
    return true;
}

}
}
