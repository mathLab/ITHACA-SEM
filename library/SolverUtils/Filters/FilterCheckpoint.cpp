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
#include <GlobalMapping/Mapping.h>

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
    it = pParams.find("OutputFrequency");
    ASSERTL0(it != pParams.end(), "Missing parameter 'OutputFrequency'.");
    LibUtilities::Equation equ(m_session, it->second);
    m_outputFrequency = floor(equ.Evaluate());

    m_fld = LibUtilities::FieldIO::CreateDefault(pSession);
}

FilterCheckpoint::~FilterCheckpoint()
{

}

void FilterCheckpoint::v_Initialise(
        const LibUtilities::FieldMetaDataMap &fieldMetaDataMap,
        const Array<OneD, const Array<OneD, NekDouble > > &coeffs,
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time)
{
    m_index = 0;
    m_outputIndex = 0;

    LibUtilities::FieldMetaDataMap tmp = fieldMetaDataMap;
    LibUtilities::FieldMetaDataMap::iterator iter;
    iter = tmp.find("ChkFileNum");
    if (iter != fieldMetaDataMap.end())
    {
        m_outputIndex = boost::lexical_cast<NekDouble>(iter->second);
    }

    v_Update(fieldMetaDataMap, coeffs, pFields, time);
}

void FilterCheckpoint::v_Update(
        const LibUtilities::FieldMetaDataMap &fieldMetaDataMap,
        const Array<OneD, const Array<OneD, NekDouble > > &coeffs,
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time)
{
    if (m_index++ % m_outputFrequency > 0)
    {
        return;
    }

    std::vector<std::string> variables;
    std::string allVars = fieldMetaDataMap["Variables"] + fieldMetaDataMap["AuxVariables"];
    ParseUtils::GenerateOrderedStringVector(allVars.c_str(), variables);

    std::string outname =  m_outputFile +  "_" +
        boost::lexical_cast<std::string>(m_outputIndex) + ".chk";

    // make sure the coeffs match the first expansion in pFields
    std::vector<Array<OneD, NekDouble> > fieldcoeffs(
        pFields.num_elements());
    for (int i = 0; i < pFields.num_elements(); ++i)
    {
        if (pFields[i]->GetNcoeffs() == pFields[0]->GetNcoeffs())
        {
            fieldcoeffs[i] = coeffs[i];
        }
        else
        {
            fieldcoeffs[i] = Array<OneD,NekDouble>(pFields[0]->GetNcoeffs());
            pFields[0]->ExtractCoeffsToCoeffs(pFields[i],
                                              coeffs[i],
                                              fieldcoeffs[i]);
        }
    }

    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
        = pFields[0]->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

    // Copy Data into FieldData and set variable
    for(int j = 0; j < fieldcoeffs.size(); ++j)
    {
        for(int i = 0; i < FieldDef.size(); ++i)
        {
            // Could do a search here to find correct variable
            FieldDef[i]->m_fields.push_back(variables[j]);
            pFields[0]->AppendFieldData(FieldDef[i],
                                        FieldData[i],
                                        fieldcoeffs[j]);
        }
    }

    // adjust the metadata before writing to file
    LibUtilities::FieldMetaDataMap metaDataMap = fieldMetaDataMap;

    // Update time in field info if required
    if(metaDataMap.find("Time") != metaDataMap.end())
    {
        metaDataMap["Time"] = boost::lexical_cast<std::string>(time);
    }

    // Update step in field info if required
    if(metaDataMap.find("ChkFileNum") != metaDataMap.end())
    {
        metaDataMap["ChkFileNum"] = boost::lexical_cast<std::string>(m_outputIndex);
    }

    // If necessary, add mapping information to metadata
    //      and output mapping coordinates
    Array<OneD, MultiRegions::ExpListSharedPtr> fields(1);
    fields[0] = pFields[0];
    GlobalMapping::MappingSharedPtr mapping =
            GlobalMapping::Mapping::Load(m_session, fields);
    mapping->Output(metaDataMap, outname);

    m_fld->Write(outname, FieldDef, FieldData, metaDataMap, true);

    m_outputIndex++;
}

void FilterCheckpoint::v_Finalise(
        const LibUtilities::FieldMetaDataMap &fieldMetaDataMap,
        const Array<OneD, const Array<OneD, NekDouble > > &coeffs,
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
