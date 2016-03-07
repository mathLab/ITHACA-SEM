///////////////////////////////////////////////////////////////////////////////
//
// File FilterSampler.cpp
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
// Description: Base clase for filters performing operations on samples
//              of the field.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Filters/FilterSampler.h>

namespace Nektar
{
namespace SolverUtils
{

FilterSampler::FilterSampler(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const ParamMap &pParams)
    : Filter(pSession)
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
        ASSERTL0(it->second.length() > 0, "Missing parameter 'OutputFile'.");
        m_outputFile = it->second;
    }

    // SampleFrequency
    it = pParams.find("SampleFrequency");
    if (it == pParams.end())
    {
        m_sampleFrequency = 1;
    }
    else
    {
        LibUtilities::Equation equ(m_session, it->second);
        m_sampleFrequency = floor(equ.Evaluate());
    }

    // OutputFrequency
    it = pParams.find("OutputFrequency");
    if (it == pParams.end())
    {
        m_outputFrequency = m_session->GetParameter("NumSteps");
    }
    else
    {
        LibUtilities::Equation equ(m_session, it->second);
        m_outputFrequency = floor(equ.Evaluate());
    }

    m_scale       = 1.0;
    m_numSamples  = 0;
    m_index       = 0;
    m_outputIndex = 0;
    m_fld = MemoryManager<LibUtilities::FieldIO>::AllocateSharedPtr(
        pSession->GetComm());
}

FilterSampler::~FilterSampler()
{
}

void FilterSampler::v_Initialise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    int ncoeff = pFields[0]->GetNcoeffs();
    // m_variables need to be filled by a derived class
    m_outFields.resize(m_variables.size());

    for (int n = 0; n < m_variables.size(); ++n)
    {
        m_outFields[n] = Array<OneD, NekDouble>(ncoeff, 0.0);
    }

    m_fieldMetaData["InitialTime"] = boost::lexical_cast<std::string>(time);
}

void FilterSampler::v_Update(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    m_index++;
    if (m_index % m_sampleFrequency > 0)
    {
        return;
    }

    m_numSamples++;
    v_ProcessSample(pFields, time);

    if (m_index % m_outputFrequency == 0)
    {
        m_fieldMetaData["FinalTime"] = boost::lexical_cast<std::string>(time);
        v_PrepareOutput(pFields, time);
        OutputField(pFields, ++m_outputIndex);
    }
}

void FilterSampler::v_Finalise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    m_fieldMetaData["FinalTime"] = boost::lexical_cast<std::string>(time);
    v_PrepareOutput(pFields, time);
    OutputField(pFields);
}

void FilterSampler::OutputField(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, int dump)
{
    for (int n = 0; n < m_outFields.size(); ++n)
    {
        Vmath::Smul(m_outFields[n].num_elements(),
                    m_scale,
                    m_outFields[n],
                    1,
                    m_outFields[n],
                    1);
    }

    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef =
        pFields[0]->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

    Array<OneD, NekDouble> fieldcoeffs;
    int ncoeffs = pFields[0]->GetNcoeffs();

    // copy Data into FieldData and set variable
    for (int j = 0; j < m_outFields.size(); ++j)
    {
        // check to see if field is same order a zeroth field
        if (m_outFields[j].num_elements() == ncoeffs)
        {
            fieldcoeffs = m_outFields[j];
        }
        else
        {
            fieldcoeffs = Array<OneD, NekDouble>(ncoeffs);
            pFields[0]->ExtractCoeffsToCoeffs(
                pFields[j], m_outFields[j], fieldcoeffs);
        }

        for (int i = 0; i < FieldDef.size(); ++i)
        {
            FieldDef[i]->m_fields.push_back(m_variables[j]);
            pFields[0]->AppendFieldData(FieldDef[i], FieldData[i], fieldcoeffs);
        }
    }

    m_fieldMetaData["NumberOfFieldDumps"] =
        boost::lexical_cast<std::string>(m_numSamples);

    std::stringstream outname;
    std::string suffix = v_GetFileSuffix();
    if (dump == -1) // final dump
    {
        outname << m_outputFile << suffix << ".fld";
    }
    else
    {
        outname << m_outputFile << "_" << dump << suffix << ".fld";
    }

    m_fld->Write(outname.str(), FieldDef, FieldData, m_fieldMetaData);

    if (dump != -1) // not final dump so rescale
    {
        for (int n = 0; n < m_outFields.size(); ++n)
        {
            Vmath::Smul(m_outFields[n].num_elements(),
                        1.0 / m_scale,
                        m_outFields[n],
                        1,
                        m_outFields[n],
                        1);
        }
    }
}
}
}
