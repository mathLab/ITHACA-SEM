///////////////////////////////////////////////////////////////////////////////
//
// File FilterAverageField.cpp
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
// Description: Average solution fields during time-stepping.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Filters/FilterAverageFields.h>

namespace Nektar
{
    namespace SolverUtils
    {
        std::string FilterAverageFields::className = GetFilterFactory().RegisterCreatorFunction("AverageFields", FilterAverageFields::create);

        FilterAverageFields::FilterAverageFields(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const std::map<std::string, std::string> &pParams) :
            Filter(pSession)
        {
            if (pParams.find("OutputFile") == pParams.end())
            {
                string outname = m_session->GetSessionName() + "_avg.fld";
                m_outputFile = outname;
            }
            else
            {
                ASSERTL0(!(pParams.find("OutputFile")->second.empty()),
                         "Missing parameter 'OutputFile'.");
                m_outputFile = pParams.find("OutputFile")->second;
            }
            ASSERTL0(pParams.find("OutputFrequency") != pParams.end(),
                     "Missing parameter 'OutputFrequency'.");
            m_outputFrequency = atoi(pParams.find("OutputFrequency")->second.c_str());

            m_numAverages  = 0;
            m_index = 0;
            m_fld = MemoryManager<LibUtilities::FieldIO>::AllocateSharedPtr(pSession->GetComm());

        }

        FilterAverageFields::~FilterAverageFields()
        {
        }

        void FilterAverageFields::v_Initialise(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time)
        {
            m_avgFields = Array<OneD,  Array<OneD, NekDouble> >(pFields.num_elements());
            for(int n =0; n < pFields.num_elements(); ++n)
            {
                m_avgFields[n] = Array<OneD, NekDouble>(pFields[n]->GetNcoeffs(),0.0);
            }
            m_avgFieldMetaData["InitialTime"] = boost::lexical_cast<std::string>(time);
        }

        void FilterAverageFields::v_Update(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time)
        {
            m_index++;
            if (m_index % m_outputFrequency > 0)
            {
                return;
            }

            for(int n = 0; n < pFields.num_elements(); ++n)
            {
                Vmath::Vadd(m_avgFields[n].num_elements(),
                            pFields[n]->GetCoeffs(),1,m_avgFields[n],1,
                            m_avgFields[n],1);
            }
            m_numAverages += 1;
            // update FinalTime here since this is when last field will be added. 
            m_avgFieldMetaData["FinalTime"] = boost::lexical_cast<std::string>(time);
        }
        
        void FilterAverageFields::v_Finalise(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time)
        {
            for(int n = 0; n < m_avgFields.num_elements(); ++n)
            {
                Vmath::Smul(m_avgFields[n].num_elements(), 1.0/m_numAverages,m_avgFields[n],1,m_avgFields[n],1);
            }

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
                    pFields[0]->AppendFieldData(FieldDef[i], FieldData[i], m_avgFields[j]);
                }
            }
            m_avgFieldMetaData["NumberOfFieldDumps"] = boost::lexical_cast<std::string>(m_numAverages);
            m_fld->Write(m_outputFile,FieldDef,FieldData,m_avgFieldMetaData);
        }

        bool FilterAverageFields::v_IsTimeDependent()
        {
            return true;
        }
    }
}
