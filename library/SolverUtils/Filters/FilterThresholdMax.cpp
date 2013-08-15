///////////////////////////////////////////////////////////////////////////////
//
// File FilterThresholdMax.cpp
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
// Description: Outputs time when solution first exceeds a threshold value.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Filters/FilterThresholdMax.h>

namespace Nektar
{
    namespace SolverUtils
    {
        std::string FilterThresholdMax::className = GetFilterFactory().RegisterCreatorFunction("ThresholdMax", FilterThresholdMax::create);

        FilterThresholdMax::FilterThresholdMax(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const std::map<std::string, std::string> &pParams) :
            Filter(pSession)
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
        }

        FilterThresholdMax::~FilterThresholdMax()
        {

        }

        void FilterThresholdMax::v_Initialise(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time)
        {
            m_threshold = Array<OneD, NekDouble> (pFields[0]->GetNpoints(), m_initialValue);
        }

        void FilterThresholdMax::v_Update(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time)
        {
            int i;
            NekDouble timestep = pFields[0]->GetSession()->GetParameter("TimeStep");

            for (i = 0; i < pFields[0]->GetNpoints(); ++i)
            {
                if (m_threshold[i] < timestep && pFields[0]->GetPhys()[i] > m_thresholdValue)
                {
                    m_threshold[i] = time;
                }
            }
        }

        void FilterThresholdMax::v_Finalise(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time)
        {
            std::stringstream vOutputFilename;
            vOutputFilename << m_outputFile;
            if (m_session->GetComm()->GetSize() > 1)
            {
                vOutputFilename << "_P" << m_session->GetComm()->GetRank();
            }
            vOutputFilename << ".fld";

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

            LibUtilities::Write(vOutputFilename.str(),FieldDef,FieldData);

        }

        bool FilterThresholdMax::v_IsTimeDependent()
        {
            return true;
        }
    }
}
