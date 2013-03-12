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

#include <CardiacEPSolver/Filters/FilterCheckpointCellModel.h>

namespace Nektar
{
    std::string FilterCheckpointCellModel::className = GetFilterFactory().RegisterCreatorFunction("CheckpointCellModel", FilterCheckpointCellModel::create);

    FilterCheckpointCellModel::FilterCheckpointCellModel(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::map<std::string, std::string> &pParams) :
        Filter(pSession)
    {
        if (pParams.find("OutputFile") == pParams.end())
        {
            m_outputFile = m_session->GetSessionName();
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
        m_outputIndex = 0;
        m_index = 0;
    }

    FilterCheckpointCellModel::~FilterCheckpointCellModel()
    {

    }

    void FilterCheckpointCellModel::v_Initialise(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time)
    {
        ASSERTL0(m_cell.get(), "Cell model has not been set by EquationSystem "
                "class. Use SetCellModel on this filter to achieve this.");

        m_index = 0;
        m_outputIndex = 0;
    }

    void FilterCheckpointCellModel::v_Update(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time)
    {
        m_index++;
        if (m_index % m_outputFrequency > 0)
        {
            return;
        }

        std::stringstream vOutputFilename;
        vOutputFilename << m_outputFile << "_" << m_outputIndex;

        if (m_session->GetComm()->GetSize() > 1)
        {
            vOutputFilename << "_P" << m_session->GetComm()->GetRank();
        }
        vOutputFilename << ".chk";

        SpatialDomains::MeshGraphSharedPtr vGraph = pFields[0]->GetGraph();

        std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
            = pFields[0]->GetFieldDefinitions();
        std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

        // copy Data into FieldData and set variable
        for(int j = 0; j < m_cell->GetNumCellVariables(); ++j)
        {
            for(int i = 0; i < FieldDef.size(); ++i)
            {
                // Retrieve data from cell model
                Array<OneD, NekDouble> data = m_cell->GetCellSolutionCoeffs(j);

                // Could do a search here to find correct variable
                FieldDef[i]->m_fields.push_back(boost::lexical_cast<std::string>(j));
                pFields[0]->AppendFieldData(FieldDef[i], FieldData[i], data);
            }
        }
        LibUtilities::Write(vOutputFilename.str(),FieldDef,FieldData);
        m_outputIndex++;
    }

    void FilterCheckpointCellModel::v_Finalise(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time)
    {

    }

    bool FilterCheckpointCellModel::v_IsTimeDependent()
    {
        return true;
    }
}
