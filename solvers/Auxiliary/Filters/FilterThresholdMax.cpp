/*
 * FilterThresholdMax.cpp
 *
 *  Created on: 24 Feb 2012
 *      Author: cc
 */

#include <Auxiliary/Filters/FilterThresholdMax.h>

namespace Nektar
{
    std::string FilterThresholdMax::className = GetFilterFactory().RegisterCreatorFunction("ThresholdMax", FilterThresholdMax::create);

    FilterThresholdMax::FilterThresholdMax(const std::map<std::string, std::string> &pParams) :
            Filter()
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
        SpatialDomains::MeshGraphSharedPtr vGraph = pFields[0]->GetGraph();

        std::vector<SpatialDomains::FieldDefinitionsSharedPtr> FieldDef
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

        vGraph->Write(m_outputFile,FieldDef,FieldData);

    }

    bool FilterThresholdMax::v_IsTimeDependent()
    {
        return true;
    }

}
