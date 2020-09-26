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

#include <boost/core/ignore_unused.hpp>

#include <SolverUtils/Filters/FilterThresholdMin.h>

using namespace std;

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
        const std::weak_ptr<EquationSystem>      &pEquation,
        const ParamMap &pParams)
    : Filter(pSession, pEquation)
{
    // ThresholdValue
    auto it = pParams.find("ThresholdValue");
    ASSERTL0(it != pParams.end(), "Missing parameter 'ThresholdValue'.");
    LibUtilities::Equation equ1(
        m_session->GetInterpreter(), it->second);
    m_thresholdValue = equ1.Evaluate();

    // InitialValue
    it = pParams.find("InitialValue");
    ASSERTL0(it != pParams.end(), "Missing parameter 'InitialValue'.");
    LibUtilities::Equation equ2(
        m_session->GetInterpreter(), it->second);
    m_initialValue = equ2.Evaluate();

    // StartTime
    it = pParams.find("StartTime");
    m_startTime = 0.0;
    if (it != pParams.end())
    {
        LibUtilities::Equation equ(
            m_session->GetInterpreter(), it->second);
        m_startTime = equ.Evaluate();
    }

    // OutputFile
    it = pParams.find("OutputFile");
    m_outputFile = pSession->GetSessionName() + "_max.fld";
    if (it != pParams.end())
    {
        m_outputFile = it->second;
    }

    // ThresholdVar
    it = pParams.find("ThresholdVar");
    m_thresholdVar = 0;
    if (it != pParams.end())
    {
        std::string var = it->second.c_str();
        std::vector<string> varlist = pSession->GetVariables();
        auto x = std::find(varlist.begin(), varlist.end(), var);
        ASSERTL0(x != varlist.end(),
                 "Specified variable " + var +
                 " in ThresholdMin filter is not available.");
        m_thresholdVar = x - varlist.begin();
    }

    m_fld = LibUtilities::FieldIO::CreateDefault(pSession);
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
    boost::ignore_unused(time);

    m_threshold = Array<OneD, NekDouble> (
            pFields[m_thresholdVar]->GetNpoints(), m_initialValue);
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
    NekDouble timestep = pFields[m_thresholdVar]->GetSession()->GetParameter("TimeStep");

    for (i = 0; i < pFields[m_thresholdVar]->GetNpoints(); ++i)
    {
        if (m_threshold[i] < timestep &&
            pFields[m_thresholdVar]->GetPhys()[i] > m_thresholdValue)
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
    boost::ignore_unused(time);

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
