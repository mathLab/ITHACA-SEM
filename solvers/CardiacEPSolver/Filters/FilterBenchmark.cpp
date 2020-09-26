///////////////////////////////////////////////////////////////////////////////
//
// File FilterBenchmark.cpp
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
// Description: Outputs times when solution crosses a threshold value.
//
///////////////////////////////////////////////////////////////////////////////

#include <CardiacEPSolver/Filters/FilterBenchmark.h>

namespace Nektar
{
std::string FilterBenchmark::className =
        SolverUtils::GetFilterFactory().RegisterCreatorFunction(
                "Benchmark",
                FilterBenchmark::create);

/**
 * @class FilterBenchmark
 *
 * This class records the sequence of activation and repolarisation times across
 * the entire domain into a two-dimensional storage structure. At each
 * timestep, the voltage at each point in the domain is examined to identify if
 * it has crossed the threshold value. If so, the time of crossing is recorded.
 * Auxiliary arrays hold the current index of each point (i.e. the number of
 * crossings of the threshold) and the type of the last crossing (activation or
 * repolarisation).
 */

/**
 * @param       pSession    Session reader for IO
 * @param       pParams     Parameters of filter
 */
FilterBenchmark::FilterBenchmark(
        const LibUtilities::SessionReaderSharedPtr         &pSession,
        const std::weak_ptr<SolverUtils::EquationSystem> &pEquation,
        const ParamMap &pParams)
    : Filter(pSession, pEquation)
{
    // ThresholdValue
    auto it = pParams.find("ThresholdValue");
    ASSERTL0(it != pParams.end(), "Missing parameter 'ThresholdValue'.");
    LibUtilities::Equation equ1(
        m_session->GetInterpreter(), it->second);
    m_thresholdValue = floor(equ1.Evaluate());

    // InitialValue
    it = pParams.find("InitialValue");
    ASSERTL0(it != pParams.end(), "Missing parameter 'InitialValue'.");
    LibUtilities::Equation equ2(
        m_session->GetInterpreter(), it->second);
    m_initialValue = floor(equ2.Evaluate());

    // OutputFile
    it = pParams.find("OutputFile");
    ASSERTL0(it->second.length() > 0, "Missing parameter 'OutputFile'.");
    m_outputFile = it->second;

    // StartTime
    m_startTime = 0.0;
    it = pParams.find("StartTime");
    if (it != pParams.end())
    {
        LibUtilities::Equation equ(m_session->GetInterpreter(), it->second);
        m_startTime = floor(equ.Evaluate());
    }

    m_fld = LibUtilities::FieldIO::CreateDefault(pSession);

}


/**
 *
 */
FilterBenchmark::~FilterBenchmark()
{

}


/*
 * Initialises the storage.
 * @param       pFields     Field storage expansion lists
 * @param       time        Current time
 */
void FilterBenchmark::v_Initialise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time)
{
    m_threshold.push_back(Array<OneD, NekDouble>(
                            pFields[0]->GetNpoints(), m_initialValue));

    m_idx      = Array<OneD, int> (pFields[0]->GetNpoints(),  0);
    m_polarity = Array<OneD, int> (pFields[0]->GetNpoints(), -1);
}


/**
 * Checks each point in the domain to determine if it has crossed the threshold.
 * The direction of crossing is determined. Additional storage is allocated if
 * needed.
 * @param       pFields     Field storage expansion lists
 * @param       time        Current time
 */
void FilterBenchmark::v_Update(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time)
{
    // Only proceed if the start time has passed
    if (time < m_startTime)
    {
        return;
    }

    // Examine each point in turn
    for (int i = 0; i < pFields[0]->GetNpoints(); ++i)
    {
        if ((m_polarity[i] == -1 &&
                pFields[0]->GetPhys()[i] > m_thresholdValue) ||
            (m_polarity[i] == 1 &&
                pFields[0]->GetPhys()[i] < m_thresholdValue))
        {
            // If APD less than 50ms, remove last activation
            if (m_polarity[i] == 1 &&
                time - m_threshold[m_idx[i]][i] < 50)
            {
                m_idx[i]--;
                m_threshold[m_idx[i]][i] = m_initialValue;
            }
            else
            {
                m_threshold[m_idx[i]][i] = time;
                m_idx[i]++;
            }
            // Update polarity of last crossing
            m_polarity[i] *= -1;
        }
    }

    // Allocate additional storage if any point has as many crossings as
    // current storage permits.
    int max_idx = Vmath::Vmax(pFields[0]->GetNpoints(), m_idx, 1);
    pFields[0]->GetSession()->GetComm()->AllReduce(max_idx,
                            LibUtilities::ReduceMax);
    if (m_threshold.size() == max_idx)
    {
        m_threshold.push_back(Array<OneD, NekDouble>(
                            pFields[0]->GetNpoints(), m_initialValue));
    }

}


/**
 * Writes out the crossings to file.
 * @param       pFields     Field storage expansion list.
 * @param       time        Current time.
 */
void FilterBenchmark::v_Finalise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time)
{
    for (int i = 0; i < m_threshold.size() - 1; ++i)
    {
        std::stringstream vOutputFilename;
        vOutputFilename << m_outputFile << "_" << i << ".fld";

        std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
            = pFields[0]->GetFieldDefinitions();
        std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

        Array<OneD, NekDouble> vCoeffs(pFields[0]->GetNcoeffs());
        pFields[0]->FwdTrans_IterPerExp(m_threshold[i], vCoeffs);

        // copy Data into FieldData and set variable
        for(int i = 0; i < FieldDef.size(); ++i)
        {
            // Could do a search here to find correct variable
            FieldDef[i]->m_fields.push_back("m");
            pFields[0]->AppendFieldData(FieldDef[i], FieldData[i], vCoeffs);
        }

        m_fld->Write(vOutputFilename.str(),FieldDef,FieldData);
    }
}


/**
 * @return This filter is time dependent.
 */
bool FilterBenchmark::v_IsTimeDependent()
{
    return true;
}

}
