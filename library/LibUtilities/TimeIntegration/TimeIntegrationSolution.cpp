///////////////////////////////////////////////////////////////////////////////
//
// File TimeIntegrationSolution.cpp
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
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/TimeIntegration/TimeIntegrationSolution.h>

namespace Nektar
{
namespace LibUtilities
{
TimeIntegrationSolution::TimeIntegrationSolution(
    const TimeIntegrationSchemeData *schemeData, const DoubleArray &y,
    const NekDouble time, const NekDouble timestep)
    : m_schemeData(schemeData), m_solVector(m_schemeData->m_numsteps),
      m_t(m_schemeData->m_numsteps)
{
    m_solVector[0] = y;
    m_t[0]         = time;

    int nsteps = m_schemeData->m_numsteps;

    int nvar           = y.size();
    int npoints        = y[0].size();
    int nMultiStepVals = m_schemeData->GetNmultiStepValues();
    const Array<OneD, const unsigned int> &timeLevels =
        m_schemeData->GetTimeLevelOffset();

    for (int i = 1; i < nsteps; i++)
    {
        m_solVector[i] = Array<OneD, Array<OneD, NekDouble>>(nvar);
        for (int j = 0; j < nvar; j++)
        {
            m_solVector[i][j] = Array<OneD, NekDouble>(npoints, 0.0);
        }
        if (i < nMultiStepVals)
        {
            m_t[i] = time - i * timestep * timeLevels[i];
        }
        else
        {
            m_t[i] = timestep;
        }
    }
}

TimeIntegrationSolution::TimeIntegrationSolution(
    const TimeIntegrationSchemeData *schemeData, const TripleArray &y,
    const Array<OneD, NekDouble> &t)
    : m_schemeData(schemeData), m_solVector(y), m_t(t)
{
    ASSERTL1(y.size() == m_schemeData->m_numsteps,
             "Amount of Entries does not match number of (multi-) steps");
}

TimeIntegrationSolution::TimeIntegrationSolution(
    const TimeIntegrationSchemeData *schemeData, const unsigned int nvar,
    const unsigned int npoints)
    : m_schemeData(schemeData),
      m_solVector(schemeData->m_numsteps),
      m_t(schemeData->m_numsteps)
{
    for (int i = 0; i < m_schemeData->m_numsteps; i++)
    {
        m_solVector[i] = Array<OneD, Array<OneD, NekDouble>>(nvar);
        for (int j = 0; j < nvar; j++)
        {
            m_solVector[i][j] = Array<OneD, NekDouble>(npoints);
        }
    }
}

TimeIntegrationSolution::TimeIntegrationSolution(
    const TimeIntegrationSchemeData *schemeData)
    : m_schemeData(schemeData), m_solVector(m_schemeData->m_numsteps),
      m_t(m_schemeData->m_numsteps)
{
}

std::string TimeIntegrationSolution::GetName() const
{
    return m_schemeData->m_parent->GetFullName();
}

} // end namespace LibUtilities
} // end namespace NekTar
