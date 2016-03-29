///////////////////////////////////////////////////////////////////////////////
//
// File FilterMovingAverage.cpp
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
// Description: Calculates exponential moving average of solution fields
//              during time-stepping.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Filters/FilterMovingAverage.h>

namespace Nektar
{
namespace SolverUtils
{
std::string FilterMovingAverage::className =
    GetFilterFactory().RegisterCreatorFunction("MovingAverage",
                                               FilterMovingAverage::create);

FilterMovingAverage::FilterMovingAverage(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const ParamMap &pParams)
    : FilterFieldConvert(pSession, pParams)
{
    ParamMap::const_iterator it;

    // Load filter parameter
    it = pParams.find("alpha");
    if (it == pParams.end())
    {
        it = pParams.find("tau");
        if (it == pParams.end())
        {
            ASSERTL0(false, "FilterMovingAverage needs either alpha or tau.");
        }
        else
        {
            // Load time constant
            LibUtilities::Equation equ(m_session, it->second);
            NekDouble tau = equ.Evaluate();
            // Load delta T between samples
            NekDouble dT;
            m_session->LoadParameter("TimeStep", dT);
            dT = dT * m_sampleFrequency;
            // Calculate alpha
            m_alpha = dT / (tau + dT);
        }
    }
    else
    {
        LibUtilities::Equation equ(m_session, it->second);
        m_alpha = equ.Evaluate();
        // Check if tau was also defined
        it = pParams.find("tau");
        if (it != pParams.end())
        {
            ASSERTL0(false,
                     "Cannot define both alpha and tau in MovingAverage.");
        }
    }
    // Check bounds of m_alpha
    ASSERTL0(m_alpha > 0 && m_alpha < 1, "Alpha out of bounds.");
}

FilterMovingAverage::~FilterMovingAverage()
{
}

void FilterMovingAverage::v_ProcessSample(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    // Take first sample as initial vector
    NekDouble alpha = m_alpha;
    if (m_numSamples == 1)
    {
        alpha = 1.0;
    }

    // \bar{u}_n = alpha * u_n + (1-alpha) * \bar{u}_{n-1}
    for (int n = 0; n < pFields.num_elements(); ++n)
    {
        Vmath::Svtsvtp(m_outFields[n].num_elements(),
                       alpha,
                       pFields[n]->GetCoeffs(),
                       1,
                       (1.0 - alpha),
                       m_outFields[n],
                       1,
                       m_outFields[n],
                       1);
    }
}

}
}
