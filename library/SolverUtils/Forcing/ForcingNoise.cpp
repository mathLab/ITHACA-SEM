///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingNoise.cpp
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
// Description: white noise forcing
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <SolverUtils/Forcing/ForcingNoise.h>
#include <MultiRegions/ExpList.h>

using namespace std;

namespace Nektar
{
namespace SolverUtils
{

    std::string ForcingNoise::className = GetForcingFactory().
                                RegisterCreatorFunction("Noise",
                                                        ForcingNoise::create,
                                                        "White Noise Forcing");

    ForcingNoise::ForcingNoise(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const std::weak_ptr<EquationSystem>      &pEquation)
        : Forcing(pSession, pEquation)
    {
    }

    void ForcingNoise::v_InitObject(
            const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
            const unsigned int& pNumForcingFields,
            const TiXmlElement* pForce)
    {
        m_NumVariable = pNumForcingFields;
        int nq         = pFields[0]->GetTotPoints();

        const TiXmlElement* noiseElmt = pForce->FirstChildElement("WHITENOISE");
        ASSERTL0(noiseElmt, "Requires WHITENOISE tag specifying "
                               "magnitude of white noise force.");

        string noiseValue = noiseElmt->GetText();

        m_noise = boost::lexical_cast<NekDouble>(noiseValue);

        // Load optional parameters
        const TiXmlElement* freqElmt = pForce->FirstChildElement("UPDATEFREQ");
        if (freqElmt)
        {
            string freqValue = freqElmt->GetText();
            m_updateFreq = boost::lexical_cast<int>(freqValue);
        }
        else
        {
            // Default is 0 (never update forcing)
            m_updateFreq = 0;
        }

        const TiXmlElement* stepsElmt = pForce->FirstChildElement("NSTEPS");
        if (stepsElmt)
        {
            string stepsValue = stepsElmt->GetText();
            m_numSteps = boost::lexical_cast<int>(stepsValue);
        }
        else
        {
            // Default is 0 (use noise in the entire simulation)
            m_numSteps = 0;
        }

        m_Forcing = Array<OneD, Array<OneD, NekDouble> > (m_NumVariable);

        // Fill forcing: use rank in seed to avoid repeated results
        int seed = - m_session->GetComm()->GetRank();
        m_Forcing[0] = Array<OneD, NekDouble> (nq, 0.0);
        Vmath::FillWhiteNoise(nq,m_noise,m_Forcing[0],1,seed);
        for (int i = 1; i < m_NumVariable; ++i)
        {
            m_Forcing[i] = Array<OneD, NekDouble> (nq, 0.0);
            Vmath::FillWhiteNoise(nq,m_noise,m_Forcing[i],1);
        }

        m_index = 0;
    }

    void ForcingNoise::v_Apply(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> > &inarray,
            Array<OneD, Array<OneD, NekDouble> > &outarray,
            const NekDouble &time)
    {
        boost::ignore_unused(fields, inarray, time);

        // Do not apply forcing if exceeded m_numSteps
        if( m_numSteps && (m_index >= m_numSteps) )
        {
            return;
        }

        // Update forcing (change seed to avoid getting same result)
        if(m_updateFreq && m_index && !((m_index)%m_updateFreq))
        {
            for (int i = 0; i < m_NumVariable; ++i)
            {
                Vmath::FillWhiteNoise(outarray[i].size(),
                                      m_noise,m_Forcing[i],1);
            }
        }

        // Apply forcing
        for (int i = 0; i < m_NumVariable; i++)
        {
            Vmath::Vadd(outarray[i].size(), outarray[i], 1,
                        m_Forcing[i], 1, outarray[i], 1);
        }

        ++m_index;
    }

}
}
