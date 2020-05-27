///////////////////////////////////////////////////////////////////////////////
//
// File StimulusPoint.cpp
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
// Description: Rectangular stimulus class
//
///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <tinyxml.h>
#include <LibUtilities/BasicUtils/VmathArray.hpp>

#include <CardiacEPSolver/Stimuli/StimulusPoint.h>

namespace Nektar
{
    std::string StimulusPoint::className
          = GetStimulusFactory().RegisterCreatorFunction(
                    "StimulusPoint",
                    StimulusPoint::create,
                     "Point stimulus.");

    /**
     * @class StimulusPoint
     *
     * The Stimulus class and derived classes implement a range of stimuli.
     * The stimulus contains input stimuli that can be applied throughout the
     * domain, on specified regions determined by the derived classes of
     * Stimulus, at specified frequencies determined by the derived classes of
     * Protocol.
     */

    /**
     * Stimulus base class constructor.
     */
    StimulusPoint::StimulusPoint(
            const LibUtilities::SessionReaderSharedPtr& pSession,
            const MultiRegions::ExpListSharedPtr& pField,
            const TiXmlElement* pXml)
            : Stimulus(pSession, pField, pXml)
    {
        m_session = pSession;
        m_field = pField;
        m_nq = pField->GetTotPoints();

        if (!pXml)
        {
            return;
        }

        const TiXmlElement *pXmlparameter;

        pXmlparameter = pXml->FirstChildElement("p_strength");
        m_strength = atof(pXmlparameter->GetText());
    }


    /**
     * Initialise the stimulus. Allocate workspace and variable storage.
     */
    void StimulusPoint::Initialise()
    {
    }


    /**
     *
     */
    void StimulusPoint::v_Update(
            Array<OneD, Array<OneD, NekDouble> >&outarray,
            const NekDouble time)
    {
        // Get the protocol amplitude
        NekDouble v_amp = m_Protocol->GetAmplitude(time) * m_strength;

        outarray[0][0] += v_amp;
    }


    /**
     *
     */
    void StimulusPoint::v_GenerateSummary(SolverUtils::SummaryList& s)
    {
    }
}
