///////////////////////////////////////////////////////////////////////////////
//
// File ProtocolS1S2.cpp
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
// Description: S1S2 protocol.
//
///////////////////////////////////////////////////////////////////////////////

#include <tinyxml.h>
#include <CardiacEPSolver/Stimuli/ProtocolS1S2.h>

namespace Nektar
{   std::string ProtocolS1S2::className
    = GetProtocolFactory().RegisterCreatorFunction(
                                                   "ProtocolS1S2",
                                                   ProtocolS1S2::create,
                                                   "S1S2 stimulus protocol.");
    /**
     * @class ProtocolS1S2
     *
     * The Stimuli class and derived classes implement a range of stimuli.
     * The stimulus contains input stimuli that can be applied throughout the
     * domain, on specified regions determined by the derived classes of Stimulus,
     * at specified frequencies determined by the derived classes of Protocol.
     *
     */
    /**
     * Protocol base class constructor.
     */
    ProtocolS1S2::ProtocolS1S2(const LibUtilities::SessionReaderSharedPtr& pSession,const TiXmlElement* pXml)
    : Protocol(pSession, pXml)
    {
        m_session = pSession;

        if (!pXml)
        {
            return;
        }

        const TiXmlElement *pXmlparameter; //Declaring variable called pxml...
        // See if we have parameters defined.  They are optional so we go on if not.

        //member variables m_* defined in ProtocolS1S2.h

        pXmlparameter = pXml->FirstChildElement("START");
        m_start = atof(pXmlparameter->GetText()); //text value within px1, convert to a floating pt and save in m_px1

        pXmlparameter = pXml->FirstChildElement("DURATION");
        m_dur = atof(pXmlparameter->GetText());

        pXmlparameter = pXml->FirstChildElement("S1CYCLELENGTH");
        m_s1cyclelength = atof(pXmlparameter->GetText());

        pXmlparameter = pXml->FirstChildElement("NUM_S1");
        m_num_s1 = atof(pXmlparameter->GetText());

        pXmlparameter = pXml->FirstChildElement("S2CYCLELENGTH");
        m_s2cyclelength = atof(pXmlparameter->GetText());

        m_s2start = m_s1cyclelength*(m_num_s1-1)+m_s2cyclelength+m_start;
    }



    /**
     * Initialise the protocol. Allocate workspace and variable storage.
     */
    void ProtocolS1S2::Initialise()
    {

    }


    NekDouble ProtocolS1S2::v_GetAmplitude(const NekDouble time)
    {
        // Number of complete S1 intervals
        NekDouble a = floor((time-m_start)/m_s1cyclelength);

        // Time since start of most recent S1 interval
        NekDouble time1 = time - a * m_s1cyclelength - m_start;

        // S1
        if( (time1 > 0) && (a < m_num_s1) && (time1 < m_dur))
        {
            return 1.0;
        }

        // S2
        if ((time > m_s2start) && (time < m_s2start + m_dur))
        {
            return 1.0;
        }

        return 0.0;
    }

    void ProtocolS1S2::v_GenerateSummary(SolverUtils::SummaryList& s)
    {

    }

    void ProtocolS1S2::v_SetInitialConditions()
    {

    }

}
