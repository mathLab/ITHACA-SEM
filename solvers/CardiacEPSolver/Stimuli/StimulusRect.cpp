///////////////////////////////////////////////////////////////////////////////
//
// File StimulusRect.cpp
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

#include <CardiacEPSolver/Stimuli/StimulusRect.h>

namespace Nektar
{
    std::string StimulusRect::className
          = GetStimulusFactory().RegisterCreatorFunction(
                    "StimulusRect",
                    StimulusRect::create,
                     "Rectangular stimulus.");

    /**
     * @class StimulusRect
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
    StimulusRect::StimulusRect(
            const LibUtilities::SessionReaderSharedPtr& pSession,
            const MultiRegions::ExpListSharedPtr& pField,
            const TiXmlElement* pXml)
            : Stimulus(pSession, pField, pXml)
    {
        m_session = pSession;
        m_field = pField;
        m_nq = pField->GetTotPoints();
        m_chiCapMembrane = m_session->GetParameter("chi")
                            * m_session->GetParameter("Cm");

        if (!pXml)
        {
            return;
        }


        const TiXmlElement *pXmlparameter;

        pXmlparameter = pXml->FirstChildElement("p_x1");
        m_px1 = atof(pXmlparameter->GetText());

        pXmlparameter = pXml->FirstChildElement("p_y1");
        m_py1 = atof(pXmlparameter->GetText());

        pXmlparameter = pXml->FirstChildElement("p_z1");
        m_pz1 = atof(pXmlparameter->GetText());

        pXmlparameter = pXml->FirstChildElement("p_x2");
        m_px2 = atof(pXmlparameter->GetText());

        pXmlparameter = pXml->FirstChildElement("p_y2");
        m_py2 = atof(pXmlparameter->GetText());

        pXmlparameter = pXml->FirstChildElement("p_z2");
        m_pz2 = atof(pXmlparameter->GetText());

        pXmlparameter = pXml->FirstChildElement("p_is");
        m_pis = atof(pXmlparameter->GetText());

        pXmlparameter = pXml->FirstChildElement("p_strength");
        m_strength = atof(pXmlparameter->GetText());
    }


    /**
     * Initialise the stimulus. Allocate workspace and variable storage.
     */
    void StimulusRect::Initialise()
    {
    }

    /**
     *
     */
    void StimulusRect::v_Update(
            Array<OneD, Array<OneD, NekDouble> >&outarray,
            const NekDouble time)
    {
        if (m_field->GetNumElmts() == 0)
        {
            return;
        }

        // Get the dimension of the expansion
        int dim = m_field->GetCoordim(0);

        // Retrieve coordinates of quadrature points
        int nq = m_field->GetNpoints();
        Array<OneD,NekDouble> x0(nq);
        Array<OneD,NekDouble> x1(nq);
        Array<OneD,NekDouble> x2(nq);
        m_field->GetCoords(x0,x1,x2);

        // Get the protocol amplitude
        NekDouble v_amp = m_Protocol->GetAmplitude(time) * m_strength
                            / m_chiCapMembrane;

        switch (dim)
        {
        case 1:
            for(int j=0; j<nq; j++)
            {
                outarray[0][j] += v_amp * ( ( tanh(m_pis*(x0[j] - m_px1))
                                              - tanh(m_pis*(x0[j] - m_px2))
                                            ) / 2.0 );
            }
            break;
        case 2:
            for(int j=0; j<nq; j++)
            {
                outarray[0][j] += v_amp * ( ( (tanh(m_pis*(x0[j] - m_px1))
                                               - tanh(m_pis*(x0[j] - m_px2)))
                                            * (tanh(m_pis*(x1[j] - m_py1))
                                               - tanh(m_pis*(x1[j] - m_py2)))
                                            ) / 2.0 );
            }
            break;
        case 3:
            for(int j=0; j<nq; j++)
            {
                outarray[0][j] += v_amp * ( ( (tanh(m_pis*(x0[j] - m_px1))
                                               - tanh(m_pis*(x0[j] - m_px2)))
                                            * (tanh(m_pis*(x1[j] - m_py1))
                                               - tanh(m_pis*(x1[j] - m_py2)))
                                            * (tanh(m_pis*(x2[j] - m_pz1))
                                               - tanh(m_pis*(x2[j] - m_pz2)))
                                            ) / 2.0 );
            }
            break;
        }
    }


    /**
     *
     */
    void StimulusRect::v_GenerateSummary(SolverUtils::SummaryList& s)
    {
    }
}
