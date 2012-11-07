///////////////////////////////////////////////////////////////////////////////
//
// File Stimulus.cpp
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
// Description: Stimulus base class.
//
///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <tinyxml/tinyxml.h>
#include <LibUtilities/BasicUtils/VmathArray.hpp>

#include <StdRegions/StdNodalTriExp.h>

#include <CardiacEPSolver/Stimuli/StimulusRect.h>



//#include <LibUtilities/LinearAlgebra/Blas.hpp>

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
     * domain, on specified regions determined by the derived classes of Stimulus,
     * at specified frequencies determined by the derived classes of Protocol.
     *
     */
    
    /**
     * Stimulus base class constructor.
     */
    StimulusRect::StimulusRect(const LibUtilities::SessionReaderSharedPtr& pSession,
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
        

        const TiXmlElement *pXmlparameter; //Declaring variable called pxml...
        // See if we have parameters defined.  They are optional so we go on if not.
    
        
        //member variables m_p defined in StimulusRect.h
        
        pXmlparameter = pXml->FirstChildElement("p_x1");
        m_px1 = atof(pXmlparameter->GetText()); //text value within px1, convert to a floating pt and save in m_px1
        
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
    }
    
    /**
     * Initialise the stimulus. Allocate workspace and variable storage.
     */
    void StimulusRect::Initialise()
    {

        
    }
    
    void StimulusRect::v_Update(Array<OneD, Array<OneD, NekDouble> >&outarray,
                          const NekDouble time)
    {

    }
    
    void StimulusRect::v_PrintSummary(std::ostream &out)
    {


    }

}
