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

#include <tinyxml/tinyxml.h>
#include <LibUtilities/BasicUtils/VmathArray.hpp>

#include <CardiacEPSolver/Stimuli/Stimulus.h>

namespace Nektar
{
    StimulusFactory& GetStimulusFactory()
    {
        typedef Loki::SingletonHolder<StimulusFactory,
        Loki::CreateUsingNew,
        Loki::NoDestroy > Type;
        return Type::Instance();
    }
    
    /**
     * @class Stimulus
     *
     * The Stimulus class and derived classes implement a range of stimuli.
     * The stimulus contains input stimuli that can be applied throughout the 
     * domain, on specified regions determined by the derived classes of
     * Stimulus, at specified frequencies determined by the derived classes of
     * Protocol.
     *
     */
    
    /**
     * Stimulus base class constructor.
     */
    Stimulus::Stimulus(const LibUtilities::SessionReaderSharedPtr& pSession,
                       const MultiRegions::ExpListSharedPtr& pField,
                       const TiXmlElement* pXml)
    {
        m_session = pSession;
        m_field = pField;
        m_nq = pField->GetTotPoints();
        
        const TiXmlElement* vProtocol = pXml->FirstChildElement("PROTOCOL");
        string vTypeP = vProtocol->Attribute("TYPE");
        
        m_Protocol = GetProtocolFactory().CreateInstance(
                                vTypeP, pSession, vProtocol);
 
    }
    
    
    /**
     * Initialise the stimulus. Allocate workspace and variable storage.
     */
    void Stimulus::Initialise()
    {
        
        
    }
    
    
    /**
     *
     */
    vector<StimulusSharedPtr> Stimulus::LoadStimuli(
                        const LibUtilities::SessionReaderSharedPtr& pSession,
                        const MultiRegions::ExpListSharedPtr& pField)
    {   
        vector<StimulusSharedPtr> vStimList;

        TiXmlElement* vStimuli = pSession->GetElement("Nektar/Stimuli");
        if (vStimuli)
        {
            
            TiXmlElement* vStimulus = vStimuli->FirstChildElement("STIMULUS");
            while (vStimulus)
            {
                
                string vType = vStimulus->Attribute("TYPE");
                //unsigned int vId = atoi(vStimulus->Attribute("ID"));

                vStimList.push_back(GetStimulusFactory().CreateInstance(
                                        vType, pSession, pField, vStimulus));
                vStimulus = vStimulus->NextSiblingElement("STIMULUS");
                
                
                
            }
        }
        return vStimList;
    }
}
