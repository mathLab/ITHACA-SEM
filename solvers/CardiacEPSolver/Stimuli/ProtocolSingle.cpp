///////////////////////////////////////////////////////////////////////////////
//
// File CellModel.cpp
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
// Description: Single impulse protocol.
//
///////////////////////////////////////////////////////////////////////////////

#include <tinyxml/tinyxml.h>
#include <CardiacEPSolver/Stimuli/ProtocolSingle.h>

namespace Nektar
{   std::string ProtocolSingle::className
    = GetProtocolFactory().RegisterCreatorFunction(
                            "ProtocolSingle",
                            ProtocolSingle::create,
                            "Single stimulus protocol.");
    /**
     * @class ProtocolSingle
     *
     * The CellModel class and derived classes implement a range of cell model
     * ODE systems. A cell model comprises a system of ion concentration
     * variables and zero or more gating variables. Gating variables are
     * time-integrated using the Rush-Larsen method and for each variable y,
     * the corresponding y_inf and tau_y value is computed by Update(). The tau
     * values are stored in separate storage to inarray/outarray, #m_gates_tau.
     */
    
    /**
     * Protocol base class constructor.
     */
    ProtocolSingle::ProtocolSingle(const LibUtilities::SessionReaderSharedPtr& pSession,const TiXmlElement* pXml)
        : Protocol(pSession, pXml)
    {
        m_session = pSession;
        
        if (!pXml)
        {
            return;
        }
        
        
        const TiXmlElement *pXmlparameter; //Declaring variable called pxml...
        // See if we have parameters defined.  They are optional so we go on if not.
        
        //member variables m_* defined in ProtocolSingle.h
        
        pXmlparameter = pXml->FirstChildElement("START");
        m_start = atof(pXmlparameter->GetText()); //text value within px1, convert to a floating pt and save in m_px1
        
        pXmlparameter = pXml->FirstChildElement("DURATION");
        m_dur = atof(pXmlparameter->GetText());
        
        
    }
    
    
    
    /**
     * Initialise the protocol. Allocate workspace and variable storage.
     */
    void ProtocolSingle::Initialise()
    {

    }

    
    NekDouble ProtocolSingle::v_GetAmplitude(
                          const NekDouble time)
    {
        if(time > m_start && time < (m_start+m_dur))
        {
            return 1.0;
        }
        else
        {
            return 0.0;
        }
        
    }

    void ProtocolSingle::v_PrintSummary(std::ostream &out)
    {

    }

    void ProtocolSingle::v_SetInitialConditions()
    {

    }

}
