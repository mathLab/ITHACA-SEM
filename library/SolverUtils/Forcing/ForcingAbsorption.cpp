///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingAbsorption.cpp
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
// Description: Absorption layer forcing
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Forcing/ForcingAbsorption.h>

namespace Nektar
{
namespace SolverUtils
{

    std::string ForcingAbsorption::className = GetForcingFactory().
                                RegisterCreatorFunction("Absorption",
                                                        ForcingAbsorption::create,
                                                        "Forcing Absorption");

    ForcingAbsorption::ForcingAbsorption(const LibUtilities::SessionReaderSharedPtr& pSession)
            : Forcing(pSession),
              m_hasRefFlow(false),	
              m_hasRefFlowTime(false)
    {
    }

    void ForcingAbsorption::v_InitObject(
            const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
            const unsigned int& pNumForcingFields,
            const TiXmlElement* pForce)
    {
        m_NumVariable = pNumForcingFields;
        int npts       = pFields[0]->GetTotPoints();

        const TiXmlElement* funcNameElmt;
        funcNameElmt = pForce->FirstChildElement("COEFF");
        ASSERTL0(funcNameElmt, "Requires COEFF tag, specifying function "
                               "name which prescribes absorption layer coefficient.");

        string funcName = funcNameElmt->GetText();
        ASSERTL0(m_session->DefinesFunction(funcName),
                 "Function '" + funcName + "' not defined.");

        std::string s_FieldStr;
        m_Absorption  = Array<OneD, Array<OneD, NekDouble> > (m_NumVariable);
        m_Forcing = Array<OneD, Array<OneD, NekDouble> > (m_NumVariable);
        for (int i = 0; i < m_NumVariable; ++i)
        {
            s_FieldStr = m_session->GetVariable(i);
            ASSERTL0(m_session->DefinesFunction(funcName, s_FieldStr),
                     "Variable '" + s_FieldStr + "' not defined.");
            m_Absorption[i]  = Array<OneD, NekDouble> (npts, 0.0);
            m_Forcing[i] = Array<OneD, NekDouble> (npts, 0.0);
            EvaluateFunction(pFields, m_session, s_FieldStr,
                             m_Absorption[i], funcName);
        }

        funcNameElmt = pForce->FirstChildElement("REFFLOW");
        if (funcNameElmt)
        {
            string funcName = funcNameElmt->GetText();
            ASSERTL0(m_session->DefinesFunction(funcName),
                     "Function '" + funcName + "' not defined.");

            m_Refflow = Array<OneD, Array<OneD, NekDouble> > (m_NumVariable);
            for (int i = 0; i < m_NumVariable; ++i)
            {
                s_FieldStr = m_session->GetVariable(i);
                ASSERTL0(m_session->DefinesFunction(funcName, s_FieldStr),
                         "Variable '" + s_FieldStr + "' not defined.");
                m_Refflow[i] = Array<OneD, NekDouble> (npts, 0.0);
                EvaluateFunction(pFields, m_session, s_FieldStr,
                                 m_Refflow[i], funcName);
            }
            m_hasRefFlow = true;
        }
   
        funcNameElmt = pForce->FirstChildElement("REFFLOWTIME");
        if (funcNameElmt)
        {
            m_funcNameTime = funcNameElmt->GetText();
            ASSERTL0(m_session->DefinesFunction(funcName),
                     "Function '" + funcName + "' not defined.");
            m_hasRefFlowTime = true;
        }

    }

    void ForcingAbsorption::v_Apply(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> > &inarray,
            Array<OneD, Array<OneD, NekDouble> > &outarray,
            const NekDouble &time)
    {
        int nq = m_Forcing[0].num_elements();
       
        std::string s_FieldStr;
        Array<OneD, NekDouble> TimeScale(1);
        Array<OneD, Array<OneD, NekDouble> > RefflowScaled(m_NumVariable);

        if (m_hasRefFlow)
        {
            for (int i = 0; i < m_NumVariable; i++)
            {
                RefflowScaled[i] = Array<OneD, NekDouble> (nq);
                if (m_hasRefFlowTime)
                {
                    s_FieldStr = m_session->GetVariable(i);
                    EvaluateTimeFunction(m_session, s_FieldStr, TimeScale, m_funcNameTime, time);
                    Vmath::Smul(nq, TimeScale[0], m_Refflow[i],1,RefflowScaled[i],1);
                }
                else
                {
                    Vmath::Vcopy(nq, m_Refflow[i],1, RefflowScaled[i],1);
                }
                

                Vmath::Vsub(nq, inarray[i], 1,
                            RefflowScaled[i], 1, m_Forcing[i], 1);
                Vmath::Vmul(nq, m_Absorption[i], 1,
                            m_Forcing[i], 1, m_Forcing[i], 1);
                Vmath::Vadd(nq, m_Forcing[i], 1,
                            outarray[i], 1, outarray[i], 1);
            }
        }
        else
        {
            for (int i = 0; i < m_NumVariable; i++)
            {
                Vmath::Vmul(nq, m_Absorption[i], 1,
                            inarray[i], 1, m_Forcing[i], 1);
                Vmath::Vadd(nq, m_Forcing[i], 1,
                            outarray[i], 1, outarray[i], 1);
            }
        }
    }
        
}
}
