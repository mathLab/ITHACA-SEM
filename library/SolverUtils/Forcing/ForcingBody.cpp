///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingBody.cpp
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
// Description: Body forcing
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Forcing/ForcingBody.h>
#include <MultiRegions/ExpList.h>

namespace Nektar
{
namespace SolverUtils
{

    std::string ForcingBody::classNameBody = GetForcingFactory().
        RegisterCreatorFunction("Body",
                                ForcingBody::create,
                                "Body Forcing");
    std::string ForcingBody::classNameField = GetForcingFactory().
        RegisterCreatorFunction("Field",
                                ForcingBody::create,
                                "Field Forcing");

    ForcingBody::ForcingBody(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        : Forcing(pSession),
          m_hasTimeFcnScaling(false)
    {
    }

    void ForcingBody::v_InitObject(
                                   const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
                                   const unsigned int& pNumForcingFields,
                                   const TiXmlElement* pForce)
    {
        m_NumVariable = pNumForcingFields;
        int nq        = pFields[0]->GetTotPoints();

        const TiXmlElement* funcNameElmt = pForce->FirstChildElement("BODYFORCE");
        if(!funcNameElmt)
        {
            funcNameElmt = pForce->FirstChildElement("FIELDFORCE");

            ASSERTL0(funcNameElmt, "Requires BODYFORCE or FIELDFORCE tag "
                     "specifying function name which prescribes body force.");
        }

        string funcName = funcNameElmt->GetText();
        ASSERTL0(m_session->DefinesFunction(funcName),
                 "Function '" + funcName + "' not defined.");

        // Time function is optional
        funcNameElmt = pForce->FirstChildElement("BODYFORCETIMEFCN");
        if(!funcNameElmt)
        {
            funcNameElmt = pForce->FirstChildElement("FIELDFORCETIMEFCN");
        }

        // Load time function if specified
        if(funcNameElmt)
        {
            std::string funcNameTime = funcNameElmt->GetText();

            ASSERTL0(!funcNameTime.empty(),
                     "Expression must be given in BODYFORCETIMEFCN or "
                     "FIELDFORCETIMEFCN.");

            m_session->SubstituteExpressions(funcNameTime);
            m_timeFcnEqn = MemoryManager<LibUtilities::Equation>
                            ::AllocateSharedPtr(m_session,funcNameTime);

            m_hasTimeFcnScaling = true;
        }

        m_Forcing = Array<OneD, Array<OneD, NekDouble> > (m_NumVariable);

        std::string  s_FieldStr;

        for (int i = 0; i < m_NumVariable; ++i)
        {
            m_Forcing[i] = Array<OneD, NekDouble> (nq, 0.0);
            s_FieldStr   = m_session->GetVariable(i);
            ASSERTL0(m_session->DefinesFunction(funcName, s_FieldStr),
                     "Variable '" + s_FieldStr + "' not defined.");
            EvaluateFunction(pFields, m_session, s_FieldStr,
                             m_Forcing[i], funcName);
        }

        bool singleMode, halfMode;
        m_session->MatchSolverInfo("ModeType","SingleMode",singleMode,false);
        m_session->MatchSolverInfo("ModeType","HalfMode",  halfMode,  false);
        bool homogeneous = pFields[0]->GetExpType() == MultiRegions::e3DH1D;

        // If singleMode or halfMode, transform the forcing term to be in
        // physical space in the plane, but Fourier space in the homogeneous
        // direction
        if (singleMode || halfMode || homogeneous)
        {
            for (int i = 0; i < m_NumVariable; ++i)
            {
                pFields[0]->HomogeneousFwdTrans(m_Forcing[i],m_Forcing[i]);
            }
        }
    }

    void ForcingBody::v_Apply(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> > &inarray,
            Array<OneD, Array<OneD, NekDouble> > &outarray,
            const NekDouble &time)
    {

        if(m_hasTimeFcnScaling)
        {
            Array<OneD, NekDouble>  TimeFcn(1);

            for (int i = 0; i < m_NumVariable; i++)
            {
                EvaluateTimeFunction(time, m_timeFcnEqn, TimeFcn);

                Vmath::Svtvp(outarray[i].num_elements(), TimeFcn[0],
                             m_Forcing[i], 1,
                             outarray[i],  1,
                             outarray[i],  1);
            }
        }
        else
        {
            for (int i = 0; i < m_NumVariable; i++)
            {
                Vmath::Vadd(outarray[i].num_elements(), outarray[i], 1,
                            m_Forcing[i], 1, outarray[i], 1);
            }
        }
    }

}
}
