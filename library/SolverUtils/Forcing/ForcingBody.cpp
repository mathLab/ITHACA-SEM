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

using namespace std;

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
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const std::weak_ptr<EquationSystem>      &pEquation)
        : Forcing(pSession, pEquation),
          m_hasTimeFcnScaling(false)
    {
    }

    void ForcingBody::v_InitObject(
                                   const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
                                   const unsigned int& pNumForcingFields,
                                   const TiXmlElement* pForce)
    {
        m_NumVariable = pNumForcingFields;

        const TiXmlElement* funcNameElmt = pForce->FirstChildElement("BODYFORCE");
        if(!funcNameElmt)
        {
            funcNameElmt = pForce->FirstChildElement("FIELDFORCE");

            ASSERTL0(funcNameElmt, "Requires BODYFORCE or FIELDFORCE tag "
                     "specifying function name which prescribes body force.");
        }

        m_funcName = funcNameElmt->GetText();
        ASSERTL0(m_session->DefinesFunction(m_funcName),
                 "Function '" + m_funcName + "' not defined.");

        bool singleMode, halfMode;
        m_session->MatchSolverInfo("ModeType","SingleMode",singleMode,false);
        m_session->MatchSolverInfo("ModeType","HalfMode",  halfMode,  false);
        bool homogeneous = pFields[0]->GetExpType() == MultiRegions::e3DH1D ||
                           pFields[0]->GetExpType() == MultiRegions::e3DH2D;
        m_transform = (singleMode || halfMode || homogeneous);

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
                ::AllocateSharedPtr(m_session->GetInterpreter(),funcNameTime);

            m_hasTimeFcnScaling = true;
        }

        m_Forcing = Array<OneD, Array<OneD, NekDouble> > (m_NumVariable);
        for (int i = 0; i < m_NumVariable; ++i)
        {
            m_Forcing[i] = Array<OneD, NekDouble> (pFields[0]->GetTotPoints(), 0.0);
        }

        Array<OneD, Array<OneD, NekDouble> > tmp(pFields.num_elements());
        for (int i = 0; i < pFields.num_elements(); ++i)
        {
            tmp[i] = pFields[i]->GetPhys();
        }

        Update(pFields, tmp, 0.0);
    }


    void ForcingBody::Update(
            const Array< OneD, MultiRegions::ExpListSharedPtr > &pFields,
            const Array<OneD, Array<OneD, NekDouble> > &inarray,
            const NekDouble &time)
    {
        LibUtilities::EquationSharedPtr eqn = m_session->GetFunction(
            m_funcName, m_session->GetVariable(0));

        if (!boost::iequals(eqn->GetVlist(), "x y z t"))
        {
            // Coupled forcing
            int nq = pFields[0]->GetNpoints();
            Array<OneD, NekDouble> xc(nq), yc(nq), zc(nq), t(nq, time);
            std::string varstr = "x y z";
            std::vector<Array<OneD, const NekDouble>> fielddata = {
                xc, yc, zc, t};

            for (int i = 0; i < pFields.num_elements(); ++i)
            {
                varstr += " " + m_session->GetVariable(i);
                fielddata.push_back(inarray[i]);
            }

            // Evaluate function
            for (int i = 0; i < m_NumVariable; ++i)
            {
                m_session->GetFunction(m_funcName, m_session->GetVariable(i))->
                    Evaluate(fielddata, m_Forcing[i]);
            }

            return;
        }

        for (int i = 0; i < m_NumVariable; ++i)
        {
            std::string  s_FieldStr   = m_session->GetVariable(i);
            ASSERTL0(m_session->DefinesFunction(m_funcName, s_FieldStr),
                     "Variable '" + s_FieldStr + "' not defined.");
            GetFunction(pFields, m_session, m_funcName, true)->Evaluate(
                s_FieldStr, m_Forcing[i], time);
        }

        // If singleMode or halfMode, transform the forcing term to be in
        // physical space in the plane, but Fourier space in the homogeneous
        // direction
        if (m_transform)
        {
            for (int i = 0; i < m_NumVariable; ++i)
            {
                pFields[0]->HomogeneousFwdTrans(m_Forcing[i], m_Forcing[i]);
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
            Update(fields, inarray, time);

            for (int i = 0; i < m_NumVariable; i++)
            {
                Vmath::Vadd(outarray[i].num_elements(), outarray[i], 1,
                            m_Forcing[i], 1, outarray[i], 1);
            }
        }
    }

    void ForcingBody::v_Apply_coeff(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> > &inarray,
            Array<OneD, Array<OneD, NekDouble> > &outarray,
            const NekDouble &time)
    {
        int ncoeff = outarray[m_NumVariable-1].num_elements();
        Array<OneD, NekDouble> tmp(ncoeff, 0.0);

        if(m_hasTimeFcnScaling)
        {
            Array<OneD, NekDouble>  TimeFcn(1);

            for (int i = 0; i < m_NumVariable; i++)
            {
                EvaluateTimeFunction(time, m_timeFcnEqn, TimeFcn);

                fields[i]->FwdTrans(m_Forcing[i],tmp);

                Vmath::Svtvp(ncoeff, TimeFcn[0],
                             tmp, 1,
                             outarray[i],  1,
                             outarray[i],  1);
            }
        }
        else
        {
            Update(fields, inarray, time);

            for (int i = 0; i < m_NumVariable; i++)
            {
                fields[i]->FwdTrans(m_Forcing[i],tmp);

                Vmath::Vadd(ncoeff, outarray[i], 1,
                            tmp, 1, outarray[i], 1);
            }
        }
    }

}
}
