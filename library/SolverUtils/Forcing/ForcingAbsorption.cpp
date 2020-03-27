///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingAbsorption.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2016 Kilian Lackhove
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
// Description: Absorption layer forcing
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <SolverUtils/Forcing/ForcingAbsorption.h>

#include <LibUtilities/BasicUtils/Equation.h>
#include <LibUtilities/BasicUtils/ParseUtils.h>

using namespace std;

namespace Nektar
{
namespace SolverUtils
{

    std::string ForcingAbsorption::className = GetForcingFactory().
                                RegisterCreatorFunction("Absorption",
                                                        ForcingAbsorption::create,
                                                        "Forcing Absorption");

    ForcingAbsorption::ForcingAbsorption(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const std::weak_ptr<EquationSystem>      &pEquation)
            : Forcing(pSession, pEquation),
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
        int npts      = pFields[0]->GetTotPoints();

        CalcAbsorption(pFields, pForce);

        m_Forcing = Array<OneD, Array<OneD, NekDouble> >(m_NumVariable);
        for (int i = 0; i < m_NumVariable; ++i)
        {
            m_Forcing[i] = Array<OneD, NekDouble>(npts, 0.0);
        }

        const TiXmlElement* funcNameElmt = pForce->FirstChildElement("REFFLOW");
        if (funcNameElmt)
        {
            string funcName = funcNameElmt->GetText();
            ASSERTL0(m_session->DefinesFunction(funcName),
                     "Function '" + funcName + "' not defined.");
            m_Refflow = Array<OneD, Array<OneD, NekDouble> > (m_NumVariable);
            for (int i = 0; i < m_NumVariable; ++i)
            {
                std::string s_FieldStr = m_session->GetVariable(i);
                ASSERTL0(m_session->DefinesFunction(funcName, s_FieldStr),
                         "Variable '" + s_FieldStr + "' not defined.");
                m_Refflow[i] = Array<OneD, NekDouble> (npts, 0.0);
                GetFunction(pFields, m_session, funcName)->Evaluate(s_FieldStr, m_Refflow[i]);
            }
            m_hasRefFlow = true;
        }

        funcNameElmt = pForce->FirstChildElement("REFFLOWTIME");
        if (funcNameElmt)
        {
            m_funcNameTime = funcNameElmt->GetText();
            m_hasRefFlowTime = true;
        }
    }

    void ForcingAbsorption::CalcAbsorption(
        const Nektar::Array<Nektar::OneD, Nektar::MultiRegions::ExpListSharedPtr>
            &pFields,
        const TiXmlElement *pForce)
    {
        const TiXmlElement *funcNameElmt = pForce->FirstChildElement("COEFF");
        ASSERTL0(funcNameElmt,
                "Requires COEFF tag, specifying function "
                "name which prescribes absorption layer coefficient.");
        string funcName = funcNameElmt->GetText();
        ASSERTL0(m_session->DefinesFunction(funcName),
                "Function '" + funcName + "' not defined.");

        int npts = pFields[0]->GetTotPoints();

        m_Absorption = Array<OneD, Array<OneD, NekDouble> >(m_NumVariable);
        for (int i = 0; i < m_NumVariable; ++i)
        {
            m_Absorption[i] = Array<OneD, NekDouble>(npts, 0.0);
        }

        funcNameElmt = pForce->FirstChildElement("BOUNDARYREGIONS");
        if (funcNameElmt)
        {
            ASSERTL0(ParseUtils::GenerateVector(funcNameElmt->GetText(),
                                                    m_bRegions),
                    "Unable to process list of BOUNDARYREGIONS in Absorption "
                    "Forcing: " +
                        std::string(funcNameElmt->GetText()));

            // alter m_bRegions so that it contains the boundaryRegions of this rank
            std::vector<unsigned int>  localBRegions;
            SpatialDomains::BoundaryConditions bcs(m_session, pFields[0]->GetGraph());
            SpatialDomains::BoundaryRegionCollection regions = bcs.GetBoundaryRegions();
            SpatialDomains::BoundaryRegionCollection::iterator it1;
            int n = 0;
            for (it1 = regions.begin(); it1 != regions.end(); ++it1)
            {
                if (std::find(m_bRegions.begin(), m_bRegions.end(), it1->first) != m_bRegions.end())
                {
                    localBRegions.push_back(n);
                }
                n++;
            }
            m_bRegions = localBRegions;

            if (m_bRegions.size() == 0)
            {
                return;
            }

            std::vector<Array<OneD, const NekDouble> > points;

            Array<OneD, Array<OneD, NekDouble> > x(3);
            for (int i = 0; i < 3; i++)
            {
                x[i] = Array<OneD, NekDouble>(npts, 0.0);
            }
            pFields[0]->GetCoords(x[0], x[1], x[2]);
            for (int i = 0; i < 3; i++)
            {
                points.push_back(x[i]);
            }

            Array<OneD, NekDouble> t(npts, 0.0);
            points.push_back(t);

            Array<OneD, NekDouble> r(npts, 0.0);
            std::vector<unsigned int>::iterator it;
            std::vector<BPointPair> inPoints;
            Array<OneD, Array<OneD, NekDouble> > b(3);
            for (it = m_bRegions.begin(); it != m_bRegions.end(); ++it)
            {
                int bpts = pFields[0]->GetBndCondExpansions()[*it]->GetNpoints();
                for (int i = 0; i < 3; i++)
                {
                    b[i] = Array<OneD, NekDouble>(bpts, 0.0);
                }
                pFields[0]->GetBndCondExpansions()[*it]->GetCoords(
                    b[0], b[1], b[2]);
                for (int i = 0;
                    i < pFields[0]->GetBndCondExpansions()[*it]->GetNpoints();
                    ++i)
                {
                    inPoints.push_back(
                        BPointPair(BPoint(b[0][i], b[1][i], b[2][i]), i));
                }
            }
            m_rtree = MemoryManager<BRTree>::AllocateSharedPtr();
            m_rtree->insert(inPoints.begin(), inPoints.end());

            for (int i = 0; i < npts; ++i)
            {
                std::vector<BPointPair> result;
                BPoint sPoint(x[0][i], x[1][i], x[2][i]);
                m_rtree->query(bgi::nearest(sPoint, 1), std::back_inserter(result));
                r[i] = bg::distance(sPoint, result[0].first);
            }
            points.push_back(r);

            std::string s_FieldStr;
            for (int i = 0; i < m_NumVariable; ++i)
            {
                s_FieldStr = m_session->GetVariable(i);
                ASSERTL0(m_session->DefinesFunction(funcName, s_FieldStr),
                        "Variable '" + s_FieldStr + "' not defined.");

                LibUtilities::EquationSharedPtr ffunc =
                    m_session->GetFunction(funcName, s_FieldStr);
                ASSERTL0(ffunc->GetVlist() == "x y z t r",
                        "EVARS of " + funcName + " must be 'r'");

                ffunc->Evaluate(points, m_Absorption[i]);
            }
        }
        else
        {
            for (int i = 0; i < m_NumVariable; ++i)
            {
                std::string s_FieldStr = m_session->GetVariable(i);
                GetFunction(pFields, m_session, funcName)->Evaluate(s_FieldStr, m_Absorption[i]);
            }
        }
    }

    void ForcingAbsorption::v_Apply(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> > &inarray,
            Array<OneD, Array<OneD, NekDouble> > &outarray,
            const NekDouble &time)
    {
        boost::ignore_unused(fields);

        int nq = m_Forcing[0].size();

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
