///////////////////////////////////////////////////////////////////////////////
//
// File: Forcing.cpp
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
// Description: Abstract base class for forcing terms.
//
///////////////////////////////////////////////////////////////////////////////

#include <FieldUtils/Interpolator.h>
#include <SolverUtils/Forcing/Forcing.h>

using namespace std;

namespace Nektar
{
    namespace SolverUtils
    {
        ForcingFactory& GetForcingFactory()
        {
            static ForcingFactory instance;
            return instance;
        }

        Forcing::Forcing(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const std::weak_ptr<EquationSystem>      &pEquation)
                : m_session(pSession), m_equ(pEquation)
        {

        }

        void Forcing::InitObject(
                const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
                const unsigned int& pNumForcingFields,
                const TiXmlElement* pForce)
        {
            v_InitObject(pFields, pNumForcingFields, pForce);
        }


        /**
         * @param   fields      Expansion lists corresponding to input arrays
         * @param   inarray     u^n from previous timestep
         * @param   outarray    output array to append forcing to
         */
        void Forcing::Apply(
                const Array<OneD, MultiRegions::ExpListSharedPtr>& fields,
                const Array<OneD, Array<OneD, NekDouble> >&        inarray,
                Array<OneD, Array<OneD, NekDouble> >&              outarray,
                const NekDouble&                                   time)
        {
            v_Apply(fields, inarray, outarray, time);
        }


        /**
         *
         */
        vector<ForcingSharedPtr> Forcing::Load(
                    const LibUtilities::SessionReaderSharedPtr &pSession,
                    const std::weak_ptr<EquationSystem>      &pEquation,
                    const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
                    const unsigned int& pNumForcingFields)
        {
            vector<ForcingSharedPtr> vForceList;

            if (!pSession->DefinesElement("Nektar/Forcing"))
            {
                return vForceList;
            }

            TiXmlElement* vForcing = pSession->GetElement("Nektar/Forcing");
            if (vForcing)
            {
                unsigned int vNumForcingFields = pNumForcingFields;
                if (!pNumForcingFields)
                {
                    vNumForcingFields = pFields.size();
                }

                TiXmlElement* vForce = vForcing->FirstChildElement("FORCE");
                while (vForce)
                {
                    string vType = vForce->Attribute("TYPE");

                    vForceList.push_back(GetForcingFactory().CreateInstance(
                                        vType, pSession, pEquation, pFields,
                                        vNumForcingFields, vForce));
                    vForce = vForce->NextSiblingElement("FORCE");
                }
            }
            return vForceList;
        }

        Nektar::Array<OneD, Array<OneD, NekDouble> > &Forcing::UpdateForces()
        {
            return m_Forcing;
        }

        const Nektar::Array<OneD, const Array<OneD, NekDouble> > &Forcing::GetForces()
        {
            return m_Forcing;
        }

        void Forcing::EvaluateTimeFunction(
                LibUtilities::SessionReaderSharedPtr              pSession,
                std::string                                       pFieldName,
                Array<OneD, NekDouble>&                           pArray,
                const std::string&                                pFunctionName,
                NekDouble                                         pTime)
        {
            ASSERTL0(pSession->DefinesFunction(pFunctionName),
                     "Function '" + pFunctionName + "' does not exist.");

            LibUtilities::EquationSharedPtr ffunc =
                pSession->GetFunction(pFunctionName, pFieldName);

            EvaluateTimeFunction(pTime,ffunc,pArray);
        }


        void Forcing::EvaluateTimeFunction(
               const NekDouble                          pTime,
               const LibUtilities::EquationSharedPtr&   pEqn,
               Array<OneD, NekDouble>&                  pArray)
        {
            // dummy array of zero pts.
            Array<OneD, NekDouble> x0(pArray.size(),0.0);

            pEqn->Evaluate(x0, x0, x0, pTime, pArray);
        }

        SessionFunctionSharedPtr Forcing::GetFunction(
                const Array<OneD, MultiRegions::ExpListSharedPtr>  &pFields,
                const LibUtilities::SessionReaderSharedPtr         &pSession,
                std::string                                         pName,
                bool                                                pCache)
        {
            if (pCache)
            {
                if ((m_sessionFunctions.find(pName) == m_sessionFunctions.end())
                    || (m_sessionFunctions[pName]->GetSession() != pSession)
                    || (m_sessionFunctions[pName]->GetExpansion() != pFields[0])
                )
                {
                    m_sessionFunctions[pName] =
                        MemoryManager<SessionFunction>::AllocateSharedPtr(
                            pSession, pFields[0], pName, pCache);
                }

                return m_sessionFunctions[pName];
            }
            else
            {
                return SessionFunctionSharedPtr(new SessionFunction(pSession, pFields[0], pName, pCache));
            }
        }

    }
}
