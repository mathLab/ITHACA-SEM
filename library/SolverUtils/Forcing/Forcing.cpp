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
// Description: Abstract base class for forcing terms.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Forcing/Forcing.h>

namespace Nektar
{
    namespace SolverUtils
    {
        ForcingFactory& GetForcingFactory()
        {
            typedef Loki::SingletonHolder<ForcingFactory,
                                          Loki::CreateUsingNew,
                                          Loki::NoDestroy,
                                          Loki::SingleThreaded> Type;
            return Type::Instance();
        }

        Forcing::Forcing(const LibUtilities::SessionReaderSharedPtr& pSession)
                : m_session(pSession)
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
                            const LibUtilities::SessionReaderSharedPtr& pSession,
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
                    vNumForcingFields = pFields.num_elements();
                }

                TiXmlElement* vForce = vForcing->FirstChildElement("FORCE");
                while (vForce)
                {
                    string vType = vForce->Attribute("TYPE");

                    vForceList.push_back(GetForcingFactory().CreateInstance(
                                            vType, pSession, pFields,
                                            vNumForcingFields, vForce));
                    vForce = vForce->NextSiblingElement("FORCE");
                }
            }
            return vForceList;
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
            
            Array<OneD, NekDouble> x0(1,0.0);
            Array<OneD, NekDouble> x1(1,0.0);
            Array<OneD, NekDouble> x2(1,0.0);
         
            ffunc->Evaluate(x0, x1, x2, pTime, pArray);
        }


        void Forcing::EvaluateFunction(
                Array<OneD, MultiRegions::ExpListSharedPtr>       pFields,
                LibUtilities::SessionReaderSharedPtr              pSession,
                std::string                                       pFieldName,
                Array<OneD, NekDouble>&                           pArray,
                const std::string&                                pFunctionName,
                NekDouble                                         pTime)
        {
            ASSERTL0(pSession->DefinesFunction(pFunctionName),
                     "Function '" + pFunctionName + "' does not exist.");

            unsigned int nq = pFields[0]->GetNpoints();
            if (pArray.num_elements() != nq)
            {
                pArray = Array<OneD, NekDouble> (nq);
            }

            LibUtilities::FunctionType vType;
            vType = pSession->GetFunctionType(pFunctionName, pFieldName);
            if (vType == LibUtilities::eFunctionTypeExpression)
            {
                Array<OneD, NekDouble> x0(nq);
                Array<OneD, NekDouble> x1(nq);
                Array<OneD, NekDouble> x2(nq);
                
                pFields[0]->GetCoords(x0, x1, x2);
                LibUtilities::EquationSharedPtr ffunc =
                    pSession->GetFunction(pFunctionName, pFieldName);
                
                ffunc->Evaluate(x0, x1, x2, pTime, pArray);
            }
            else if (vType == LibUtilities::eFunctionTypeFile)
            {
                std::string filename = pSession->GetFunctionFilename(
                                                    pFunctionName,
                                                    pFieldName);

                std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef;
                std::vector<std::vector<NekDouble> > FieldData;
                Array<OneD, NekDouble> vCoeffs(pFields[0]->GetNcoeffs());
                Vmath::Zero(vCoeffs.num_elements(), vCoeffs, 1);

                LibUtilities::FieldIOSharedPtr fld =
                    MemoryManager<LibUtilities::FieldIO>::AllocateSharedPtr(m_session->GetComm());
                fld->Import(filename, FieldDef, FieldData);

                int idx = -1;
                for (int i = 0; i < FieldDef.size(); ++i)
                {
                    for (int j = 0; j < FieldDef[i]->m_fields.size(); ++j)
                    {
                        if (FieldDef[i]->m_fields[j] == pFieldName)
                        {
                            idx = j;
                        }
                    }

                    if (idx >= 0)
                    {
                        pFields[0]->ExtractDataToCoeffs(
                                                    FieldDef[i],
                                                    FieldData[i],
                                                    FieldDef[i]->m_fields[idx],
                                                    vCoeffs);
                    }
                    else
                    {
                        cout << "Field " + pFieldName + " not found." << endl;
                    }
                }
                pFields[0]->BwdTrans_IterPerExp(vCoeffs, pArray);
            }
        }

    }
}
