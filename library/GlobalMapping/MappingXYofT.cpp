///////////////////////////////////////////////////////////////////////////////
//
// File: MappingXYofT.cpp
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
// Description: Mapping of the type X = x + f(t), Y = y + g(t)
//
///////////////////////////////////////////////////////////////////////////////

#include <GlobalMapping/MappingXYofT.h>
#include <MultiRegions/ExpList.h>

namespace Nektar
{
namespace GlobalMapping
{

    std::string MappingXYofT::className =
            GetMappingFactory().RegisterCreatorFunction("XYofT",
                    MappingXYofT::create, "X = x + f(t), Y = y +g(t)");

    /**
     *
     */
    MappingXYofT::MappingXYofT(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const Array<OneD, MultiRegions::ExpListSharedPtr>&   pFields)
        : MappingIdentity(pSession, pFields)
    {
    }


    /**
     *
     */
    void MappingXYofT::v_InitObject(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
            const TiXmlElement                                *pMapping)
    {
        MappingIdentity::v_InitObject(pFields, pMapping);
        
        int phystot         = pFields[0]->GetTotPoints();
        
        ASSERTL0(m_nConvectiveFields>=2,
               "Mapping X = x + f(t), Y = y+g(t) needs 2 velocity components.");
       
        // Allocation of geometry memory
        m_GeometricInfo =  Array<OneD, Array<OneD, NekDouble> >(12);
        for (int i = 0; i < m_GeometricInfo.num_elements(); i++)
        {
            m_GeometricInfo[i] = Array<OneD, NekDouble>(phystot, 0.0);
        }

        // Read and evaluate function
        const TiXmlElement* funcNameElmt;
        funcNameElmt = pMapping->FirstChildElement("COORDS");
        ASSERTL0(funcNameElmt, "Requires COORDS tag, specifying function "
                "name which prescribes mapping.");

        m_funcName = funcNameElmt->GetText();
        ASSERTL0(m_session->DefinesFunction(m_funcName),
                "Function '" + m_funcName + "' not defined.");

        std::string s_XFieldStr = m_session->GetVariable(0);
        ASSERTL0(m_session->DefinesFunction(m_funcName, s_XFieldStr),
                "Variable '" + s_XFieldStr + "' not defined.");
        
        std::string s_YFieldStr = m_session->GetVariable(1);
        ASSERTL0(m_session->DefinesFunction(m_funcName, s_YFieldStr),
                "Variable '" + s_YFieldStr + "' not defined.");     
        
        funcNameElmt = pMapping->FirstChildElement("VEL");
        ASSERTL0(funcNameElmt, "Requires VEL tag, specifying function "
                "name which prescribes mapping velocity.");

        m_velFuncName = funcNameElmt->GetText();
        ASSERTL0(m_session->DefinesFunction(m_velFuncName),
                "Function '" + m_velFuncName + "' not defined.");

        ASSERTL0(m_session->DefinesFunction(m_velFuncName, s_XFieldStr),
                "Variable '" + s_XFieldStr + "' not defined.");
        
        ASSERTL0(m_session->DefinesFunction(m_velFuncName, s_YFieldStr),
                "Variable '" + s_YFieldStr + "' not defined."); 
        
        // Evaluate x-functions --> GeometricInfo 0-1
        EvaluateFunction(pFields, m_session, s_XFieldStr, m_GeometricInfo[0],
                m_funcName); 
        EvaluateFunction(pFields, m_session, s_XFieldStr, m_GeometricInfo[1],
                m_velFuncName);
        
        // Evaluate y-functions --> GeometricInfo 2-3
        EvaluateFunction(pFields, m_session, s_YFieldStr, m_GeometricInfo[2],
                m_funcName);
        EvaluateFunction(pFields, m_session, s_YFieldStr, m_GeometricInfo[3],
                m_velFuncName);

    }

    void MappingXYofT::v_GetCartesianCoordinates(
                Array<OneD, NekDouble>               &out0,
                Array<OneD, NekDouble>               &out1,
                Array<OneD, NekDouble>               &out2)
    {
        int physTot = m_fields[0]->GetTotPoints();
        
        Array<OneD, NekDouble> x0(physTot);
        Array<OneD, NekDouble> x1(physTot);
        Array<OneD, NekDouble> x2(physTot);

        m_fields[0]->GetCoords(x0, x1, x2);
        
        // x' = m_GeometricInfo[0]
        Vmath::Vcopy(physTot, m_GeometricInfo[0], 1, out0, 1);
        
        // y' = m_GeometricInfo[2]
        Vmath::Vcopy(physTot, m_GeometricInfo[2], 1, out1, 1);
        
        // z' = z
        Vmath::Vcopy(physTot, x2, 1, out2, 1);        
    }
    
    void MappingXYofT::v_GetCoordVelocity(
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        int nvel = m_nConvectiveFields;

        for (int i=0; i<nvel; i++)
        {
            outarray[i] = Array<OneD, NekDouble> (physTot, 0.0); 
        }

        Vmath::Vcopy(physTot, m_GeometricInfo[1], 1, outarray[0], 1);
        Vmath::Vcopy(physTot, m_GeometricInfo[3], 1, outarray[1], 1);
    }

    bool MappingXYofT::v_IsTimeDependent()
    {
        return true;
    }

    void MappingXYofT::v_UpdateMapping(const NekDouble time)
    {
        std::string s_XFieldStr = m_session->GetVariable(0);
        std::string s_YFieldStr = m_session->GetVariable(1);
        // Evaluate x-functions --> GeometricInfo 0-1
        EvaluateFunction(m_fields, m_session, s_XFieldStr, m_GeometricInfo[0],
                m_funcName, time); 
        EvaluateFunction(m_fields, m_session, s_XFieldStr, m_GeometricInfo[1],
                m_velFuncName, time);
        
        // Evaluate y-functions --> GeometricInfo 2-3
        EvaluateFunction(m_fields, m_session, s_YFieldStr, m_GeometricInfo[2],
                m_funcName, time);
        EvaluateFunction(m_fields, m_session, s_YFieldStr, m_GeometricInfo[3],
                m_velFuncName, time);
    }

}
}
