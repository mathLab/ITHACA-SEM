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

#include <SolverUtils/Mapping/MappingXYofT.h>
#include <MultiRegions/ExpList.h>

namespace Nektar
{
namespace SolverUtils
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
        : Mapping(pSession, pFields)
    {
    }


    /**
     *
     */
    void MappingXYofT::v_InitObject(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
            const TiXmlElement                                *pMapping)
    {
        Mapping::v_InitObject(pFields, pMapping);
        
        int phystot         = pFields[0]->GetTotPoints();
        
        ASSERTL0(m_nConvectiveFields>=2,"Mapping X = x + f(t), Y = y+g(t) needs 2 velocity components.");
       
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

    void MappingXYofT::v_ContravarToCartesian(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        
        // U1 = u1
        Vmath::Vcopy(physTot, inarray[0], 1, outarray[0], 1);
        
        // U2 = u2
        Vmath::Vcopy(physTot, inarray[1], 1, outarray[1], 1);
        
        // U3 = u3
        if (m_nConvectiveFields ==3)
        {
            Vmath::Vcopy(physTot, inarray[2], 1, outarray[2], 1);
        }        
    }

    void MappingXYofT::v_CovarToCartesian(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        
        // U1 = u1
        Vmath::Vcopy(physTot, inarray[0], 1, outarray[0], 1);
        
        // U2 = u2
        Vmath::Vcopy(physTot, inarray[1], 1, outarray[1], 1);
        
        // U3 = u3
        if (m_nConvectiveFields ==3)
        {
            Vmath::Vcopy(physTot, inarray[2], 1, outarray[2], 1);
        } 
    }

    void MappingXYofT::v_ContravarFromCartesian(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        
        // U1 = u1
        Vmath::Vcopy(physTot, inarray[0], 1, outarray[0], 1);
        
        // U2 = u2
        Vmath::Vcopy(physTot, inarray[1], 1, outarray[1], 1);
        
        // U3 = u3
        if (m_nConvectiveFields ==3)
        {
            Vmath::Vcopy(physTot, inarray[2], 1, outarray[2], 1);
        }         
    }

    void MappingXYofT::v_CovarFromCartesian(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        
        // U1 = u1
        Vmath::Vcopy(physTot, inarray[0], 1, outarray[0], 1);
        
        // U2 = u2
        Vmath::Vcopy(physTot, inarray[1], 1, outarray[1], 1);
        
        // U3 = u3
        if (m_nConvectiveFields ==3)
        {
            Vmath::Vcopy(physTot, inarray[2], 1, outarray[2], 1);
        } 
    }

    void MappingXYofT::v_CoordinatesToCartesian(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray)
    {
        int npoints = outarray[0].num_elements();
        
        // x' = f(x,z)
        LibUtilities::EquationSharedPtr ffunc =
                    m_session->GetFunction(m_funcName, m_session->GetVariable(0));
                
        ffunc->Evaluate(inarray[0], inarray[1], inarray[2], 0.0, outarray[0]);       
        // y' = g(y,z)
        ffunc = m_session->GetFunction(m_funcName, m_session->GetVariable(1));
                
        ffunc->Evaluate(inarray[0], inarray[1], inarray[2], 0.0, outarray[1]);
        // z' = z
        // U3 = u3
        if (m_nConvectiveFields ==3)
        {
            Vmath::Vcopy(npoints, inarray[2], 1, outarray[2], 1);
        }        
    }

    void MappingXYofT::v_GetJacobian(
        Array<OneD, NekDouble>               &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        Vmath::Fill(physTot, 1.0, outarray, 1);
    }
    
    void MappingXYofT::v_DotGradJacobian(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, NekDouble>               &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        
        Vmath::Zero(physTot, outarray, 1);   
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

    void MappingXYofT::v_GetMetricTensor(
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        int nvel = m_nConvectiveFields;

        for (int i=0; i<nvel*nvel; i++)
        {
            outarray[i] = Array<OneD, NekDouble> (physTot, 0.0); 
        }
        // Fill diagonal with 1.0
        for (int i=0; i<nvel; i++)
        {
            Vmath::Sadd(physTot, 1.0, outarray[i+nvel*i], 1, outarray[i+nvel*i], 1); 
        }            
    }

    void MappingXYofT::v_GetInvMetricTensor(
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        int nvel = m_nConvectiveFields;

        for (int i=0; i<nvel*nvel; i++)
        {
            outarray[i] = Array<OneD, NekDouble> (physTot, 0.0); 
        }
        // Fill diagonal with 1.0
        for (int i=0; i<nvel; i++)
        {
            Vmath::Sadd(physTot, 1.0, outarray[i+nvel*i], 1, outarray[i+nvel*i], 1); 
        }            
    }
    
    void MappingXYofT::v_RaiseIndex(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        int nvel = m_nConvectiveFields;

        for (int i=0; i<nvel*nvel; i++)
        {
            outarray[i] = Array<OneD, NekDouble> (physTot, 0.0); 
        }
        // Copy
        for (int i=0; i<nvel; i++)
        {
            Vmath::Vcopy(physTot, inarray[i], 1, outarray[i], 1); 
        }            
    }
    
    void MappingXYofT::v_LowerIndex(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        int nvel = m_nConvectiveFields;

        for (int i=0; i<nvel*nvel; i++)
        {
            outarray[i] = Array<OneD, NekDouble> (physTot, 0.0); 
        }
        // Copy
        for (int i=0; i<nvel; i++)
        {
            Vmath::Vcopy(physTot, inarray[i], 1, outarray[i], 1); 
        }            
    }

    void MappingXYofT::v_ApplyChristoffelContravar(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        int nvel = m_nConvectiveFields;
        
        for (int i = 0; i< nvel; i++)
        {
            for (int j = 0; j< nvel; j++)
            {
                outarray[i*nvel+j] = Array<OneD, NekDouble>(physTot,0.0);
            }            
        }        
    }

    void MappingXYofT::v_ApplyChristoffelCovar(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        int nvel = m_nConvectiveFields;
        
        for (int i = 0; i< nvel; i++)
        {
            for (int j = 0; j< nvel; j++)
            {
                outarray[i*nvel+j] = Array<OneD, NekDouble>(physTot,0.0);
            }            
        }
    }

    void MappingXYofT::v_IncNSAdvectionCorrection(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        int nvel = m_nConvectiveFields;

        for (int i=0; i<nvel; i++)
        {
            outarray[i] = Array<OneD, NekDouble> (physTot, 0.0); 
        }
    }

    void MappingXYofT::v_IncNSPressureCorrection(
        const Array<OneD, NekDouble>                      &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        int nvel = m_nConvectiveFields;

        for (int i=0; i<nvel; i++)
        {
            outarray[i] = Array<OneD, NekDouble> (physTot, 0.0); 
        }
    }

    void MappingXYofT::v_IncNSViscousCorrection(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        int nvel = m_nConvectiveFields;

        for (int i=0; i<nvel; i++)
        {
            outarray[i] = Array<OneD, NekDouble> (physTot, 0.0); 
        }           
    }    

    bool MappingXYofT::v_IsTimeDependent()
    {
        return true;
    }

    bool MappingXYofT::v_HasConstantJacobian()
    {
        return true;
    }

    void MappingXYofT::v_UpdateMapping()
    {
        std::string s_XFieldStr = m_session->GetVariable(0);
        std::string s_YFieldStr = m_session->GetVariable(1);
        // Evaluate x-functions --> GeometricInfo 0-1
        EvaluateFunction(m_fields, m_session, s_XFieldStr, m_GeometricInfo[0],
                m_funcName); 
        EvaluateFunction(m_fields, m_session, s_XFieldStr, m_GeometricInfo[1],
                m_velFuncName);
        
        // Evaluate y-functions --> GeometricInfo 2-3
        EvaluateFunction(m_fields, m_session, s_YFieldStr, m_GeometricInfo[2],
                m_funcName);
        EvaluateFunction(m_fields, m_session, s_YFieldStr, m_GeometricInfo[3],
                m_velFuncName);
    }

}
}
