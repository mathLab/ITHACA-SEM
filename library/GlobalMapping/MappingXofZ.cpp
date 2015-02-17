///////////////////////////////////////////////////////////////////////////////
//
// File: MappingXofZ.cpp
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
// Description: Mapping of the type X = x + f(z)
//
///////////////////////////////////////////////////////////////////////////////

#include <GlobalMapping/MappingXofZ.h>
#include <MultiRegions/ExpList.h>

namespace Nektar
{
namespace GlobalMapping
{

    std::string MappingXofZ::className =
            GetMappingFactory().RegisterCreatorFunction("XofZ",
                    MappingXofZ::create, "X = x + f(z)");

    /**
     *
     */
    MappingXofZ::MappingXofZ(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const Array<OneD, MultiRegions::ExpListSharedPtr>&   pFields)
        : Mapping(pSession, pFields)
    {
    }


    /**
     *
     */
    void MappingXofZ::v_InitObject(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
            const TiXmlElement                                *pMapping)
    {
        Mapping::v_InitObject(pFields, pMapping);     
                
        int phystot         = pFields[0]->GetTotPoints();
        
        ASSERTL0(m_nConvectiveFields==3,
                    "Mapping X = x + f(z) needs 3 velocity components.");
        
        // Allocation of geometry memory
        m_GeometricInfo =  Array<OneD, Array<OneD, NekDouble> >(6);
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

        std::string s_FieldStr = m_session->GetVariable(0);
        ASSERTL0(m_session->DefinesFunction(m_funcName, s_FieldStr),
                "Variable '" + s_FieldStr + "' not defined.");
        
        bool waveSpace = pFields[0]->GetWaveSpace();
        pFields[0]->SetWaveSpace(false);
        
        EvaluateFunction(pFields, m_session, s_FieldStr, m_GeometricInfo[0],
                m_funcName);

        // Calculate derivatives of transformation
        for(int i = 1; i < 4; i++)
        {
            pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],m_GeometricInfo[i-1],m_GeometricInfo[i]);
        }
        // m_GeometricInfo[4] = fz^2
        Vmath::Vmul(phystot,m_GeometricInfo[1],1,m_GeometricInfo[1],1,m_GeometricInfo[4],1);
        // m_GeometricInfo[5] = fz*fzz
        Vmath::Vmul(phystot,m_GeometricInfo[1],1,m_GeometricInfo[2],1,m_GeometricInfo[5],1);

        pFields[0]->SetWaveSpace(waveSpace);

    }

    void MappingXofZ::v_ContravarToCartesian(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        Array<OneD, NekDouble> wk(physTot, 0.0);
        
        // U1 = u1 + fz*u3
        Vmath::Vmul(physTot, inarray[2], 1, m_GeometricInfo[1], 1, wk, 1);
        Vmath::Vadd(physTot, wk, 1, inarray[0], 1, outarray[0],1);
        
        // U2 = u2
        Vmath::Vcopy(physTot, inarray[1], 1, outarray[1], 1);
        
        // U3 = u3
        Vmath::Vcopy(physTot, inarray[2], 1, outarray[2], 1);
    }

    void MappingXofZ::v_CovarToCartesian(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        Array<OneD, NekDouble> wk(physTot, 0.0);
        
        // U1 = u1
        Vmath::Vcopy(physTot, inarray[0], 1, outarray[0], 1);
        
        // U2 = u2
        Vmath::Vcopy(physTot, inarray[1], 1, outarray[1], 1);
        
        // U3 = u3 - fz*u1
        Vmath::Vmul(physTot, m_GeometricInfo[1], 1, inarray[0], 1, wk, 1);
        Vmath::Vsub(physTot, inarray[2], 1, wk, 1, outarray[2], 1);
    }

    void MappingXofZ::v_ContravarFromCartesian(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        Array<OneD, NekDouble> wk(physTot, 0.0);
        
        // U1 = u1 - fz * u3
        Vmath::Vmul(physTot, m_GeometricInfo[1], 1, inarray[2], 1, wk, 1);
        Vmath::Vsub(physTot, inarray[0], 1, wk, 1, outarray[0], 1);        
        
        // U2 = u2
        Vmath::Vcopy(physTot, inarray[1], 1, outarray[1], 1);
        
        // U3 = u3
        Vmath::Vcopy(physTot, inarray[2], 1, outarray[2], 1);        
    }

    void MappingXofZ::v_CovarFromCartesian(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        
        // U1 = u1
        Vmath::Vcopy(physTot, inarray[0], 1, outarray[0], 1);
        
        // U2 = u2
        Vmath::Vcopy(physTot, inarray[1], 1, outarray[1], 1);
        
        // U3 = u3 + fz*u1
        Vmath::Vmul(physTot, m_GeometricInfo[1], 1, 
                             inarray[0], 1, outarray[2], 1);
        Vmath::Vadd(physTot, inarray[2], 1, outarray[2], 1, outarray[2], 1);
    }
    
    void MappingXofZ::v_GetCartesianCoordinates(
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
        
        // y' = y
        Vmath::Vcopy(physTot, x1, 1, out1, 1);    
        
        // z' = z
        Vmath::Vcopy(physTot, x2, 1, out2, 1);        
    }

    void MappingXofZ::v_GetJacobian(
        Array<OneD, NekDouble>               &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        Vmath::Fill(physTot, 1.0, outarray, 1);
    }
    
    void MappingXofZ::v_DotGradJacobian(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, NekDouble>               &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        
        Vmath::Zero(physTot, outarray, 1);   
    }

    void MappingXofZ::v_GetMetricTensor(
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
    {
            int physTot = m_fields[0]->GetTotPoints();
            int nvel = m_nConvectiveFields;
            Array<OneD, NekDouble> wk(physTot, 0.0);
            
            for (int i=0; i<nvel*nvel; i++)
            {
                outarray[i] = Array<OneD, NekDouble> (physTot, 0.0); 
            }
            // Fill diagonal with 1.0
            for (int i=0; i<nvel; i++)
            {
                Vmath::Sadd(physTot, 1.0, outarray[i+nvel*i], 1, 
                                        outarray[i+nvel*i], 1); 
            }            
            
            // G_{13} and G_{31} = fz
            Vmath::Vcopy(physTot, m_GeometricInfo[1], 1, outarray[0*nvel+2], 1);
            Vmath::Vcopy(physTot, m_GeometricInfo[1], 1, outarray[2*nvel+0], 1);
            
            // G^{33} = (1+fz^2)
            Vmath::Vmul(physTot, m_GeometricInfo[1], 1,
                                m_GeometricInfo[1], 1, wk, 1); // fz^2
            Vmath::Vadd(physTot, wk, 1, outarray[2*nvel+2], 1, 
                                        outarray[2*nvel+2], 1);
    }

    void MappingXofZ::v_GetInvMetricTensor(
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
    {
            int physTot = m_fields[0]->GetTotPoints();
            int nvel = m_nConvectiveFields;
            Array<OneD, NekDouble> wk(physTot, 0.0);
            
            for (int i=0; i<nvel*nvel; i++)
            {
                outarray[i] = Array<OneD, NekDouble> (physTot, 0.0); 
            }
            // Fill diagonal with 1.0
            for (int i=0; i<nvel; i++)
            {
                Vmath::Sadd(physTot, 1.0, outarray[i+nvel*i], 1, 
                                            outarray[i+nvel*i], 1); 
            }            
            
            // G^{13} and G^{31} = -fz
            Vmath::Vcopy(physTot, m_GeometricInfo[1],1,wk,1); // fz
            Vmath::Neg(physTot, wk, 1);
            Vmath::Vcopy(physTot, wk, 1, outarray[0*nvel+2], 1);
            Vmath::Vcopy(physTot, wk, 1, outarray[2*nvel+0], 1);

            // G^{11} = (1+fz^2)
            Vmath::Vmul(physTot, m_GeometricInfo[1], 1,
                                m_GeometricInfo[1], 1, wk, 1); // fz^2
            Vmath::Sadd(physTot, 1.0, wk, 1, outarray[0*nvel+0], 1);
    }

    void MappingXofZ::v_LowerIndex(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        Array<OneD, NekDouble> wk(physTot, 0.0);
        
        Vmath::Vmul(physTot, m_GeometricInfo[1], 1, inarray[2], 1, 
                                                    outarray[0], 1); //in[2]*fz
        Vmath::Vmul(physTot, m_GeometricInfo[1], 1, inarray[0], 1, 
                                                    outarray[2], 1); //in[0]*fz
        
        Vmath::Vadd(physTot, outarray[0], 1, inarray[0], 1, 
                                                    outarray[0], 1); // out[0] = in[0] + in[2] * fz
        
        Vmath::Vcopy(physTot, inarray[1], 1, outarray[1], 1); // out[1] = in[1]]

        Vmath::Vmul(physTot, m_GeometricInfo[1], 1, m_GeometricInfo[1], 1, 
                                                        wk, 1); // fz^2
        Vmath::Sadd(physTot, 1.0, wk, 1, wk, 1); // 1+fz^2
        Vmath::Vmul(physTot, wk, 1, inarray[2],1, wk, 1); // (1+fz^2)*in[2]
        Vmath::Vadd(physTot, wk, 1, outarray[2],1, outarray[2], 1); // out[2] = fz*in[0] + (1+fz^2)*in[2]
    }

    void MappingXofZ::v_RaiseIndex(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        Array<OneD, NekDouble> wk(physTot, 0.0);
        Array<OneD, NekDouble> wk_2(physTot, 0.0);
        
        Vmath::Vmul(physTot, m_GeometricInfo[1], 1, inarray[2], 1, 
                                                    outarray[0], 1); // in[2]*fz
        Vmath::Vmul(physTot, m_GeometricInfo[1], 1, inarray[0], 1, 
                                                    outarray[2], 1); //in[0]*fz       
        Vmath::Vsub(physTot, inarray[2], 1, outarray[2], 1, 
                                                    outarray[2], 1); // out[2] = in[2] - in[0] * fz
        
        Vmath::Vcopy(physTot, inarray[1], 1, outarray[1], 1); // out[1] = in[1]]

        Vmath::Vmul(physTot, m_GeometricInfo[1], 1, m_GeometricInfo[1], 1, 
                                                        wk, 1); // fz^2
        Vmath::Sadd(physTot, 1.0, wk, 1, wk, 1); // 1+fz^2
        Vmath::Vmul(physTot, wk, 1, inarray[0],1, wk, 1); // in[0]*(1+fz^2)
        Vmath::Vsub(physTot, wk, 1, outarray[0], 1, outarray[0], 1); // out[0] = in[0]*(1+fz^2)- in[2] * fz
    }

    void MappingXofZ::v_ApplyChristoffelContravar(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        int nvel = m_nConvectiveFields;
        Array<OneD, NekDouble> wk(physTot, 0.0);
        
        for (int i = 0; i< nvel; i++)
        {
            for (int j = 0; j< nvel; j++)
            {
                outarray[i*nvel+j] = Array<OneD, NekDouble>(physTot,0.0);
            }            
        }
        
        // Calculate non-zero terms  
        
        // outarray(0,2) = U3 * fzz
        Vmath::Vmul(physTot,m_GeometricInfo[2],1,inarray[2],1,
                                            outarray[0*nvel+2],1); // U1 * fzz
        
    }

    void MappingXofZ::v_ApplyChristoffelCovar(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        int nvel = m_nConvectiveFields;
        Array<OneD, NekDouble> wk(physTot, 0.0);
        
        for (int i = 0; i< nvel; i++)
        {
            for (int j = 0; j< nvel; j++)
            {
                outarray[i*nvel+j] = Array<OneD, NekDouble>(physTot,0.0);
            }            
        }
        
        // Calculate non-zero terms
        
        // outarray(2,2) = U1 * fzz
        Vmath::Vmul(physTot,m_GeometricInfo[2],1,inarray[0],1,
                                        outarray[2*nvel+2],1); // U1 * fzz 
    }   

    bool MappingXofZ::v_IsTimeDependent()
    {
        return false;
    }

    bool MappingXofZ::v_HasConstantJacobian()
    {
        return true;
    }

    void MappingXofZ::v_UpdateMapping(const NekDouble time)
    {

    }


}
}
