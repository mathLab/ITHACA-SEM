///////////////////////////////////////////////////////////////////////////////
//
// File: MappingIdentity.cpp
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
// Description: Empty mapping for identity transformation
//
///////////////////////////////////////////////////////////////////////////////

#include <GlobalMapping/MappingIdentity.h>
#include <MultiRegions/ExpList.h>

namespace Nektar
{
namespace GlobalMapping
{

    std::string MappingIdentity::className =
            GetMappingFactory().RegisterCreatorFunction("Identity",
                    MappingIdentity::create, "Identity mapping");

    /**
     *
     */
    MappingIdentity::MappingIdentity(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const Array<OneD, MultiRegions::ExpListSharedPtr>&   pFields)
        : Mapping(pSession, pFields)
    {
    }


    /**
     *
     */
    void MappingIdentity::v_InitObject(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
            const TiXmlElement                                *pMapping)
    {
        Mapping::v_InitObject(pFields, pMapping);  
        
        m_constantJacobian = true;
        m_timeDependent    = false;
        
    }

    void MappingIdentity::v_ContravarToCartesian(
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

    void MappingIdentity::v_CovarToCartesian(
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

    void MappingIdentity::v_ContravarFromCartesian(
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

    void MappingIdentity::v_CovarFromCartesian(
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

    void MappingIdentity::v_GetCartesianCoordinates(
                Array<OneD, NekDouble>               &out0,
                Array<OneD, NekDouble>               &out1,
                Array<OneD, NekDouble>               &out2)
    {
        m_fields[0]->GetCoords(out0, out1, out2);    
    }

    void MappingIdentity::v_GetJacobian(
        Array<OneD, NekDouble>               &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        Vmath::Fill(physTot, 1.0, outarray, 1);
    }
    
    void MappingIdentity::v_DotGradJacobian(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, NekDouble>               &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        
        Vmath::Zero(physTot, outarray, 1);   
    }

    void MappingIdentity::v_GetMetricTensor(
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
            Vmath::Sadd(physTot, 1.0, outarray[i+nvel*i], 1, 
                                        outarray[i+nvel*i], 1); 
        }            
    }

    void MappingIdentity::v_GetInvMetricTensor(
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
            Vmath::Sadd(physTot, 1.0, outarray[i+nvel*i], 1, 
                                        outarray[i+nvel*i], 1); 
        }            
    }
    
    void MappingIdentity::v_RaiseIndex(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        int nvel = m_nConvectiveFields;

        for (int i=0; i<nvel; i++)
        {
            outarray[i] = Array<OneD, NekDouble> (physTot, 0.0); 
        }
        // Copy
        for (int i=0; i<nvel; i++)
        {
            Vmath::Vcopy(physTot, inarray[i], 1, outarray[i], 1); 
        }            
    }
    
    void MappingIdentity::v_LowerIndex(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        int nvel = m_nConvectiveFields;

        for (int i=0; i<nvel; i++)
        {
            outarray[i] = Array<OneD, NekDouble> (physTot, 0.0); 
        }
        // Copy
        for (int i=0; i<nvel; i++)
        {
            Vmath::Vcopy(physTot, inarray[i], 1, outarray[i], 1); 
        }            
    }

    void MappingIdentity::v_ApplyChristoffelContravar(
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

    void MappingIdentity::v_ApplyChristoffelCovar(
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

    void MappingIdentity::v_UpdateMapping(const NekDouble time)
    {

    }

}
}
