///////////////////////////////////////////////////////////////////////////////
//
// File: MappingTranslation.cpp
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
// Description: Trivial mapping for translation transformation
//
///////////////////////////////////////////////////////////////////////////////

#include <GlobalMapping/MappingTranslation.h>
#include <MultiRegions/ExpList.h>

namespace Nektar
{
namespace GlobalMapping
{

    std::string MappingTranslation::className =
            GetMappingFactory().RegisterCreatorFunction("Translation",
                    MappingTranslation::create, "Translation mapping (X_i = x_i + constant)");

    /**
     *
     */
    MappingTranslation::MappingTranslation(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const Array<OneD, MultiRegions::ExpListSharedPtr>&   pFields)
        : Mapping(pSession, pFields)
    {
    }


    /**
     *
     */
    void MappingTranslation::v_InitObject(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
            const TiXmlElement                                *pMapping)
    {   
        m_constantJacobian = true;
        // When there is no Mapping object defined, use the identity
        //      transformation as a default
        if (m_session->DefinesElement("Nektar/Mapping"))
        {
            Mapping::v_InitObject(pFields, pMapping);
        }
        else
        {            
            int phystot         = pFields[0]->GetTotPoints();
            
            m_timeDependent    = false;

            m_coords    = Array<OneD, Array<OneD, NekDouble> > (3);
            m_coordsVel = Array<OneD, Array<OneD, NekDouble> > (3);
            for (int i = 0; i < 3; i++)
            {
                m_coords[i]    = Array<OneD, NekDouble> (phystot);
                m_coordsVel[i] = Array<OneD, NekDouble> (phystot, 0.0);
            }

            m_fields[0]->GetCoords(m_coords[0], m_coords[1], m_coords[2]);
        }
        
    }

    void MappingTranslation::v_ContravarToCartesian(
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

    void MappingTranslation::v_CovarToCartesian(
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

    void MappingTranslation::v_ContravarFromCartesian(
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

    void MappingTranslation::v_CovarFromCartesian(
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

    void MappingTranslation::v_GetJacobian(
        Array<OneD, NekDouble>               &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        Vmath::Fill(physTot, 1.0, outarray, 1);
    }
    
    void MappingTranslation::v_DotGradJacobian(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, NekDouble>               &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        
        Vmath::Zero(physTot, outarray, 1);   
    }

    void MappingTranslation::v_GetMetricTensor(
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

    void MappingTranslation::v_GetInvMetricTensor(
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
    
    void MappingTranslation::v_RaiseIndex(
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
    
    void MappingTranslation::v_LowerIndex(
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

    void MappingTranslation::v_ApplyChristoffelContravar(
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

    void MappingTranslation::v_ApplyChristoffelCovar(
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

    void MappingTranslation::v_UpdateGeomInfo()
    {

    }

}
}
