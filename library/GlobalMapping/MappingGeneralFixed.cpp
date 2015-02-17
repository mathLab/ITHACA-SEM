///////////////////////////////////////////////////////////////////////////////
//
// File: MappingGeneralFixed.cpp
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
// Description: Mapping of the type X = X(x,y), Y = Y(x,y)
//
///////////////////////////////////////////////////////////////////////////////

#include <GlobalMapping/MappingGeneralFixed.h>
#include <MultiRegions/ExpList.h>
#include <iomanip> 

namespace Nektar
{
namespace GlobalMapping
{

    std::string MappingGeneralFixed::className =
            GetMappingFactory().RegisterCreatorFunction("GeneralFixed",
                    MappingGeneralFixed::create, "X = X(x,y,z), Y = Y(x,y,z), Z=Z(x,y,z)");

    /**
     *
     */
    MappingGeneralFixed::MappingGeneralFixed(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const Array<OneD, MultiRegions::ExpListSharedPtr>&   pFields)
        : Mapping(pSession, pFields)
    {
    }


    /**
     *
     */
    void MappingGeneralFixed::v_InitObject(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
            const TiXmlElement                                *pMapping)
    {
        Mapping::v_InitObject(pFields, pMapping); 
        
        int phystot         = pFields[0]->GetTotPoints();
        
        ASSERTL0(m_nConvectiveFields>=2,
                "General Mapping needs at least 2 velocity components.");   

        // Allocation of geometry memory
        m_GeometricInfo =  Array<OneD, Array<OneD, NekDouble> >(3);
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

        // Get coordinates in the domain
        Array<OneD, Array<OneD, NekDouble> > coords(3);
        for (int i = 0; i < 3; i++)
        {
            coords[i] = Array<OneD, NekDouble> (phystot);
        }
        m_fields[0]->GetCoords(coords[0], coords[1], coords[2]);
        
        std::string s_FieldStr; 
        for(int i = 0; i < m_nConvectiveFields; i++)
        {
            s_FieldStr = m_session->GetVariable(i);
            if ( m_session->DefinesFunction(m_funcName, s_FieldStr))
            {
                EvaluateFunction(pFields, m_session, s_FieldStr, m_GeometricInfo[i],
                                        m_funcName);
                if ( i==2 && m_fields[0]->GetExpType() == MultiRegions::e3DH1D)
                {
                    ASSERTL0 (false,
                            "3DH1D does not support mapping in the z-direction.");
                }
            }
            else
            {
                // This coordinate is not defined, so use (x^i)' = x^i
                Vmath::Vcopy(phystot, coords[i], 1, m_GeometricInfo[i], 1);
            }
        }
       
        CalculateMetricTerms();
        CalculateChristoffel();

    }

    void MappingGeneralFixed::v_ContravarToCartesian(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        int nvel = m_nConvectiveFields;

        for (int i=0; i<nvel; i++)
        {
            outarray[i] = Array<OneD, NekDouble> (physTot, 0.0); 
            for (int j=0; j<nvel; j++)
            {
                Vmath::Vvtvp(physTot, inarray[j], 1,
                                        m_deriv[i*nvel+j], 1,
                                        outarray[i], 1,
                                        outarray[i], 1);
            }
            
        }    
    }

    void MappingGeneralFixed::v_CovarToCartesian(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        int nvel = m_nConvectiveFields;

        for (int i=0; i<nvel; i++)
        {
            outarray[i] = Array<OneD, NekDouble> (physTot, 0.0); 
            for (int j=0; j<nvel; j++)
            {
                Vmath::Vvtvp(physTot, inarray[j], 1,
                                        m_invDeriv[i*nvel+j], 1,
                                        outarray[i], 1,
                                        outarray[i], 1);
            }
            
        }           
    }

    void MappingGeneralFixed::v_ContravarFromCartesian(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        int nvel = m_nConvectiveFields;
      
        for (int i=0; i<nvel; i++)
        {
            outarray[i] = Array<OneD, NekDouble> (physTot, 0.0); 
            for (int j=0; j<nvel; j++)
            {
                Vmath::Vvtvp(physTot, inarray[j], 1,
                                        m_invDeriv[j*nvel+i], 1,
                                        outarray[i], 1,
                                        outarray[i], 1);
            }
            
        }
    }

    void MappingGeneralFixed::v_CovarFromCartesian(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        int nvel = m_nConvectiveFields;
        
        for (int i=0; i<nvel; i++)
        {
            outarray[i] = Array<OneD, NekDouble> (physTot, 0.0); 
            for (int j=0; j<nvel; j++)
            {
                Vmath::Vvtvp(physTot, inarray[j], 1,
                                        m_deriv[j*nvel+i], 1,
                                        outarray[i], 1,
                                        outarray[i], 1);
            }
            
        } 
    }
    
    void MappingGeneralFixed::v_GetCartesianCoordinates(
                Array<OneD, NekDouble>               &out0,
                Array<OneD, NekDouble>               &out1,
                Array<OneD, NekDouble>               &out2)
    {
        int physTot = m_fields[0]->GetTotPoints();
        
        // x' = m_GeometricInfo[0]
        Vmath::Vcopy(physTot, m_GeometricInfo[0], 1, out0, 1);
        
        // y' = m_GeometricInfo[1]
        Vmath::Vcopy(physTot, m_GeometricInfo[1], 1, out1, 1);
        
        // z' = m_GeometricInfo[2]
        Vmath::Vcopy(physTot, m_GeometricInfo[2], 1, out2, 1);      
    }

    void MappingGeneralFixed::v_GetJacobian(
        Array<OneD, NekDouble>               &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        Vmath::Vcopy(physTot, m_jac, 1, outarray, 1);
    }

    void MappingGeneralFixed::v_GetMetricTensor(
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        int nvel = m_nConvectiveFields;

        for (int i=0; i<nvel; i++)
        {
            for (int j=0; j<nvel; j++)
            {
                outarray[i*nvel+j] = Array<OneD, NekDouble> (physTot, 0.0); 
                Vmath::Vcopy(physTot, m_metricTensor[i*nvel+j], 1,
                                        outarray[i*nvel+j], 1);
            }
            
        }
    }

    void MappingGeneralFixed::v_GetInvMetricTensor(
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        int nvel = m_nConvectiveFields;

        for (int i=0; i<nvel; i++)
        {
            for (int j=0; j<nvel; j++)
            {
                outarray[i*nvel+j] = Array<OneD, NekDouble> (physTot, 0.0); 
                Vmath::Vcopy(physTot, m_invMetricTensor[i*nvel+j], 1,
                                        outarray[i*nvel+j], 1);
            }
            
        }            
    }

    void MappingGeneralFixed::v_ApplyChristoffelContravar(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        int nvel = m_nConvectiveFields;
        
        // Calculate {i,jk} u^j
        for (int i = 0; i <nvel; i++)
        {
            for (int k = 0; k < nvel; k++)            
            {
                outarray[i*nvel+k] = Array<OneD, NekDouble>(physTot,0.0);
                for (int j = 0; j < nvel; j++)
                {
                    Vmath::Vvtvp(physTot, inarray[j], 1, 
                                        m_Christoffel[i*nvel*nvel+j*nvel+k], 1,
                                        outarray[i*nvel+k], 1,
                                        outarray[i*nvel+k], 1);
                }                
            }
        } 
        
    }

    void MappingGeneralFixed::v_ApplyChristoffelCovar(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        int nvel = m_nConvectiveFields;
        
        // Calculate {i,jk} u_i
        for (int j = 0; j < nvel; j++)
        {
            for (int k = 0; k < nvel; k++)
            {
                outarray[j*nvel+k] = Array<OneD, NekDouble>(physTot,0.0);
                for (int i = 0; i <nvel; i++)
                {
                    Vmath::Vvtvp(physTot, inarray[i], 1, 
                                        m_Christoffel[i*nvel*nvel+j*nvel+k], 1,
                                        outarray[j*nvel+k], 1,
                                        outarray[j*nvel+k], 1);
                }                
            }
        } 
    }

    bool MappingGeneralFixed::v_IsTimeDependent()
    {
        return false;
    }

    bool MappingGeneralFixed::v_HasConstantJacobian()
    {
        return false;
    }

    void MappingGeneralFixed::v_UpdateMapping(const NekDouble time)
    {

    }
    
    void MappingGeneralFixed::CalculateMetricTerms()
    {
        int physTot = m_fields[0]->GetTotPoints();
        int nvel = m_nConvectiveFields;
        
        // Set wavespace to false and store current value
        bool wavespace = m_fields[0]->GetWaveSpace();
        m_fields[0]->SetWaveSpace(false);
        
        // Allocate memory
        m_metricTensor =  Array<OneD, Array<OneD, NekDouble> >(nvel*nvel);
        m_invMetricTensor =  Array<OneD, Array<OneD, NekDouble> >(nvel*nvel);
        m_deriv =  Array<OneD, Array<OneD, NekDouble> >(nvel*nvel);
        m_invDeriv =  Array<OneD, Array<OneD, NekDouble> >(nvel*nvel);
        for (int i = 0; i < m_metricTensor.num_elements(); i++)
        {
            m_metricTensor[i] = Array<OneD, NekDouble>(physTot, 0.0);
            m_invMetricTensor[i] = Array<OneD, NekDouble>(physTot, 0.0);
            m_deriv[i] = Array<OneD, NekDouble>(physTot, 0.0);
            m_invDeriv[i] = Array<OneD, NekDouble>(physTot, 0.0);
        }
        m_jac = Array<OneD, NekDouble>(physTot, 0.0);
        
        // First, calculate derivatives of the mapping ->  dX^i/dx^j = c^i_j
        for( int i = 0; i<nvel; i++)
        {
            for( int j = 0; j<nvel; j++)
            {
                m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[j],
                                            m_GeometricInfo[i],
                                            m_deriv[i*nvel+j]);
            }
        }
        // In Homogeneous case, m_deriv(2,2) needs to be set to 1
        //    because differentiation in wavespace is not accurate for mappings
        if (m_fields[0]->GetExpType() == MultiRegions::e3DH1D)
        {
            Vmath::Fill(physTot, 1.0, m_deriv[2*nvel+2], 1);
        }
        
        // Now calculate the metric tensor -->   g_ij = sum_k  { c^k_i c^k_j }
        for( int i = 0; i<nvel; i++)
        {
            for( int j = 0; j<nvel; j++)
            {
                for( int k = 0; k<nvel; k++)
                {
                    Vmath::Vvtvp(physTot, m_deriv[k*nvel+i], 1,
                                            m_deriv[k*nvel+j], 1,
                                            m_metricTensor[i*nvel+j], 1,
                                            m_metricTensor[i*nvel+j], 1);
                }
            }
        } 
        
        // Put the adjoint of g in m_invMetricTensor
        switch (nvel)
        {
            case 1:
                Vmath::Fill (physTot,  1.0, m_invMetricTensor[0], 1);
                break;
            case 2:
                Vmath::Vcopy(physTot,       m_metricTensor[1*nvel+1],    1,
                                            m_invMetricTensor[0*nvel+0], 1);
                Vmath::Smul (physTot, -1.0, m_metricTensor[0*nvel+1],    1, 
                                            m_invMetricTensor[1*nvel+0], 1);
                Vmath::Smul (physTot, -1.0, m_metricTensor[1*nvel+0],    1, 
                                            m_invMetricTensor[0*nvel+1], 1);
                Vmath::Vcopy(physTot,       m_metricTensor[0*nvel+0],    1,
                                            m_invMetricTensor[1*nvel+1], 1);
                break;
            case 3:
            {
                int a, b, c, d, e, i, j;

                // Compute g^{ij} by computing Cofactors(g_ij)^T
                for (i = 0; i < nvel; ++i)
                {
                    for (j = 0; j < nvel; ++j)
                    {
                        a = ((i+1)%nvel) * nvel + ((j+1)%nvel);
                        b = ((i+1)%nvel) * nvel + ((j+2)%nvel);
                        c = ((i+2)%nvel) * nvel + ((j+1)%nvel);
                        d = ((i+2)%nvel) * nvel + ((j+2)%nvel);
                        e = i*nvel + j;
                        // a*d - b*c
                        Vmath::Vmul(physTot,    m_metricTensor[b], 1,
                                                m_metricTensor[c], 1,
                                                m_invMetricTensor[e], 1);
                        Vmath::Vvtvm(physTot,   m_metricTensor[a], 1,
                                                m_metricTensor[d], 1,
                                                m_invMetricTensor[e], 1,
                                                m_invMetricTensor[e], 1);
                    }
                }
                break;
            }
        }
        
        // Compute g = det(g_{ij}) (= Jacobian squared) and store
        // temporarily in m_jac.
        for (int i = 0; i < nvel; ++i)
        {
            Vmath::Vvtvp(physTot, m_metricTensor[i], 1, 
                                 m_invMetricTensor[i*nvel], 1,
                                 m_jac, 1, m_jac, 1);
        }

        // Calculate g^ij (the inverse of g_ij) by dividing by jac
        for (int i = 0; i < nvel*nvel; ++i)
        {
            Vmath::Vdiv(physTot, m_invMetricTensor[i], 1, m_jac, 1,
                                    m_invMetricTensor[i], 1);
        }

        // Compute the Jacobian = sqrt(g)
        Vmath::Vsqrt(physTot, m_jac, 1, m_jac, 1);        
        
        // Calculate the derivatives of the inverse transformation
        //          c'^j_i = dx^j/dX^i = sum_k {g^jk c^i_k}
        for (int i = 0; i < nvel; ++i)
        {
            for (int j = 0; j < nvel; ++j)
            {
                for (int k = 0; k < nvel; ++k)
                {
                    Vmath::Vvtvp(physTot, m_deriv[i*nvel+k],           1,
                                          m_invMetricTensor[j*nvel+k], 1,
                                          m_invDeriv[i*nvel+j],        1,
                                          m_invDeriv[i*nvel+j],        1);
                }
            }
        }

        // Restore value of wavespace
        m_fields[0]->SetWaveSpace(wavespace);
    }
    
    void MappingGeneralFixed::CalculateChristoffel()
    {
        int physTot = m_fields[0]->GetTotPoints();
        int nvel = m_nConvectiveFields;

        Array<OneD, Array<OneD, NekDouble> >   gradG(nvel*nvel*nvel);
        Array<OneD, Array<OneD, NekDouble> >   tmp(nvel*nvel*nvel);
        m_Christoffel = Array<OneD, Array<OneD, NekDouble> > (nvel*nvel*nvel);
        // Allocate memory
        for (int i = 0; i < gradG.num_elements(); i++)
        {
            gradG[i] = Array<OneD, NekDouble>(physTot, 0.0);
            tmp[i] = Array<OneD, NekDouble>(physTot, 0.0);
            m_Christoffel[i] = Array<OneD, NekDouble>(physTot, 0.0);
        }

        // Set wavespace to false and store current value
        bool waveSpace = m_fields[0]->GetWaveSpace();
        m_fields[0]->SetWaveSpace(false);
        
        //Calculate gradients of g_ij
        for (int i = 0; i <nvel; i++)
        {
            for(int j=0; j<nvel; j++)
            {
                for(int k=0; k<nvel; k++)
                {
                    m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[k],
                                            m_metricTensor[i*nvel+j],
                                            gradG[i*nvel*nvel + j*nvel + k]);
                }
            }
        }
        
        // Calculate tmp[p,j,k] = 1/2( gradG[pj,k]+ gradG[pk,j]-gradG[jk,p])
        for (int p = 0; p <nvel; p++)
        {
            for (int j = 0; j < nvel; j++)
            {
                for (int k = 0; k < nvel; k++)
                {
                    Vmath::Vadd(physTot, gradG[p*nvel*nvel + j*nvel + k], 1,
                                         gradG[p*nvel*nvel + k*nvel + j], 1,
                                         tmp[p*nvel*nvel + j*nvel + k], 1);
                    Vmath::Vsub(physTot, tmp[p*nvel*nvel + j*nvel + k], 1,
                                         gradG[j*nvel*nvel + k*nvel + p], 1,
                                         tmp[p*nvel*nvel + j*nvel + k], 1);
                    Vmath::Smul(physTot, 0.5, tmp[p*nvel*nvel + j*nvel + k], 1,
                                               tmp[p*nvel*nvel + j*nvel + k], 1);
                }
            }
        }
        
        // Calculate Christoffel symbols = g^ip tmp[p,j,k]
        for (int i = 0; i <nvel; i++)
        {
            for (int j = 0; j < nvel; j++)
            {
                for (int k = 0; k < nvel; k++)
                {
                    for (int p = 0; p < nvel; p++)
                    {
                        Vmath::Vvtvp(physTot, m_invMetricTensor[i*nvel+p], 1,
                                        tmp[p*nvel*nvel + j*nvel + k], 1,
                                        m_Christoffel[i*nvel*nvel+j*nvel+k], 1,
                                        m_Christoffel[i*nvel*nvel+j*nvel+k], 1);
                    }
                }                
            }
        }        
        // Restore wavespace
        m_fields[0]->SetWaveSpace(waveSpace);
        
    }    


}
}
