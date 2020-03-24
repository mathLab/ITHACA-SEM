///////////////////////////////////////////////////////////////////////////////
//
// File: MappingXofXZ.cpp
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
// Description: Mapping of the type X = X(x,z)
//
///////////////////////////////////////////////////////////////////////////////

#include <GlobalMapping/MappingXofXZ.h>
#include <MultiRegions/ExpList.h>

namespace Nektar
{
namespace GlobalMapping
{

std::string MappingXofXZ::className =
        GetMappingFactory().RegisterCreatorFunction("XofXZ",
                MappingXofXZ::create, "X = X(x,z)");

/**
 * @class MappingXofXZ
 * This class implements a mapping defined by a transformation of the type
 * \f[ \bar{x} = \bar{x}(x,z) \f]
 * \f[ \bar{y} = y \f]
 * \f[ \bar{z} = z \f]
 * where \f$(\bar{x},\bar{y},\bar{z})\f$ are the Cartesian (physical)
 * coordinates and \f$(x,y,z)\f$ are the transformed (computational)
 *  coordinates.
 */
MappingXofXZ::MappingXofXZ(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const Array<OneD, MultiRegions::ExpListSharedPtr>&   pFields)
    : Mapping(pSession, pFields)
{
}


/**
 *
 */
void MappingXofXZ::v_InitObject(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const TiXmlElement                                *pMapping)
{
    Mapping::v_InitObject(pFields, pMapping);

    m_constantJacobian = false;

    ASSERTL0(m_nConvectiveFields==3,
            "Mapping X = X(x,z) needs 3 velocity components.");
}

void MappingXofXZ::v_ContravarToCartesian(
    const Array<OneD, Array<OneD, NekDouble> >        &inarray,
    Array<OneD, Array<OneD, NekDouble> >              &outarray)
{
    int physTot = m_fields[0]->GetTotPoints();

    // U1 = fx*u1 + fz*u3
    Vmath::Vmul(physTot, m_GeometricInfo[0], 1, inarray[0], 1,
                                                outarray[0], 1);
    Vmath::Vvtvp(physTot, m_GeometricInfo[1], 1, inarray[2], 1,
                            outarray[0], 1, outarray[0],1);

    // U2 = u2
    Vmath::Vcopy(physTot, inarray[1], 1, outarray[1], 1);

    // U3 = u3
    Vmath::Vcopy(physTot, inarray[2], 1, outarray[2], 1);
}

void MappingXofXZ::v_CovarToCartesian(
    const Array<OneD, Array<OneD, NekDouble> >        &inarray,
    Array<OneD, Array<OneD, NekDouble> >              &outarray)
{
    int physTot = m_fields[0]->GetTotPoints();
    Array<OneD, NekDouble> wk(physTot, 0.0);

    // U1 = u1/fx
    Vmath::Vdiv(physTot, inarray[0], 1, m_GeometricInfo[0], 1,
                                                outarray[0], 1);

    // U2 = u2
    Vmath::Vcopy(physTot, inarray[1], 1, outarray[1], 1);

    // U3 = u3 - fz/fx*u1
    Vmath::Vdiv(physTot, m_GeometricInfo[1], 1,
                         m_GeometricInfo[0], 1, wk, 1);
    Vmath::Vmul(physTot, wk, 1, inarray[0], 1, wk, 1);
    Vmath::Vsub(physTot, inarray[2], 1, wk, 1, outarray[2], 1);
}

void MappingXofXZ::v_ContravarFromCartesian(
    const Array<OneD, Array<OneD, NekDouble> >        &inarray,
    Array<OneD, Array<OneD, NekDouble> >              &outarray)
{
    int physTot = m_fields[0]->GetTotPoints();
    Array<OneD, NekDouble> wk(physTot, 0.0);

    // U1 = u1/fx - fz/fx * u3
    Vmath::Vdiv(physTot, inarray[0], 1,
                         m_GeometricInfo[0], 1, outarray[0], 1);
    Vmath::Vdiv(physTot, m_GeometricInfo[1], 1,
                         m_GeometricInfo[0], 1, wk, 1);
    Vmath::Vmul(physTot, wk, 1, inarray[2], 1, wk, 1);
    Vmath::Vsub(physTot, outarray[0], 1, wk, 1, outarray[0], 1);

    // U2 = u2
    Vmath::Vcopy(physTot, inarray[1], 1, outarray[1], 1);

    // U3 = u3
    Vmath::Vcopy(physTot, inarray[2], 1, outarray[2], 1);
}

void MappingXofXZ::v_CovarFromCartesian(
    const Array<OneD, Array<OneD, NekDouble> >        &inarray,
    Array<OneD, Array<OneD, NekDouble> >              &outarray)
{
    int physTot = m_fields[0]->GetTotPoints();

    // U1 = u1*fx
    Vmath::Vmul(physTot, inarray[0], 1, m_GeometricInfo[0], 1,
                                                    outarray[0], 1);

    // U2 = u2
    Vmath::Vcopy(physTot, inarray[1], 1, outarray[1], 1);

    // U3 = u3 + fz*u1
    Vmath::Vmul(physTot, m_GeometricInfo[1], 1,
                         inarray[0], 1, outarray[2], 1);
    Vmath::Vadd(physTot, inarray[2], 1, outarray[2], 1, outarray[2], 1);
}

void MappingXofXZ::v_GetJacobian(
    Array<OneD, NekDouble>               &outarray)
{
    int physTot = m_fields[0]->GetTotPoints();
    Vmath::Vcopy(physTot, m_GeometricInfo[0], 1, outarray, 1);
}

void MappingXofXZ::v_DotGradJacobian(
    const Array<OneD, Array<OneD, NekDouble> >        &inarray,
    Array<OneD, NekDouble>                            &outarray)
{
    int physTot = m_fields[0]->GetTotPoints();

    Vmath::Vmul(physTot, m_GeometricInfo[2], 1, inarray[0], 1, outarray, 1);
    Vmath::Vvtvp(physTot, m_GeometricInfo[3], 1, inarray[2], 1,
                            outarray, 1, outarray,1);
}

void MappingXofXZ::v_GetMetricTensor(
    Array<OneD, Array<OneD, NekDouble> >              &outarray)
{
    int physTot = m_fields[0]->GetTotPoints();
    int nvel = m_nConvectiveFields;
    Array<OneD, NekDouble> wk(physTot, 0.0);

    for (int i=0; i<nvel*nvel; i++)
    {
        outarray[i] = Array<OneD, NekDouble> (physTot, 0.0);
    }
    // Fill G^{22} and G^{33} with 1.0
    for (int i=1; i<nvel; i++)
    {
        Vmath::Sadd(physTot, 1.0, outarray[i+nvel*i], 1,
                                    outarray[i+nvel*i], 1);
    }

    // G_{13} and G_{31} = fz*fx
    Vmath::Vmul(physTot,m_GeometricInfo[1],1,
                        m_GeometricInfo[0],1,wk,1); // fz*fx
    Vmath::Vcopy(physTot, wk, 1, outarray[0*nvel+2], 1);
    Vmath::Vcopy(physTot, wk, 1, outarray[2*nvel+0], 1);

    // G^{11} = (fx^2)
    Vmath::Vmul(physTot, m_GeometricInfo[0], 1,
                         m_GeometricInfo[0], 1, outarray[0*nvel+0], 1);

    // G^{33} = (1+fz^2)
    Vmath::Vmul(physTot, m_GeometricInfo[1], 1,
                        m_GeometricInfo[1], 1, wk, 1); // fz^2
    Vmath::Vadd(physTot, wk, 1, outarray[2*nvel+2], 1,
                                    outarray[2*nvel+2], 1);
}

void MappingXofXZ::v_GetInvMetricTensor(
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

    // G^{13} and G^{31} = -fz/fx
    Vmath::Vdiv(physTot,m_GeometricInfo[1],1,
                        m_GeometricInfo[0],1,wk,1); // fz/fx
    Vmath::Neg(physTot, wk, 1);
    Vmath::Vcopy(physTot, wk, 1, outarray[0*nvel+2], 1);
    Vmath::Vcopy(physTot, wk, 1, outarray[2*nvel+0], 1);

    // G^{11} = (1+fz^2)/(fx^2)
    Vmath::Vmul(physTot, m_GeometricInfo[1], 1,
                        m_GeometricInfo[1], 1, wk, 1); // fz^2
    Vmath::Vadd(physTot, wk, 1, outarray[0*nvel+0], 1,
                                    outarray[0*nvel+0], 1);

    Vmath::Vmul(physTot, m_GeometricInfo[0], 1,
                         m_GeometricInfo[0], 1, wk, 1); // fx^2
    Vmath::Vdiv(physTot, outarray[0*nvel+0], 1, wk,1,
                                    outarray[0*nvel+0], 1);
}

void MappingXofXZ::v_LowerIndex(
    const Array<OneD, Array<OneD, NekDouble> >        &inarray,
    Array<OneD, Array<OneD, NekDouble> >              &outarray)
{
    int physTot = m_fields[0]->GetTotPoints();
    Array<OneD, NekDouble> wk(physTot, 0.0);

    // out[0] = in[0]*fx^2 + in[2] * fz*fx
    Vmath::Vmul(physTot,m_GeometricInfo[1],1,m_GeometricInfo[0],1,
                                                wk,1); // fz*fx
    Vmath::Vmul(physTot, wk, 1, inarray[2], 1, outarray[0], 1); //in[2]*fz*fx
    Vmath::Vmul(physTot, wk, 1, inarray[0], 1, outarray[2], 1); //in[0]*fz*fx

    Vmath::Vmul(physTot, m_GeometricInfo[0], 1, m_GeometricInfo[0], 1,
                                                wk, 1); //fx^2
    Vmath::Vmul(physTot, wk, 1, inarray[0], 1, wk, 1);  //in[0]*fx^2

    Vmath::Vadd(physTot, outarray[0], 1, wk, 1, outarray[0], 1);

    // out[1] = in[1]
    Vmath::Vcopy(physTot, inarray[1], 1, outarray[1], 1);

    // out[2] = fx*fz*in[0] + (1+fz^2)*in[2]
    Vmath::Vmul(physTot, m_GeometricInfo[1], 1, m_GeometricInfo[1], 1,
                                                wk, 1); // fz^2
    Vmath::Sadd(physTot, 1.0, wk, 1, wk, 1); // 1+fz^2
    Vmath::Vmul(physTot, wk, 1, inarray[2],1, wk, 1); // (1+fz^2)*in[2]

    Vmath::Vadd(physTot, wk, 1, outarray[2],1, outarray[2], 1);
}

void MappingXofXZ::v_RaiseIndex(
    const Array<OneD, Array<OneD, NekDouble> >        &inarray,
    Array<OneD, Array<OneD, NekDouble> >              &outarray)
{
    int physTot = m_fields[0]->GetTotPoints();
    Array<OneD, NekDouble> wk(physTot, 0.0);
    Array<OneD, NekDouble> wk_2(physTot, 0.0);

    // out[2] = in[2] - in[0] * fz/fx
    Vmath::Vdiv(physTot,m_GeometricInfo[1],1,m_GeometricInfo[0],1,
                                                    wk,1);
    Vmath::Vmul(physTot, wk, 1, inarray[0], 1, outarray[2], 1);
    Vmath::Vsub(physTot, inarray[2], 1, outarray[2], 1, outarray[2], 1);

    // out[0] = in[0]*(1+fz^2)/(fx^2) - in[2] * fz/fx
    Vmath::Vmul(physTot, wk, 1, inarray[2], 1, outarray[0], 1);
    Vmath::Vmul(physTot, m_GeometricInfo[1], 1, m_GeometricInfo[1], 1,
                                                        wk, 1);
    Vmath::Sadd(physTot, 1.0, wk, 1, wk, 1);
    Vmath::Vmul(physTot, m_GeometricInfo[0], 1, m_GeometricInfo[0], 1,
                                                        wk_2, 1);
    Vmath::Vdiv(physTot, wk, 1, wk_2,1, wk, 1);
    Vmath::Vmul(physTot, wk, 1, inarray[0],1, wk, 1);
    Vmath::Vsub(physTot, wk, 1, outarray[0], 1, outarray[0], 1);

    // out[1] = in[1]
    Vmath::Vcopy(physTot, inarray[1], 1, outarray[1], 1);


}

void MappingXofXZ::v_ApplyChristoffelContravar(
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

    // outarray(0,0) = U1 * fxx/fx + U3 * fxz/fx
    Vmath::Vdiv(physTot,m_GeometricInfo[2],1,m_GeometricInfo[0],1,wk,1);
    Vmath::Vmul(physTot,wk,1,inarray[0],1,outarray[0*nvel+0],1);
    Vmath::Vdiv(physTot,m_GeometricInfo[3],1,m_GeometricInfo[0],1,wk,1);
    Vmath::Vvtvp(physTot,wk,1,inarray[2],1,outarray[0*nvel+0],1,
                                            outarray[0*nvel+0],1);

    // outarray(0,2) = U1 * fxz/fx + U3 * fzz/fx
    Vmath::Vmul(physTot,wk,1,inarray[0],1,outarray[0*nvel+2],1);
    Vmath::Vdiv(physTot,m_GeometricInfo[4],1,m_GeometricInfo[0],1,wk,1);
    Vmath::Vvtvp(physTot,wk,1,inarray[2],1,outarray[0*nvel+2],1,
                                            outarray[0*nvel+2],1);

}

void MappingXofXZ::v_ApplyChristoffelCovar(
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

    // outarray(0,0) = U1 * fxx/fx
    Vmath::Vdiv(physTot,m_GeometricInfo[2],1,m_GeometricInfo[0],1,wk,1);
    Vmath::Vmul(physTot,wk,1,inarray[0],1,outarray[0*nvel+0],1);

    //outarray(0,2) = outarray(2,0) = U1 * fxz/fx
    Vmath::Vdiv(physTot,m_GeometricInfo[3],1,m_GeometricInfo[0],1,wk,1);
    Vmath::Vmul(physTot,wk,1,inarray[0],1,outarray[0*nvel+2],1);
    Vmath::Vcopy(physTot,outarray[0*nvel+2],1,outarray[2*nvel+0],1);

    // outarray(2,2) = U1 * fzz/fx
    Vmath::Vdiv(physTot,m_GeometricInfo[4],1,m_GeometricInfo[0],1,wk,1);
    Vmath::Vmul(physTot,wk,1,inarray[0],1,outarray[2*nvel+2],1);
}

void MappingXofXZ::v_UpdateGeomInfo()
{
    int phystot         = m_fields[0]->GetTotPoints();
    // Allocation of geometry memory
    m_GeometricInfo =  Array<OneD, Array<OneD, NekDouble> >(5);
    for (int i = 0; i < m_GeometricInfo.size(); i++)
    {
        m_GeometricInfo[i] = Array<OneD, NekDouble>(phystot, 0.0);
    }

    bool waveSpace = m_fields[0]->GetWaveSpace();
    m_fields[0]->SetWaveSpace(false);

    // Calculate derivatives of transformation
    m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0],
                m_coords[0], m_GeometricInfo[0]);  //f_x
    m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                m_coords[0], m_GeometricInfo[1]);  //f_z

    m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0],
                m_GeometricInfo[0], m_GeometricInfo[2]);  //f_xx
    m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                m_GeometricInfo[0], m_GeometricInfo[3]);  //f_xz
    m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                m_GeometricInfo[1], m_GeometricInfo[4]);  //f_zz

    m_fields[0]->SetWaveSpace(waveSpace);
}


}
}
