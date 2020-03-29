///////////////////////////////////////////////////////////////////////////////
//
// File: MappingXYofXY.cpp
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
// Description: Mapping of the type X = X(x,y), Y = Y(x,y)
//
///////////////////////////////////////////////////////////////////////////////

#include <GlobalMapping/MappingXYofXY.h>
#include <MultiRegions/ExpList.h>

namespace Nektar
{
namespace GlobalMapping
{

std::string MappingXYofXY::className =
        GetMappingFactory().RegisterCreatorFunction("XYofXY",
                MappingXYofXY::create, "X = X(x,y), Y = Y(x,y)");

/**
 * @class MappingXYofXY
 * This class implements a mapping defined by the transformation
 * \f[ \bar{x} = \bar{x}(x,y) \f]
 * \f[ \bar{y} = \bar{y}(x,y) \f]
 * \f[ \bar{z} = z \f]
 * where \f$(\bar{x},\bar{y},\bar{z})\f$ are the Cartesian (physical)
 * coordinates and \f$(x,y,z)\f$ are the transformed (computational)
 *  coordinates.
 */
MappingXYofXY::MappingXYofXY(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const Array<OneD, MultiRegions::ExpListSharedPtr>&   pFields)
    : Mapping(pSession, pFields)
{
}


/**
 *
 */
void MappingXYofXY::v_InitObject(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const TiXmlElement                                *pMapping)
{
    Mapping::v_InitObject(pFields, pMapping);

    m_constantJacobian = false;

    ASSERTL0(m_nConvectiveFields>=2,
            "Mapping X = X(x,y), Y = Y(x,y) needs 2 velocity components.");
}

void MappingXYofXY::v_ContravarToCartesian(
    const Array<OneD, Array<OneD, NekDouble> >        &inarray,
    Array<OneD, Array<OneD, NekDouble> >              &outarray)
{
    int physTot = m_fields[0]->GetTotPoints();

    // U1 = fx*u1 + fy*u2
    Vmath::Vmul(physTot, m_GeometricInfo[0], 1, inarray[0], 1,
                                                outarray[0], 1);
    Vmath::Vvtvp(physTot, m_GeometricInfo[1], 1, inarray[1], 1,
                            outarray[0], 1, outarray[0],1);

    // U2 = gx*u1+gy*u2
    Vmath::Vmul(physTot, m_GeometricInfo[2], 1, inarray[0], 1,
                                                outarray[1], 1);
    Vmath::Vvtvp(physTot, m_GeometricInfo[3], 1, inarray[1], 1,
                            outarray[1], 1, outarray[1],1);

    // U3 = u3
    if (m_nConvectiveFields ==3)
    {
        Vmath::Vcopy(physTot, inarray[2], 1, outarray[2], 1);
    }
}

void MappingXYofXY::v_CovarToCartesian(
    const Array<OneD, Array<OneD, NekDouble> >        &inarray,
    Array<OneD, Array<OneD, NekDouble> >              &outarray)
{
    int physTot = m_fields[0]->GetTotPoints();
    Array<OneD, NekDouble> wk(physTot, 0.0);

    // U1 = [gy*u1-gx*u2]/(fx*gy-gx*fy)
    Vmath::Vmul(physTot, inarray[1], 1, m_GeometricInfo[2], 1,
                                        outarray[0], 1);
    Vmath::Vvtvm(physTot, inarray[0], 1, m_GeometricInfo[3], 1,
                                        outarray[0], 1,
                                        outarray[0], 1);
    Vmath::Vdiv(physTot, outarray[0], 1, m_GeometricInfo[4], 1,
                                        outarray[0], 1);

    // U2 = [fx*u2 - fy*u1]/(fx*gy-gx*fy)
    Vmath::Vmul(physTot, inarray[0], 1, m_GeometricInfo[1], 1,
                                        outarray[1], 1);
    Vmath::Vvtvm(physTot, inarray[1], 1, m_GeometricInfo[0], 1,
                                        outarray[1], 1,
                                        outarray[1], 1);
    Vmath::Vdiv(physTot, outarray[1], 1, m_GeometricInfo[4], 1,
                                        outarray[1], 1);

    // U3 = u3
    if (m_nConvectiveFields ==3)
    {
        Vmath::Vcopy(physTot, inarray[2], 1, outarray[2], 1);
    }
}

void MappingXYofXY::v_ContravarFromCartesian(
    const Array<OneD, Array<OneD, NekDouble> >        &inarray,
    Array<OneD, Array<OneD, NekDouble> >              &outarray)
{
    int physTot = m_fields[0]->GetTotPoints();
    Array<OneD, NekDouble> wk(physTot, 0.0);

    // U1 = [gy*u1-fy*u2]/(fx*gy-gx*fy)
    Vmath::Vmul(physTot, inarray[1], 1, m_GeometricInfo[1], 1,
                                        outarray[0], 1);
    Vmath::Vvtvm(physTot, inarray[0], 1, m_GeometricInfo[3], 1,
                                        outarray[0], 1,
                                        outarray[0], 1);
    Vmath::Vdiv(physTot, outarray[0], 1, m_GeometricInfo[4], 1,
                                        outarray[0], 1);

    // U2 = [fx*u2-gx*u1]/(fx*gy-gx*fy)
    Vmath::Vmul(physTot, inarray[0], 1, m_GeometricInfo[2], 1,
                                        outarray[1], 1);
    Vmath::Vvtvm(physTot, inarray[1], 1, m_GeometricInfo[0], 1,
                                        outarray[1], 1,
                                        outarray[1], 1);
    Vmath::Vdiv(physTot, outarray[1], 1, m_GeometricInfo[4], 1,
                                        outarray[1], 1);

    // U3 = u3
    if (m_nConvectiveFields ==3)
    {
        Vmath::Vcopy(physTot, inarray[2], 1, outarray[2], 1);
    }
}

void MappingXYofXY::v_CovarFromCartesian(
    const Array<OneD, Array<OneD, NekDouble> >        &inarray,
    Array<OneD, Array<OneD, NekDouble> >              &outarray)
{
    int physTot = m_fields[0]->GetTotPoints();

    // U1 = u1*fx +gx*u2
    Vmath::Vmul(physTot, m_GeometricInfo[0], 1, inarray[0], 1,
                                                outarray[0], 1);
    Vmath::Vvtvp(physTot, m_GeometricInfo[2], 1, inarray[1], 1,
                            outarray[0], 1, outarray[0],1);

    // U2 = fy*u1 + gy*u2
    Vmath::Vmul(physTot, m_GeometricInfo[1], 1, inarray[0], 1,
                                                outarray[1], 1);
    Vmath::Vvtvp(physTot, m_GeometricInfo[3], 1, inarray[1], 1,
                            outarray[1], 1, outarray[1],1);

    // U3 = u3
    if (m_nConvectiveFields ==3)
    {
        Vmath::Vcopy(physTot, inarray[2], 1, outarray[2], 1);
    }
}

void MappingXYofXY::v_GetJacobian(
    Array<OneD, NekDouble>               &outarray)
{
    int physTot = m_fields[0]->GetTotPoints();
    Vmath::Vabs(physTot, m_GeometricInfo[4], 1, outarray, 1);
}

void MappingXYofXY::v_GetMetricTensor(
    Array<OneD, Array<OneD, NekDouble> >              &outarray)
{
    int physTot = m_fields[0]->GetTotPoints();
    int nvel = m_nConvectiveFields;

    for (int i=0; i<nvel*nvel; i++)
    {
        outarray[i] = Array<OneD, NekDouble> (physTot, 0.0);
    }

    // g_{1,1} = m_metricTensor[0]
    Vmath::Vcopy(physTot, m_metricTensor[0], 1, outarray[0*nvel+0], 1);

    // g_{2,2} = m_metricTensor[1]
    Vmath::Vcopy(physTot, m_metricTensor[1], 1, outarray[1*nvel+1], 1);

    // g_{1,2}=g{2,1} = m_metricTensor[2]
    Vmath::Vcopy(physTot, m_metricTensor[2], 1, outarray[0*nvel+1], 1);
    Vmath::Vcopy(physTot, m_metricTensor[2], 1, outarray[1*nvel+0], 1);

    // g_{3,3}  = 1
    if (m_nConvectiveFields ==3)
    {
        Vmath::Sadd(physTot, 1.0, outarray[2*nvel+2], 1,
                                    outarray[2*nvel+2], 1);
    }
}

void MappingXYofXY::v_GetInvMetricTensor(
    Array<OneD, Array<OneD, NekDouble> >              &outarray)
{
    int physTot = m_fields[0]->GetTotPoints();
    int nvel = m_nConvectiveFields;

    for (int i=0; i<nvel*nvel; i++)
    {
        outarray[i] = Array<OneD, NekDouble> (physTot, 0.0);
    }

    // Get Jacobian
    Array<OneD, NekDouble> Jac(physTot, 0.0);
    GetJacobian(Jac);

    // Get Jacobian squared
    Array<OneD, NekDouble> wk(physTot, 0.0);
    Vmath::Vmul(physTot, Jac, 1, Jac, 1, wk, 1);
    // G^{1,1} = m_metricTensor[1]/Jac^2
    Vmath::Vcopy(physTot, m_metricTensor[1], 1, outarray[0*nvel+0], 1);
    Vmath::Vdiv(physTot, outarray[0*nvel+0], 1, wk,1,
                                            outarray[0*nvel+0], 1);

    // G^{2,2} = m_metricTensor[0]/Jac^2
    Vmath::Vcopy(physTot, m_metricTensor[0], 1, outarray[1*nvel+1], 1);
    Vmath::Vdiv(physTot, outarray[1*nvel+1], 1, wk,1,
                                            outarray[1*nvel+1], 1);

    // G^{1,2} = G^{2,1} = -m_metricTensor[2]/Jac^2
    Vmath::Vcopy(physTot, m_metricTensor[2], 1, outarray[0*nvel+1], 1);
    Vmath::Neg(physTot, outarray[0*nvel+1], 1);
    Vmath::Vdiv(physTot, outarray[0*nvel+1], 1, wk,1,
                                            outarray[0*nvel+1], 1);
    Vmath::Vcopy(physTot, outarray[0*nvel+1], 1, outarray[1*nvel+0], 1);

    // G^{3,3}  = 1
    if (m_nConvectiveFields ==3)
    {
        Vmath::Sadd(physTot, 1.0, outarray[2*nvel+2], 1,
                                        outarray[2*nvel+2], 1);
    }
}

void MappingXYofXY::v_ApplyChristoffelContravar(
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

    // Calculate non-zero terms

    // outarray(0,0) = U1 * m_Christoffel[0] + U2 * m_Christoffel[1]
    Vmath::Vmul(physTot,m_Christoffel[0],1,inarray[0],1,
                                            outarray[0*nvel+0],1);
    Vmath::Vvtvp(physTot,m_Christoffel[1],1,inarray[1],1,
                        outarray[0*nvel+0],1,outarray[0*nvel+0],1);

    // outarray(0,1) = U1 * m_Christoffel[1] + U2 * m_Christoffel[2]
    Vmath::Vmul(physTot,m_Christoffel[1],1,inarray[0],1,
                                            outarray[0*nvel+1],1);
    Vmath::Vvtvp(physTot,m_Christoffel[2],1,inarray[1],1,
                        outarray[0*nvel+1],1,outarray[0*nvel+1],1);

    // outarray(1,0) = U1 * m_Christoffel[3] + U2 * m_Christoffel[4]
    Vmath::Vmul(physTot,m_Christoffel[3],1,inarray[0],1,
                                            outarray[1*nvel+0],1);
    Vmath::Vvtvp(physTot,m_Christoffel[4],1,inarray[1],1,
                        outarray[1*nvel+0],1,outarray[1*nvel+0],1);

    // outarray(1,1) = U1 * m_Christoffel[4] + U2 * m_Christoffel[5]
    Vmath::Vmul(physTot,m_Christoffel[4],1,inarray[0],1,
                                            outarray[1*nvel+1],1);
    Vmath::Vvtvp(physTot,m_Christoffel[5],1,inarray[1],1,
                        outarray[1*nvel+1],1,outarray[1*nvel+1],1);

}

void MappingXYofXY::v_ApplyChristoffelCovar(
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

    // Calculate non-zero terms

    // outarray(0,0) = U1 * m_Christoffel[0] + U2 * m_Christoffel[3]
    Vmath::Vmul(physTot,m_Christoffel[0],1,inarray[0],1,
                                            outarray[0*nvel+0],1);
    Vmath::Vvtvp(physTot,m_Christoffel[3],1,inarray[1],1,
                        outarray[0*nvel+0],1,outarray[0*nvel+0],1);

    // outarray(0,1) = U1 * m_Christoffel[1] + U2 * m_Christoffel[4]
    Vmath::Vmul(physTot,m_Christoffel[1],1,inarray[0],1,
                                            outarray[0*nvel+1],1);
    Vmath::Vvtvp(physTot,m_Christoffel[4],1,inarray[1],1,
                        outarray[0*nvel+1],1,outarray[0*nvel+1],1);

    // outarray(1,0) = U1 * m_Christoffel[1] + U2 * m_Christoffel[4]
    Vmath::Vmul(physTot,m_Christoffel[1],1,inarray[0],1,
                                            outarray[1*nvel+0],1);
    Vmath::Vvtvp(physTot,m_Christoffel[4],1,inarray[1],1,
                        outarray[1*nvel+0],1,outarray[1*nvel+0],1);

    // outarray(1,1) = U1 * m_Christoffel[2] + U2 * m_Christoffel[5]
    Vmath::Vmul(physTot,m_Christoffel[2],1,inarray[0],1,
                                            outarray[1*nvel+1],1);
    Vmath::Vvtvp(physTot,m_Christoffel[5],1,inarray[1],1,
                        outarray[1*nvel+1],1,outarray[1*nvel+1],1);
}

void MappingXYofXY::v_UpdateGeomInfo()
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

    // Calculate derivatives of x transformation --> m_GeometricInfo 0-1
    m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0],m_coords[0],m_GeometricInfo[0]);
    m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1],m_coords[0],m_GeometricInfo[1]);

    // Calculate derivatives of y transformation m_GeometricInfo 2-3
    m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0],m_coords[1],m_GeometricInfo[2]);
    m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1],m_coords[1],m_GeometricInfo[3]);

    // Calculate fx*gy-gx*fy --> m_GeometricInfo4
    Vmath::Vmul(phystot, m_GeometricInfo[1], 1, m_GeometricInfo[2], 1, m_GeometricInfo[4], 1);
    Vmath::Vvtvm(phystot, m_GeometricInfo[0], 1, m_GeometricInfo[3], 1,
                                                m_GeometricInfo[4], 1,
                                                m_GeometricInfo[4], 1);
    //
    CalculateMetricTensor();
    CalculateChristoffel();

    m_fields[0]->SetWaveSpace(waveSpace);
}

void MappingXYofXY::CalculateMetricTensor()
{
    int physTot = m_fields[0]->GetTotPoints();
    // Allocate memory
    m_metricTensor =  Array<OneD, Array<OneD, NekDouble> >(3);
    for (int i = 0; i < m_metricTensor.size(); i++)
    {
        m_metricTensor[i] = Array<OneD, NekDouble>(physTot, 0.0);
    }
    // g_{1,1} = fx^2+gx^2
    Vmath::Vmul(physTot, m_GeometricInfo[0], 1, m_GeometricInfo[0], 1,
                                                m_metricTensor[0], 1);
    Vmath::Vvtvp(physTot, m_GeometricInfo[2], 1, m_GeometricInfo[2], 1,
                                                m_metricTensor[0], 1,
                                                m_metricTensor[0], 1);
    //g_{2,2} = fy^2+gy^2
    Vmath::Vmul(physTot, m_GeometricInfo[1], 1, m_GeometricInfo[1], 1,
                                                m_metricTensor[1], 1);
    Vmath::Vvtvp(physTot, m_GeometricInfo[3], 1, m_GeometricInfo[3], 1,
                                                m_metricTensor[1], 1,
                                                m_metricTensor[1], 1);
    //g_{1,2} = g_{2,1} = fy*fx+gx*gy
    Vmath::Vmul(physTot, m_GeometricInfo[0], 1, m_GeometricInfo[1], 1,
                                                m_metricTensor[2], 1);
    Vmath::Vvtvp(physTot, m_GeometricInfo[2], 1, m_GeometricInfo[3], 1,
                                                m_metricTensor[2], 1,
                                                m_metricTensor[2], 1);
}

void MappingXYofXY::CalculateChristoffel()
{
    int physTot = m_fields[0]->GetTotPoints();
    int nvel = m_nConvectiveFields;

    Array<OneD, Array<OneD, NekDouble> >   G(nvel*nvel);
    Array<OneD, Array<OneD, NekDouble> >   G_inv(nvel*nvel);
    Array<OneD, Array<OneD, NekDouble> >   gradG(2*2*2);
    Array<OneD, Array<OneD, NekDouble> >   tmp(2*2*2);
    m_Christoffel = Array<OneD, Array<OneD, NekDouble> > (6);
    // Allocate memory
    for (int i = 0; i < gradG.size(); i++)
    {
        gradG[i] = Array<OneD, NekDouble>(physTot, 0.0);
        tmp[i] = Array<OneD, NekDouble>(physTot, 0.0);
    }
    for (int i = 0; i < G.size(); i++)
    {
        G[i] = Array<OneD, NekDouble>(physTot, 0.0);
        G_inv[i] = Array<OneD, NekDouble>(physTot, 0.0);
    }

    // Get the metric tensor and its inverse
    GetMetricTensor(G);
    GetInvMetricTensor(G_inv);

    bool waveSpace = m_fields[0]->GetWaveSpace();
    m_fields[0]->SetWaveSpace(false);
    //Calculate gradients of g
    //   consider only 2 dimensions, since the 3rd is trivial
    for (int i = 0; i <2; i++)
    {
        for(int j=0; j<2; j++)
        {
            for(int k=0; k<2; k++)
            {
                m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[k],
                                        G[i*nvel+j],gradG[i*2*2 + j*2 + k]);
            }
        }
    }

    // Calculate tmp[p,j,k] = 1/2( gradG[pj,k]+ gradG[pk,j]-gradG[jk,p])
    for (int p = 0; p <2; p++)
    {
        for (int j = 0; j < 2; j++)
        {
            for (int k = 0; k < 2; k++)
            {
                Vmath::Vadd(physTot, gradG[p*2*2 + j*2 + k], 1,
                                     gradG[p*2*2 + k*2 + j], 1,
                                     tmp[p*2*2 + j*2 + k], 1);
                Vmath::Vsub(physTot, tmp[p*2*2 + j*2 + k], 1,
                                     gradG[j*2*2 + k*2 + p], 1,
                                     tmp[p*2*2 + j*2 + k], 1);
                Vmath::Smul(physTot, 0.5, tmp[p*2*2 + j*2 + k], 1,
                                           tmp[p*2*2 + j*2 + k], 1);
            }
        }
    }

    // Calculate Christoffel symbols = g^ip tmp[p,j,k]
    int n=0;
    for (int i = 0; i <2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            for (int k = 0; k <= j; k++)
            {
                m_Christoffel[n] = Array<OneD, NekDouble>(physTot, 0.0);
                for (int p = 0; p < 2; p++)
                {
                    Vmath::Vvtvp(physTot, G_inv[i*nvel+p], 1,
                                          tmp[p*2*2 + j*2 + k], 1,
                                          m_Christoffel[n], 1,
                                          m_Christoffel[n], 1);
                }
                n = n+1;
            }
        }
    }

    m_fields[0]->SetWaveSpace(waveSpace);
}

}
}
