////////////////////////////////////////////////////////////////////////////////
//
//  File: octree.h
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: cad object methods.
//
////////////////////////////////////////////////////////////////////////////////

#include <MeshUtils/SurfaceMeshing/SurfaceMeshing.h>

using namespace std;
namespace Nektar{
namespace MeshUtils{

Array<OneD, NekDouble> SurfaceMeshing::EdgeGrad(Array<OneD, NekDouble> uv1,
                                                Array<OneD, NekDouble> uv2,
                                                Array<OneD, NekDouble> uvx,
                                                int surf, bool &valid)
{
    NekDouble sig = m_order/2.0;

    Array<OneD, NekDouble> df(2);

    Array<OneD, NekDouble> ra = m_cad->GetSurf(surf)->P(uv1);
    Array<OneD, NekDouble> rb = m_cad->GetSurf(surf)->P(uv2);
    Array<OneD, NekDouble> rm = m_cad->GetSurf(surf)->D1(uvx);

    NekDouble dfdu,dfdv;

    dfdu     = ((rm[0] - rb[0])*rm[3] +
                (rm[1] - rb[1])*rm[4] +
                (rm[2] - rb[2])*rm[5]
                +
                (rm[0] - ra[0])*rm[3] +
                (rm[1] - ra[1])*rm[4] +
                (rm[2] - ra[2])*rm[5]) * 2.0*sig;

    dfdv     = ((rm[0] - rb[0])*rm[6] +
                (rm[1] - rb[1])*rm[7] +
                (rm[2] - rb[2])*rm[8]
                +
                (rm[0] - ra[0])*rm[6] +
                (rm[1] - ra[1])*rm[7] +
                (rm[2] - ra[2])*rm[8]) * 2.0*sig;

    df[0] = dfdu; df[1] = dfdv;

    NekDouble dfmag = sqrt(df[0]*df[0] + df[1]*df[1]);
    df[0] = df[0]/dfmag; df[1] = df[1]/dfmag;

    if(dfmag < 1E-30)
    {
        valid = false;
    }
    else
    {
        valid = true;
    }
    return df;
}

NekDouble SurfaceMeshing::EdgeF(NekDouble ux, NekDouble vx,
                                vector<Array<OneD,NekDouble> > bcs, int surf)
{
    ASSERTL0(bcs.size()==2,"need two bc for edge optmisation");
    NekDouble sig = m_order/2.0;

    NekDouble F;

    Array<OneD, NekDouble> uvx(2); uvx[0] = ux; uvx[1] = vx;

    Array<OneD, NekDouble> ra = m_cad->GetSurf(surf)->P(bcs[0]);
    Array<OneD, NekDouble> rb = m_cad->GetSurf(surf)->P(bcs[1]);
    Array<OneD, NekDouble> rm = m_cad->GetSurf(surf)->P(uvx);

    F        = ((rb[0] - rm[0])*(rb[0] - rm[0]) +
                (rb[1] - rm[1])*(rb[1] - rm[1]) +
                (rb[2] - rm[2])*(rb[2] - rm[2])
                +
                (rm[0] - ra[0])*(rm[0] - ra[0]) +
                (rm[1] - ra[1])*(rm[1] - ra[1]) +
                (rm[2] - ra[2])*(rm[2] - ra[2])) * sig;

    return F;

}

Array<OneD, NekDouble> SurfaceMeshing::FaceGrad(Array<OneD, NekDouble> uvx,
                                                vector<Array<OneD, NekDouble> > bcs,
                                                int surf, bool &valid)
{
    Array<OneD, NekDouble> df(2);

    Array<OneD, NekDouble> rx = m_cad->GetSurf(surf)->D1(uvx);
    Array<OneD, NekDouble> ra = m_cad->GetSurf(surf)->P(bcs[0]);
    Array<OneD, NekDouble> rb = m_cad->GetSurf(surf)->P(bcs[1]);
    Array<OneD, NekDouble> rc = m_cad->GetSurf(surf)->P(bcs[2]);
    Array<OneD, NekDouble> rd = m_cad->GetSurf(surf)->P(bcs[3]);
    Array<OneD, NekDouble> re = m_cad->GetSurf(surf)->P(bcs[4]);
    Array<OneD, NekDouble> rf = m_cad->GetSurf(surf)->P(bcs[5]);


    NekDouble dfdu,dfdv;

    dfdu = ((rx[0] - 0.5*(ra[0] + rb[0]))*(rx[3] - 0.5*(ra[3] + rb[3])) +
            (rx[1] - 0.5*(ra[1] + rb[1]))*(rx[4] - 0.5*(ra[4] + rb[4])) +
            (rx[2] - 0.5*(ra[2] + rb[2]))*(rx[5] - 0.5*(ra[5] + rb[5])))* 2.0*W[0]
           +
           ((rx[0] - 0.5*(rc[0] + rd[0]))*(rx[3] - 0.5*(rc[3] + rd[3])) +
            (rx[1] - 0.5*(rc[1] + rd[1]))*(rx[4] - 0.5*(rc[4] + rd[4])) +
            (rx[2] - 0.5*(rc[2] + rd[2]))*(rx[5] - 0.5*(rc[5] + rd[5])))* 2.0*W[1]
           +
           ((rx[0] - 0.5*(re[0] + rf[0]))*(rx[3] - 0.5*(re[3] + rf[3])) +
            (rx[1] - 0.5*(re[1] + rf[1]))*(rx[4] - 0.5*(re[4] + rf[4])) +
            (rx[2] - 0.5*(re[2] + rf[2]))*(rx[5] - 0.5*(re[5] + rf[5])))* 2.0*W[2];

    dfdv = ((rx[0] - 0.5*(ra[0] + rb[0]))*(rx[3] - 0.5*(ra[6] + rb[6])) +
            (rx[1] - 0.5*(ra[1] + rb[1]))*(rx[4] - 0.5*(ra[7] + rb[7])) +
            (rx[2] - 0.5*(ra[2] + rb[2]))*(rx[5] - 0.5*(ra[8] + rb[8])))* 2.0*W[0]
           +
           ((rx[0] - 0.5*(rc[0] + rd[0]))*(rx[3] - 0.5*(rc[6] + rd[6])) +
            (rx[1] - 0.5*(rc[1] + rd[1]))*(rx[4] - 0.5*(rc[7] + rd[7])) +
            (rx[2] - 0.5*(rc[2] + rd[2]))*(rx[5] - 0.5*(rc[8] + rd[8])))* 2.0*W[1]
           +
           ((rx[0] - 0.5*(re[0] + rf[0]))*(rx[3] - 0.5*(re[6] + rf[6])) +
            (rx[1] - 0.5*(re[1] + rf[1]))*(rx[4] - 0.5*(re[7] + rf[7])) +
            (rx[2] - 0.5*(re[2] + rf[2]))*(rx[5] - 0.5*(re[8] + rf[8])))* 2.0*W[2];

    df[0] = dfdu; df[1] = dfdv;
    NekDouble dfmag = sqrt(df[0]*df[0] + df[1]*df[1]);
    df[0] = df[0]/dfmag; df[1] = df[1]/dfmag;
    if(dfmag < 1E-30)
    {
        valid = false;
    }
    else
    {
        valid = true;
    }
    return df;
}

NekDouble SurfaceMeshing::FaceF(NekDouble ux, NekDouble vx,
                                vector<Array<OneD,NekDouble> > bcs, int surf)
{
    ASSERTL0(bcs.size()==6,"face optimsation needs 6 boundary springs");
    NekDouble F;

    Array<OneD, NekDouble> uvx(2); uvx[0] = ux; uvx[1] = vx;

    Array<OneD, NekDouble> rx = m_cad->GetSurf(surf)->P(uvx);
    Array<OneD, NekDouble> ra = m_cad->GetSurf(surf)->P(bcs[0]);
    Array<OneD, NekDouble> rb = m_cad->GetSurf(surf)->P(bcs[1]);
    Array<OneD, NekDouble> rc = m_cad->GetSurf(surf)->P(bcs[2]);
    Array<OneD, NekDouble> rd = m_cad->GetSurf(surf)->P(bcs[3]);
    Array<OneD, NekDouble> re = m_cad->GetSurf(surf)->P(bcs[4]);
    Array<OneD, NekDouble> rf = m_cad->GetSurf(surf)->P(bcs[5]);

    F = ((rx[0] - 0.5*(ra[0] + rb[0]))*(rx[0] - 0.5*(ra[0] + rb[0])) +
         (rx[1] - 0.5*(ra[1] + rb[1]))*(rx[1] - 0.5*(ra[1] + rb[1])) +
         (rx[2] - 0.5*(ra[2] + rb[2]))*(rx[2] - 0.5*(ra[2] + rb[2])))*W[0]
        +
        ((rx[0] - 0.5*(rc[0] + rd[0]))*(rx[0] - 0.5*(rc[0] + rd[0])) +
         (rx[1] - 0.5*(rc[1] + rd[1]))*(rx[1] - 0.5*(rc[1] + rd[1])) +
         (rx[2] - 0.5*(rc[2] + rd[2]))*(rx[2] - 0.5*(rc[2] + rd[2])))*W[1]
        +
        ((rx[0] - 0.5*(re[0] + rf[0]))*(rx[0] - 0.5*(re[0] + rf[0])) +
         (rx[1] - 0.5*(re[1] + rf[1]))*(rx[1] - 0.5*(re[1] + rf[1])) +
         (rx[2] - 0.5*(re[2] + rf[2]))*(rx[2] - 0.5*(re[2] + rf[2])))*W[2];

    return F;

}

}
}
