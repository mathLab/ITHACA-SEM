////////////////////////////////////////////////////////////////////////////////
//
//  File: CADSurfCFI.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
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
//  Description: cad object surface methods.
//
////////////////////////////////////////////////////////////////////////////////

#include "CADSurfCFI.h"

using namespace std;

namespace Nektar
{
namespace NekMeshUtils
{

std::string CADSurfCFI::key = GetCADSurfFactory().RegisterCreatorFunction(
    "cfi", CADSurfCFI::create, "CADSurfCFI");

void CADSurfCFI::Initialise(int i, cfi::Face *in, NekDouble s)
{
    m_cfiSurface = in;
    m_id         = i;
    m_scal       = s;
}

Array<OneD, NekDouble> CADSurfCFI::GetBounds()
{
    Array<OneD, NekDouble> b(4);

    cfi::UVBox bx = m_cfiSurface->calcUVBox();
    b[0]          = bx.uLower;
    b[1]          = bx.uUpper;
    b[2]          = bx.vLower;
    b[3]          = bx.vUpper;

    return b;
}

void CADSurfCFI::GetBounds(NekDouble &umin, NekDouble &umax,
                           NekDouble &vmin, NekDouble &vmax)
{
    cfi::UVBox bx = m_cfiSurface->calcUVBox();
    umin          = bx.uLower;
    umax          = bx.uUpper;
    vmin          = bx.vLower;
    vmax          = bx.vUpper;
}

Array<OneD, NekDouble> CADSurfCFI::locuv(Array<OneD, NekDouble> p,
                                         NekDouble &dist)
{
    Array<OneD, NekDouble> uv(2);
    cfi::Position px;
    px.x = p[0] / m_scal;
    px.y = p[1] / m_scal;
    px.z = p[2] / m_scal;

    boost::optional<cfi::UVPosition> r = m_cfiSurface->calcUVFromXYZ(px);

    uv[0] = r.value().u;
    uv[1] = r.value().v;

    Array<OneD, NekDouble> p2 = P(uv);

    dist =
        sqrt((p[0] - p2[0]) * (p[0] - p2[0]) + (p[1] - p2[1]) * (p[1] - p2[1]) +
             (p[2] - p2[2]) * (p[2] - p2[2])) *
        m_scal;

    return uv;
}

NekDouble CADSurfCFI::Curvature(Array<OneD, NekDouble> uv)
{
    cfi::UVPosition uvp(uv[0], uv[1]);
    cfi::MaxMinCurvaturePair mxp = m_cfiSurface->calcCurvAtUV(uvp);

    return mxp.maxCurv.curvature * m_scal;
}

Array<OneD, NekDouble> CADSurfCFI::P(Array<OneD, NekDouble> uv)
{
    cfi::UVPosition uvp(uv[0], uv[1]);
    cfi::Position p = m_cfiSurface->calcXYZAtUV(uvp);
    Array<OneD, NekDouble> out(3);
    out[0] = p.x * m_scal;
    out[1] = p.y * m_scal;
    out[2] = p.z * m_scal;

    return out;
}

void CADSurfCFI::P(Array<OneD, NekDouble> uv, NekDouble &x, NekDouble &y, NekDouble &z)
{
    cfi::UVPosition uvp(uv[0], uv[1]);
    cfi::Position p = m_cfiSurface->calcXYZAtUV(uvp);
    x               = p.x * m_scal;
    y               = p.y * m_scal;
    z               = p.z * m_scal;
}

Array<OneD, NekDouble> CADSurfCFI::N(Array<OneD, NekDouble> uv)
{
    cfi::UVPosition uvp(uv[0], uv[1]);
    cfi::Direction d = m_cfiSurface->calcFaceNormalAtUV(uvp);

    Array<OneD, NekDouble> normal(3);
    normal[0] = d.x;
    normal[1] = d.y;
    normal[2] = d.z;

    return normal;
}

Array<OneD, NekDouble> CADSurfCFI::D1(Array<OneD, NekDouble> uv)
{
    Array<OneD, NekDouble> p = P(uv);
    cfi::UVPosition uvp(uv[0], uv[1]);
    vector<cfi::DerivativeList> *l = m_cfiSurface->calcDerivAtUV(uvp);
    Array<OneD, NekDouble> r(9);
    r[0] = p[0] * m_scal;                    // x
    r[1] = p[1] * m_scal;                    // y
    r[2] = p[2] * m_scal;                    // z
    r[3] = l->at(0).derivatives[0] * m_scal; // dx/dx
    r[4] = l->at(0).derivatives[1] * m_scal; // dy/dy
    r[5] = l->at(0).derivatives[2] * m_scal; // dz/dz
    r[6] = l->at(0).derivatives[3] * m_scal; // dx/dx
    r[7] = l->at(0).derivatives[4] * m_scal; // dy/dy
    r[8] = l->at(0).derivatives[5] * m_scal; // dz/dz

    return r;
}

Array<OneD, NekDouble> CADSurfCFI::D2(Array<OneD, NekDouble> uv)
{
    Array<OneD, NekDouble> p = P(uv);
    cfi::UVPosition uvp(uv[0], uv[1]);
    vector<cfi::DerivativeList> *l = m_cfiSurface->calcDerivAtUV(uvp);
    Array<OneD, NekDouble> r(18);
    r[0]  = p[0] * m_scal;                    // x
    r[1]  = p[1] * m_scal;                    // y
    r[2]  = p[2] * m_scal;                    // z
    r[3]  = l->at(0).derivatives[0] * m_scal; // dx/dx
    r[4]  = l->at(0).derivatives[1] * m_scal; // dy/dy
    r[5]  = l->at(0).derivatives[2] * m_scal; // dz/dz
    r[6]  = l->at(0).derivatives[3] * m_scal; // dx/dx
    r[7]  = l->at(0).derivatives[4] * m_scal; // dy/dy
    r[8]  = l->at(0).derivatives[5] * m_scal; // dz/dz
    r[9]  = l->at(1).derivatives[0] * m_scal; // d2x/du2
    r[10] = l->at(1).derivatives[1] * m_scal; // d2y/du2
    r[11] = l->at(1).derivatives[2] * m_scal; // d2z/du2
    r[12] = l->at(1).derivatives[6] * m_scal; // d2x/dv2
    r[13] = l->at(1).derivatives[7] * m_scal; // d2y/dv2
    r[14] = l->at(1).derivatives[8] * m_scal; // d2z/dv2
    r[15] = l->at(1).derivatives[3] * m_scal; // d2x/dudv
    r[16] = l->at(1).derivatives[4] * m_scal; // d2y/dudv
    r[17] = l->at(1).derivatives[5] * m_scal; // d2z/dudv

    return r;
}
}
}
