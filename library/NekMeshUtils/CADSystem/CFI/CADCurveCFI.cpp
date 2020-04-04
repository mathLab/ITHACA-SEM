////////////////////////////////////////////////////////////////////////////////
//
//  File: CADCurveCFI.cpp
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
//  Description: CAD object curve methods.
//
////////////////////////////////////////////////////////////////////////////////

#include "CADCurveCFI.h"

using namespace std;

namespace Nektar
{
namespace NekMeshUtils
{

std::string CADCurveCFI::key = GetCADCurveFactory().RegisterCreatorFunction(
    "cfi", CADCurveCFI::create, "CADCurveCFI");

void CADCurveCFI::Initialise(int i, cfi::Line *in, NekDouble s)
{
    m_cfiEdge = in;
    m_scal    = s;
    m_length  = m_cfiEdge->calcLength() * m_scal;

    m_id = i;
}

NekDouble CADCurveCFI::tAtArcLength(NekDouble s)
{
    s /= m_scal;
    Array<OneD, NekDouble> bds = GetBounds();
    NekDouble dt = (bds[1] - bds[0]) / 1000;

    NekDouble t   = bds[0];
    NekDouble len = 0.0;

    while (len <= s)
    {
        Array<OneD, NekDouble> drdt1, drdt2;
        drdt1 = D2(t);
        t += dt;
        drdt2 = D2(t);

        NekDouble mag1 = sqrt(drdt1[3] * drdt1[3] + drdt1[4] * drdt1[4] +
                              drdt1[5] * drdt1[5]);
        NekDouble mag2 = sqrt(drdt2[3] * drdt2[3] + drdt2[4] * drdt2[4] +
                              drdt2[5] * drdt2[5]);

        len += (mag1 + mag2) / 2.0 * dt;
    }

    return t - dt;
}

NekDouble CADCurveCFI::loct(Array<OneD, NekDouble> xyz, NekDouble &t)
{
    cfi::Position p;
    p.x = xyz[0] / m_scal;
    p.y = xyz[1] / m_scal;
    p.z = xyz[2] / m_scal;

    boost::optional<cfi::Projected<double>> pj = m_cfiEdge->calcTFromXYZ(p, -1);

    t = pj.value().parameters;

    return pj.value().distance * m_scal;
}

NekDouble CADCurveCFI::Length(NekDouble ti, NekDouble tf)
{
    Array<OneD, NekDouble> bds = GetBounds();
    NekDouble dt = (bds[1] - bds[0]) / 1000;

    NekDouble t   = ti;
    NekDouble len = 0.0;

    while (t <= tf)
    {
        Array<OneD, NekDouble> drdt1, drdt2;
        drdt1 = D2(t);
        t += dt;
        drdt2 = D2(t);

        NekDouble mag1 = sqrt(drdt1[3] * drdt1[3] + drdt1[4] * drdt1[4] +
                              drdt1[5] * drdt1[5]);
        NekDouble mag2 = sqrt(drdt2[3] * drdt2[3] + drdt2[4] * drdt2[4] +
                              drdt2[5] * drdt2[5]);

        len += (mag1 + mag2) / 2.0 * dt;
    }

    return len * m_scal;
}

Array<OneD, NekDouble> CADCurveCFI::P(NekDouble t)
{
    cfi::Position p = m_cfiEdge->calcXYZAtT(t);

    Array<OneD, NekDouble> out(3);

    out[0] = p.x * m_scal;
    out[1] = p.y * m_scal;
    out[2] = p.z * m_scal;

    return out;
}

void CADCurveCFI::P(NekDouble t, NekDouble &x, NekDouble &y, NekDouble &z)
{
    cfi::Position p = m_cfiEdge->calcXYZAtT(t);

    x = p.x * m_scal;
    y = p.y * m_scal;
    z = p.z * m_scal;
}

Array<OneD, NekDouble> CADCurveCFI::D2(NekDouble t)
{
    vector<cfi::DerivativeList> *d = m_cfiEdge->calcDerivAtT(t);
    cfi::Position p                = m_cfiEdge->calcXYZAtT(t);

    Array<OneD, NekDouble> out(9);

    out[0] = p.x * m_scal;
    out[1] = p.y * m_scal;
    out[2] = p.z * m_scal;

    cfi::DerivativeList d1 = d->at(0);
    cfi::DerivativeList d2 = d->at(1);

    out[3] = d1.getDeriv(0) * m_scal;
    out[4] = d1.getDeriv(1) * m_scal;
    out[5] = d1.getDeriv(2) * m_scal;
    out[6] = d2.getDeriv(0) * m_scal;
    out[7] = d2.getDeriv(1) * m_scal;
    out[8] = d2.getDeriv(2) * m_scal;

    return out;
}

Array<OneD, NekDouble> CADCurveCFI::GetBounds()
{
    Array<OneD, NekDouble> t(2);
    t[0] = 0.0;
    t[1] = 1.0;

    return t;
}

void CADCurveCFI::GetBounds(NekDouble &tmin, NekDouble &tmax)
{
    tmin = 0.0;
    tmax = 1.0;
}

Array<OneD, NekDouble> CADCurveCFI::GetMinMax()
{
    Array<OneD, NekDouble> bds = GetBounds();

    cfi::Position x1 = m_cfiEdge->calcXYZAtT(bds[0]);
    cfi::Position x2 = m_cfiEdge->calcXYZAtT(bds[1]);

    Array<OneD, NekDouble> locs(6);

    locs[0] = x1.x * m_scal;
    locs[1] = x1.y * m_scal;
    locs[2] = x1.z * m_scal;
    locs[3] = x2.x * m_scal;
    locs[4] = x2.y * m_scal;
    locs[5] = x2.z * m_scal;

    return locs;
}
}
}
