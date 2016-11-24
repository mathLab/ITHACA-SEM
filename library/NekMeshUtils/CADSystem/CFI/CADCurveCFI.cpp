////////////////////////////////////////////////////////////////////////////////
//
//  File: CADCurve.cpp
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


void CADCurveCFI::Initialise(int i, cfi::Line* in)
{
    m_cfiEdge = in;
    m_length = m_cfiEdge->calcLength();

    m_id   = i;
}

NekDouble CADCurveCFI::tAtArcLength(NekDouble s)
{
    Array<OneD, NekDouble> bds = Bounds();
    NekDouble dt = (bds[1] - bds[0]) / 5000;

    NekDouble t = bds[0];
    NekDouble len = 0.0;

    while (len <= s)
    {
        Array<OneD, NekDouble> drdt1, drdt2;
        drdt1 = D2(t);
        t += dt;
        drdt2 = D2(t);

        NekDouble mag1 = sqrt(drdt1[3]*drdt1[3] + drdt1[4]*drdt1[4] + drdt1[5]*drdt1[5]);
        NekDouble mag2 = sqrt(drdt2[3]*drdt2[3] + drdt2[4]*drdt2[4] + drdt2[5]*drdt2[5]);

        len += (mag1 + mag2) / 2.0 * dt;
    }

    return t - dt;
}

NekDouble CADCurveCFI::loct(Array<OneD, NekDouble> xyz)
{
    cfi::Position p;
    p.x = xyz[0];
    p.y = xyz[1];
    p.z = xyz[2];

    boost::optional<cfi::Projected<double> > pj = m_cfiEdge->calcTFromXYZ(p,-1);

    if(pj.value().distance > 1e-5)
    {
        cerr << "large loct distance" << endl;
    }

    return pj.value().parameters;
}

NekDouble CADCurveCFI::Length(NekDouble ti, NekDouble tf)
{
    Array<OneD, NekDouble> bds = Bounds();
    NekDouble dt = (bds[1] - bds[0]) / 5000;

    NekDouble t = ti;
    NekDouble len = 0.0;

    while (t <= tf)
    {
        Array<OneD, NekDouble> drdt1, drdt2;
        drdt1 = D2(t);
        t += dt;
        drdt2 = D2(t);

        NekDouble mag1 = sqrt(drdt1[3]*drdt1[3] + drdt1[4]*drdt1[4] + drdt1[5]*drdt1[5]);
        NekDouble mag2 = sqrt(drdt2[3]*drdt2[3] + drdt2[4]*drdt2[4] + drdt2[5]*drdt2[5]);

        len += (mag1 + mag2) / 2.0 * dt;
    }

    return len;
}

Array<OneD, NekDouble> CADCurveCFI::P(NekDouble t)
{
    cfi::Position p = m_cfiEdge->calcXYZAtT(t);

    Array<OneD, NekDouble> out(3);

    out[0] = p.x;
    out[1] = p.y;
    out[2] = p.z;

    return out;
}

Array<OneD, NekDouble> CADCurveCFI::D2(NekDouble t)
{
    vector<cfi::DerivativeList>* d = m_cfiEdge->calcDerivAtT(t);
    cfi::Position p = m_cfiEdge->calcXYZAtT(t);

    Array<OneD, NekDouble> out(9);

    out[0] = p.x;
    out[1] = p.y;
    out[2] = p.z;

    cfi::DerivativeList d1 = d->at(0);
    cfi::DerivativeList d2 = d->at(1);

    out[3] = d1.getDeriv(0);
    out[4] = d1.getDeriv(1);
    out[5] = d1.getDeriv(2);
    out[6] = d2.getDeriv(0);
    out[7] = d2.getDeriv(1);
    out[8] = d2.getDeriv(2);

    return out;
}

Array<OneD, NekDouble> CADCurveCFI::Bounds()
{
    Array<OneD, NekDouble> t(2);
    t[0] = 0.0;
    t[1] = 1.0;

    return t;
}

Array<OneD, NekDouble> CADCurveCFI::GetMinMax()
{
    Array<OneD, NekDouble> bds = Bounds();

    cfi::Position x1 = m_cfiEdge->calcXYZAtT(bds[0]);
    cfi::Position x2 = m_cfiEdge->calcXYZAtT(bds[1]);

    Array<OneD, NekDouble> locs(6);

    locs[0] = x1.x;
    locs[1] = x1.y;
    locs[2] = x1.z;
    locs[3] = x2.x;
    locs[4] = x2.y;
    locs[5] = x2.z;

    return locs;
}

}
}
