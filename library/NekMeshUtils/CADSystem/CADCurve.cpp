////////////////////////////////////////////////////////////////////////////////
//
//  File: CADSystem.cpp
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
//  Description: cad object methods.
//
////////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include "CADCurve.h"
#include "CADSurf.h"

using namespace std;

namespace Nektar
{
namespace NekMeshUtils
{

Array<OneD, NekDouble> CADCurve::NormalWRT(NekDouble t, int surf)
{
    boost::ignore_unused(surf);

    Array<OneD, NekDouble> p = P(t);
    pair<weak_ptr<CADSurf>, CADOrientation::Orientation> surface;
    ASSERTL0(m_adjSurfs.size() == 1,
             "This will only work in 2D for one surface at the moment");
    surface = m_adjSurfs[0];

    Array<OneD, NekDouble> uv = surface.first.lock()->locuv(p);
    Array<OneD, NekDouble> d1 = surface.first.lock()->D1(uv);

    NekDouble t1 = t - 1e-8;
    NekDouble t2 = t + 1e-8;

    if (surface.second == CADOrientation::eBackwards)
    {
        swap(t1, t2);
    }

    Array<OneD, NekDouble> uv1 = surface.first.lock()->locuv(P(t1));
    Array<OneD, NekDouble> uv2 = surface.first.lock()->locuv(P(t2));

    NekDouble du = uv2[1] - uv1[1];
    NekDouble dv = -1.0 * (uv2[0] - uv1[0]);

    Array<OneD, NekDouble> N(3, 0.0);
    N[0] = (d1[3] * du + d1[6] * dv) / 2.0;
    N[1] = (d1[4] * du + d1[7] * dv) / 2.0;
    N[2] = (d1[5] * du + d1[8] * dv) / 2.0;

    NekDouble mag = sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
    N[0] /= mag;
    N[1] /= mag;
    N[2] /= mag;

    return N;
}

CADOrientation::Orientation CADCurve::GetOrienationWRT(int surf)
{
    for (int i = 0; i < m_adjSurfs.size(); i++)
    {
        if (m_adjSurfs[i].first.lock()->GetId() == surf)
        {
            return m_adjSurfs[i].second;
        }
    }

    ASSERTL0(false, "surf not in adjecency list");
    return CADOrientation::eUnknown;
}
}
}
