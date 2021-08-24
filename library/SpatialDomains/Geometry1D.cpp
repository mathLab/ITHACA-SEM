////////////////////////////////////////////////////////////////////////////////
//
//  File: Geometry1D.cpp
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
//  Description:  1D geometry information
//
//
////////////////////////////////////////////////////////////////////////////////
#include <SpatialDomains/Geometry1D.h>

namespace Nektar
{
namespace SpatialDomains
{

Geometry1D::Geometry1D()
{
}

Geometry1D::Geometry1D(const int coordim) : Geometry(coordim)
{
}

Geometry1D::~Geometry1D()
{
}

int Geometry1D::v_GetShapeDim() const
{
    return 1;
}

NekDouble Geometry1D::v_GetLocCoords(const Array<OneD, const NekDouble> &coords,
                                  Array<OneD, NekDouble> &Lcoords)
{
    NekDouble dist = 0.;
    v_FillGeom();

    // calculate local coordinate for coord
    if (GetMetricInfo()->GetGtype() == eRegular)
    {
        NekDouble len = 0.0;
        NekDouble xi = 0.0;

        const int npts = m_xmap->GetTotPoints();
        Array<OneD, NekDouble> pts(npts);

        for (int i = 0; i < m_coordim; ++i)
        {
            m_xmap->BwdTrans(m_coeffs[i], pts);
            len += (pts[npts - 1] - pts[0]) * (pts[npts - 1] - pts[0]);
            xi += (coords[i] - pts[0]) * (pts[npts-1] - pts[0]);
        }
        xi = xi / len;
        if(xi<0.)
        {
            dist = - xi * sqrt(len);
        }
        else if(xi>1.)
        {
            dist = (xi - 1.) * sqrt(len);
        }

        Lcoords[0] = 2. * xi - 1.0;
    }
    else
    {
        NEKERROR(ErrorUtil::efatal,
                 "inverse mapping must be set up to use this call");
    }
    return dist;
}
}
}
