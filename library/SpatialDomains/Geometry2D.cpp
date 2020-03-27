////////////////////////////////////////////////////////////////////////////////
//
//  File:  Geometry2D.cpp
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
//  Description:  2D geometry information
//
//
////////////////////////////////////////////////////////////////////////////////

#include <iomanip>

#include <SpatialDomains/Geometry2D.h>
#include <SpatialDomains/SegGeom.h>

#include <iomanip>

using namespace std;

namespace Nektar
{
namespace SpatialDomains
{

Geometry2D::Geometry2D()
{
}

Geometry2D::Geometry2D(const int coordim, CurveSharedPtr curve)
    : Geometry(coordim), m_curve(curve)
{
    ASSERTL0(m_coordim > 1,
             "Coordinate dimension should be at least 2 for a 2D geometry");
}

Geometry2D::~Geometry2D()
{
}

void Geometry2D::NewtonIterationForLocCoord(
    const Array<OneD, const NekDouble> &coords,
    const Array<OneD, const NekDouble> &ptsx,
    const Array<OneD, const NekDouble> &ptsy,
    Array<OneD, NekDouble> &Lcoords,
    NekDouble &resid)
{
    // Maximum iterations for convergence
    const int MaxIterations = 51;
    // |x-xp|^2 < EPSILON  error    tolerance
    const NekDouble Tol = 1.e-8;
    // |r,s|    > LcoordDIV stop   the search
    const NekDouble LcoordDiv = 15.0;

    Array<OneD, const NekDouble> Jac =
        m_geomFactors->GetJac(m_xmap->GetPointsKeys());

    NekDouble ScaledTol = Vmath::Vsum(Jac.size(), Jac, 1) /
                          ((NekDouble)Jac.size());
    ScaledTol *= Tol;

    NekDouble xmap, ymap, F1, F2;
    NekDouble derx_1, derx_2, dery_1, dery_2, jac;

    // save intiial guess for later reference if required.
    NekDouble init0 = Lcoords[0], init1 = Lcoords[1];

    Array<OneD, NekDouble> DxD1(ptsx.size());
    Array<OneD, NekDouble> DxD2(ptsx.size());
    Array<OneD, NekDouble> DyD1(ptsx.size());
    Array<OneD, NekDouble> DyD2(ptsx.size());

    // Ideally this will be stored in m_geomfactors
    m_xmap->PhysDeriv(ptsx, DxD1, DxD2);
    m_xmap->PhysDeriv(ptsy, DyD1, DyD2);

    int cnt = 0;
    Array<OneD, DNekMatSharedPtr> I(2);
    Array<OneD, NekDouble> eta(2);

    F1 = F2 = 2000; // Starting value of Function

    while (cnt++ < MaxIterations)
    {
        //  evaluate lagrange interpolant at Lcoords
        m_xmap->LocCoordToLocCollapsed(Lcoords, eta);
        I[0] = m_xmap->GetBasis(0)->GetI(eta);
        I[1] = m_xmap->GetBasis(1)->GetI(eta + 1);

        // calculate the global point `corresponding to Lcoords
        xmap = m_xmap->PhysEvaluate(I, ptsx);
        ymap = m_xmap->PhysEvaluate(I, ptsy);

        F1 = coords[0] - xmap;
        F2 = coords[1] - ymap;

        if (F1 * F1 + F2 * F2 < ScaledTol)
        {
            resid = sqrt(F1 * F1 + F2 * F2);
            break;
        }

        // Interpolate derivative metric at Lcoords
        derx_1 = m_xmap->PhysEvaluate(I, DxD1);
        derx_2 = m_xmap->PhysEvaluate(I, DxD2);
        dery_1 = m_xmap->PhysEvaluate(I, DyD1);
        dery_2 = m_xmap->PhysEvaluate(I, DyD2);

        jac = dery_2 * derx_1 - dery_1 * derx_2;

        // use analytical inverse of derivitives which are
        // also similar to those of metric factors.
        Lcoords[0] =
            Lcoords[0] +
            (dery_2 * (coords[0] - xmap) - derx_2 * (coords[1] - ymap)) / jac;

        Lcoords[1] =
            Lcoords[1] +
            (-dery_1 * (coords[0] - xmap) + derx_1 * (coords[1] - ymap)) / jac;

        if (fabs(Lcoords[0]) > LcoordDiv || fabs(Lcoords[1]) > LcoordDiv)
        {
            break; // lcoords have diverged so stop iteration
        }
    }

    resid = sqrt(F1 * F1 + F2 * F2);

    if (cnt >= MaxIterations)
    {
        Array<OneD, NekDouble> collCoords(2);
        m_xmap->LocCoordToLocCollapsed(Lcoords, collCoords);

        // if coordinate is inside element dump error!
        if ((collCoords[0] >= -1.0 && collCoords[0] <= 1.0) &&
            (collCoords[1] >= -1.0 && collCoords[1] <= 1.0))
        {
            std::ostringstream ss;

            ss << "Reached MaxIterations (" << MaxIterations
               << ") in Newton iteration ";
            ss << "Init value (" << setprecision(4) << init0 << "," << init1
               << ","
               << ") ";
            ss << "Fin  value (" << Lcoords[0] << "," << Lcoords[1] << ","
               << ") ";
            ss << "Resid = " << resid << " Tolerance = " << sqrt(ScaledTol);

            WARNINGL1(cnt < MaxIterations, ss.str());
        }
    }
}

int Geometry2D::v_GetNumVerts() const
{
    return m_verts.size();
}

int Geometry2D::v_GetNumEdges() const
{
    return m_edges.size();
}

PointGeomSharedPtr Geometry2D::v_GetVertex(int i) const
{
    ASSERTL2(i >= 0 && i < m_verts.size(), "Index out of range");
    return m_verts[i];
}

Geometry1DSharedPtr Geometry2D::v_GetEdge(int i) const
{
    ASSERTL2(i >= 0 && i < m_edges.size(), "Index out of range");
    return m_edges[i];
}

StdRegions::Orientation Geometry2D::v_GetEorient(const int i) const
{
    ASSERTL2(i >= 0 && i < m_eorient.size(), "Index out of range");
    return m_eorient[i];
}

int Geometry2D::v_GetShapeDim() const
{
    return 2;
}

}
}
