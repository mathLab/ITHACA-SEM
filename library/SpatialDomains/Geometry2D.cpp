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
    NekDouble &dist)
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
    NekDouble resid;
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

    m_xmap->LocCoordToLocCollapsed(Lcoords, eta);
    if(ClampLocCoords(eta, 0.))
    {
        I[0] = m_xmap->GetBasis(0)->GetI(eta);
        I[1] = m_xmap->GetBasis(1)->GetI(eta + 1);
        // calculate the global point corresponding to Lcoords
        xmap = m_xmap->PhysEvaluate(I, ptsx);
        ymap = m_xmap->PhysEvaluate(I, ptsy);
        F1 = coords[0] - xmap;
        F2 = coords[1] - ymap;
        dist = sqrt(F1 * F1 + F2 * F2);
    }
    else
    {
        dist = 0.;
    }

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

NekDouble Geometry2D::v_GetLocCoords(const Array<OneD, const NekDouble> &coords,
                                   Array<OneD, NekDouble> &Lcoords)
{
    NekDouble dist = 0.;
    if (GetMetricInfo()->GetGtype() == eRegular)
    {
        int v2;
        if(m_shapeType == LibUtilities::eTriangle)
        {
            v2 = 2;
        }
        else if(m_shapeType == LibUtilities::eQuadrilateral)
        {
            v2 = 3;
        }
        else
        {
            v2 = 2;
            ASSERTL0(false, "unrecognized 2D element type");
        }

        NekDouble coords2 = (m_coordim == 3) ? coords[2] : 0.0;
        PointGeom r(m_coordim, 0, coords[0], coords[1], coords2);

        // Edges
        PointGeom er0, e10, e20;
        PointGeom norm, orth1, orth2;
        er0.Sub(r, *m_verts[0]);
        e10.Sub(*m_verts[1], *m_verts[0]);
        e20.Sub(*m_verts[v2], *m_verts[0]);

        // Obtain normal to element plane
        norm.Mult(e10, e20);
        // Obtain vector which are proportional to normal of e10 and e20.
        orth1.Mult(norm, e10);
        orth2.Mult(norm, e20);

        // Calculate length using L/|dv1| = (x-v0).n1/(dv1.n1) for coordiante 1
        // Then rescale to [-1,1].
        Lcoords[0] = er0.dot(orth2) / e10.dot(orth2);
        Lcoords[0] = 2. * Lcoords[0] - 1.;
        Lcoords[1] = er0.dot(orth1) / e20.dot(orth1);
        Lcoords[1] = 2. * Lcoords[1] - 1.;

        // Set distance
        Array<OneD, NekDouble> eta(2, 0.);
        m_xmap->LocCoordToLocCollapsed(Lcoords, eta);
        if(ClampLocCoords(eta, 0.))
        {
            Array<OneD, NekDouble> xi(2, 0.);
            m_xmap->LocCollapsedToLocCoord(eta, xi);
            xi[0] = (xi[0] + 1.) * 0.5; //re-scaled to ratio [0, 1]
            xi[1] = (xi[1] + 1.) * 0.5;
            for (int i = 0; i < m_coordim; ++i)
            {
                NekDouble tmp = xi[0]*e10[i] + xi[1]*e20[i] - er0[i];
                dist += tmp * tmp;
            }
            dist = sqrt(dist);
        }
    }
    else
    {
        v_FillGeom();
        // Determine nearest point of coords  to values in m_xmap
        int npts = m_xmap->GetTotPoints();
        Array<OneD, NekDouble> ptsx(npts), ptsy(npts);
        Array<OneD, NekDouble> tmpx(npts), tmpy(npts);

        m_xmap->BwdTrans(m_coeffs[0], ptsx);
        m_xmap->BwdTrans(m_coeffs[1], ptsy);

        const Array<OneD, const NekDouble> za = m_xmap->GetPoints(0);
        const Array<OneD, const NekDouble> zb = m_xmap->GetPoints(1);

        // guess the first local coords based on nearest point
        Vmath::Sadd(npts, -coords[0], ptsx, 1, tmpx, 1);
        Vmath::Sadd(npts, -coords[1], ptsy, 1, tmpy, 1);
        Vmath::Vmul(npts, tmpx, 1, tmpx, 1, tmpx, 1);
        Vmath::Vvtvp(npts, tmpy, 1, tmpy, 1, tmpx, 1, tmpx, 1);

        int min_i = Vmath::Imin(npts, tmpx, 1);

        Array<OneD, NekDouble> eta(2, 0.);
        eta[0] = za[min_i % za.size()];
        eta[1] = zb[min_i / za.size()];
        m_xmap->LocCollapsedToLocCoord(eta, Lcoords);

        // Perform newton iteration to find local coordinates
        NewtonIterationForLocCoord(coords, ptsx, ptsy, Lcoords, dist);
    }
    return dist;
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
