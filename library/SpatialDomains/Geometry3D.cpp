////////////////////////////////////////////////////////////////////////////////
//
//  File:  Geometry3D.cpp
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
//  Description: 3D geometry information.
//
////////////////////////////////////////////////////////////////////////////////

#include <iomanip>

#include <SpatialDomains/Geometry3D.h>
#include <SpatialDomains/Geometry2D.h>
#include <SpatialDomains/GeomFactors.h>
#include <SpatialDomains/PointGeom.h>
#include <SpatialDomains/SegGeom.h>
#include <iomanip>

using namespace std;

namespace Nektar
{
namespace SpatialDomains
{

Geometry3D::Geometry3D()
{
}

Geometry3D::Geometry3D(const int coordim) : Geometry(coordim)
{
    ASSERTL0(m_coordim > 2,
             "Coordinate dimension should be at least 3 for a 3D geometry.");
}

Geometry3D::~Geometry3D()
{
}

//---------------------------------------
// Helper functions
//---------------------------------------

/**
 * @brief Returns the element coordinate direction corresponding to a given face
 * coordinate direction
 */
int Geometry3D::GetDir(const int faceidx, const int facedir) const
{
    return v_GetDir(faceidx, facedir);
}

//---------------------------------------
// 3D Geometry Methods
//---------------------------------------
void Geometry3D::NewtonIterationForLocCoord(
    const Array<OneD, const NekDouble> &coords,
    const Array<OneD, const NekDouble> &ptsx,
    const Array<OneD, const NekDouble> &ptsy,
    const Array<OneD, const NekDouble> &ptsz,
    Array<OneD, NekDouble> &Lcoords,
    NekDouble &resid)
{
    // maximum iterations for convergence
    const int MaxIterations = 51;
    // |x-xp|^2 < EPSILON  error tolerance
    const NekDouble Tol = 1.e-8;
    // |r,s|    > LcoordDIV stop the search
    const NekDouble LcoordDiv = 15.0;

    Array<OneD, const NekDouble> Jac =
        m_geomFactors->GetJac(m_xmap->GetPointsKeys());

    NekDouble ScaledTol = Vmath::Vsum(Jac.size(), Jac, 1) /
                          ((NekDouble)Jac.size());
    ScaledTol *= Tol;

    NekDouble xmap, ymap, zmap, F1, F2, F3;

    NekDouble derx_1, derx_2, derx_3, dery_1, dery_2, dery_3, derz_1, derz_2,
        derz_3, jac;

    // save intiial guess for later reference if required.
    NekDouble init0 = Lcoords[0], init1 = Lcoords[1], init2 = Lcoords[2];

    Array<OneD, NekDouble> DxD1(ptsx.size());
    Array<OneD, NekDouble> DxD2(ptsx.size());
    Array<OneD, NekDouble> DxD3(ptsx.size());
    Array<OneD, NekDouble> DyD1(ptsx.size());
    Array<OneD, NekDouble> DyD2(ptsx.size());
    Array<OneD, NekDouble> DyD3(ptsx.size());
    Array<OneD, NekDouble> DzD1(ptsx.size());
    Array<OneD, NekDouble> DzD2(ptsx.size());
    Array<OneD, NekDouble> DzD3(ptsx.size());

    // Ideally this will be stored in m_geomfactors
    m_xmap->PhysDeriv(ptsx, DxD1, DxD2, DxD3);
    m_xmap->PhysDeriv(ptsy, DyD1, DyD2, DyD3);
    m_xmap->PhysDeriv(ptsz, DzD1, DzD2, DzD3);

    int cnt = 0;
    Array<OneD, DNekMatSharedPtr> I(3);
    Array<OneD, NekDouble> eta(3);

    F1 = F2 = F3 = 2000; // Starting value of Function

    while (cnt++ < MaxIterations)
    {
        //  evaluate lagrange interpolant at Lcoords
        m_xmap->LocCoordToLocCollapsed(Lcoords, eta);
        I[0] = m_xmap->GetBasis(0)->GetI(eta);
        I[1] = m_xmap->GetBasis(1)->GetI(eta + 1);
        I[2] = m_xmap->GetBasis(2)->GetI(eta + 2);

        // calculate the global point `corresponding to Lcoords
        xmap = m_xmap->PhysEvaluate(I, ptsx);
        ymap = m_xmap->PhysEvaluate(I, ptsy);
        zmap = m_xmap->PhysEvaluate(I, ptsz);

        F1 = coords[0] - xmap;
        F2 = coords[1] - ymap;
        F3 = coords[2] - zmap;

        if (F1 * F1 + F2 * F2 + F3 * F3 < ScaledTol)
        {
            resid = sqrt(F1 * F1 + F2 * F2 + F3 * F3);
            break;
        }

        // Interpolate derivative metric at Lcoords
        derx_1 = m_xmap->PhysEvaluate(I, DxD1);
        derx_2 = m_xmap->PhysEvaluate(I, DxD2);
        derx_3 = m_xmap->PhysEvaluate(I, DxD3);
        dery_1 = m_xmap->PhysEvaluate(I, DyD1);
        dery_2 = m_xmap->PhysEvaluate(I, DyD2);
        dery_3 = m_xmap->PhysEvaluate(I, DyD3);
        derz_1 = m_xmap->PhysEvaluate(I, DzD1);
        derz_2 = m_xmap->PhysEvaluate(I, DzD2);
        derz_3 = m_xmap->PhysEvaluate(I, DzD3);

        jac = derx_1 * (dery_2 * derz_3 - dery_3 * derz_2) -
              derx_2 * (dery_1 * derz_3 - dery_3 * derz_1) +
              derx_3 * (dery_1 * derz_2 - dery_2 * derz_1);

        // use analytical inverse of derivitives which are also similar to
        // those of metric factors.
        Lcoords[0] =
            Lcoords[0] +
            ((dery_2 * derz_3 - dery_3 * derz_2) * (coords[0] - xmap) -
             (derx_2 * derz_3 - derx_3 * derz_2) * (coords[1] - ymap) +
             (derx_2 * dery_3 - derx_3 * dery_2) * (coords[2] - zmap)) /
                jac;

        Lcoords[1] =
            Lcoords[1] -
            ((dery_1 * derz_3 - dery_3 * derz_1) * (coords[0] - xmap) -
             (derx_1 * derz_3 - derx_3 * derz_1) * (coords[1] - ymap) +
             (derx_1 * dery_3 - derx_3 * dery_1) * (coords[2] - zmap)) /
                jac;

        Lcoords[2] =
            Lcoords[2] +
            ((dery_1 * derz_2 - dery_2 * derz_1) * (coords[0] - xmap) -
             (derx_1 * derz_2 - derx_2 * derz_1) * (coords[1] - ymap) +
             (derx_1 * dery_2 - derx_2 * dery_1) * (coords[2] - zmap)) /
                jac;

        if (fabs(Lcoords[0]) > LcoordDiv || fabs(Lcoords[1]) > LcoordDiv ||
            fabs(Lcoords[0]) > LcoordDiv)
        {
            break; // lcoords have diverged so stop iteration
        }
    }

    resid = sqrt(F1 * F1 + F2 * F2 + F3 * F3);

    if (cnt >= MaxIterations)
    {
        Array<OneD, NekDouble> collCoords(3);
        m_xmap->LocCoordToLocCollapsed(Lcoords, collCoords);

        // if coordinate is inside element dump error!
        if ((collCoords[0] >= -1.0 && collCoords[0] <= 1.0) &&
            (collCoords[1] >= -1.0 && collCoords[1] <= 1.0) &&
            (collCoords[2] >= -1.0 && collCoords[2] <= 1.0))
        {
            std::ostringstream ss;

            ss << "Reached MaxIterations (" << MaxIterations
               << ") in Newton iteration ";
            ss << "Init value (" << setprecision(4) << init0 << "," << init1
               << "," << init2 << ") ";
            ss << "Fin  value (" << Lcoords[0] << "," << Lcoords[1] << ","
               << Lcoords[2] << ") ";
            ss << "Resid = " << resid << " Tolerance = " << sqrt(ScaledTol);

            WARNINGL1(cnt < MaxIterations, ss.str());
        }
    }
}

/**
 * @brief Put all quadrature information into face/edge structure and
 * backward transform.
 *
 * Note verts, edges, and faces are listed according to anticlockwise
 * convention but points in _coeffs have to be in array format from left
 * to right.
 */
void Geometry3D::v_FillGeom()
{
    if (m_state == ePtsFilled)
    {
        return;
    }

    int i, j, k;

    for (i = 0; i < m_forient.size(); i++)
    {
        m_faces[i]->FillGeom();

        int nFaceCoeffs = m_faces[i]->GetXmap()->GetNcoeffs();

        Array<OneD, unsigned int> mapArray(nFaceCoeffs);
        Array<OneD, int> signArray(nFaceCoeffs);

        if (m_forient[i] < 9)
        {
            m_xmap->GetFaceToElementMap(
                i,
                m_forient[i],
                mapArray,
                signArray,
                m_faces[i]->GetXmap()->GetEdgeNcoeffs(0),
                m_faces[i]->GetXmap()->GetEdgeNcoeffs(1));
        }
        else
        {
            m_xmap->GetFaceToElementMap(
                i,
                m_forient[i],
                mapArray,
                signArray,
                m_faces[i]->GetXmap()->GetEdgeNcoeffs(1),
                m_faces[i]->GetXmap()->GetEdgeNcoeffs(0));
        }

        for (j = 0; j < m_coordim; j++)
        {
            const Array<OneD, const NekDouble> &coeffs =
                m_faces[i]->GetCoeffs(j);

            for (k = 0; k < nFaceCoeffs; k++)
            {
                NekDouble v = signArray[k] * coeffs[k];
                m_coeffs[j][mapArray[k]] = v;
            }
        }
    }

    m_state = ePtsFilled;
}

/**
* @brief Given local collapsed coordinate Lcoord return the value of
* physical coordinate in direction i.
*/
NekDouble Geometry3D::v_GetCoord(const int i,
                                 const Array<OneD, const NekDouble> &Lcoord)
{
    ASSERTL1(m_state == ePtsFilled, "Geometry is not in physical space");

    Array<OneD, NekDouble> tmp(m_xmap->GetTotPoints());
    m_xmap->BwdTrans(m_coeffs[i], tmp);

    return m_xmap->PhysEvaluate(Lcoord, tmp);
}

//---------------------------------------
// Helper functions
//---------------------------------------

int Geometry3D::v_GetShapeDim() const
{
    return 3;
}

int Geometry3D::v_GetNumVerts() const
{
    return m_verts.size();
}

int Geometry3D::v_GetNumEdges() const
{
    return m_edges.size();
}

int Geometry3D::v_GetNumFaces() const
{
    return m_faces.size();
}

PointGeomSharedPtr Geometry3D::v_GetVertex(int i) const
{
    return m_verts[i];
}

Geometry1DSharedPtr Geometry3D::v_GetEdge(int i) const
{
    ASSERTL2(i >= 0 && i <= m_edges.size() - 1,
             "Edge ID must be between 0 and " +
                 boost::lexical_cast<string>(m_edges.size() - 1));
    return m_edges[i];
}

Geometry2DSharedPtr Geometry3D::v_GetFace(int i) const
{
    ASSERTL2((i >= 0) && (i <= 5), "Edge id must be between 0 and 4");
    return m_faces[i];
}

inline StdRegions::Orientation Geometry3D::v_GetEorient(const int i) const
{
    ASSERTL2(i >= 0 && i <= m_edges.size() - 1,
             "Edge ID must be between 0 and " +
                 boost::lexical_cast<string>(m_edges.size() - 1));
    return m_eorient[i];
}

StdRegions::Orientation Geometry3D::v_GetForient(const int i) const
{
    ASSERTL2(i >= 0 && i <= m_faces.size() - 1,
             "Face ID must be between 0 and " +
                 boost::lexical_cast<string>(m_faces.size() - 1));
    return m_forient[i];
}

}
}
