////////////////////////////////////////////////////////////////////////////////
//
//  File:  SegGeom.cpp
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
//  Description:
//
//
////////////////////////////////////////////////////////////////////////////////

#include <SpatialDomains/SegGeom.h>
#include <SpatialDomains/GeomFactors.h>

#include <StdRegions/StdRegions.hpp>
#include <StdRegions/StdSegExp.h>
#include <LibUtilities/Foundations/ManagerAccess.h> // for PointsManager, etc

namespace Nektar
{
namespace SpatialDomains
{
SegGeom::SegGeom()
{
    m_shapeType = LibUtilities::eSegment;
}

SegGeom::SegGeom(int id,
                 const int coordim,
                 const PointGeomSharedPtr vertex[],
                 const CurveSharedPtr curve)
    : Geometry1D(coordim)
{
    m_shapeType = LibUtilities::eSegment;
    m_globalID = id;
    m_state = eNotFilled;
    m_curve = curve;

    m_verts[0] = vertex[0];
    m_verts[1] = vertex[1];
}

SegGeom::SegGeom(const SegGeom &in)
    : Geometry1D(in)
{
    // From Geometry class
    m_shapeType = in.m_shapeType;

    // info from EdgeComponent class
    m_globalID = in.m_globalID;
    m_xmap = in.m_xmap;
    SetUpCoeffs(m_xmap->GetNcoeffs());

    // info from SegGeom class
    m_coordim = in.m_coordim;
    m_verts[0] = in.m_verts[0];
    m_verts[1] = in.m_verts[1];

    m_state = in.m_state;
}

void SegGeom::SetUpXmap()
{
    if (m_curve)
    {
        int npts = m_curve->m_points.size();
        LibUtilities::PointsKey pkey(npts + 1,
                                     LibUtilities::eGaussLobattoLegendre);
        const LibUtilities::BasisKey B(LibUtilities::eModified_A, npts, pkey);
        m_xmap = MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(B);
    }
    else
    {
        const LibUtilities::BasisKey B(
            LibUtilities::eModified_A,
            2,
            LibUtilities::PointsKey(2, LibUtilities::eGaussLobattoLegendre));
        m_xmap = MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(B);
    }
}

/**
 * \brief Generate a one dimensional space segment geometry where the vert[0]
 * has the same x value and vert[1] is set to vert[0] plus the length of the
 * original segment
 **/
SegGeomSharedPtr SegGeom::GenerateOneSpaceDimGeom(void)
{
    SegGeomSharedPtr returnval = MemoryManager<SegGeom>::AllocateSharedPtr();

    // info about numbering
    returnval->m_globalID = m_globalID;

    // geometric information.
    returnval->m_coordim = 1;
    NekDouble x0 = (*m_verts[0])[0];
    PointGeomSharedPtr vert0 = MemoryManager<PointGeom>::AllocateSharedPtr(
        1, m_verts[0]->GetGlobalID(), x0, 0.0, 0.0);
    vert0->SetGlobalID(vert0->GetGlobalID());
    returnval->m_verts[0] = vert0;

    // Get information to calculate length.
    const Array<OneD, const LibUtilities::BasisSharedPtr> base =
        m_xmap->GetBase();
    LibUtilities::PointsKeyVector v;
    v.push_back(base[0]->GetPointsKey());
    v_GenGeomFactors();

    const Array<OneD, const NekDouble> jac = m_geomFactors->GetJac(v);

    NekDouble len = 0.0;
    if (jac.size() == 1)
    {
        len = jac[0] * 2.0;
    }
    else
    {
        Array<OneD, const NekDouble> w0 = base[0]->GetW();
        len = 0.0;

        for (int i = 0; i < jac.size(); ++i)
        {
            len += jac[i] * w0[i];
        }
    }
    // Set up second vertex.
    PointGeomSharedPtr vert1 = MemoryManager<PointGeom>::AllocateSharedPtr(
        1, m_verts[1]->GetGlobalID(), x0 + len, 0.0, 0.0);
    vert1->SetGlobalID(vert1->GetGlobalID());

    returnval->m_verts[1] = vert1;

    // at present just use previous m_xmap[0];
    returnval->m_xmap = m_xmap;
    returnval->SetUpCoeffs(m_xmap->GetNcoeffs());
    returnval->m_state = eNotFilled;

    return returnval;
}

SegGeom::~SegGeom()
{
}

LibUtilities::ShapeType SegGeom::v_GetShapeType() const
{
    return LibUtilities::eSegment;
}

NekDouble SegGeom::v_GetCoord(const int i,
                              const Array<OneD, const NekDouble> &Lcoord)
{
    ASSERTL1(m_state == ePtsFilled, "Geometry is not in physical space");

    Array<OneD, NekDouble> tmp(m_xmap->GetTotPoints());
    m_xmap->BwdTrans(m_coeffs[i], tmp);

    return m_xmap->PhysEvaluate(Lcoord, tmp);
}

/**
 * @brief Get the orientation of @p edge1.
 *
 * If @p edge1 is connected to @p edge2 in the same direction as the points
 * comprising @p edge1 then it is forward, otherwise it is backward.
 *
 * For example, assume @p edge1 is comprised of points 1 and 2, and @p edge2 is
 * comprised of points 2 and 3, then @p edge1 is forward.
 *
 * If @p edge1 is comprised of points 2 and 1 and @p edge2 is comprised of
 * points 3 and 2, then @p edge1 is backward.
 *
 * Since both edges are passed, it does not need any information from the
 * EdgeComponent instance.
 */
StdRegions::Orientation SegGeom::GetEdgeOrientation(const SegGeom &edge1,
                                                    const SegGeom &edge2)
{
    StdRegions::Orientation returnval = StdRegions::eForwards;

    if ((*edge1.GetVertex(0) == *edge2.GetVertex(0)) ||
        (*edge1.GetVertex(0) == *edge2.GetVertex(1)))
    {
        // Backward direction.  Vertex 0 is connected to edge 2.
        returnval = StdRegions::eBackwards;
    }
    else if ((*edge1.GetVertex(1) != *edge2.GetVertex(0)) &&
             (*edge1.GetVertex(1) != *edge2.GetVertex(1)))
    {
        // Not forward either, then we have a problem.
        std::ostringstream errstrm;
        errstrm << "Connected edges do not share a vertex. Edges ";
        errstrm << edge1.GetGlobalID() << ", " << edge2.GetGlobalID();
        ASSERTL0(false, errstrm.str());
    }

    return returnval;
}

void SegGeom::v_GenGeomFactors()
{
    if(!m_setupState)
    {
        SegGeom::v_Setup();
    }

    if (m_geomFactorsState != ePtsFilled)
    {
        SpatialDomains::GeomType gType = eRegular;
        SegGeom::v_FillGeom();

        if (m_xmap->GetBasisNumModes(0) != 2)
        {
            gType = eDeformed;
        }

        m_geomFactors = MemoryManager<GeomFactors>::AllocateSharedPtr(
            gType, m_coordim, m_xmap, m_coeffs);
        m_geomFactorsState = ePtsFilled;
    }
}

void SegGeom::v_FillGeom()
{
    if (m_state != ePtsFilled)
    {
        int i;

        if (m_coordim > 0 && m_curve)
        {
            int npts = m_curve->m_points.size();
            LibUtilities::PointsKey pkey(npts + 1,
                                         LibUtilities::eGaussLobattoLegendre);
            Array<OneD, NekDouble> tmp(npts);

            if (m_verts[0]->dist(*(m_curve->m_points[0])) >
                NekConstants::kVertexTheSameDouble)
            {
                std::string err =
                    "Vertex 0 is separated from first point by more than ";
                std::stringstream strstrm;
                strstrm << NekConstants::kVertexTheSameDouble << " in edge "
                        << m_globalID;
                err += strstrm.str();
                NEKERROR(ErrorUtil::ewarning, err.c_str());
            }

            if (m_verts[1]->dist(*(m_curve->m_points[npts - 1])) >
                NekConstants::kVertexTheSameDouble)
            {
                std::string err =
                    "Vertex 1 is separated from last point by more than ";
                std::stringstream strstrm;
                strstrm << NekConstants::kVertexTheSameDouble << " in edge "
                        << m_globalID;
                err += strstrm.str();
                NEKERROR(ErrorUtil::ewarning, err.c_str());
            }

            LibUtilities::PointsKey fkey(npts, m_curve->m_ptype);
            DNekMatSharedPtr I0 =
                LibUtilities::PointsManager()[fkey]->GetI(pkey);
            NekVector<NekDouble> out(npts + 1);

            for (int i = 0; i < m_coordim; ++i)
            {
                // Load up coordinate values into tmp
                for (int j = 0; j < npts; ++j)
                {
                    tmp[j] = (m_curve->m_points[j]->GetPtr())[i];
                }

                // Interpolate to GLL points
                NekVector<NekDouble> in(npts, tmp, eWrapper);
                out = (*I0) * in;

                m_xmap->FwdTrans(out.GetPtr(), m_coeffs[i]);
            }
        }

        for (i = 0; i < m_coordim; ++i)
        {
            m_coeffs[i][0] = (*m_verts[0])[i];
            m_coeffs[i][1] = (*m_verts[1])[i];
        }

        m_state = ePtsFilled;
    }
}

void SegGeom::v_Reset(CurveMap &curvedEdges, CurveMap &curvedFaces)
{
    Geometry::v_Reset(curvedEdges, curvedFaces);
    CurveMap::iterator it = curvedEdges.find(m_globalID);

    if (it != curvedEdges.end())
    {
        m_curve = it->second;
    }

    SetUpXmap();
    SetUpCoeffs(m_xmap->GetNcoeffs());
}

void SegGeom::v_Setup()
{
    if(!m_setupState)
    {
        SetUpXmap();
        SetUpCoeffs(m_xmap->GetNcoeffs());
        m_setupState = true;
    }
}

NekDouble SegGeom::v_GetLocCoords(const Array<OneD, const NekDouble> &coords,
                                  Array<OneD, NekDouble> &Lcoords)
{
    int i;
    NekDouble resid = 0.0;
    SegGeom::v_FillGeom();

    // calculate local coordinate for coord
    if (GetMetricInfo()->GetGtype() == eRegular)
    {
        NekDouble len = 0.0;
        NekDouble xi = 0.0;

        const int npts = m_xmap->GetTotPoints();
        Array<OneD, NekDouble> pts(npts);

        for (i = 0; i < m_coordim; ++i)
        {
            m_xmap->BwdTrans(m_coeffs[i], pts);
            len += (pts[npts - 1] - pts[0]) * (pts[npts - 1] - pts[0]);
            xi += (coords[i] - pts[0]) * (pts[npts-1] - pts[0]);
        }

        Lcoords[0] = 2 * xi / len - 1.0;
    }
    else
    {
        NEKERROR(ErrorUtil::efatal,
                 "inverse mapping must be set up to use this call");
    }
    return resid;
}

bool SegGeom::v_ContainsPoint(const Array<OneD, const NekDouble> &gloCoord,
                              Array<OneD, NekDouble> &stdCoord,
                              NekDouble tol,
                              NekDouble &resid)
{
    //Rough check if within twice min/max point
    if (GetMetricInfo()->GetGtype() != eRegular)
    {
        if (!MinMaxCheck(gloCoord))
        {
            return false;
        }
    }

    // Convert to the local (eta) coordinates.
    resid = GetLocCoords(gloCoord, stdCoord);

    if (stdCoord[0] >= -(1 + tol) && stdCoord[0] <= 1 + tol)
    {
        return true;
    }

    //Clamp local coords
    ClampLocCoords(stdCoord, tol);

    return false;
}

PointGeomSharedPtr SegGeom::v_GetVertex(const int i) const
{
    PointGeomSharedPtr returnval;

    if (i >= 0 && i < kNverts)
    {
        returnval = m_verts[i];
    }

    return returnval;
}

int SegGeom::v_GetNumVerts() const
{
    return kNverts;
}

int SegGeom::v_GetNumEdges() const
{
    return kNedges;
}

}
}
