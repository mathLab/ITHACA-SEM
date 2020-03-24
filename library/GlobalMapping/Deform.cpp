///////////////////////////////////////////////////////////////////////////////
//
// File: Deform.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Deformation of mesh from fields.
//
///////////////////////////////////////////////////////////////////////////////

#include <string>

#include <LibUtilities/Foundations/Interp.h>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <GlobalMapping/Deform.h>
#include <StdRegions/StdNodalTriExp.h>
#include <StdRegions/StdSegExp.h>
#include <StdRegions/StdQuadExp.h>
#include <SpatialDomains/MeshGraph.h>
#include <MultiRegions/ExpList.h>

using namespace std;

namespace Nektar {
namespace GlobalMapping {

    /**
     * @brief Update geometry according to displacement that is in current
     * fields.
     *
     * @param graph   The MeshGraph of the current geometry.
     * @param fields  The fields containing the displacement.
     */
    void UpdateGeometry(
        SpatialDomains::MeshGraphSharedPtr           graph,
        Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        bool modal)
    {
        // Clear existing curvature.
        SpatialDomains::CurveMap &curvedEdges = graph->GetCurvedEdges();
        SpatialDomains::CurveMap &curvedFaces = graph->GetCurvedFaces();
        curvedEdges.clear();
        curvedFaces.clear();

        int i, j, k, l, dim;

        // Sets to hold IDs of updated vertices to avoid duplicating effort.
        set<int> updatedVerts, updatedEdges, updatedFaces;

        dim = graph->GetSpaceDimension();
        Array<OneD, Array<OneD, NekDouble> > phys (dim);
        Array<OneD, Array<OneD, NekDouble> > coord(dim);

        for (i = 0; i < fields[0]->GetExpSize(); ++i)
        {
            LocalRegions::ExpansionSharedPtr exp = fields[0]->GetExp(i);
            int offset = fields[0]->GetPhys_Offset(i);
            int nquad  = exp->GetTotPoints();

            // Extract displacement for this element, allocate storage for
            // elemental coordinates.
            for (j = 0; j < dim; ++j)
            {
                phys[j] = Array<OneD, NekDouble>(
                    nquad, fields[j]->UpdatePhys() + offset);
                coord[j] = Array<OneD, NekDouble>(nquad);
            }

            // In 2D loop over edges.
            if (dim == 2)
            {
                exp->GetCoords(coord[0], coord[1]);

                SpatialDomains::Geometry2DSharedPtr geom =
                    std::dynamic_pointer_cast<SpatialDomains::Geometry2D>(
                        exp->GetGeom());

                for (j = 0; j < exp->GetNedges(); ++j)
                {
                    SpatialDomains::Geometry1DSharedPtr edge = geom->GetEdge(j);

                    // This edge has already been processed.
                    if (updatedEdges.find(edge->GetGlobalID()) !=
                            updatedEdges.end())
                    {
                        continue;
                    }

                    // Extract edge displacement.
                    int nEdgePts = exp->GetEdgeNumPoints(j);
                    Array<OneD, Array<OneD, NekDouble> > edgePhys (dim);
                    Array<OneD, Array<OneD, NekDouble> > edgeCoord(dim);

                    const LibUtilities::BasisKey B(
                        LibUtilities::eModified_A, nEdgePts,
                        LibUtilities::PointsKey(
                            nEdgePts, LibUtilities::eGaussLobattoLegendre));
                    StdRegions::StdExpansion1DSharedPtr seg = MemoryManager<
                        StdRegions::StdSegExp>::AllocateSharedPtr(B);

                    for (k = 0; k < dim; ++k)
                    {
                        edgePhys [k] = Array<OneD, NekDouble>(nEdgePts);
                        edgeCoord[k] = Array<OneD, NekDouble>(nEdgePts);
                        exp->GetEdgePhysVals(j, seg, phys [k], edgePhys [k]);
                        exp->GetEdgePhysVals(j, seg, coord[k], edgeCoord[k]);
                    }

                    // Update verts
                    for (k = 0; k < 2; ++k)
                    {
                        int id = edge->GetVid(k);
                        if (updatedVerts.find(id) != updatedVerts.end())
                        {
                            continue;
                        }

                        SpatialDomains::PointGeomSharedPtr pt =
                            edge->GetVertex(k);

                        pt->UpdatePosition(
                            (*pt)(0) + edgePhys[0][k*(nEdgePts-1)],
                            (*pt)(1) + edgePhys[1][k*(nEdgePts-1)],
                            (*pt)(2));

                        updatedVerts.insert(id);
                    }

                    // Update curve
                    SpatialDomains::CurveSharedPtr curve = MemoryManager<
                        SpatialDomains::Curve>::AllocateSharedPtr(
                            edge->GetGlobalID(),
                            LibUtilities::eGaussLobattoLegendre);

                    for (k = 0; k < nEdgePts; ++k)
                    {
                        SpatialDomains::PointGeomSharedPtr vert =
                            MemoryManager<SpatialDomains::PointGeom>
                            ::AllocateSharedPtr(
                                dim, edge->GetGlobalID(),
                                edgeCoord[0][k] + edgePhys[0][k],
                                edgeCoord[1][k] + edgePhys[1][k], 0.0);

                        curve->m_points.push_back(vert);
                    }

                    curvedEdges[edge->GetGlobalID()] = curve;
                    updatedEdges.insert(edge->GetGlobalID());
                }
            }
            else if (dim == 3)
            {
                exp->GetCoords(coord[0], coord[1], coord[2]);

                SpatialDomains::Geometry3DSharedPtr geom =
                    std::dynamic_pointer_cast<SpatialDomains::Geometry3D>(
                        exp->GetGeom());

                for (j = 0; j < exp->GetNfaces(); ++j)
                {
                    SpatialDomains::Geometry2DSharedPtr face = geom->GetFace(j);

                    // This edge has already been processed.
                    if (updatedFaces.find(face->GetGlobalID()) !=
                            updatedFaces.end())
                    {
                        continue;
                    }

                    // Extract face displacement.
                    LibUtilities::BasisKey B0 = exp->DetFaceBasisKey(j,0);
                    LibUtilities::BasisKey B1 = exp->DetFaceBasisKey(j,1);
                    int nq0 = B0.GetNumPoints();
                    int nq1 = B1.GetNumPoints();

                    ASSERTL1(B0.GetPointsType()
                                 == LibUtilities::eGaussLobattoLegendre &&
                             B1.GetPointsType()
                                 == LibUtilities::eGaussLobattoLegendre,
                             "Deformation requires GLL points in both "
                             "directions on a face.");

                    Array<OneD, Array<OneD, NekDouble> > newPos(dim);

                    StdRegions::StdExpansion2DSharedPtr faceexp;
                    StdRegions::Orientation orient = exp->GetForient(j);

                    if (face->GetShapeType() == LibUtilities::eTriangle)
                    {
                        faceexp = MemoryManager<StdRegions::StdTriExp>::
                            AllocateSharedPtr(B0, B1);
                    }
                    else
                    {
                        faceexp = MemoryManager<StdRegions::StdQuadExp>::
                            AllocateSharedPtr(B0, B1);
                    }

                    for (k = 0; k < dim; ++k)
                    {
                        Array<OneD, NekDouble> tmp(nq0*nq1);
                        newPos[k] = Array<OneD, NekDouble>(nq0*nq1);
                        exp->GetFacePhysVals(
                            j, faceexp, phys [k], tmp,       orient);
                        exp->GetFacePhysVals(
                            j, faceexp, coord[k], newPos[k], orient);
                        Vmath::Vadd(
                            nq0*nq1, tmp, 1, newPos[k], 1, newPos[k], 1);
                    }

                    // Now interpolate face onto a more reasonable set of
                    // points.
                    int nq = max(nq0, nq1);
                    if(!modal)
                        nq--;

                    LibUtilities::PointsKey edgePts(
                        nq, LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::PointsKey triPts(
                        nq, LibUtilities::eNodalTriElec);

                    Array<OneD, Array<OneD, NekDouble> > intPos(dim);

                    for (k = 0; k < dim; ++k)
                    {
                        intPos[k] = Array<OneD, NekDouble>(nq*nq);
                        LibUtilities::Interp2D(
                            faceexp->GetPointsKeys()[0],
                            faceexp->GetPointsKeys()[1],
                            newPos[k], edgePts, edgePts, intPos[k]);
                    }

                    int edgeOff[2][4][2] = {
                        {
                            {0,           1},
                            {nq-1,       nq},
                            {nq*(nq-1), -nq},
                            {-1,-1}
                        },
                        {
                            {0,           1},
                            {nq-1,       nq},
                            {nq*nq-1,    -1},
                            {nq*(nq-1), -nq}
                        }
                    };

                    for (k = 0; k < face->GetNumVerts(); ++k)
                    {
                        // Update verts
                        int id = face->GetVid(k);
                        const int o =
                            face->GetShapeType() - LibUtilities::eTriangle;

                        if (updatedVerts.find(id) == updatedVerts.end())
                        {
                            SpatialDomains::PointGeomSharedPtr pt =
                                face->GetVertex(k);
                            pt->UpdatePosition(
                                intPos[0][edgeOff[o][k][0]],
                                intPos[1][edgeOff[o][k][0]],
                                intPos[2][edgeOff[o][k][0]]);
                            updatedVerts.insert(id);
                        }

                        // Update edges
                        id = face->GetEid(k);
                        if (updatedEdges.find(id) == updatedEdges.end())
                        {
                            SpatialDomains::Geometry1DSharedPtr edge
                                = face->GetEdge(k);
                            SpatialDomains::CurveSharedPtr curve =
                                MemoryManager<SpatialDomains::Curve>
                                ::AllocateSharedPtr(
                                    edge->GetGlobalID(),
                                    LibUtilities::eGaussLobattoLegendre);

                            const int offset = edgeOff[o][k][0];
                            const int pos    = edgeOff[o][k][1];

                            if (face->GetEorient(k) == StdRegions::eBackwards)
                            {
                                for (l = nq-1; l >= 0; --l)
                                {
                                    int m = offset + pos*l;
                                    SpatialDomains::PointGeomSharedPtr vert =
                                        MemoryManager<SpatialDomains::PointGeom>
                                        ::AllocateSharedPtr(
                                            dim, edge->GetGlobalID(),
                                            intPos[0][m], intPos[1][m],
                                            intPos[2][m]);
                                    curve->m_points.push_back(vert);
                                }
                            }
                            else
                            {
                                for (l = 0; l < nq; ++l)
                                {
                                    int m = offset + pos*l;
                                    SpatialDomains::PointGeomSharedPtr vert =
                                        MemoryManager<SpatialDomains::PointGeom>
                                        ::AllocateSharedPtr(
                                            dim, edge->GetGlobalID(),
                                            intPos[0][m], intPos[1][m],
                                            intPos[2][m]);
                                    curve->m_points.push_back(vert);
                                }
                            }

                            curvedEdges[edge->GetGlobalID()] = curve;
                            updatedEdges.insert(edge->GetGlobalID());
                        }
                    }

                    // Update face-interior curvature
                    LibUtilities::PointsType pType =
                        face->GetShapeType() == LibUtilities::eTriangle ?
                        LibUtilities::eNodalTriElec :
                        LibUtilities::eGaussLobattoLegendre;

                    SpatialDomains::CurveSharedPtr curve = MemoryManager<
                        SpatialDomains::Curve>::AllocateSharedPtr(
                            face->GetGlobalID(),
                            pType);

                    if (face->GetShapeType() == LibUtilities::eTriangle)
                    {
                        // This code is probably pretty crappy. Have to go from
                        // GLL-GLL points -> GLL-Gauss-Radau -> nodal triangle
                        // points.
                        const LibUtilities::BasisKey B0(
                            LibUtilities::eOrtho_A, nq,
                            LibUtilities::PointsKey(
                                nq, LibUtilities::eGaussLobattoLegendre));
                        const LibUtilities::BasisKey B1(
                            LibUtilities::eOrtho_B, nq,
                            LibUtilities::PointsKey(
                                nq, LibUtilities::eGaussRadauMAlpha1Beta0));
                        StdRegions::StdNodalTriExp nodalTri(B0, B1, pType);
                        StdRegions::StdTriExp      tri     (B0, B1);

                        for (k = 0; k < dim; ++k)
                        {
                            Array<OneD, NekDouble> nodal(nq*nq);

                            LibUtilities::Interp2D(
                                faceexp->GetBasis(0)->GetBasisKey(),
                                faceexp->GetBasis(1)->GetBasisKey(),
                                newPos[k], B0, B1, nodal);

                            Array<OneD, NekDouble> tmp1(nq*(nq+1)/2);
                            Array<OneD, NekDouble> tmp2(nq*(nq+1)/2);

                            tri.FwdTrans(nodal, tmp1);
                            nodalTri.ModalToNodal(tmp1, tmp2);
                            newPos[k] = tmp2;
                        }

                        for (l = 0; l < nq*(nq+1)/2; ++l)
                        {
                            SpatialDomains::PointGeomSharedPtr vert =
                                MemoryManager<SpatialDomains::PointGeom>
                                ::AllocateSharedPtr(
                                    dim, face->GetGlobalID(),
                                    newPos[0][l], newPos[1][l], newPos[2][l]);
                            curve->m_points.push_back(vert);
                        }
                    }
                    else
                    {
                        for (l = 0; l < nq*nq; ++l)
                        {
                            SpatialDomains::PointGeomSharedPtr vert =
                                MemoryManager<SpatialDomains::PointGeom>
                                ::AllocateSharedPtr(
                                    dim, face->GetGlobalID(),
                                    intPos[0][l], intPos[1][l], intPos[2][l]);
                            curve->m_points.push_back(vert);
                        }
                    }

                    curvedFaces[face->GetGlobalID()] = curve;
                    updatedFaces.insert(face->GetGlobalID());
                }
            }
        }

        // Reset geometry information
        for (i = 0; i < fields.size(); ++i)
        {
            fields[i]->Reset();
        }
    }
}
}
