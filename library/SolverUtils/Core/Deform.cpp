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
// License for the specific language governing rights and limitations under
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

#ifndef NEKTAR_SOLVERUTILS_CORE_DEFORM_H
#define NEKTAR_SOLVERUTILS_CORE_DEFORM_H

#include <string>

#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <SolverUtils/SolverUtilsDeclspec.h>

namespace Nektar {
namespace SolverUtils {

    /**
     * @brief Update geometry according to displacement that is in current
     * fields.
     */
    void UpdateGeometry(
        SpatialDomains::MeshGraphSharedPtr           graph,
        Array<OneD, MultiRegions::ExpListSharedPtr> &exp)
    {
        SpatialDomains::CurveMap &curvedEdges = graph->GetCurvedEdges();
        SpatialDomains::CurveMap &curvedFaces = graph->GetCurvedFaces();
        curvedEdges.clear();
        curvedFaces.clear();

        int i, j, k, l, dim;
        set<int> updatedVerts, updatedEdges, updatedFaces;

        dim = graph->GetSpaceDimension();
        Array<OneD, Array<OneD, NekDouble> > phys (dim);
        Array<OneD, Array<OneD, NekDouble> > coord(dim);

        for (i = 0; i < fields[0]->GetExpSize(); ++i)
        {
            LocalRegions::ExpansionSharedPtr exp = fields[0]->GetExp(i);
            int offset = fields[0]->GetPhys_Offset(i);
            int nquad  = exp->GetTotPoints();

            for (j = 0; j < dim; ++j)
            {
                phys[j] = Array<OneD, NekDouble>(
                    nquad, fields[j]->UpdatePhys() + offset);
                coord[j] = Array<OneD, NekDouble>(nquad);
            }

            // In 2D loop over edges. 3D TODO
            if (dim == 2)
            {
                exp->GetCoords(coord[0], coord[1]);

                SpatialDomains::Geometry2DSharedPtr geom =
                    boost::dynamic_pointer_cast<SpatialDomains::Geometry2D>(
                        exp->GetGeom());

                for (j = 0; j < exp->GetNedges(); ++j)
                {
                    SpatialDomains::Geometry1DSharedPtr edge = geom->GetEdge(j);

                    // This edge has already been processed.
                    if (updatedEdges.find(edge->GetGlobalID()) != updatedEdges.end())
                    {
                        continue;
                    }

                    // Extract edge displacement.
                    int nEdgePts = exp->GetEdgeNumPoints(j);
                    Array<OneD, Array<OneD, NekDouble> > edgePhys(dim);
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
                    boost::dynamic_pointer_cast<SpatialDomains::Geometry3D>(
                        exp->GetGeom());

                for (j = 0; j < exp->GetNfaces(); ++j)
                {
                    SpatialDomains::Geometry2DSharedPtr face = geom->GetFace(j);

                    // This edge has already been processed.
                    if (updatedFaces.find(face->GetGlobalID()) != updatedFaces.end())
                    {
                        continue;
                    }

                    // Extract face displacement.
                    int nFacePts = exp->GetFaceNumPoints(j);
                    int nq0 = exp->GetNumPoints(0);
                    int nq1 = exp->GetNumPoints(1);
                    Array<OneD, Array<OneD, NekDouble> > facePhys(dim);
                    Array<OneD, Array<OneD, NekDouble> > faceCoord(dim);

                    const LibUtilities::BasisKey B0(
                        LibUtilities::eModified_A, nq0,
                        LibUtilities::PointsKey(
                            nq0, LibUtilities::eGaussLobattoLegendre));
                    const LibUtilities::BasisKey B1(
                        LibUtilities::eModified_A, nq1,
                        LibUtilities::PointsKey(
                            nq0, LibUtilities::eGaussLobattoLegendre));
                    StdRegions::StdExpansion2DSharedPtr faceexp;

                    if (geom->GetShapeType() == LibUtilities::eTriangle)
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
                        facePhys [k] = Array<OneD, NekDouble>(nFacePts);
                        faceCoord[k] = Array<OneD, NekDouble>(nFacePts);
                        exp->GetFacePhysVals(j, faceexp, phys [k], facePhys [k], exp->GetForient(j));
                        exp->GetFacePhysVals(j, faceexp, coord[k], faceCoord[k], exp->GetForient(j));
                    }

                    int edgeOff[2][4][2] = {
                        {
                            {0,           1},
                            {nq0-1,       nq0},
                            {nq0*(nq1-1), -nq0},
                            {-1,-1}
                        },
                        {
                            {0,           1},
                            {nq0-1,       nq0},
                            {nq0*nq1-1,   -1},
                            {nq0*(nq1-1), -nq0}
                        }
                    };

                    for (k = 0; k < face->GetNumVerts(); ++k)
                    {
                        // Update verts
                        int id = face->GetVid(k);
                        const int o = face->GetShapeType() - LibUtilities::eTriangle;

                        if (updatedVerts.find(id) == updatedVerts.end())
                        {
                            SpatialDomains::PointGeomSharedPtr pt =
                                face->GetVertex(k);

                            pt->UpdatePosition(
                                (*pt)(0) + facePhys[0][edgeOff[o][k][0]],
                                (*pt)(1) + facePhys[1][edgeOff[o][k][0]],
                                (*pt)(2) + facePhys[2][edgeOff[o][k][0]]);

                            updatedVerts.insert(id);
                        }

                        // Update edges
                        id = face->GetEid(k);
                        if (updatedEdges.find(id) == updatedEdges.end())
                        {
                            SpatialDomains::Geometry1DSharedPtr edge = face->GetEdge(k);
                            SpatialDomains::CurveSharedPtr curve = MemoryManager<
                                SpatialDomains::Curve>::AllocateSharedPtr(
                                    edge->GetGlobalID(),
                                    LibUtilities::eGaussLobattoLegendre);

                            int nEdgePts;
                            if (face->GetNumVerts() == 3)
                            {
                                nEdgePts = k > 0 ? nq1 : nq0;
                            }
                            else
                            {
                                nEdgePts = k % 2 ? nq1 : nq0;
                            }

                            const int offset = edgeOff[o][k][0];
                            const int pos    = edgeOff[o][k][1];

                            if (face->GetEorient(k) == StdRegions::eBackwards)
                            {
                                for (l = nEdgePts-1; l >= 0; --l)
                                {
                                    int m = offset + pos*l;
                                    SpatialDomains::PointGeomSharedPtr vert =
                                        MemoryManager<SpatialDomains::PointGeom>
                                        ::AllocateSharedPtr(
                                            dim, edge->GetGlobalID(),
                                            faceCoord[0][m] + facePhys[0][m],
                                            faceCoord[1][m] + facePhys[1][m],
                                            faceCoord[2][m] + facePhys[2][m]);
                                    curve->m_points.push_back(vert);
                                }
                            }
                            else
                            {
                                for (l = 0; l < nEdgePts; ++l)
                                {
                                    int m = offset + pos*l;
                                    SpatialDomains::PointGeomSharedPtr vert =
                                        MemoryManager<SpatialDomains::PointGeom>
                                        ::AllocateSharedPtr(
                                            dim, edge->GetGlobalID(),
                                            faceCoord[0][m] + facePhys[0][m],
                                            faceCoord[1][m] + facePhys[1][m],
                                            faceCoord[2][m] + facePhys[2][m]);
                                    curve->m_points.push_back(vert);
                                }
                            }

                            curvedEdges[edge->GetGlobalID()] = curve;
                            updatedEdges.insert(edge->GetGlobalID());
                        }
                    }

                    // Update face-interior curvature
                    SpatialDomains::CurveSharedPtr curve = MemoryManager<
                        SpatialDomains::Curve>::AllocateSharedPtr(
                            face->GetGlobalID(),
                            LibUtilities::eGaussLobattoLegendre);

                    for (l = 0; l < nFacePts; ++l)
                    {
                        SpatialDomains::PointGeomSharedPtr vert =
                            MemoryManager<SpatialDomains::PointGeom>
                            ::AllocateSharedPtr(
                                dim, face->GetGlobalID(),
                                faceCoord[0][l] + facePhys[0][l],
                                faceCoord[1][l] + facePhys[1][l],
                                faceCoord[2][l] + facePhys[2][l]);
                        curve->m_points.push_back(vert);
                    }

                    curvedFaces[face->GetGlobalID()] = curve;
                    updatedFaces.insert(face->GetGlobalID());
                }
            }
        }

        // Reset geometry information
        for (i = 0; i < fields.num_elements(); ++i)
        {
            fields[i]->Reset();
        }
    }

}
}

#endif
