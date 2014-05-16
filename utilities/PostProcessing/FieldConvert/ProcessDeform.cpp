////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessDeform.cpp
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
//  Description: Computes Q Criterion field.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
using namespace std;

#include "ProcessDeform.h"

#include <StdRegions/StdSegExp.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>

namespace Nektar
{
    namespace Utilities
    {
        ModuleKey ProcessDeform::className =
        GetModuleFactory().RegisterCreatorFunction(
            ModuleKey(eProcessModule, "deform"), ProcessDeform::create,
            "Deform a mesh given an input field defining displacement");
        
        ProcessDeform::ProcessDeform(FieldSharedPtr f) :
            ProcessModule(f)
        {
        }
        
        ProcessDeform::~ProcessDeform()
        {
        }
        
        void ProcessDeform::Process(po::variables_map &vm)
        {
            if (m_f->m_verbose)
            {
                cout << "ProcessDeform: Deforming grid..." << endl;
            }

            // Maybe create a new copy of MeshGraph to work on?
            SpatialDomains::MeshGraphSharedPtr graph = m_f->m_graph;
            SpatialDomains::CurveVector &curvedEdges = m_f->m_graph->GetCurvedEdges();
            curvedEdges.clear();

            int i, j, k, dim;
            set<int> updatedVerts, updatedEdges, updatedFaces;

            dim = graph->GetSpaceDimension();
            Array<OneD, Array<OneD, NekDouble> > phys(dim);
            Array<OneD, Array<OneD, NekDouble> > coord(dim);

            for (i = 0; i < m_f->m_exp[0]->GetExpSize(); ++i)
            {
                LocalRegions::ExpansionSharedPtr exp = m_f->m_exp[0]->GetExp(i);
                int offset = m_f->m_exp[0]->GetPhys_Offset(i);
                int nquad  = exp->GetTotPoints();

                for (j = 0; j < dim; ++j)
                {
                    phys[j] = Array<OneD, NekDouble>(
                        nquad, m_f->m_exp[j]->UpdatePhys() + offset);
                    coord[j] = Array<OneD, NekDouble>(nquad);
                }

                exp->GetCoords(coord[0], coord[1]);

                // In 2D loop over edges. 3D TODO
                if (dim == 2)
                {
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

                        curvedEdges.push_back(curve);

                        updatedEdges.insert(edge->GetGlobalID());
                    }
                }
#if 0
                else if (dim == 3)
                {
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
                        Array<OneD, Array<OneD, NekDouble> > facePhys(dim);
                        Array<OneD, Array<OneD, NekDouble> > faceCoord(dim);

                        for (k = 0; k < dim; ++k)
                        {
                            facePhys [k] = Array<OneD, NekDouble>(nFacePts);
                            faceCoord[k] = Array<OneD, NekDouble>(nFacePts);
                            exp->GetFacePhysVals(j, StdRegions::StdExpansionSharedPtr(), phys [k], facePhys [k], exp->GetFaceOrient(j));
                            exp->GetFacePhysVals(j, StdRegions::StdExpansionSharedPtr(), coord[k], faceCoord[k], exp->GetFaceOrient(j));
                        }

                        // Update edges

                        // Update verts
                        for (k = 0; k < geom->GetNverts(); ++k)
                        {
                            int id = face->GetVid(k);
                            if (updatedVerts.find(id) != updatedVerts.end())
                            {
                                continue;
                            }

                            SpatialDomains::PointGeomSharedPtr pt =
                                face->GetVertex(k);

                            pt->UpdatePosition(
                                (*pt)(0) + facePhys[0][k*(nFacePts-1)],
                                (*pt)(1) + facePhys[1][k*(nFacePts-1)],
                                (*pt)(2));

                            updatedVerts.insert(id);
                        }

                        int dir;
                        if (exp->GetNverts() == 3)
                        {
                            dir = j == 0 ? 0 : 1;
                        }
                        else
                        {
                            dir = j % 2;
                        }

                        // Update curve
                        SpatialDomains::CurveSharedPtr curve = MemoryManager<
                            SpatialDomains::Curve>::AllocateSharedPtr(
                                edge->GetGlobalID(),
                                exp->GetBasis(dir)->GetPointsType());

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

                        curvedEdges.push_back(curve);

                        updatedEdges.insert(edge->GetGlobalID());
                    }
                }
#endif
            }
        }
    }
}
