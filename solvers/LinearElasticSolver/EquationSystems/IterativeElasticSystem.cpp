///////////////////////////////////////////////////////////////////////////////
//
// File IterativeElasticSystem.cpp
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
// Description: IterativeElasticSystem solve routines 
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/format.hpp>

#include <LibUtilities/Foundations/Interp.h>
#include <LibUtilities/BasicUtils/FileSystem.h>
#include <StdRegions/StdQuadExp.h>
#include <StdRegions/StdTriExp.h>
#include <LocalRegions/MatrixKey.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/GlobalLinSysDirectStaticCond.h>

#include <LinearElasticSolver/EquationSystems/IterativeElasticSystem.h>

namespace Nektar
{
    string IterativeElasticSystem::className = GetEquationSystemFactory().
        RegisterCreatorFunction("IterativeElasticSystem",
                                IterativeElasticSystem::create);

    IterativeElasticSystem::IterativeElasticSystem(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        : LinearElasticSystem(pSession)
    {
    }

    void IterativeElasticSystem::v_InitObject()
    {
        LinearElasticSystem::v_InitObject();

        const int nVel = m_fields[0]->GetCoordim(0);

        // Read in number of steps to take.
        m_session->LoadParameter("NumSteps", m_numSteps, 0);
        ASSERTL0(m_numSteps > 0, "You must specify at least one step");

        // Read in whether to repeatedly apply boundary conditions (for e.g.
        // rotation purposes).
        string bcType;
        m_session->LoadSolverInfo("BCType", bcType, "Normal");
        m_repeatBCs = bcType != "Normal";

        if (!m_repeatBCs)
        {
            // Loop over BCs, identify which ones we need to deform.
            const Array<OneD, const SpatialDomains::BoundaryConditionShPtr>
                &bndCond = m_fields[0]->GetBndConditions();

            for (int i = 0; i < bndCond.num_elements(); ++i)
            {
                if (bndCond[i]->GetUserDefined() == SpatialDomains::eWall)
                {
                    m_toDeform.push_back(i);
                }
            }

            ASSERTL0(m_toDeform.size() > 0, "Must tag at least one boundary "
                                            "with the WALL user-defined type");

            m_savedBCs  = Array<OneD, Array<OneD, Array<OneD, NekDouble> > >(
                m_toDeform.size());

            for (int i = 0; i < m_toDeform.size(); ++i)
            {
                m_savedBCs[i] = Array<OneD, Array<OneD, NekDouble> >(nVel);
                for (int j = 0; j < nVel; ++j)
                {
                    const int id = m_toDeform[i];
                    MultiRegions::ExpListSharedPtr bndCondExp =
                        m_fields[j]->GetBndCondExpansions()[id];
                    int nCoeffs = bndCondExp->GetNcoeffs();

                    m_savedBCs[i][j] = Array<OneD, NekDouble>(nCoeffs);
                    Vmath::Smul(nCoeffs, 1.0/m_numSteps,
                                bndCondExp->GetCoeffs(),    1,
                                bndCondExp->UpdateCoeffs(), 1);
                    Vmath::Vcopy(nCoeffs, bndCondExp->GetCoeffs(), 1,
                                 m_savedBCs[i][j], 1);
                }
            }
        }
    }

    void IterativeElasticSystem::v_GenerateSummary(SolverUtils::SummaryList& s)
    {
        LinearElasticSystem::v_GenerateSummary(s);
    }

    void IterativeElasticSystem::v_DoSolve()
    {
        int i, j, k;

        // Write initial geometry for consistency/script purposes
        WriteGeometry(0);

        // Now loop over desired number of steps
        for (i = 1; i <= m_numSteps; ++i)
        {
            int invalidElmtId = -1;

            // Perform solve for this iteration and update geometry accordingly.
            LinearElasticSystem::v_DoSolve();
            UpdateGeometry();

            // Check for invalid elements.
            for (j = 0; j < m_fields[0]->GetExpSize(); ++j)
            {
                SpatialDomains::GeomFactorsSharedPtr geomFac =
                    m_fields[0]->GetExp(j)->GetGeom()->GetGeomFactors();

                if (!geomFac->IsValid())
                {
                    invalidElmtId =
                        m_fields[0]->GetExp(j)->GetGeom()->GetGlobalID();
                    break;
                }
            }

            m_session->GetComm()->AllReduce(invalidElmtId, LibUtilities::ReduceMax);

            // If we found an invalid element, exit loop without writing output.
            if (invalidElmtId >= 0)
            {
                if (m_session->GetComm()->GetRank() == 0)
                {
                    cout << "- Detected negative Jacobian in element "
                         << invalidElmtId << "; terminating" << endl;
                }

                break;
            }

            WriteGeometry(i);

            if (m_session->GetComm()->GetRank() == 0)
            {
                cout << "Step: " << i << endl;
            }

            // Update boundary conditions
            if (m_repeatBCs)
            {
                for (j = 0; j < m_fields.num_elements(); ++j)
                {
                    string varName = m_session->GetVariable(j);
                    m_fields[j]->EvaluateBoundaryConditions(m_time, varName);
                }
            }
            else
            {
                for (j = 0; j < m_fields.num_elements(); ++j)
                {
                    const Array<OneD, const MultiRegions::ExpListSharedPtr> &bndCondExp =
                        m_fields[j]->GetBndCondExpansions();

                    for (k = 0; k < m_toDeform.size(); ++k)
                    {
                        const int id = m_toDeform[k];
                        const int nCoeffs = bndCondExp[id]->GetNcoeffs();
                        Vmath::Vcopy(nCoeffs,
                                     m_savedBCs[k][j],               1,
                                     bndCondExp[id]->UpdateCoeffs(), 1);
                    }
                }
            }
        }
    }

    /**
     * @brief Write out a file in serial or directory in parallel containing new
     * mesh geometry.
     */
    void IterativeElasticSystem::WriteGeometry(const int i)
    {
        fs::path filename;
        stringstream s;
        s << m_session->GetSessionName() << "-" << i;

        if (m_session->GetComm()->GetSize() > 1)
        {
            s << "_xml";

            if(!fs::is_directory(s.str()))
            {
                fs::create_directory(s.str());
            }

            boost::format pad("P%1$07d.xml");
            pad % m_session->GetComm()->GetRank();
            filename = fs::path(s.str()) / fs::path(pad.str());
        }
        else
        {
            s << ".xml";
            filename = fs::path(s.str());
        }

        string fname = LibUtilities::PortablePath(filename);
        m_fields[0]->GetGraph()->WriteGeometry(fname);
    }

    /**
     * @brief Update geometry according to displacement that is in current
     * fields.
     */
    void IterativeElasticSystem::UpdateGeometry()
    {
        SpatialDomains::MeshGraphSharedPtr graph = m_fields[0]->GetGraph();
        SpatialDomains::CurveVector &curvedEdges = graph->GetCurvedEdges();
        SpatialDomains::CurveVector &curvedFaces = graph->GetCurvedFaces();

        int i, j, k, l, dim;
        set<int> updatedVerts, updatedEdges, updatedFaces;

        dim = graph->GetSpaceDimension();
        Array<OneD, Array<OneD, NekDouble> > phys(dim);
        Array<OneD, Array<OneD, NekDouble> > coord(dim);
        Array<OneD, Array<OneD, NekDouble> > coordtmp(dim);

        map<int, SpatialDomains::CurveSharedPtr> curveEdgeMap, curveFaceMap;
        map<int, SpatialDomains::CurveSharedPtr>::iterator it;

        for (i = 0; i < curvedEdges.size(); ++i)
        {
            curveEdgeMap[curvedEdges[i]->m_curveID] = curvedEdges[i];
        }

        for (i = 0; i < curvedFaces.size(); ++i)
        {
            curveFaceMap[curvedFaces[i]->m_curveID] = curvedFaces[i];
        }

        for (i = 0; i < m_fields[0]->GetExpSize(); ++i)
        {
            LocalRegions::ExpansionSharedPtr exp = m_fields[0]->GetExp(i);
            LibUtilities::PointsKeyVector pKeyFrom = exp->GetPointsKeys();
            LibUtilities::PointsKeyVector pKeyTo(dim);

            int offset = m_fields[0]->GetPhys_Offset(i);
            int nquad  = exp->GetTotPoints();

            for (j = 0; j < dim; ++j)
            {
                pKeyTo[j] = LibUtilities::PointsKey(
                    pKeyFrom[j].GetNumPoints(),
                    LibUtilities::eGaussLobattoLegendre);
            }

            for (j = 0; j < dim; ++j)
            {
                Array<OneD, NekDouble> tmp(
                    nquad, m_fields[j]->UpdatePhys() + offset);
                phys[j] = Array<OneD, NekDouble>(nquad);

                if (dim == 2)
                {
                    LibUtilities::Interp2D(
                        pKeyFrom[0], pKeyFrom[1], tmp,
                        pKeyTo  [0], pKeyTo  [1], phys[j]);
                }

                coord   [j] = Array<OneD, NekDouble>(nquad);
                coordtmp[j] = Array<OneD, NekDouble>(nquad);
            }

            exp->GetCoords(coordtmp[0], coordtmp[1]);

            // In 2D loop over edges. 3D TODO
            if (dim == 2)
            {
                SpatialDomains::Geometry2DSharedPtr geom =
                    boost::dynamic_pointer_cast<SpatialDomains::Geometry2D>(
                        exp->GetGeom());

                // Do an interpolation to GLL points to avoid issues with
                // triangles and curved edges which may use Gauss-Radau points.
                LibUtilities::Interp2D(
                    pKeyFrom[0], pKeyFrom[1], coordtmp[0],
                    pKeyTo  [0], pKeyTo  [1], coord   [0]);
                LibUtilities::Interp2D(
                    pKeyFrom[0], pKeyFrom[1], coordtmp[1],
                    pKeyTo  [0], pKeyTo  [1], coord   [1]);

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

                    for (k = 0; k < dim; ++k)
                    {
                        edgePhys [k] = Array<OneD, NekDouble>(nEdgePts);
                        edgeCoord[k] = Array<OneD, NekDouble>(nEdgePts);
                        exp->GetEdgePhysVals(j, phys [k], edgePhys [k]);
                        exp->GetEdgePhysVals(j, coord[k], edgeCoord[k]);
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

                    // Update curve: do interpolation to GaussLobattoLegendre
                    // points to avoid warnings; probably not necessary...
                    SpatialDomains::CurveSharedPtr curve;
                    it = curveEdgeMap.find(edge->GetGlobalID());
                    if (it == curveEdgeMap.end())
                    {
                        curve = MemoryManager<
                            SpatialDomains::Curve>::AllocateSharedPtr(
                                edge->GetGlobalID(),
                                LibUtilities::eGaussLobattoLegendre);
                    }
                    else
                    {
                        curve = it->second;
                        curve->m_ptype = LibUtilities::eGaussLobattoLegendre;
                        curve->m_points.clear();
                    }

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

                    if (it == curveEdgeMap.end())
                    {
                        curvedEdges.push_back(curve);
                    }

                    updatedEdges.insert(edge->GetGlobalID());
                }
            }
            else if (dim == 3)
            {
                SpatialDomains::Geometry3DSharedPtr geom =
                    boost::dynamic_pointer_cast<SpatialDomains::Geometry3D>(
                        exp->GetGeom());

                // Do an interpolation to GLL points to avoid issues with
                // triangles and curved edges which may use Gauss-Radau points.
                LibUtilities::Interp3D(
                    pKeyFrom[0], pKeyFrom[1], pKeyFrom[2], coordtmp[0],
                    pKeyTo  [0], pKeyTo  [1], pKeyFrom[2], coord   [0]);
                LibUtilities::Interp3D(
                    pKeyFrom[0], pKeyFrom[1], pKeyFrom[2], coordtmp[1],
                    pKeyTo  [0], pKeyTo  [1], pKeyFrom[2], coord   [1]);
                LibUtilities::Interp3D(
                    pKeyFrom[0], pKeyFrom[1], pKeyFrom[2], coordtmp[2],
                    pKeyTo  [0], pKeyTo  [1], pKeyFrom[2], coord   [2]);

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
                        exp->GetFacePhysVals(j, faceexp, phys [k], facePhys [k], exp->GetFaceOrient(j));
                        exp->GetFacePhysVals(j, faceexp, coord[k], faceCoord[k], exp->GetFaceOrient(j));
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

                            SpatialDomains::CurveSharedPtr curve;
                            it = curveEdgeMap.find(edge->GetGlobalID());
                            if (it == curveEdgeMap.end())
                            {
                                curve = MemoryManager<
                                    SpatialDomains::Curve>::AllocateSharedPtr(
                                        edge->GetGlobalID(),
                                        LibUtilities::eGaussLobattoLegendre);
                            }
                            else
                            {
                                curve = it->second;
                                curve->m_ptype = LibUtilities::eGaussLobattoLegendre;
                                curve->m_points.clear();
                            }

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

                            curvedEdges.push_back(curve);
                            updatedEdges.insert(edge->GetGlobalID());
                        }
                    }

                    // Update face-interior curvature
                    SpatialDomains::CurveSharedPtr curve;
                    it = curveFaceMap.find(face->GetGlobalID());
                    if (it == curveFaceMap.end())
                    {
                        curve = MemoryManager<
                            SpatialDomains::Curve>::AllocateSharedPtr(
                                face->GetGlobalID(),
                                LibUtilities::eGaussLobattoLegendre);
                    }
                    else
                    {
                        curve = it->second;
                        curve->m_ptype = LibUtilities::eGaussLobattoLegendre;
                        curve->m_points.clear();
                    }

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

                    curvedFaces.push_back(curve);
                    updatedFaces.insert(face->GetGlobalID());
                }
            }
        }

        // Reset geometry information
        for (i = 0; i < m_fields.num_elements(); ++i)
        {
            m_fields[i]->Reset();
        }
    }
}
