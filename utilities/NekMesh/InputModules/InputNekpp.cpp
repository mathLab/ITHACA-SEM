////////////////////////////////////////////////////////////////////////////////
//
//  File: InputNekpp.cpp
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
//  Description: GMSH converter.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
using namespace std;

#include <SpatialDomains/MeshGraph.h>
#include <NekMeshUtils/MeshElements/Element.h>
#include "InputNekpp.h"

namespace Nektar
{
namespace Utilities
{

using namespace Nektar::NekMeshUtils;

ModuleKey InputNekpp::className =
    GetModuleFactory().RegisterCreatorFunction(ModuleKey(eInputModule, "xml"),
                                               InputNekpp::create,
                                               "Reads Nektar++ xml file.");
/**
 * @brief Set up InputNekpp object.
 */
InputNekpp::InputNekpp(MeshSharedPtr m) : InputModule(m)
{
}

InputNekpp::~InputNekpp()
{
}

/**
 *
 */
void InputNekpp::Process()
{
    vector<string> filename;
    filename.push_back(m_config["infile"].as<string>());

    char *prgname = const_cast<char *>("NekMesh");
    LibUtilities::SessionReaderSharedPtr pSession =
        LibUtilities::SessionReader::CreateInstance(1, &prgname, filename);
    SpatialDomains::MeshGraphSharedPtr graph =
        SpatialDomains::MeshGraph::Read(pSession);

    auto comm = pSession->GetComm();
    if (comm->GetType().find("MPI") != std::string::npos)
    {
        m_mesh->m_comm = comm;
    }

    m_mesh->m_expDim   = graph->GetMeshDimension();
    m_mesh->m_spaceDim = graph->GetSpaceDimension();

    // Copy vertices.
    map<int, NodeSharedPtr> vIdMap;
    for (auto &vit : graph->GetAllPointGeoms())
    {
        SpatialDomains::PointGeomSharedPtr vert = vit.second;
        NodeSharedPtr n = std::make_shared<Node>(
            vert->GetGlobalID(), (*vert)(0), (*vert)(1), (*vert)(2));
        m_mesh->m_vertexSet.insert(n);
        vIdMap[vert->GetVid()] = n;
    }

    std::unordered_map<int, EdgeSharedPtr> eIdMap;
    std::unordered_map<int, FaceSharedPtr> fIdMap;

    // Load up all edges from graph
    {
        for (auto &it : graph->GetAllSegGeoms())
        {
            pair<int, SpatialDomains::GeometrySharedPtr> tmp2(
                it.first,
                std::dynamic_pointer_cast<SpatialDomains::Geometry>(
                    it.second));

            // load up edge set in order of SegGeomMap;
            vector<NodeSharedPtr> curve; // curved nodes if deformed
            int id0 = it.second->GetVid(0);
            int id1 = it.second->GetVid(1);
            LibUtilities::PointsType ptype =
                it.second->GetXmap()->GetPointsKeys()[0].GetPointsType();
            EdgeSharedPtr ed = std::make_shared<Edge>(
                vIdMap[id0], vIdMap[id1], curve, ptype);

            auto testIns = m_mesh->m_edgeSet.insert(ed);
            (*(testIns.first))->m_id = it.second->GetGlobalID();
            eIdMap[it.second->GetGlobalID()] = ed;
        }
    }

    // load up all faces from graph
    {
        for (auto &it : graph->GetAllTriGeoms())
        {
            pair<int, SpatialDomains::GeometrySharedPtr> tmp2(
                it.first,
                std::dynamic_pointer_cast<SpatialDomains::Geometry>(
                    it.second));

            vector<NodeSharedPtr> faceVertices;
            vector<EdgeSharedPtr> faceEdges;
            vector<NodeSharedPtr> faceNodes;

            for (int i = 0; i < 3; ++i)
            {
                faceVertices.push_back(vIdMap[it.second->GetVid(i)]);
                faceEdges.push_back(eIdMap[it.second->GetEid(i)]);
            }

            FaceSharedPtr fac = std::make_shared<Face>(
                faceVertices, faceNodes, faceEdges,
                LibUtilities::ePolyEvenlySpaced);
            auto testIns = m_mesh->m_faceSet.insert(fac);
            (*(testIns.first))->m_id = it.second->GetGlobalID();
            fIdMap[it.second->GetGlobalID()] = fac;
        }

        for (auto &it : graph->GetAllQuadGeoms())
        {
            pair<int, SpatialDomains::GeometrySharedPtr> tmp2(
                it.first,
                std::dynamic_pointer_cast<SpatialDomains::Geometry>(
                    it.second));

            vector<NodeSharedPtr> faceVertices;
            vector<EdgeSharedPtr> faceEdges;
            vector<NodeSharedPtr> faceNodes;

            for (int i = 0; i < 4; ++i)
            {
                faceVertices.push_back(vIdMap[it.second->GetVid(i)]);
                faceEdges.push_back(eIdMap[it.second->GetEid(i)]);
            }

            FaceSharedPtr fac = std::make_shared<Face>(
                faceVertices, faceNodes, faceEdges,
                LibUtilities::ePolyEvenlySpaced);
            auto testIns = m_mesh->m_faceSet.insert(fac);
            (*(testIns.first))->m_id = it.second->GetGlobalID();
            fIdMap[it.second->GetGlobalID()] = fac;
        }
    }

    // Set up curved information

    // Curved Edges
    for (auto &it : graph->GetCurvedEdges())
    {
        SpatialDomains::CurveSharedPtr curve = it.second;
        int id = curve->m_curveID;
        ASSERTL1(eIdMap.find(id) != eIdMap.end(), "Failed to find curved edge");
        EdgeSharedPtr edg = eIdMap[id];
        edg->m_curveType = curve->m_ptype;
        for (int j = 0; j < curve->m_points.size() - 2; ++j)
        {
            edg->m_edgeNodes.push_back(
                std::make_shared<Node>(
                    j,
                    (*curve->m_points[j + 1])(0),
                    (*curve->m_points[j + 1])(1),
                    (*curve->m_points[j + 1])(2)));
        }
    }

    // Curved Faces
    for (auto &it : graph->GetCurvedFaces())
    {
        SpatialDomains::CurveSharedPtr curve = it.second;
        int id = curve->m_curveID;
        ASSERTL1(fIdMap.find(id) != fIdMap.end(), "Failed to find curved edge");
        FaceSharedPtr fac = fIdMap[id];
        fac->m_curveType  = curve->m_ptype;
        int Ntot          = curve->m_points.size();

        if (fac->m_curveType == LibUtilities::eNodalTriFekete ||
            fac->m_curveType == LibUtilities::eNodalTriEvenlySpaced ||
            fac->m_curveType == LibUtilities::eNodalTriElec)
        {
            int N = ((int)sqrt(8.0 * Ntot + 1.0) - 1) / 2;
            for (int j = 3 + 3 * (N - 2); j < Ntot; ++j)
            {
                fac->m_faceNodes.push_back(
                    std::make_shared<Node>(
                        j,
                        (*curve->m_points[j])(0),
                        (*curve->m_points[j])(1),
                        (*curve->m_points[j])(2)));
            }
        }
        else // quad face.
        {
            int N = (int)sqrt((double)Ntot);

            for (int j = 1; j < N - 1; ++j)
            {
                for (int k = 1; k < N - 1; ++k)
                {
                    fac->m_faceNodes.push_back(
                        std::make_shared<Node>(
                            (j - 1) * (N - 2) + k - 1,
                            (*curve->m_points[j * N + k])(0),
                            (*curve->m_points[j * N + k])(1),
                            (*curve->m_points[j * N + k])(2)));
                }
            }
        }
    }

    // Get hold of mesh composites and set up m_mesh->m_elements. Loop over all
    // composites and set up elements with edges and faces from the maps above.
    for (auto &compIt : graph->GetComposites())
    {
        // Get hold of dimension
        int dim = compIt.second->m_geomVec[0]->GetShapeDim();

        // compIt->second is a GeometryVector
        for (auto &geomIt : compIt.second->m_geomVec)
        {
            ElmtConfig conf(geomIt->GetShapeType(), 1, true, true, false);

            // Get hold of geometry
            vector<NodeSharedPtr> nodeList;
            for (int i = 0; i < geomIt->GetNumVerts(); ++i)
            {
                nodeList.push_back(vIdMap[geomIt->GetVid(i)]);
            }

            vector<int> tags;
            tags.push_back(compIt.first);

            ElementSharedPtr E = GetElementFactory().CreateInstance(
                geomIt->GetShapeType(), conf, nodeList, tags);

            E->SetId(geomIt->GetGlobalID());
            m_mesh->m_element[dim].push_back(E);

            if (dim == 1)
            {
                EdgeSharedPtr edg = eIdMap[geomIt->GetGlobalID()];
                E->SetVolumeNodes(edg->m_edgeNodes);
                E->SetCurveType(edg->m_curveType);
            }

            if (dim > 1)
            {
                // reset edges
                for (int i = 0; i < geomIt->GetNumEdges(); ++i)
                {
                    EdgeSharedPtr edg = eIdMap[geomIt->GetEid(i)];
                    E->SetEdge(i, edg);
                    // set up link back to this element
                    edg->m_elLink.push_back(pair<ElementSharedPtr, int>(E, i));
                }

                if (dim == 2)
                {
                    FaceSharedPtr fac = fIdMap[geomIt->GetGlobalID()];
                    E->SetVolumeNodes(fac->m_faceNodes);
                    E->SetCurveType(fac->m_curveType);
                }
            }

            if (dim == 3)
            {
                // reset faces
                for (int i = 0; i < geomIt->GetNumFaces(); ++i)
                {
                    FaceSharedPtr fac = fIdMap[geomIt->GetFid(i)];
                    E->SetFace(i, fac);
                    // set up link back to this slement
                    fac->m_elLink.push_back(pair<ElementSharedPtr, int>(E, i));
                }
            }
        }
    }

    // set up composite labels if they exist
    m_mesh->m_faceLabels = graph->GetCompositesLabels();

    ProcessEdges(false);
    ProcessFaces(false);
    ProcessComposites();
}
}
}
