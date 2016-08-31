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
//  Description: GMSH converter.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
using namespace std;

#include <tinyxml.h>
#include <boost/algorithm/string.hpp>

#ifdef NEKTAR_USE_MESHGEN
#include <NekMeshUtils/CADSystem/CADSystem.h>
#endif

#include <SpatialDomains/MeshGraph.h>
#include <NekMeshUtils/MeshElements/Element.h>
#include "InputNekpp.h"

namespace Nektar
{
namespace Utilities
{
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

struct cadVar
{
    string type;
    int id;
    NekDouble u,v;
};

/**
 *
 */
void InputNekpp::Process()
{
    vector<string> filename;
    filename.push_back(m_config["infile"].as<string>());

    LibUtilities::SessionReaderSharedPtr pSession =
        LibUtilities::SessionReader::CreateInstance(0, NULL, filename);
    SpatialDomains::MeshGraphSharedPtr graph =
        SpatialDomains::MeshGraph::Read(pSession);

    m_mesh->m_expDim   = graph->GetMeshDimension();
    m_mesh->m_spaceDim = graph->GetSpaceDimension();

    // Copy vertices.
    map<int, NodeSharedPtr> vIdMap;
    const SpatialDomains::PointGeomMap vertset = graph->GetVertSet();
    SpatialDomains::PointGeomMap::const_iterator vit;

    for (vit = vertset.begin(); vit != vertset.end(); ++vit)
    {
        SpatialDomains::PointGeomSharedPtr vert = vit->second;
        NodeSharedPtr n(
            new Node(vert->GetVid(), (*vert)(0), (*vert)(1), (*vert)(2)));
        m_mesh->m_vertexSet.insert(n);
        vIdMap[vert->GetVid()] = n;
    }

    map<int, EdgeSharedPtr> eIdMap;
    map<int, EdgeSharedPtr>::iterator itEmap;
    map<int, FaceSharedPtr> fIdMap;
    map<int, FaceSharedPtr>::iterator itFmap;

    // Load up all edges from graph
    {
        int nel                        = 0;
        SpatialDomains::SegGeomMap tmp = graph->GetAllSegGeoms();
        SpatialDomains::SegGeomMap::iterator it;
        pair<EdgeSet::iterator, bool> testIns;
        nel += tmp.size();

        for (it = tmp.begin(); it != tmp.end(); ++it)
        {
            pair<int, SpatialDomains::GeometrySharedPtr> tmp2(
                it->first,
                boost::dynamic_pointer_cast<SpatialDomains::Geometry>(
                    it->second));

            // load up edge set in order of SegGeomMap;
            vector<NodeSharedPtr> curve; // curved nodes if deformed
            int id0 = it->second->GetVid(0);
            int id1 = it->second->GetVid(1);
            LibUtilities::PointsType ptype =
                it->second->GetPointsKeys()[0].GetPointsType();
            EdgeSharedPtr ed =
                EdgeSharedPtr(new Edge(vIdMap[id0], vIdMap[id1], curve, ptype));

            testIns = m_mesh->m_edgeSet.insert(ed);
            (*(testIns.first))->m_id = it->second->GetEid();
            eIdMap[it->second->GetEid()] = ed;
        }
    }

    // load up all faces from graph
    {
        int nel                        = 0;
        SpatialDomains::TriGeomMap tmp = graph->GetAllTriGeoms();
        SpatialDomains::TriGeomMap::iterator it;
        pair<FaceSet::iterator, bool> testIns;
        nel += tmp.size();

        for (it = tmp.begin(); it != tmp.end(); ++it)
        {
            pair<int, SpatialDomains::GeometrySharedPtr> tmp2(
                it->first,
                boost::dynamic_pointer_cast<SpatialDomains::Geometry>(
                    it->second));

            vector<NodeSharedPtr> faceVertices;
            vector<EdgeSharedPtr> faceEdges;
            vector<NodeSharedPtr> faceNodes;

            for (int i = 0; i < 3; ++i)
            {
                faceVertices.push_back(vIdMap[it->second->GetVid(i)]);
                faceEdges.push_back(eIdMap[it->second->GetEid(i)]);
            }

            FaceSharedPtr fac =
                FaceSharedPtr(new Face(faceVertices,
                                       faceNodes,
                                       faceEdges,
                                       LibUtilities::ePolyEvenlySpaced));
            testIns = m_mesh->m_faceSet.insert(fac);
            (*(testIns.first))->m_id = it->second->GetFid();
            fIdMap[it->second->GetFid()] = fac;
        }

        SpatialDomains::QuadGeomMap tmp3 = graph->GetAllQuadGeoms();
        SpatialDomains::QuadGeomMap::iterator it2;

        for (it2 = tmp3.begin(); it2 != tmp3.end(); ++it2)
        {
            pair<int, SpatialDomains::GeometrySharedPtr> tmp2(
                it2->first,
                boost::dynamic_pointer_cast<SpatialDomains::Geometry>(
                    it2->second));

            vector<NodeSharedPtr> faceVertices;
            vector<EdgeSharedPtr> faceEdges;
            vector<NodeSharedPtr> faceNodes;

            for (int i = 0; i < 4; ++i)
            {
                faceVertices.push_back(vIdMap[it2->second->GetVid(i)]);
                faceEdges.push_back(eIdMap[it2->second->GetEid(i)]);
            }

            FaceSharedPtr fac =
                FaceSharedPtr(new Face(faceVertices,
                                       faceNodes,
                                       faceEdges,
                                       LibUtilities::ePolyEvenlySpaced));
            testIns = m_mesh->m_faceSet.insert(fac);
            (*(testIns.first))->m_id = it2->second->GetFid();
            fIdMap[it2->second->GetFid()] = fac;
        }
    }

    // Set up curved information

    // Curved Edges
    SpatialDomains::CurveMap &curvedEdges = graph->GetCurvedEdges();
    SpatialDomains::CurveMap::iterator it;

    for (it = curvedEdges.begin(); it != curvedEdges.end(); ++it)
    {
        SpatialDomains::CurveSharedPtr curve = it->second;
        int id = curve->m_curveID;
        ASSERTL1(eIdMap.find(id) != eIdMap.end(), "Failed to find curved edge");
        EdgeSharedPtr edg = eIdMap[id];
        edg->m_curveType = curve->m_ptype;
        for (int j = 0; j < curve->m_points.size() - 2; ++j)
        {
            NodeSharedPtr n(new Node(j,
                                     (*curve->m_points[j + 1])(0),
                                     (*curve->m_points[j + 1])(1),
                                     (*curve->m_points[j + 1])(2)));
            edg->m_edgeNodes.push_back(n);
        }
    }

    // Curved Faces
    SpatialDomains::CurveMap &curvedFaces = graph->GetCurvedFaces();
    for (it = curvedFaces.begin(); it != curvedFaces.end(); ++it)
    {
        SpatialDomains::CurveSharedPtr curve = it->second;
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
                NodeSharedPtr n(new Node(j,
                                         (*curve->m_points[j])(0),
                                         (*curve->m_points[j])(1),
                                         (*curve->m_points[j])(2)));
                fac->m_faceNodes.push_back(n);
            }
        }
        else // quad face.
        {
            int N = (int)sqrt((double)Ntot);

            for (int j = 1; j < N - 1; ++j)
            {
                for (int k = 1; k < N - 1; ++k)
                {
                    NodeSharedPtr n(new Node((j - 1) * (N - 2) + k - 1,
                                             (*curve->m_points[j * N + k])(0),
                                             (*curve->m_points[j * N + k])(1),
                                             (*curve->m_points[j * N + k])(2)));
                    fac->m_faceNodes.push_back(n);
                }
            }
        }
    }

    // Get hold of mesh composites and set up m_mesh->m_elements
    SpatialDomains::CompositeMap GraphComps = graph->GetComposites();
    SpatialDomains::CompositeMapIter compIt;
    SpatialDomains::GeometryVectorIter geomIt;

    // loop over all composites and set up elements with edges
    // and faces from the maps above.
    for (compIt = GraphComps.begin(); compIt != GraphComps.end(); ++compIt)
    {
        // Get hold of dimension
        int dim = (*compIt->second)[0]->GetShapeDim();

        // compIt->second is a GeometryVector
        for (geomIt = (*compIt->second).begin();
             geomIt != (*compIt->second).end();
             ++geomIt)
        {
            ElmtConfig conf((*geomIt)->GetShapeType(), 1, true, true, false);

            // Get hold of geometry
            vector<NodeSharedPtr> nodeList;
            for (int i = 0; i < (*geomIt)->GetNumVerts(); ++i)
            {
                nodeList.push_back(vIdMap[(*geomIt)->GetVid(i)]);
            }

            vector<int> tags;
            tags.push_back(compIt->first);

            ElementSharedPtr E = GetElementFactory().CreateInstance(
                (*geomIt)->GetShapeType(), conf, nodeList, tags);

            E->SetId((*geomIt)->GetGlobalID());
            m_mesh->m_element[dim].push_back(E);

            if (dim == 1)
            {
                EdgeSharedPtr edg = eIdMap[(*geomIt)->GetGlobalID()];
                E->SetVolumeNodes(edg->m_edgeNodes);
            }

            if (dim > 1)
            {
                // reset edges
                for (int i = 0; i < (*geomIt)->GetNumEdges(); ++i)
                {
                    EdgeSharedPtr edg = eIdMap[(*geomIt)->GetEid(i)];
                    E->SetEdge(i, edg);
                    // set up link back to this element
                    edg->m_elLink.push_back(pair<ElementSharedPtr, int>(E, i));
                }

                if (dim == 2)
                {
                    FaceSharedPtr fac = fIdMap[(*geomIt)->GetGlobalID()];
                    E->SetVolumeNodes(fac->m_faceNodes);
                }
            }

            if (dim == 3)
            {
                // reset faces
                for (int i = 0; i < (*geomIt)->GetNumFaces(); ++i)
                {
                    FaceSharedPtr fac = fIdMap[(*geomIt)->GetFid(i)];
                    E->SetFace(i, fac);
                    // set up link back to this slement
                    fac->m_elLink.push_back(pair<ElementSharedPtr, int>(E, i));
                }
            }
        }
    }

    // set up composite labels if they exist
    m_mesh->m_faceLabels = graph->GetCompositesLabels();

    m_mesh->m_hasCAD = false;

#ifdef NEKTAR_USE_MESHGEN

    if(pSession->DefinesElement("NEKTAR/GEOMETRY/CADID"))
    {
        m_mesh->m_hasCAD = true;
        TiXmlElement* id = pSession->GetElement("NEKTAR/GEOMETRY/CADID");

        id->QueryStringAttribute("NAME",&m_mesh->m_CADId);

        m_mesh->m_cad =
            MemoryManager<CADSystem>::AllocateSharedPtr(m_mesh->m_CADId);
        ASSERTL0(m_mesh->m_cad->LoadCAD(), "Failed to load CAD");
    }

    if(m_mesh->m_hasCAD && !pSession->DefinesElement("NEKTAR/GEOMETRY/CAD"))
    {
        ASSERTL0(false, "CAD file but no data");
    }

    map<int, vector<cadVar> > vertToString;
    map<pair<int,int>, pair<int,vector<cadVar> > > edgeToString;
    map<pair<int,int>, vector<cadVar> > faceToString;

    if(pSession->DefinesElement("NEKTAR/GEOMETRY/CAD"))
    {
        int err;

        TiXmlElement* cad = pSession->GetElement("NEKTAR/GEOMETRY/CAD");

        TiXmlElement *vertex = cad->FirstChildElement("V");

        int vid;

        while(vertex)
        {
            std::string vert(vertex->ValueStr());
            ASSERTL0(vert == "V", (std::string("Unknown CAD type:") + vert).c_str());

            /// Read edge id attribute.
            err = vertex->QueryIntAttribute("VERTID", &vid);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read curve attribute VERTID.");

            vector<cadVar> items;

            TiXmlElement *curve = vertex->FirstChildElement("C");

            while(curve)
            {
                cadVar cv;
                cv.type = "C";
                curve->QueryIntAttribute("CADID", &cv.id);
                TiXmlNode* elementChild = curve->FirstChild();
                string str = elementChild->ToText()->ValueStr();
                istringstream ss(str);
                ss >> cv.u;
                items.push_back(cv);
                curve = curve->NextSiblingElement("C");
            }

            TiXmlElement *surf = vertex->FirstChildElement("S");

            while(surf)
            {
                cadVar cv;
                cv.type = "S";
                surf->QueryIntAttribute("CADID", &cv.id);
                TiXmlNode* elementChild = surf->FirstChild();
                string str = elementChild->ToText()->ValueStr();
                istringstream ss(str);
                ss >> cv.u >> cv.v;
                items.push_back(cv);
                surf = surf->NextSiblingElement("S");
            }

            vertToString[vid] = items;

            vertex = vertex->NextSiblingElement("V");
        }

        TiXmlElement *edge = cad->FirstChildElement("E");

        int eid;
        int node;

        while(edge)
        {
            std::string e(edge->ValueStr());
            ASSERTL0(e == "E", (std::string("Unknown CAD type:") + e).c_str());

            /// Read edge id attribute.
            err = edge->QueryIntAttribute("EDGEID", &eid);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read curve attribute EDGEID.");
            int cu;
            /// Read edge node attribute.
            err = edge->QueryIntAttribute("ONCURVE", &cu);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read curve attribute EDGEID.");

            TiXmlElement *node = edge->FirstChildElement("N");

            while(node)
            {
                int nid;
                node->QueryIntAttribute("NODEID",&nid);
                vector<cadVar> items;

                TiXmlElement *curve = node->FirstChildElement("C");

                while(curve)
                {
                    cadVar cv;
                    cv.type = "C";
                    curve->QueryIntAttribute("CADID", &cv.id);
                    TiXmlNode* elementChild = curve->FirstChild();
                    string str = elementChild->ToText()->ValueStr();
                    istringstream ss(str);
                    ss >> cv.u;
                    items.push_back(cv);
                    curve = curve->NextSiblingElement("C");
                }

                TiXmlElement *surf = node->FirstChildElement("S");

                while(surf)
                {
                    cadVar cv;
                    cv.type = "S";
                    surf->QueryIntAttribute("CADID", &cv.id);
                    TiXmlNode* elementChild = surf->FirstChild();
                    string str = elementChild->ToText()->ValueStr();
                    istringstream ss(str);
                    ss >> cv.u >> cv.v;
                    items.push_back(cv);
                    surf = surf->NextSiblingElement("S");
                }

                node = node->NextSiblingElement("N");
                edgeToString[pair<int,int>(eid,nid)] = pair<int,vector<cadVar>(cu,items);
            }

            edge = edge->NextSiblingElement("E");
        }

        TiXmlElement *face = cad->FirstChildElement("F");

        int fid;

        while(face)
        {
            std::string f(face->ValueStr());
            ASSERTL0(f == "F", (std::string("Unknown CAD type:") + f).c_str());

            /// Read edge id attribute.
            err = face->QueryIntAttribute("FACEID", &fid);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read curve attribute FACEID.");
            int su;
            err = edge->QueryIntAttribute("ONSURF", &su);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read curve attribute EDGEID.");

            TiXmlElement *node = face->FirstChildElement("N");

            while(node)
            {
                int nid;
                node->QueryIntAttribute("NODEID",&nid);
                vector<cadVar> items;

                TiXmlElement *surf = node->FirstChildElement("S");

                while(surf)
                {
                    cadVar cv;
                    cv.type = "S";
                    surf->QueryIntAttribute("CADID", &cv.id);
                    TiXmlNode* elementChild = surf->FirstChild();
                    string str = elementChild->ToText()->ValueStr();
                    istringstream ss(str);
                    ss >> cv.u >> cv.v;
                    items.push_back(cv);
                    surf = surf->NextSiblingElement("S");
                }

                node = node->NextSiblingElement("N");
                edgeToString[pair<int,int>(fid,nid)] = pair<int,vector<cadVar>(su,items);
            }

            face = face->NextSiblingElement("F");
        }
    }

#endif

    ProcessEdges(false);
    ProcessFaces(false);
    ProcessComposites();

#ifdef NEKTAR_USE_MESHGEN

    map<int, vector<string> >::iterator vsit;
    map<pair<int,int>, vector<string> >::iterator esit;
    map<pair<int,int>, vector<string> >::iterator fsit;

    {
        int ct= 0;
        NodeSet::iterator it;
        for(it = m_mesh->m_vertexSet.begin(); it != m_mesh->m_vertexSet.end();
            it++)
        {
            vsit = vertToString.find((*it)->m_id);
            if(vsit != vertToString.end())
            {
                ct++;
                for(int i = 0; i < vsit->second.size(); i++)
                {
                    istringstream iss(vsit->second[i]);
                    string t;
                    int id;
                    NekDouble u, v;
                    iss >> t >> id >> u >> v;
                    if(t == "C")
                    {
                        (*it)->SetCADCurve(id,m_mesh->m_cad->GetCurve(id),u);
                    }
                    else if(t == "S")
                    {
                        Array<OneD,NekDouble> uv(2);
                        uv[0] = u;
                        uv[1] = v;
                        (*it)->SetCADSurf(id,m_mesh->m_cad->GetSurf(id),uv);
                    }
                    else
                    {
                        ASSERTL0(false,"unsure on type");
                    }
                }
            }
        }
        ASSERTL0(ct == vertToString.size(), "did not find all CAD information");
    }

    {
        int ct = 0;
        EdgeSet::iterator it;
        for(it = m_mesh->m_edgeSet.begin(); it != m_mesh->m_edgeSet.end();
            it++)
        {
            for(int j = 0; j < (*it)->m_edgeNodes.size(); j++)
            {
                esit = edgeToString.find(pair<int,int>((*it)->m_id,j));
                if(esit != edgeToString.end())
                {
                    ct++;
                    for(int i = 0; i < esit->second.size(); i++)
                    {
                        istringstream iss(esit->second[i]);
                        string t;
                        int id;
                        NekDouble u, v;
                        iss >> t >> id >> u >> v;
                        if(t == "C")
                        {
                            (*it)->m_edgeNodes[j]->SetCADCurve(id,m_mesh->m_cad->GetCurve(id),u);
                        }
                        else if(t == "S")
                        {
                            Array<OneD,NekDouble> uv(2);
                            uv[0] = u;
                            uv[1] = v;
                            (*it)->m_edgeNodes[j]->SetCADSurf(id,m_mesh->m_cad->GetSurf(id),uv);
                        }
                        else
                        {
                            ASSERTL0(false,"unsure on type");
                        }
                    }
                }
            }
        }
        ASSERTL0(ct == edgeToString.size(), "did not find all CAD information");
    }

    {
        int ct = 0;
        FaceSet::iterator it;
        for(it = m_mesh->m_faceSet.begin(); it != m_mesh->m_faceSet.end();
            it++)
        {
            for(int j = 0; j < (*it)->m_faceNodes.size(); j++)
            {
                fsit = faceToString.find(pair<int,int>((*it)->m_id,j));
                if(fsit != faceToString.end())
                {
                    ct++;
                    for(int i = 0; i < fsit->second.size(); i++)
                    {
                        istringstream iss(fsit->second[i]);
                        string t;
                        int id;
                        NekDouble u, v;
                        iss >> t >> id >> u >> v;
                        if(t == "C")
                        {
                            (*it)->m_faceNodes[j]->SetCADCurve(id,m_mesh->m_cad->GetCurve(id),u);
                        }
                        else if(t == "S")
                        {
                            Array<OneD,NekDouble> uv(2);
                            uv[0] = u;
                            uv[1] = v;
                            (*it)->m_faceNodes[j]->SetCADSurf(id,m_mesh->m_cad->GetSurf(id),uv);
                        }
                        else
                        {
                            ASSERTL0(false,"unsure on type");
                        }
                    }
                }
            }
        }
        ASSERTL0(ct == faceToString.size(), "did not find all CAD information");
    }
#endif

}
}
}
