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

#include <SpatialDomains/MeshGraph.h>

#include "MeshElements.h"
#include "InputNekpp.h"

namespace Nektar
{
    namespace Utilities
    {
        ModuleKey InputNekpp::className = 
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eInputModule, "xml"), InputNekpp::create,
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
         * Gmsh file contains a list of nodes and their coordinates, along with
         * a list of elements and those nodes which define them. We read in and
         * store the list of nodes in #m_node and store the list of elements in
         * #m_element. Each new element is supplied with a list of entries from
         * #m_node which defines the element. Finally some mesh statistics are
         * printed.
         *
         * @param   pFilename           Filename of Gmsh file to read.
         */
        void InputNekpp::Process()
        {
            vector<string> filename;
            filename.push_back(m_config["infile"].as<string>());
            
            LibUtilities::SessionReaderSharedPtr pSession =
                LibUtilities::SessionReader::CreateInstance(0, NULL, filename);
            SpatialDomains::MeshGraphSharedPtr graph =
                SpatialDomains::MeshGraph::Read(pSession);

            m_mesh->m_expDim   = graph->GetMeshDimension ();
            m_mesh->m_spaceDim = graph->GetSpaceDimension();

            // Copy vertices.
            map<int, NodeSharedPtr> vIdMap;
            int nVerts = graph->GetNvertices();
            for (int i = 0; i < nVerts; ++i)
            {
                SpatialDomains::PointGeomSharedPtr vert =
                    graph->GetVertex(i);
                NodeSharedPtr n(new Node(vert->GetVid(),
                    (*vert)(0), (*vert)(1), (*vert)(2)));
                m_mesh->m_vertexSet.insert(n);
                vIdMap[vert->GetVid()] = n;
            }

            map<int, EdgeSharedPtr> eIdMap;
            map<int, EdgeSharedPtr>::iterator itEmap;
            map<int, FaceSharedPtr> fIdMap;
            map<int, FaceSharedPtr>::iterator itFmap;

            // Load up all edges from graph
            {
                int nel = 0;
                SpatialDomains::SegGeomMap tmp = graph->GetAllSegGeoms();
                SpatialDomains::SegGeomMap::iterator it;
                pair<EdgeSet::iterator,bool> testIns;
                nel += tmp.size();

                for (it = tmp.begin(); it != tmp.end(); ++it)
                {
                    pair<int, SpatialDomains::GeometrySharedPtr> tmp2(
                            it->first, boost::dynamic_pointer_cast<
                            SpatialDomains::Geometry>(it->second));
                    
                    // load up edge set in order of SegGeomMap; 
                    vector<NodeSharedPtr> curve; // curved nodes if deformed
                    int id0 = it->second->GetVid(0);
                    int id1 = it->second->GetVid(1);
                    LibUtilities::PointsType ptype = it->second->GetPointsKeys()[0].GetPointsType();
                    EdgeSharedPtr ed = EdgeSharedPtr(new Edge(vIdMap[id0], vIdMap[id1],
                                                              curve, ptype));
                    
                    testIns = m_mesh->m_edgeSet.insert(ed);
                    (*(testIns.first))->m_id = it->second->GetEid();
                    eIdMap[it->second->GetEid()] = ed;
                }
            }


            // load up all faces from graph
            {
                int nel = 0;
                SpatialDomains::TriGeomMap tmp = graph->GetAllTriGeoms();
                SpatialDomains::TriGeomMap::iterator it;
                pair<FaceSet::iterator,bool> testIns;
                nel += tmp.size();

                for (it = tmp.begin(); it != tmp.end(); ++it)
                {
                    pair<int, SpatialDomains::GeometrySharedPtr> tmp2(
                            it->first, boost::dynamic_pointer_cast<
                            SpatialDomains::Geometry>(it->second));
                            
                    vector<NodeSharedPtr> faceVertices;
                    vector<EdgeSharedPtr> faceEdges;
                    vector<NodeSharedPtr> faceNodes;
                    
                    for(int i = 0; i < 3; ++i)
                    {
                        faceVertices.push_back(vIdMap[it->second->GetVid(i)]);
                        faceEdges.push_back(eIdMap[it->second->GetEid(i)]);
                    }
                    
                    FaceSharedPtr fac = FaceSharedPtr( new Face(faceVertices,faceNodes,faceEdges,
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
                            it2->first, boost::dynamic_pointer_cast<
                            SpatialDomains::Geometry>(it2->second));
                            
                    vector<NodeSharedPtr> faceVertices;
                    vector<EdgeSharedPtr> faceEdges;
                    vector<NodeSharedPtr> faceNodes;
                    
                    for(int i = 0; i < 4; ++i)
                    {
                        faceVertices.push_back(vIdMap[it2->second->GetVid(i)]);
                        faceEdges.push_back(eIdMap[it2->second->GetEid(i)]);
                    }
                    
                    FaceSharedPtr fac = FaceSharedPtr( new Face(faceVertices,faceNodes,faceEdges,
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
                ASSERTL1(eIdMap.find(id) != eIdMap.end(),
                         "Failed to find curved edge");
                EdgeSharedPtr edg = eIdMap[id];
                edg->m_curveType = curve->m_ptype;
                for(int j = 0; j < curve->m_points.size()-2; ++j)
                {
                    NodeSharedPtr n(new Node(j, (*curve->m_points[j+1])(0),
                                             (*curve->m_points[j+1])(1),
                                             (*curve->m_points[j+1])(2)));
                    edg->m_edgeNodes.push_back(n);
                }
            }

            // Curved Faces
            SpatialDomains::CurveMap &curvedFaces = graph->GetCurvedFaces();
            for (it = curvedFaces.begin(); it != curvedFaces.end(); ++it)
            {
                SpatialDomains::CurveSharedPtr curve = it->second;
                int id = curve->m_curveID;
                ASSERTL1(fIdMap.find(id) != fIdMap.end(),
                         "Failed to find curved edge");
                FaceSharedPtr fac = fIdMap[id];
                fac->m_curveType = curve->m_ptype;
                int Ntot = curve->m_points.size();

                if (fac->m_curveType == LibUtilities::eNodalTriFekete       ||
                    fac->m_curveType == LibUtilities::eNodalTriEvenlySpaced ||
                    fac->m_curveType == LibUtilities::eNodalTriElec)
                {
                    int N    = ((int)sqrt(8.0*Ntot+1.0)-1)/2;
                    for(int j = 3+3*(N-2); j < Ntot; ++j)
                    {
                        NodeSharedPtr n(new Node(j, (*curve->m_points[j])(0),
                                                 (*curve->m_points[j])(1),
                                                 (*curve->m_points[j])(2)));
                        fac->m_faceNodes.push_back(n);
                    }
                }
                else // quad face.
                {
                    int N    = (int)sqrt((double)Ntot);
                    for(int j = 1; j < N-1; ++j)
                    {
                        for(int k = 1; k < N-1; ++k)
                        {
                            NodeSharedPtr n(new Node((j-1)*(N-2)+k-1,
                                                     (*curve->m_points[j*N+k])(0),
                                                     (*curve->m_points[j*N+k])(1),
                                                     (*curve->m_points[j*N+k])(2)));
                            fac->m_faceNodes.push_back(n);
                        }
                    }
                }
            }

            // Get hold of mesh composites and set up m_mesh->m_elements

            SpatialDomains::CompositeMap       GraphComps= graph->GetComposites();
            SpatialDomains::CompositeMapIter   compIt;
            SpatialDomains::GeometryVectorIter geomIt;


            // calculate the number of element of dimension
            // m_mesh->m_expDim in composite list so we can set up
            // element vector of this size to allow for
            // non-consecutive insertion to list (Might consider
            // setting element up as a map)?
            int nel = 0; 
            for(compIt = GraphComps.begin(); compIt != GraphComps.end(); ++compIt)
            {
                // Get hold of dimension
                int dim = (*compIt->second)[0]->GetShapeDim();
                
                if(dim == m_mesh->m_expDim) 
                {
                    nel += (*compIt->second).size();
                }
            }
            m_mesh->m_element[m_mesh->m_expDim].resize(nel);

            // loop over all composites and set up elements with edges and faces from the maps above. 
            for(compIt = GraphComps.begin(); compIt != GraphComps.end(); ++compIt)
            {
                // Get hold of dimension
                int dim = (*compIt->second)[0]->GetShapeDim();

                // compIt->second is a GeometryVector
                for(geomIt  = (*compIt->second).begin(); 
                    geomIt != (*compIt->second).end();
                        ++geomIt)
                {
                    ElmtConfig conf((*geomIt)->GetShapeType(),1,true,true);
                    
                    // Get hold of geometry
                    vector<NodeSharedPtr> nodeList;
                    for (int i = 0; i < (*geomIt)->GetNumVerts(); ++i)
                    {
                        nodeList.push_back(vIdMap[(*geomIt)->GetVid(i)]);
                    }
                
                    vector<int> tags;
                    tags.push_back(compIt->first);
                    
                    ElementSharedPtr E = GetElementFactory().
                        CreateInstance((*geomIt)->GetShapeType(),conf,nodeList,tags);
                    
                    E->SetId((*geomIt)->GetGlobalID());
                    
                    if(dim == m_mesh->m_expDim) // load mesh into location baded on globalID
                    {
                        m_mesh->m_element[dim][(*geomIt)->GetGlobalID()] = E;
                    }
                    else // push onto vector for later usage as composite region
                    {
                        m_mesh->m_element[dim].push_back(E);
                    }
                    
                    if(dim > 1)
                    {
                        // reset edges 
                        for (int i = 0; i < (*geomIt)->GetNumEdges(); ++i)
                        {
                            EdgeSharedPtr edg = eIdMap[(*geomIt)->GetEid(i)];
                            E->SetEdge(i,edg);
                            // set up link back to this element
                            edg->m_elLink.push_back(pair<ElementSharedPtr,int>(E,i));
                        }
                    }
                    
                    if(dim  == 3)
                    {
                        // reset faces 
                        for (int i = 0; i < (*geomIt)->GetNumFaces(); ++i)
                        {
                            FaceSharedPtr fac = fIdMap[(*geomIt)->GetFid(i)];
                            E->SetFace(i,fac);
                            // set up link back to this slement
                            fac->m_elLink.push_back(pair<ElementSharedPtr,int>(E,i));
                        }
                    }
                }                             
            }
            ProcessEdges(false); 
            ProcessFaces(false);
            ProcessComposites();
        }
    }
}
