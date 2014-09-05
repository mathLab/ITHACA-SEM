////////////////////////////////////////////////////////////////////////////////
//
//  File: InputNek.cpp
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
//  Description: Nektar file format converter.
//
////////////////////////////////////////////////////////////////////////////////


#include "MeshElements.h"
#include "InputNek.h"

#include <LibUtilities/Foundations/ManagerAccess.h>
#include <boost/algorithm/string.hpp>

#include <map>
#include <vector>
#include <sstream>
#include <string>
using namespace std;

namespace Nektar
{
    namespace Utilities
    {
        ModuleKey InputNek::className = 
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eInputModule, "rea"), InputNek::create,
                "Reads Nektar rea file.");

        InputNek::InputNek(MeshSharedPtr m) : InputModule(m)
        {
            
        }

        InputNek::~InputNek()
        {
            
        }

        /**
         * @brief Processes Nektar file format.
         * 
         * Nektar sessions are defined by rea files, and contain sections
         * defining a DNS simulation in a specific order. The converter only
         * reads mesh information, curve information if it exists and boundary
         * information.
         * 
         * @param pFilename Filename of Nektar session file to read.
         */
        void InputNek::Process()
        {
            // Open the file stream.
            OpenStream();
            
            string      line, word;
            int         nParam, nElements, nCurves;
            int         i, j, k, nodeCounter = 0;
            int         nComposite = 0;
            LibUtilities::ShapeType elType;
            double      vertex[3][6];
            map<LibUtilities::ShapeType,int> domainComposite;
            map<LibUtilities::ShapeType,vector< vector<NodeSharedPtr> > > elNodes;
            map<LibUtilities::ShapeType,vector<int> > elIds;
            boost::unordered_map<int,int> elMap;
            vector<LibUtilities::ShapeType> elmOrder;

            // Set up vector of processing orders.
            elmOrder.push_back(LibUtilities::eSegment);
            elmOrder.push_back(LibUtilities::eTriangle);
            elmOrder.push_back(LibUtilities::eQuadrilateral);
            elmOrder.push_back(LibUtilities::ePrism);
            elmOrder.push_back(LibUtilities::ePyramid);
            elmOrder.push_back(LibUtilities::eTetrahedron);
            elmOrder.push_back(LibUtilities::eHexahedron);
             
            m_mesh->m_expDim   = 0;
            m_mesh->m_spaceDim = 0;
            
            if (m_mesh->m_verbose)
            {
                cout << "InputNek: Start reading file..." << endl;
            }
            
            // -- Read in parameters.

            // Ignore first 3 lines. 4th line contains number of parameters.
            for (i = 0; i < 4; ++i)
            {
                getline(m_mshFile, line);
            }

            stringstream s(line);
            s >> nParam;
            
            for (i = 0; i < nParam; ++i)
            {
                string tmp1, tmp2;
                getline(m_mshFile, line);
                s.str(line);
                s >> tmp1 >> tmp2;
            }
            
            // -- Read in passive scalars (ignore)
            getline(m_mshFile, line);
            s.clear();
            s.str(line);
            s >> j;
            for (i = 0; i < j; ++i)
            {
                getline(m_mshFile, line);
            }

            // -- Read in logical switches (ignore)
            getline(m_mshFile, line);
            s.clear();
            s.str(line);
            s >> j;
            for (i = 0; i < j; ++i)
            {
                getline(m_mshFile, line);
            }

            // -- Read in mesh data.
            
            // First hunt for MESH tag
            bool foundMesh = false;
            while (!m_mshFile.eof())
            {
                getline(m_mshFile, line);
                if (line.find("MESH") != string::npos)
                {
                    foundMesh = true;
                    break;
                }
            }
            
            if (!foundMesh)
            {
                cerr << "Couldn't find MESH tag inside file." << endl;
                abort();
            }
            
            // Now read in number of elements and space dimension.
            getline(m_mshFile, line);
            s.clear(); s.str(line);
            s >> nElements >> m_mesh->m_expDim;
            m_mesh->m_spaceDim = m_mesh->m_expDim;
            
            // Set up field names.
            m_mesh->m_fields.push_back("u");
            m_mesh->m_fields.push_back("v");
            if (m_mesh->m_spaceDim > 2)
            {
                m_mesh->m_fields.push_back("w");
            }
            m_mesh->m_fields.push_back("p");
            
            // Loop over and create elements.
            for (i = 0; i < nElements; ++i)
            {
                getline(m_mshFile, line);
                
                if (m_mesh->m_expDim == 2)
                {
                    if (line.find("Qua") != string::npos || 
                        line.find("qua") != string::npos)
                    {
                        elType = LibUtilities::eQuadrilateral;
                    }
                    else
                    {
                        // Default element type in 2D is triangle.
                        elType = LibUtilities::eTriangle;
                    }
                } 
                else
                {
                    if (line.find("Tet") != string::npos || 
                        line.find("tet") != string::npos) 
                    {
                        elType = LibUtilities::eTetrahedron;
                    }
                    else if (line.find("Hex") != string::npos || 
                             line.find("hex") != string::npos) 
                    {
                        elType = LibUtilities::eHexahedron;
                    }
                    else if (line.find("Prism") != string::npos || 
                             line.find("prism") != string::npos) 
                    {
                        elType = LibUtilities::ePrism;
                    }
                    else if (line.find("Pyr") != string::npos || 
                             line.find("pyr") != string::npos) 
                    {
                        cerr << "Pyramid elements not yet supported." << endl;
                        abort();
                    }
                    else if (line.find("Qua") != string::npos || 
                             line.find("qua") != string::npos) 
                    {
                        elType = LibUtilities::eQuadrilateral;
                    }
                    else
                    {
                        // Default element type in 2D is tetrahedron.
                        elType = LibUtilities::eTetrahedron;
                    }
                }
                
                // Read in number of vertices for element type.
                const int nNodes = GetNnodes(elType);

                for (j = 0; j < m_mesh->m_expDim; ++j)
                {
                    getline(m_mshFile,line);
                    s.clear(); s.str(line);
                    for (k = 0; k < nNodes; ++k)
                    {
                        s >> vertex[j][k];
                    }
                }
                
                // Zero co-ordinates bigger than expansion dimension.
                for (j = m_mesh->m_expDim; j < 3; ++j)
                {
                    for (k = 0; k < nNodes; ++k)
                    {
                        vertex[j][k] = 0.0;
                    }
                }

                // Nektar meshes do not contain a unique list of nodes, so this
                // block constructs a unique set so that elements can be created
                // with unique nodes.
                vector<NodeSharedPtr> nodeList;
                for (k = 0; k < nNodes; ++k)
                {
                    NodeSharedPtr n = boost::shared_ptr<Node>(
                        new Node(nodeCounter++, vertex[0][k], 
                                 vertex[1][k],  vertex[2][k]));
                    nodeList.push_back(n);
                }
                
                elNodes[elType].push_back(nodeList);
                elIds  [elType].push_back(i);
            }
            
            int reorderedId = 0;
            nodeCounter = 0;
            
            for (i = 0; i < elmOrder.size(); ++i)
            {
                LibUtilities::ShapeType elType = elmOrder[i];
                vector<vector<NodeSharedPtr> > &tmp = elNodes[elType];
                
                for (j = 0; j < tmp.size(); ++j)
                {
                    vector<int> tags;
                    map<LibUtilities::ShapeType,int>::iterator compIt = 
                        domainComposite.find(elType);
                    if (compIt == domainComposite.end())
                    {
                        tags.push_back(nComposite);
                        domainComposite[elType] = nComposite;
                        nComposite++;
                    }
                    else
                    {
                        tags.push_back(compIt->second);
                    }
                    
                    elMap[elIds[elType][j]] = reorderedId++;
                    
                    vector<NodeSharedPtr> nodeList = tmp[j];
                    
                    for (k = 0; k < nodeList.size(); ++k)
                    {
                        pair<NodeSet::iterator, bool> testIns = 
                            m_mesh->m_vertexSet.insert(nodeList[k]);
                        
                        if (!testIns.second)
                        {
                            nodeList[k] = *(testIns.first);
                        }
                        else
                        {
                            nodeList[k]->m_id = nodeCounter++;
                        }
                    }
                    
                    // Create linear element
                    ElmtConfig conf(elType,1,false,false);
                    ElementSharedPtr E = GetElementFactory().
                        CreateInstance(elType,conf,nodeList,tags);
                    m_mesh->m_element[E->GetDim()].push_back(E);
                }
            }

            // -- Read in curved data.
            getline(m_mshFile, line);
            if (line.find("CURVE") == string::npos)
            {
                cerr << "Cannot find curved side data." << endl;
                abort();
            }
            
            // Read number of curves.
            getline(m_mshFile, line);
            s.clear(); s.str(line);
            s >> nCurves; 
            
            if (nCurves > 0)
            {
                string curveTag;
                
                for (i = 0; i < nCurves; ++i)
                {
                    getline(m_mshFile, line);
                    s.clear(); s.str(line);
                    s >> word;
                    
                    if (word == "File")
                    {
                        // Next line contains filename and curve tag.
                        getline(m_mshFile, line);
                        s.clear(); s.str(line);
                        s >> word >> curveTag;
                        curveTags[curveTag] = make_pair(eFile, word);
                    }
                    else if (word == "Recon")
                    {
                        // Next line contains curve tag.
                        getline(m_mshFile, line);
                        s.clear(); s.str(line);
                        s >> word >> curveTag;
                        curveTags[curveTag] = make_pair(eRecon, word);
                    }
                    else
                    {
                        cerr << "Unsupported curve type " << word << endl;
                        abort();
                    }
                }
                
                // Load high order surface information.
                LoadHOSurfaces();
                
                // Read in curve information. First line should contain number
                // of curved sides.
                getline(m_mshFile,line);
                
                if (line.find("side") == string::npos)
                {
                    cerr << "Unable to read number of curved sides" << endl;
                    abort();
                }
                
                int nCurvedSides;
                int faceId, elId;
                map<string,pair<NekCurve, string> >::iterator it;
                HOSurfSet::iterator hoIt;

                s.clear(); s.str(line);
                s >> nCurvedSides;
                
                // Iterate over curved sides, and look up high-order surface
                // information in the HOSurfSet, then map this onto faces.
                for (i = 0; i < nCurvedSides; ++i)
                {
                    getline(m_mshFile, line);
                    s.clear(); s.str(line);
                    s >> faceId >> elId >> word;
                    faceId--;
                    elId = elMap[elId-1];
                    ElementSharedPtr el = m_mesh->m_element[m_mesh->m_expDim][elId];
                    
                    if (el->GetConf().m_e == LibUtilities::ePrism && faceId % 2 == 0)
                    {
                        boost::shared_ptr<Prism> p = 
                            boost::dynamic_pointer_cast<Prism>(el);
                        if (p->m_orientation == 1)
                        {
                            faceId = (faceId+2) % 6;
                        }
                        else if (p->m_orientation == 2)
                        {
                            faceId = (faceId+4) % 6;
                        }
                    }
                    else if (el->GetConf().m_e == LibUtilities::eTetrahedron)
                    {
                        boost::shared_ptr<Tetrahedron> t =
                            boost::dynamic_pointer_cast<Tetrahedron>(el);
                        faceId = t->m_orientationMap[faceId];
                    }
                    
                    it = curveTags.find(word);
                    if (it == curveTags.end())
                    {
                        cerr << "Unrecognised curve tag " << word 
                             << " in curved lines" << endl;
                        abort();
                    }
                    
                    if (it->second.first == eRecon)
                    {
                        // Spherigon information: read in vertex normals.
                        vector<NodeSharedPtr> &tmp = 
                            el->GetFace(faceId)->m_vertexList;
                        vector<Node> n(tmp.size());

                        int offset = 0;
                        if (el->GetConf().m_e == LibUtilities::ePrism && faceId % 2 == 1)
                        {
                            offset = boost::dynamic_pointer_cast<Prism>(
                                el)->m_orientation;
                        }
                        
                        // Read x/y/z coordinates.
                        getline(m_mshFile, line);
                        s.clear(); s.str(line);
                        for (j = 0; j < tmp.size(); ++j)
                        {
                            s >> n[j].m_x;
                        }

                        getline(m_mshFile, line);
                        s.clear(); s.str(line);
                        for (j = 0; j < tmp.size(); ++j)
                        {
                            s >> n[j].m_y;
                        }

                        getline(m_mshFile, line);
                        s.clear(); s.str(line);
                        for (j = 0; j < tmp.size(); ++j)
                        {
                            s >> n[j].m_z;
                        }
                            
                        for (j = 0; j < tmp.size(); ++j)
                        {
                            int id = tmp[(j+offset) % tmp.size()]->m_id;
                            boost::unordered_map<int, Node>::iterator vIt =
                                m_mesh->m_vertexNormals.find(id);
                            
                            if (vIt == m_mesh->m_vertexNormals.end())
                            {
                                m_mesh->m_vertexNormals[id] = n[j];
                            }
                        }
                        
                        // Add edge/face to list of faces to apply spherigons
                        // to.
                        m_mesh->m_spherigonSurfs.insert(make_pair(elId, faceId));
                    }
                    else if (it->second.first == eFile)
                    {
                        vector<unsigned int> vertId(3);
                        s >> vertId[0] >> vertId[1] >> vertId[2];
                        
                        // Find vertex combination in hoData.
                        hoIt = hoData[word].find(HOSurfSharedPtr(
                            new HOSurf(vertId)));
                        
                        if (hoIt == hoData[word].end())
                        {
                            cerr << "Unable to find high-order surface data "
                                 << "for element id " << elId+1 << endl;
                            abort();
                        }
                        
                        // Depending on order of vertices in rea file, surface
                        // information may need to be rotated or
                        // reflected. These procedures are taken from
                        // nektar/Hlib/src/HOSurf.C
                        HOSurfSharedPtr surf = *hoIt;
                        
                        if (vertId[0] == surf->vertId[0]) 
                        {
                            if (vertId[1] == surf->vertId[1] || 
                                vertId[1] == surf->vertId[2])
                            {
                                if (vertId[1] == surf->vertId[2])
                                {
                                    surf->Rotate(1);
                                    surf->Reflect();
                                }
                            }
                        }
                        else if (vertId[0] == surf->vertId[1])
                        {
                            if (vertId[1] == surf->vertId[0] ||
                                vertId[1] == surf->vertId[2]) 
                            {
                                if (vertId[1] == surf->vertId[0])
                                {
                                    surf->Reflect();
                                }
                                else
                                {
                                    surf->Rotate(2);
                                }
                            }
                        }
                        else if (vertId[0] == surf->vertId[2])
                        {
                            if (vertId[1] == surf->vertId[0] ||
                                vertId[1] == surf->vertId[1])
                            {
                                if (vertId[1] == surf->vertId[1])
                                {
                                    surf->Rotate(2);
                                    surf->Reflect();
                                }
                                else
                                {
                                    surf->Rotate(1);
                                }
                            }
                        }
                        
                        // If the element is a prism, check to see if
                        // orientation has changed and update order of surface
                        // vertices.
                        int reverseSide = 2;
                        
                        // Prisms may have been rotated by OrientPrism routine
                        // and break curved faces. This block rotates faces
                        // accordingly. TODO: Add similar routine for tets.
                        if (el->GetConf().m_e == LibUtilities::ePrism)
                        {
                            boost::shared_ptr<Prism> pr = 
                                boost::static_pointer_cast<Prism>(el);
                            if (pr->m_orientation == 1)
                            {
                                // Prism has been rotated counter-clockwise;
                                // rotate face, reverse what was the last edge
                                // (now located at edge 0).
                                (*hoIt)->Rotate(1);
                                reverseSide = 0;
                            }
                            else if (pr->m_orientation == 2)
                            {
                                // Prism has been rotated clockwise; rotate
                                // face, reverse what was the last edge (now
                                // located at edge 1).
                                (*hoIt)->Rotate(2);
                                reverseSide = 1;
                            }
                        }
                        
                        // Finally, add high order data to appropriate
                        // face. NOTE: this is a bit of a hack since the
                        // elements are technically linear, but should work just
                        // fine.
                        FaceSharedPtr f    = el->GetFace(faceId);
                        int           Ntot = (*hoIt)->surfVerts.size();
                        int           N    = ((int)sqrt(8.0*Ntot+1.0)-1)/2;
                        EdgeSharedPtr edge;
                        
                        // Apply high-order map to convert face data to Nektar++
                        // ordering (vertices->edges->internal).
                        vector<NodeSharedPtr> tmpVerts = (*hoIt)->surfVerts;
                        for (j = 0; j < tmpVerts.size(); ++j)
                        {
                            (*hoIt)->surfVerts[hoMap[j]] = tmpVerts[j];
                        }
                        
                        for (j = 0; j < tmpVerts.size(); ++j)
                        {
                            NodeSharedPtr a = (*hoIt)->surfVerts[j];
                        }
                        
                        for (j = 0; j < f->m_edgeList.size(); ++j)
                        {
                            edge = f->m_edgeList[j];
                            
                            // Skip over edges which have already been
                            // populated, apart from those which need to be
                            // reoriented.
                            if (edge->m_edgeNodes.size() > 0 && reverseSide == 2)
                            {
                                continue;
                            }
                            
                            edge->m_edgeNodes.clear();
                            edge->m_curveType = LibUtilities::eGaussLobattoLegendre;
                            
                            for (int k = 0; k < N-2; ++k)
                            {
                                edge->m_edgeNodes.push_back(
                                    (*hoIt)->surfVerts[3+j*(N-2)+k]);
                            }
                            
                            // Reverse order of modes along correct side.
                            if (j == reverseSide)
                            {
                                reverse(edge->m_edgeNodes.begin(), 
                                        edge->m_edgeNodes.end());
                            }
                        }

                        f->m_curveType = LibUtilities::eNodalTriElec;
                        for (int j = 3+3*(N-2); j < Ntot; ++j)
                        {
                            f->m_faceNodes.push_back((*hoIt)->surfVerts[j]);
                        }
                    }
                }
            }
            
            // -- Process fluid boundary conditions.

            // Define a fairly horrendous map: key is the condition ID, the
            // value is a vector of pairs of composites and element
            // types. Essentially this map takes conditions -> composites for
            // each element type.
            map<int,vector<pair<int,LibUtilities::ShapeType> > > surfaceCompMap;
            
            // Skip boundary conditions line.
            getline(m_mshFile, line);
            getline(m_mshFile, line);
            
            int nSurfaces = 0;

            while (true)
            {
                getline(m_mshFile, line);
                
                // Break out of loop at end of boundary conditions section.
                if (line.find("*") != string::npos || m_mshFile.eof() ||
                    line.length() == 0)
                {
                    break;
                }
                
                // Read boundary type, element ID and face ID.
                char bcType;
                int elId, faceId;
                s.clear(); s.str(line);
                s >> bcType >> elId >> faceId;
                faceId--;
                elId = elMap[elId-1];
                
                vector<string>        vals;
                vector<ConditionType> type;
                ConditionSharedPtr    c = 
                    MemoryManager<Condition>::AllocateSharedPtr();
                
                // First character on each line describes type of BC. Currently
                // only support V, W, and O. In this switch statement we
                // construct the quantities needed to search for the condition.
                switch(bcType)
                {
                    // Wall boundary.
                    case 'W':
                    {
                        for (i = 0; i < m_mesh->m_fields.size()-1; ++i)
                        {
                            vals.push_back("0");
                            type.push_back(eDirichlet);
                        }
                        // Set high-order boundary condition for wall.
                        vals.push_back("0");
                        type.push_back(eHOPCondition);
                        break;
                    }

                    // Velocity boundary condition (either constant or dependent
                    // upon x,y,z).
                    case 'V':
                    case 'v':
                    {
                        for (i = 0; i < m_mesh->m_fields.size()-1; ++i)
                        {
                            getline(m_mshFile, line);
                            size_t p = line.find_first_of('=');
                            vals.push_back(boost::algorithm::trim_copy(
                                               line.substr(p+1)));
                            type.push_back(eDirichlet);
                        }
                        // Set high-order boundary condition for Dirichlet
                        // condition.
                        vals.push_back("0");
                        type.push_back(eHOPCondition);
                        break;
                    }

                    // Natural outflow condition (default value = 0.0?)
                    case 'O':
                    {
                        for (i = 0; i < m_mesh->m_fields.size(); ++i)
                        {
                            vals.push_back("0");
                            type.push_back(eNeumann);
                        }
                        // Set zero Dirichlet condition for outflow.
                        type[m_mesh->m_fields.size()-1] = eDirichlet;
                        break;
                    }
                    
                    // Ignore unsupported BCs
                    case 'P':
                    case 'E':
                    case 'F':
                    case 'f':
                    case 'N':
                        continue;
                        break;
                        
                    default:
                        cerr << "Unknown boundary condition type " 
                             << line[0] << endl;
                        abort();
                }
                
                // Populate condition information.
                c->field = m_mesh->m_fields;
                c->type  = type;
                c->value = vals;
                
                // Now attempt to find this boundary condition inside
                // m_mesh->condition. This is currently a linear search and should
                // probably be made faster!
                ConditionMap::iterator it;
                bool found = false;
                for (it = m_mesh->m_condition.begin(); it != m_mesh->m_condition.end(); ++it)
                {
                    if (c == it->second)
                    {
                        found = true;
                        break;
                    }
                }
                
                int compTag, conditionId;
                ElementSharedPtr elm = m_mesh->m_element[m_mesh->m_spaceDim][elId];
                ElementSharedPtr surfEl;

                // Create element for face (3D) or segment (2D). At the moment
                // this is a bit of a hack since high-order nodes are not
                // copied, so some output modules (e.g. Gmsh) will not output
                // correctly.
                if (elm->GetDim() == 3)
                {
                    // 3D elements may have been reoriented, so face IDs will
                    // change.
                    if (elm->GetConf().m_e == LibUtilities::ePrism && faceId % 2 == 0)
                    {
                        boost::shared_ptr<Prism> p = 
                            boost::dynamic_pointer_cast<Prism>(elm);
                        if (p->m_orientation == 1)
                        {
                            faceId = (faceId+2) % 6;
                        }
                        else if (p->m_orientation == 2)
                        {
                            faceId = (faceId+4) % 6;
                        }
                    }
                    else if (elm->GetConf().m_e == LibUtilities::eTetrahedron)
                    {
                        boost::shared_ptr<Tetrahedron> t =
                            boost::dynamic_pointer_cast<Tetrahedron>(elm);
                        faceId = t->m_orientationMap[faceId];
                    }
                    
                    FaceSharedPtr f = elm->GetFace(faceId);
                    bool tri = f->m_vertexList.size() == 3;
                    
                    vector<NodeSharedPtr> nodeList;
                    nodeList.insert(nodeList.begin(), 
                                    f->m_vertexList.begin(), 
                                    f->m_vertexList.end());
                                        
                    vector<int> tags;
                    
                    LibUtilities::ShapeType seg = tri ? LibUtilities::eTriangle : 
                        LibUtilities::eQuadrilateral;
                    ElmtConfig conf(seg,1,true,true,false,
                                    LibUtilities::eGaussLobattoLegendre);
                    surfEl = GetElementFactory().
                        CreateInstance(seg,conf,nodeList,tags);
                    
                    // Copy high-order surface information from edges.
                    for (int i = 0; i < f->m_vertexList.size(); ++i)
                    {
                        surfEl->GetEdge(i)->m_edgeNodes = 
                            f->m_edgeList[i]->m_edgeNodes;
                        surfEl->GetEdge(i)->m_curveType = 
                            f->m_edgeList[i]->m_curveType;
                    }
                }
                else
                {
                    EdgeSharedPtr f = elm->GetEdge(faceId);
                    
                    vector<NodeSharedPtr> nodeList;
                    nodeList.push_back(f->m_n1);
                    nodeList.push_back(f->m_n2);
                    
                    vector<int> tags;
                    
                    ElmtConfig conf(LibUtilities::eSegment,1,true,true,false,
                                    LibUtilities::eGaussLobattoLegendre);
                    surfEl = GetElementFactory().
                        CreateInstance(LibUtilities::eSegment,conf,nodeList,tags);
                }
                
                if (!found)
                {
                    // If condition does not already exist, add to condition
                    // list, create new composite tag and put inside
                    // surfaceCompMap.
                    conditionId = m_mesh->m_condition.size();
                    compTag     = nComposite;
                    c->m_composite.push_back(compTag);
                    m_mesh->m_condition[conditionId] = c;
                    
                    surfaceCompMap[conditionId].push_back(
                        pair<int,LibUtilities::ShapeType>(nComposite,surfEl->GetConf().m_e));
                    
                    nComposite++;
                }
                else
                {
                    // Otherwise find existing composite inside surfaceCompMap.
                    map<int,vector<pair<int,LibUtilities::ShapeType> > >::iterator it2;
                    it2 = surfaceCompMap.find(it->first);
                    
                    found = false;
                    if (it2 == surfaceCompMap.end())
                    {
                        // This should never happen!
                        cerr << "Unable to find condition!" << endl;
                        abort();
                    }
                    
                    for (j = 0; j < it2->second.size(); ++j)
                    {
                        pair<int,LibUtilities::ShapeType> tmp = it2->second[j];
                        if (tmp.second == surfEl->GetConf().m_e)
                        {
                            found   = true;
                            compTag = tmp.first;
                            break;
                        }
                    }
                    
                    // If no pairs were found, then this condition contains
                    // multiple element types (i.e. both triangles and
                    // quads). Create another composite for the new shape type
                    // and insert into the map.
                    if (!found)
                    {
                        it2->second.push_back(pair<int,LibUtilities::ShapeType>(
                            nComposite,surfEl->GetConf().m_e));
                        compTag = nComposite;
                        m_mesh->m_condition[it->first]->m_composite.push_back(compTag);
                        nComposite++;
                    }
                    
                    conditionId = it->first;
                }
                
                // Insert composite tag into element and insert element into
                // mesh.
                vector<int> existingTags = surfEl->GetTagList();
                
                existingTags.insert(existingTags.begin(), compTag);
                surfEl->SetTagList (existingTags);
                surfEl->SetId      (nSurfaces);
                
                m_mesh->m_element[surfEl->GetDim()].push_back(surfEl);
                nSurfaces++;
            }

            m_mshFile.close();
            
            // -- Process rest of mesh.
            ProcessEdges     ();
            ProcessFaces     ();
            ProcessElements  ();
            ProcessComposites();
        }

        /**
         * Load high order surface information from hsf file.
         */
        void InputNek::LoadHOSurfaces()
        {
            map<string, pair<NekCurve, string> >::iterator it;
            int nodeId = m_mesh->GetNumEntities();
            
            for (it = curveTags.begin(); it != curveTags.end(); ++it)
            {
                ifstream hsf;
                string   line, fileName = it->second.second;
                size_t   pos;
                int      N, Nface, dot;

                if (it->second.first != eFile)
                {
                    continue;
                }
                
                // Replace fro extension with hsf.
                dot = fileName.find_last_of('.');
                fileName = fileName.substr(0,dot);
                fileName += ".hsf";
                
                // Open hsf file.
                hsf.open(fileName.c_str());
                if (!hsf.is_open())
                {
                    cerr << "Could not open surface file " << fileName << endl;
                    abort();
                }	       
                
                // Read in header line; determine element order, number of faces
                // from this line.
                getline(hsf, line);
                pos = line.find("=");
                if (pos == string::npos)
                {
                    cerr << "hsf header error: cannot read number of "
                         << "nodal points." << endl;
                    abort();
                }
                line = line.substr(pos+1);
                stringstream ss(line);
                ss >> N;
                
                pos = line.find("=");
                if (pos == string::npos)
                {
                    cerr << "hsf header error: cannot read number of "
                         << "faces." << endl;
                    abort();
                }
                line = line.substr(pos+1);
                ss.clear(); ss.str(line);
                ss >> Nface;
                
                int Ntot = N*(N+1)/2;

                // Skip a line, then read in r,s positions inside the next
                // comments.
                Array<OneD, NekDouble> r(Ntot), s(Ntot);
                getline(hsf, line);
                
                for (int i = 0; i < 2; ++i)
                {
                    string word;
                    
                    getline(hsf, line);
                    ss.clear(); ss.str(line);
                    ss >> word;
                    
                    if (word != "#")
                    {
                        cerr << "hsf header error: cannot read in "
                             << "r/s points" << endl;
                        abort();
                    }
                    
                    for (int j = 0; j < Ntot; ++j)
                    {
                        ss >> (i == 0 ? r[j] : s[j]);
                    }
                }

                // Generate electrostatic points so that re-mapping array can
                // be constructed.
                Array<OneD, NekDouble> rp(Ntot), sp(Ntot);
                LibUtilities::PointsKey elec(N, LibUtilities::eNodalTriElec);
                LibUtilities::PointsManager()[elec]->GetPoints(rp,sp);

                // Expensively construct remapping array nodemap. This will
                // map nodal ordering from hsf order to Nektar++ ordering
                // (i.e. vertices followed by edges followed by interior
                // points.)
                for (int i = 0; i < Ntot; ++i)
                {
                    for (int j = 0; j < Ntot; ++j)
                    {
                        if (fabs(r[i]-rp[j]) < 1e-5 && fabs(s[i]-sp[j]) < 1e-5)
                        {
                            hoMap[i] = j;
                            break;
                        }
                    }
                }
               
                // Skip variables line
                getline(hsf, line);

                // Read in nodal points for each face.
                map<unsigned int, vector<NodeSharedPtr> > faceMap;
                for (int i = 0; i < Nface; ++i)
                {
                    getline(hsf, line);
                    vector<NodeSharedPtr> faceNodes(Ntot);
                    for (int j = 0; j < Ntot; ++j, ++nodeId)
                    {
                        double x, y, z;
                        getline(hsf, line);
                        ss.clear(); ss.str(line);
                        ss >> x >> y >> z;
                        faceNodes[j] = NodeSharedPtr(new Node(nodeId, x, y, z));
                    }
                    // Skip over tecplot connectivity information.
                    for (int j = 0; j < (N-1)*(N-1); ++j)
                    {
                        getline(hsf, line);
                    }
                    faceMap[i] = faceNodes;
                }
                
                // Finally, read in connectivity information to set up after
                // reading rea file.
                getline(hsf, line);
                for (int i = 0; i < Nface; ++i)
                {
                    string               tmp;
                    int                  fid;
                    vector<unsigned int> nodeIds(3);

                    getline(hsf, line);
                    ss.clear(); ss.str(line);
                    ss >> tmp >> fid >> nodeIds[0] >> nodeIds[1] >> nodeIds[2];
                    
                    if (tmp != "#")
                    {
                        cerr << "Unable to read hsf connectivity information." 
                             << endl;
                        abort();
                    }
                    
                    hoData[it->first].insert(
                        HOSurfSharedPtr(new HOSurf(nodeIds, faceMap[i])));
                }
                
                hsf.close();
            }
        }
        
        /**
         * This routine aids the reading of Nektar files only; it returns the
         * number of nodes for a given entity typw.
         */
        int InputNek::GetNnodes(LibUtilities::ShapeType InputNekEntity)
        {
            int nNodes = 0;

            switch(InputNekEntity)
            {
            case LibUtilities::ePoint:         nNodes = 1;  break;
            case LibUtilities::eSegment:       nNodes = 2;  break;
            case LibUtilities::eTriangle:      nNodes = 3;  break;
            case LibUtilities::eQuadrilateral: nNodes = 4;  break;
            case LibUtilities::eTetrahedron:   nNodes = 4;  break;
            case LibUtilities::ePrism:         nNodes = 6;  break;
            case LibUtilities::eHexahedron:    nNodes = 8;  break;
            default:
                cerr << "unknown Nektar element type" << endl;
            }

            return nNodes;
        }

        /** 
         * @brief Compares two %HOSurf objects referred to as shared pointers.
         *
         * Two %HOSurf objects are defined to be equal if they contain identical
         * vertex ids contained in HOSurf::vertId.
         */
        bool operator==(HOSurfSharedPtr const &p1, HOSurfSharedPtr const &p2)
        {
            if (p1->vertId.size() != p2->vertId.size())
            {
                return false;
            }
            
            vector<unsigned int> ids1 = p1->vertId;
            vector<unsigned int> ids2 = p2->vertId;
            sort(ids1.begin(), ids1.end());
            sort(ids2.begin(), ids2.end());
            
            for (int i = 0; i < ids1.size(); ++i)
            {
                if (ids1[i] != ids2[i])
                    return false;
            }
            
            return true;
        }
        
        /** 
         * Rotates the triangle of data points inside #surfVerts
         * counter-clockwise nrot times.
         *
         * @param nrot Number of times to rotate triangle.
         */
        void HOSurf::Rotate(int nrot)
        {
            int n, i, j, cnt;
            int np = ((int)sqrt(8.0*surfVerts.size()+1.0)-1)/2;
			NodeSharedPtr* tmp = new NodeSharedPtr[np*np];
            //NodeSharedPtr tmp[np][np];
            
            for (n = 0; n < nrot; ++n) 
            {
                for (cnt = i = 0; i < np; ++i)
                {
                    for (j = 0; j < np-i; ++j, cnt++)
                    {
                        tmp[i*np+j] = surfVerts[cnt];
                    }
                }
                for (cnt = i = 0; i < np; ++i)
                {
                    for (j = 0; j < np-i; ++j,cnt++)
                    {
                        surfVerts[cnt] = tmp[(np-1-i-j)*np+i];
                    }
                }
            }
            
            delete[] tmp;
        }
        
        void HOSurf::Reflect()
        {
            int i, j, cnt;
            int np = ((int)sqrt(8.0*surfVerts.size()+1.0)-1)/2;
            NodeSharedPtr* tmp = new NodeSharedPtr[np*np];
            
            for (cnt = i = 0; i < np; ++i)
            {
                for (j = 0; j < np-i; ++j,cnt++)
                {
                    tmp[i*np+np-i-1-j] = surfVerts[cnt];
                }
            }
            
            for(cnt = i = 0; i < np; ++i)
            {
                for(j = 0; j < np-i; ++j,cnt++)
                {
                    surfVerts[cnt] = tmp[i*np+j];
                }
            }

            delete[] tmp;
        }
    }
}
