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
            int         nParam, nModes, nElements, nCurves, nCurveTypes;
            int         i, j, k, nodeCounter = 0;
            int         nComposite = 0;
            ElementType elType;
            double      vertex[3][6];
            map<ElementType,int> domainComposite;
            map<ElementType,vector< vector<NodeSharedPtr> > > elNodes;
            map<ElementType,vector<int> > elIds;
            boost::unordered_map<int,int> elMap;
            vector<ElementType> elmOrder;

            // Set up vector of processing orders.
            elmOrder.push_back(eLine);
            elmOrder.push_back(eTriangle);
            elmOrder.push_back(eQuadrilateral);
            elmOrder.push_back(ePrism);
            elmOrder.push_back(ePyramid);
            elmOrder.push_back(eTetrahedron);
            elmOrder.push_back(eHexahedron);
             
            m->expDim   = 0;
            m->spaceDim = 0;
            
            if (m->verbose)
            {
                cout << "InputNek: Start reading file..." << endl;
            }
            
            // -- Read in parameters.

            // Ignore first 3 lines. 4th line contains number of parameters.
            for (i = 0; i < 4; ++i)
            {
                getline(mshFile, line);
            }

            stringstream s(line);
            s >> nParam;
            
            for (i = 0; i < nParam; ++i)
            {
                string tmp1, tmp2;
                getline(mshFile, line);
                s.str(line);
                s >> tmp1 >> tmp2;
                if (tmp2 == "MODES")
                {
                    nModes = (int)boost::lexical_cast<double>(tmp1);
                }
            }
            
            // -- Read in passive scalars (ignore)
            getline(mshFile, line);
            s.clear();
            s.str(line);
            s >> j;
            for (i = 0; i < j; ++i)
            {
                getline(mshFile, line);
            }

            // -- Read in logical switches (ignore)
            getline(mshFile, line);
            s.clear();
            s.str(line);
            s >> j;
            for (i = 0; i < j; ++i)
            {
                getline(mshFile, line);
            }

            // -- Read in mesh data.
            
            // First hunt for MESH tag
            bool foundMesh = false;
            while (!mshFile.eof())
            {
                getline(mshFile, line);
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
            getline(mshFile, line);
            s.clear(); s.str(line);
            s >> nElements >> m->expDim;
            m->spaceDim = m->expDim;
            
            // Set up field names.
            m->fields.push_back("u");
            m->fields.push_back("v");
            if (m->spaceDim > 2)
            {
                m->fields.push_back("w");
            }
            m->fields.push_back("p");
            
            // Loop over and create elements.
            for (i = 0; i < nElements; ++i)
            {
                getline(mshFile, line);
                
                if (m->expDim == 2)
                {
                    if (line.find("Qua") != string::npos || 
                        line.find("qua") != string::npos)
                    {
                        elType = eQuadrilateral;
                    }
                    else
                    {
                        // Default element type in 2D is triangle.
                        elType = eTriangle;
                    }
                } 
                else
                {
                    if (line.find("Tet") != string::npos || 
                        line.find("tet") != string::npos) 
                    {
                        elType = eTetrahedron;
                    }
                    else if (line.find("Hex") != string::npos || 
                             line.find("hex") != string::npos) 
                    {
                        elType = eHexahedron;
                    }
                    else if (line.find("Prism") != string::npos || 
                             line.find("prism") != string::npos) 
                    {
                        elType = ePrism;
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
                        elType = eQuadrilateral;
                    }
                    else
                    {
                        // Default element type in 2D is tetrahedron.
                        elType = eTetrahedron;
                    }
                }
                
                // Read in number of vertices for element type.
                const int nNodes = GetNnodes(elType);

                for (j = 0; j < m->expDim; ++j)
                {
                    getline(mshFile,line);
                    s.clear(); s.str(line);
                    for (k = 0; k < nNodes; ++k)
                    {
                        s >> vertex[j][k];
                    }
                }
                
                // Zero co-ordinates bigger than expansion dimension.
                for (j = m->expDim; j < 3; ++j)
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
                ElementType elType = elmOrder[i];
                vector<vector<NodeSharedPtr> > &tmp = elNodes[elType];
                
                for (j = 0; j < tmp.size(); ++j)
                {
                    vector<int> tags;
                    map<ElementType,int>::iterator compIt = 
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
                            m->vertexSet.insert(nodeList[k]);
                        
                        if (!testIns.second)
                        {
                            nodeList[k] = *(testIns.first);
                        }
                        else
                        {
                            nodeList[k]->id = nodeCounter++;
                        }
                    }
                    
                    // Create linear element
                    ElmtConfig conf(elType,1,false,false);
                    ElementSharedPtr E = GetElementFactory().
                        CreateInstance(elType,conf,nodeList,tags);
                    m->element[E->GetDim()].push_back(E);
                }
            }

            // -- Read in curved data.
            getline(mshFile, line);
            if (line.find("CURVE") == string::npos)
            {
                cerr << "Cannot find curved side data." << endl;
                abort();
            }
            
            // Read number of curves.
            getline(mshFile, line);
            s.clear(); s.str(line);
            s >> nCurves; 
            
            if (nCurves > 0)
            {
                string curveTag;
                
                for (i = 0; i < nCurves; ++i)
                {
                    getline(mshFile, line);
                    s.clear(); s.str(line);
                    s >> word;
                    
                    if (word == "File")
                    {
                        // Next line contains filename and curve tag.
                        getline(mshFile, line);
                        s.clear(); s.str(line);
                        s >> word >> curveTag;
                        curveTags[curveTag] = make_pair(eFile, word);
                    }
                    else if (word == "Recon")
                    {
                        // Next line contains curve tag.
                        getline(mshFile, line);
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
                getline(mshFile,line);
                
                if (line.find("side") == string::npos)
                {
                    cerr << "Unable to read number of curved sides" << endl;
                    abort();
                }
                
                int nCurvedSides;
                int faceId, elId, vid1, vid2, vid3;
                map<string,pair<NekCurve, string> >::iterator it;
                HOSurfSet::iterator hoIt;

                s.clear(); s.str(line);
                s >> nCurvedSides;
                int skip = 0;
                
                // Iterate over curved sides, and look up high-order surface
                // information in the HOSurfSet, then map this onto faces.
                for (i = 0; i < nCurvedSides; ++i)
                {
                    getline(mshFile, line);
                    s.clear(); s.str(line);
                    s >> faceId >> elId >> word;
                    faceId--;
                    elId = elMap[elId-1];
                    ElementSharedPtr el = m->element[m->expDim][elId];
                    
                    if (el->GetConf().e == ePrism && faceId % 2 == 0)
                    {
                        boost::shared_ptr<Prism> p = 
                            boost::dynamic_pointer_cast<Prism>(el);
                        if (p->orientation == 1)
                        {
                            faceId = (faceId+2) % 6;
                        }
                        else if (p->orientation == 2)
                        {
                            faceId = (faceId+4) % 6;
                        }
                    }
                    else if (el->GetConf().e == eTetrahedron)
                    {
                        boost::shared_ptr<Tetrahedron> t =
                            boost::dynamic_pointer_cast<Tetrahedron>(el);
                        faceId = t->orientationMap[faceId];
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
                            el->GetFace(faceId)->vertexList;
                        vector<Node> n(tmp.size());

                        int offset = 0;
                        if (el->GetConf().e == ePrism && faceId % 2 == 1)
                        {
                            offset = boost::dynamic_pointer_cast<Prism>(
                                el)->orientation;
                        }
                        
                        // Read x/y/z coordinates.
                        getline(mshFile, line);
                        s.clear(); s.str(line);
                        for (j = 0; j < tmp.size(); ++j)
                        {
                            s >> n[j].x;
                        }

                        getline(mshFile, line);
                        s.clear(); s.str(line);
                        for (j = 0; j < tmp.size(); ++j)
                        {
                            s >> n[j].y;
                        }

                        getline(mshFile, line);
                        s.clear(); s.str(line);
                        for (j = 0; j < tmp.size(); ++j)
                        {
                            s >> n[j].z;
                        }
                            
                        for (j = 0; j < tmp.size(); ++j)
                        {
                            int id = tmp[(j+offset) % tmp.size()]->id;
                            boost::unordered_map<int, Node>::iterator vIt =
                                m->vertexNormals.find(id);
                            
                            if (vIt == m->vertexNormals.end())
                            {
                                m->vertexNormals[id] = n[j];
                            }
                        }
                        
                        // Add edge/face to list of faces to apply spherigons
                        // to.
                        m->spherigonFaces.insert(make_pair(elId, faceId));
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
                        // accordingly.
                        if (el->GetConf().e == ePrism)
                        {
                            boost::shared_ptr<Prism> pr = 
                                boost::static_pointer_cast<Prism>(el);
                            if (pr->orientation == 1)
                            {
                                // Prism has been rotated clockwise; rotate
                                // face, reverse what was the last edge (now
                                // located at edge 0).
                                (*hoIt)->Rotate(1);
                                reverseSide = 0;
                            }
                            else if (pr->orientation == 2)
                            {
                                // Prism has been rotated counter-clockwise;
                                // rotate face, reverse what was the last edge
                                // (now located at edge 1).
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
                        
                        for (j = 0; j < f->edgeList.size(); ++j)
                        {
                            edge = f->edgeList[j];
                            
                            // Skip over edges which have already been populated,
                            // apart from those which need to be reoriented.
                            if (edge->edgeNodes.size() > 0 && reverseSide == 2)
                            {
                                continue;
                            }
                            
                            edge->edgeNodes.clear();
                            edge->curveType = LibUtilities::eGaussLobattoLegendre;
                            
                            for (int k = 0; k < N-2; ++k)
                            {
                                edge->edgeNodes.push_back(
                                    (*hoIt)->surfVerts[3+j*(N-2)+k]);
                            }
                            
                            // Reverse order of modes along correct side.
                            if (j == reverseSide)
                            {
                                reverse(edge->edgeNodes.begin(), 
                                        edge->edgeNodes.end());
                            }
                            
                            for (int k = 3+3*(N-2); k < Ntot; ++k)
                            {
                                f->faceNodes.push_back((*hoIt)->surfVerts[k]);
                            }
                        }
                    }
                }
            }
            
            // -- Process fluid boundary conditions.

            // Define a fairly horrendous map: key is the condition ID, the
            // value is a vector of pairs of composites and element
            // types. Essentially this map takes conditions -> composites for
            // each element type.
            map<int,vector<pair<int,ElementType> > > surfaceCompMap;
            
            // Skip boundary conditions line.
            getline(mshFile, line);
            getline(mshFile, line);
            
            int nSurfaces = 0;

            while (true)
            {
                getline(mshFile, line);
                
                // Break out of loop at end of boundary conditions section.
                if (line.find("*") != string::npos || mshFile.eof() ||
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
                        for (i = 0; i < m->fields.size()-1; ++i)
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
                        for (i = 0; i < m->fields.size()-1; ++i)
                        {
                            getline(mshFile, line);
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
                        for (i = 0; i < m->fields.size(); ++i)
                        {
                            vals.push_back("0");
                            type.push_back(eNeumann);
                        }
                        // Set zero Dirichlet condition for outflow.
                        type[m->fields.size()-1] = eDirichlet;
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
                c->field = m->fields;
                c->type  = type;
                c->value = vals;
                
                // Now attempt to find this boundary condition inside
                // m->condition. This is currently a linear search and should
                // probably be made faster!
                ConditionMap::iterator it;
                bool found = false;
                for (it = m->condition.begin(); it != m->condition.end(); ++it)
                {
                    if (c == it->second)
                    {
                        found = true;
                        break;
                    }
                }
                
                int compTag, conditionId;
                ElementSharedPtr elm = m->element[m->spaceDim][elId];
                ElementSharedPtr surfEl;

                // Create element for face (3D) or segment (2D). At the moment
                // this is a bit of a hack since high-order nodes are not
                // copied, so some output modules (e.g. Gmsh) will not output
                // correctly.
                if (elm->GetDim() == 3)
                {
                    // 3D elements may have been reoriented, so face IDs will
                    // change.
                    if (elm->GetConf().e == ePrism && faceId % 2 == 0)
                    {
                        boost::shared_ptr<Prism> p = 
                            boost::dynamic_pointer_cast<Prism>(elm);
                        if (p->orientation == 1)
                        {
                            faceId = (faceId+2) % 6;
                        }
                        else if (p->orientation == 2)
                        {
                            faceId = (faceId+4) % 6;
                        }
                    }
                    else if (elm->GetConf().e == eTetrahedron)
                    {
                        boost::shared_ptr<Tetrahedron> t =
                            boost::dynamic_pointer_cast<Tetrahedron>(elm);
                        faceId = t->orientationMap[faceId];
                    }
                    
                    FaceSharedPtr f = elm->GetFace(faceId);
                    bool tri = f->vertexList.size() == 3;
                    
                    vector<NodeSharedPtr> nodeList;
                    nodeList.insert(nodeList.begin(), 
                                    f->vertexList.begin(), 
                                    f->vertexList.end());
                                        
                    vector<int> tags;
                    
                    ElementType seg = tri ? eTriangle : eQuadrilateral;
                    ElmtConfig conf(seg,1,true,true,false,
                                    LibUtilities::eGaussLobattoLegendre);
                    surfEl = GetElementFactory().
                        CreateInstance(seg,conf,nodeList,tags);
                    
                    // Copy high-order surface information from edges.
                    for (int i = 0; i < f->vertexList.size(); ++i)
                    {
                        surfEl->GetEdge(i)->edgeNodes = 
                            f->edgeList[i]->edgeNodes;
                        surfEl->GetEdge(i)->curveType = 
                            f->edgeList[i]->curveType;
                    }
                }
                else
                {
                    EdgeSharedPtr f = elm->GetEdge(faceId);
                    
                    vector<NodeSharedPtr> nodeList;
                    nodeList.push_back(f->n1);
                    nodeList.push_back(f->n2);
                    
                    vector<int> tags;
                    
                    ElmtConfig conf(eLine,1,true,true,false,
                                    LibUtilities::eGaussLobattoLegendre);
                    surfEl = GetElementFactory().
                        CreateInstance(eLine,conf,nodeList,tags);
                }
                
                if (!found)
                {
                    // If condition does not already exist, add to condition
                    // list, create new composite tag and put inside
                    // surfaceCompMap.
                    conditionId = m->condition.size();
                    compTag     = nComposite;
                    c->composite.push_back(compTag);
                    m->condition[conditionId] = c;
                    
                    surfaceCompMap[conditionId].push_back(
                        pair<int,ElementType>(nComposite,surfEl->GetConf().e));
                    
                    nComposite++;
                }
                else
                {
                    // Otherwise find existing composite inside surfaceCompMap.
                    map<int,vector<pair<int,ElementType> > >::iterator it2;
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
                        pair<int,ElementType> tmp = it2->second[j];
                        if (tmp.second == surfEl->GetConf().e)
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
                        it2->second.push_back(pair<int,ElementType>(
                            nComposite,surfEl->GetConf().e));
                        compTag = nComposite;
                        m->condition[it->first]->composite.push_back(compTag);
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
                
                m->element[surfEl->GetDim()].push_back(surfEl);
                nSurfaces++;
            }

            mshFile.close();
            
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
            int nodeId = m->GetNumEntities();
            
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
        int InputNek::GetNnodes(ElementType InputNekEntity)
        {
            int nNodes;

            switch(InputNekEntity)
            {
            case ePoint:         nNodes = 1;  break;
            case eLine:          nNodes = 2;  break;
            case eTriangle:      nNodes = 3;  break;
            case eQuadrilateral: nNodes = 4;  break;
            case eTetrahedron:   nNodes = 4;  break;
            case ePrism:         nNodes = 6;  break;
            case eHexahedron:    nNodes = 8;  break;
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
            int n, i, j, cnt;
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
