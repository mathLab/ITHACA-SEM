////////////////////////////////////////////////////////////////////////////////
//
//  File: InputSem.cpp
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
//  Description: Semtex session converter.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <fstream>
#include <iostream>
using namespace std;

#include "MeshElements.h"
#include "InputSem.h"

namespace Nektar
{
    namespace Utilities
    {
        ModuleKey InputSem::className = 
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey("sem",eInputModule), InputSem::create);

        /**
         * @brief Initialises the InputSem class.
         * 
         * This function populated the #sectionMap map, which stores the
         * position of sections in the input file.
         */
        InputSem::InputSem(MeshSharedPtr m) : InputModule(m)
        {
            // Read through input file and populate the section map.
            map<string,streampos>::iterator it;
            string                          line, word;
            stringstream                    ss;
            streampos                       linePos;
            
            sectionMap["NODES"]    = 0;
            sectionMap["ELEMENTS"] = 0;
            sectionMap["CURVES"]   = 0;
            sectionMap["SURFACES"] = 0;
            
            while (!mshFile.eof())
            {
                linePos = mshFile.tellg();
                getline(mshFile, line);
                ss.clear();
                ss.str(line);
                ss >> word;
                
                // Iterate over all tokens and see if section exists on this
                // line.
                for (it = sectionMap.begin(); it != sectionMap.end(); ++it)
                {
                    if (word == "<"+it->first)
                    {
                        sectionMap[it->first] = linePos;
                    }
                }
            }
            
            // Clear eofbit and go back to the beginning of the file.
            mshFile.clear();
            mshFile.seekg(0);

            // Check that required sections exist in the file.
            if (sectionMap["NODES"] == 0)
            {
                cerr << "Unable to locate NODES section in session file." << endl;
                abort();
            }
            
            if (sectionMap["ELEMENTS"] == 0)
            {
                cerr << "Unable to locate ELEMENTS section in session file." << endl;
                abort();
            }
        }

        InputSem::~InputSem()
        {

        }


        /**
         * @brief Process a Semtex session file.
         * 
         * Semtex files are defined by a tokenized markup format. These
         * sections have already been located in the file by the constructor
         * and their positions are stored in #sectionMap. The converter only
         * requires the NODES and ELEMENTS sections to exist, but can also
         * read CURVES and SURFACES. High-order curves rely on the meshfile
         * session.msh to be created with the Semtex utility meshpr first.
         * 
         * @param pFilename Filename of Semtex session to read.
         */
        void InputSem::Process()
        {
            m->expDim = 0;
            string line, word, tag;
            int start, end, nVertices, nEntities, nCurves, nSurfaces;
            int id, i, j, k;
            vector<double> hoXData, hoYData;
            ElementType elType = eQuadrilateral;
            ifstream homeshFile;
            stringstream ss;

            cout << "Start reading InputSem..." << endl;
            
            // Begin by reading in list of nodes which define the linear
            // elements.
            mshFile.seekg(sectionMap["NODES"]);
            getline(mshFile, line);
            ss.clear(); ss.str(line);
            ss >> word;
            
            tag       = ss.str();
            start     = tag.find_first_of('=');
            end       = tag.find_first_of('>');
            nVertices = atoi(tag.substr(start+1,end).c_str());
            
            i = id = 0;
            while (i < nVertices)
            {
                getline(mshFile, line);
                if (line.length() < 7) continue;
                ss.clear(); ss.str(line);
                double x = 0, y = 0, z = 0;
                ss >> id >> x >> y >> z;
                
                if ((y * y) > 0.000001 && m->spaceDim != 3)
                {
                    m->spaceDim = 2;
                }
                if ((z * z) > 0.000001)
                {
                    m->spaceDim = 3;
                }
                id -= 1; // counter starts at 0
                m->node.push_back(boost::shared_ptr<Node>(new Node(id, x, y, z)));
                ++i;
            }

            // Now read in elements
            mshFile.seekg(sectionMap["ELEMENTS"]);
            getline(mshFile, line);
            ss.clear(); ss.str(line);
            ss >> word;

            tag       = ss.str();
            start     = tag.find_first_of('=');
            end       = tag.find_first_of('>');
            nEntities = atoi(tag.substr(start+1,end).c_str());

            i = id = 0;
            while (i < nEntities)
            {
                getline(mshFile, line);
                if (line.length() < 18)
                {
                    continue;
                }
                int num_tag = 0, num_nodes = 0;
                
                // Create element tags
                vector<int> tags;
                tags.push_back(0); // composite
                tags.push_back(elType); // element type
                
                // Read element node list
                ss.clear(); ss.str(line);
                ss >> id >> word;
                vector<NodeSharedPtr> nodeList;
                for (j = 0; j < 4; ++j)
                {
                    int node = 0;
                    ss >> node;
                    nodeList.push_back(m->node[node-1]);
                }
                
                // Create element
                ElmtConfig conf(elType,1,false,false);
                ElementSharedPtr E = GetElementFactory().
                    CreateInstance(elType,conf,nodeList,tags);
                
                // Determine mesh expansion dimension
                if (E->GetDim() > m->expDim) {
                    m->expDim = E->GetDim();
                }
                m->element[E->GetDim()].push_back(E);
                ++i;
            }
        
            // Finally, process curves.
            if (sectionMap["CURVES"] != 0)
            {
                int np, nel, nodeId = m->node.size();
                
                mshFile.seekg(sectionMap["CURVES"]);
                getline(mshFile, line);
                ss.clear(); ss.str(line);
                ss >> word;
                
                tag     = ss.str();
                start   = tag.find_first_of('=');
                end     = tag.find_first_of('>');
                nCurves = atoi(tag.substr(start+1,end).c_str());

                // Some session files have empty curves sections; if nCurves
                // is 0, no nead to load high order mesh file.
                if (nCurves > 0)
                {
                    int    ext      = m->inFilename.find_last_of('.');
                    string meshfile = m->inFilename.substr(0,ext) + ".msh";
                    
                    homeshFile.open(meshfile.c_str());
                    if (!homeshFile.is_open())
                    {
                        cerr << "Cannot open or find mesh file: " << meshfile << endl;
                        cerr << "Make sure to run meshpr on your session file first." << endl;
                        abort();
                    }

                    // Make sure we have matching header.
                    getline(homeshFile, line);
                    ss.clear(); ss.str(line);
                    ss >> np >> nel >> nel >> nel;
                    
                    if (nel != m->element[m->expDim].size())
                    {
                        cerr << "Number of elements mismatch in mesh file." << endl;
                        abort();
                    }
                    
                    // Now read in all mesh data. This is horribly inefficient
                    // since not all elements are curved, but it is the
                    // easiest way of finding element data.
                    hoXData.resize(nel*np*np);
                    hoYData.resize(nel*np*np);
                    
                    for (j = 0; j < nel*np*np; ++j)
                    {
                        getline(homeshFile, line);
                        ss.clear(); ss.str(line);
                        ss >> hoXData[j] >> hoYData[j];
                    }
                    
                    homeshFile.close();
                }
                
                i = id = 0;
                while (i < nCurves)
                {
                    getline(mshFile, line);
                    if (line.length() < 18)
                    {
                        continue;
                    }
                    int elmt = 0, side = 0;
                    ss.clear(); ss.str(line);
                    ss >> id >> elmt >> side >> word;
                    id--;
                    elmt--;
                    
                    vector<NodeSharedPtr> edgeNodes;
                    
                    if (word != "<SPLINE>" && word != "<ARC>")
                    {
                        cerr << "Unknown curve tag: " << word << endl;
                        abort();
                    }
                        
                    // See if we have already retrieved high-order data
                    // for this elements; prevents unnecessary computation
                    // for elements with multiple curves.
                    if (m->element[2][elmt]->GetConf().order > 1)
                    {
                        ++i;
                        continue;
                    }
                    
                    // Now set high order data for requested element.
                    for (int side = 0; side < 4; ++side) 
                    {
                        int offset = elmt*np*np;
                        int stride = 0;
                        
                        switch(side)
                        {
                            case 0: // Bottom edge
                                offset += 0;
                                stride  = 1;
                                break;
                            case 1: // Right edge
                                offset += np-1;
                                stride  = np;
                                break;
                            case 2: // Top edge
                                offset += np*np-1;
                                stride  = -1;
                                break;
                            case 3: // Left edge
                                offset += np*(np-1);
                                stride  = -np;
                                break;
                            default:
                                cerr << "Unknown side for curve id " << id << endl;
                                abort();
                        }
                        
                        for (j = 1; j < np-1; ++j, ++nodeId)
                        {
                            double x = hoXData[offset+j*stride];
                            double y = hoYData[offset+j*stride];
                            edgeNodes.push_back(boost::shared_ptr<Node>(
                                new Node(nodeId, x, y, 0.0)));
                        }
                    }
                    
                    // Grab existing element from list and retrieve tags and
                    // vertices; insert these into existing edge nodes.
                    ElementSharedPtr      e      = m->element[2][elmt];
                    vector<NodeSharedPtr> elvert = e->GetVertexList();
                    vector<int>           tags   = e->GetTagList();
                    edgeNodes.insert(edgeNodes.begin(), elvert.begin(), elvert.end());
                    
                    // Create new element and replace with an incomplete
                    // quadrilateral of the correct order.
                    ElmtConfig conf(elType,np-1,false,false,
                                    LibUtilities::eGaussLobattoLegendre);
                    m->element[2][elmt] = GetElementFactory().
                        CreateInstance(elType,conf,edgeNodes,tags);
                    
                    ++i;
                }
            }
            
            // Process surfaces if they exist. Surface support is fairly
            // rudimentary and does not distinguish between different types of
            // boundary condition yet. This is deliberately done after curves
            // to ensure high-order points are preserved.
            if (sectionMap["SURFACES"] != 0)
            {
                mshFile.seekg(sectionMap["SURFACES"]);
                getline(mshFile, line);
                ss.clear(); ss.str(line);
                ss >> word;
                
                tag       = ss.str();
                start     = tag.find_first_of('=');
                end       = tag.find_first_of('>');
                nSurfaces = atoi(tag.substr(start+1,end).c_str());
                
                i = id = 0;
                int elmt, side;
                
                while (i < nSurfaces)
                {
                    getline(mshFile, line);
                    ss.clear(); ss.str(line);
                    ss >> id >> elmt >> side >> word >> tag;
                    elmt--;
                    side--;
                    
                    EdgeSharedPtr edge = m->element[2][elmt]->GetEdge(side);
                    vector<NodeSharedPtr> edgeNodes = edge->edgeNodes;
                    edgeNodes.insert(edgeNodes.begin(),edge->n2);
                    edgeNodes.insert(edgeNodes.begin(),edge->n1);
                    int order = edgeNodes.size()-1;

                    vector<int> tags;
                    tags.push_back(1);
                    tags.push_back(eLine);
                    
                    ElementType seg = eLine;
                    ElmtConfig conf(eLine,order,true,false,
                                    LibUtilities::eGaussLobattoLegendre);
                    ElementSharedPtr E = GetElementFactory().
                        CreateInstance(eLine,conf,edgeNodes,tags);
                    m->element[1].push_back(E);
                    
                    ++i;
                }
            }

            PrintSummary();
            mshFile.close();

            // Process rest of mesh.
            ProcessVertices();
            ProcessEdges();
            ProcessFaces();
            ProcessElements();
            ProcessComposites();
        }
    }
}
