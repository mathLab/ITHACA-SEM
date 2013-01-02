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
                ModuleKey(eInputModule, "sem"), InputSem::create,
                "Reads Semtex session files.");

        /**
         * @brief Initialises the InputSem class.
         */
        InputSem::InputSem(MeshSharedPtr m) : InputModule(m)
        {
            
        }

        InputSem::~InputSem()
        {
            
        }

        /**
         * @brief Process a Semtex session file.
         * 
         * Semtex files are defined by a tokenized markup format. We first
         * populate #sectionMap which stores the location of the various tags in
         * the session file so that they can be jumped to, since no ordering is
         * defined. The converter only requires the NODES and ELEMENTS sections
         * to exist, but can also read CURVES and SURFACES. High-order curves
         * rely on the meshfile session.msh to be created with the Semtex
         * utility meshpr first.
         * 
         * @param pFilename Filename of Semtex session to read.
         */
        void InputSem::Process()
        {
            // Open the file stream.
            OpenStream();

            if (m->verbose)
            {
                cout << "InputSem: Start reading file..." << endl;
            }

            // Read through input file and populate the section map.
            map<string,streampos>::iterator it;
            string                          line, word;
            stringstream                    ss;
            streampos                       linePos;
            
            sectionMap["NODES"]    = -1;
            sectionMap["ELEMENTS"] = -1;
            sectionMap["CURVES"]   = -1;
            sectionMap["SURFACES"] = -1;
            sectionMap["GROUPS"]   = -1;
            sectionMap["BCS"]      = -1;
            sectionMap["FIELDS"]   = -1;
            
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
                    if (word == "<"+it->first || word == "<"+it->first+">")
                    {
                        sectionMap[it->first] = linePos;
                    }
                }
            }
            
            // Clear eofbit and go back to the beginning of the file.
            mshFile.clear();
            mshFile.seekg(0);

            // Check that required sections exist in the file.
            if (sectionMap["NODES"] == std::streampos(-1))
            {
                cerr << "Unable to locate NODES section in session file." 
                     << endl;
                abort();
            }
            
            if (sectionMap["ELEMENTS"] == std::streampos(-1))
            {
                cerr << "Unable to locate ELEMENTS section in session file." 
                     << endl;
                abort();
            }

            if (sectionMap["SURFACES"] != std::streampos(-1))
            {
                if (sectionMap["BCS"] == std::streampos(-1))
                {
                    cerr << "SURFACES section defined but BCS section not "
                         << "found." << endl;
                    abort();
                }
                
                if (sectionMap["GROUPS"] == std::streampos(-1))
                {
                    cerr << "SURFACES section defined but GROUPS section not "
                         << "found." << endl;
                    abort();
                }

                if (sectionMap["FIELDS"] == std::streampos(-1))
                {
                    cerr << "SURFACES section defined but FIELDS section not "
                         << "found." << endl;
                    abort();
                }
            }

            m->expDim = 0;
            string tag;
            int start, end, nVertices, nEntities, nCurves, nSurf, nGroups, nBCs;
            int id, i, j, k;
            vector<double> hoXData, hoYData;
            ElementType elType = eQuadrilateral;
            ifstream homeshFile;

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
            if (sectionMap["CURVES"] != std::streampos(-1))
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
                    string fname    = config["infile"].as<string>();
                    int    ext      = fname.find_last_of('.');
                    string meshfile = fname.substr(0,ext) + ".msh";
                    
                    homeshFile.open(meshfile.c_str());
                    if (!homeshFile.is_open())
                    {
                        cerr << "Cannot open or find mesh file: " 
                             << meshfile << endl
                             << "Make sure to run meshpr on your session "
                             << "file first." << endl;
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
                            NodeSharedPtr n = boost::shared_ptr<Node>(
                                new Node(nodeId, x, y, 0.0));
                            edgeNodes.push_back(n);
                        }
                    }
                    
                    // Add internal points.
                    for (j = 1; j < np-1; ++j)
                    {
                        int offset = j*np+1;
                        for (k = 1; k < np-1; ++k, ++nodeId)
                        {
                            double x = hoXData[offset+k];
                            double y = hoYData[offset+k];
                            NodeSharedPtr n = boost::shared_ptr<Node>(
                                new Node(nodeId, x, y, 0.0));
                            edgeNodes.push_back(n);
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
                    ElmtConfig conf(elType,np-1,true,false,true,
                                    LibUtilities::eGaussLobattoLegendre);
                    m->element[2][elmt] = GetElementFactory().
                        CreateInstance(elType,conf,edgeNodes,tags);
                    
                    ++i;
                }
            }

            // Process field names
            if (sectionMap["FIELDS"] != std::streampos(-1))
            {
                mshFile.seekg(sectionMap["FIELDS"]);
                getline(mshFile, line);
                getline(mshFile, line);
                ss.clear(); ss.str(line);
                
                while (ss >> tag)
                {
                    m->fields.push_back(tag);
                }
            }
            
            // Process surfaces if they exist. This is deliberately done after
            // curves to ensure high-order points are preserved.
            if (sectionMap["SURFACES"] != std::streampos(-1))
            {
                map<string,int> conditionMap;
                int             maxTag = -1;
                
                // First read in list of groups, which defines each condition tag.
                mshFile.seekg(sectionMap["GROUPS"]);
                getline(mshFile, line);
                ss.clear(); ss.str(line);
                ss >> word;
                
                tag     = ss.str();
                start   = tag.find_first_of('=');
                end     = tag.find_first_of('>');
                nGroups = atoi(tag.substr(start+1,end).c_str());
                
                i = id = 0;
                while (i < nGroups)
                {
                    getline(mshFile, line);
                    ss.clear(); ss.str(line);
                    ss >> id >> tag;
                    conditionMap[tag] = i++;
                }
                
                maxTag = i;

                // Now read in actual values for boundary conditions from BCS
                // section.
                mshFile.seekg(sectionMap["BCS"]);
                getline(mshFile, line);
                ss.clear(); ss.str(line);
                ss >> word;
                
                tag   = ss.str();
                start = tag.find_first_of('=');
                end   = tag.find_first_of('>');
                nBCs  = atoi(tag.substr(start+1,end).c_str());
                
                i = id = 0;
                while (i < nBCs)
                {
                    int                nF;
                    string             tmp;
                    ConditionSharedPtr p;
                    getline(mshFile, line);
                    ss.clear(); ss.str(line);
                    ss >> id >> tag >> nF;

                    p = ConditionSharedPtr(new Condition());
                    m->condition[conditionMap[tag]] = p;
                    
                    // Read boundary condition.
                    j = 0;
                    while (j < nF)
                    {
                        getline(mshFile, line);
                        ss.clear(); ss.str(line);
                        ss >> tmp;
                        
                        // First string should be condition type.
                        if (tmp == "<D>")
                        {
                            p->type.push_back(eDirichlet);
                        }
                        else if (tmp == "<N>")
                        {
                            p->type.push_back(eNeumann);
                        }
                        else if (tmp == "<H>")
                        {
                            p->type. push_back(eHOPCondition);
                            p->value.push_back("0");
                            p->field.push_back("p");
                            ++j;
                            continue;
                        }
                        else
                        {
                            cerr << "Unsupported boundary condition type " << tmp << endl;
                            abort();
                        }
                        
                        // Second string should be field.
                        ss >> tmp;
                        p->field.push_back(tmp);
                        
                        // Third string should be equals sign.
                        ss >> tmp;
                        if (tmp != "=")
                        {
                            cerr << "Couldn't read boundary condition type " << tag << endl;
                            abort();
                        }
                        
                        // Fourth string should be value. CAUTION: Assumes
                        // expression is defined without any spaces in it!
                        ss >> tmp;
                        p->value.push_back(tmp);

                        ++j;
                    }
                    
                    // Finally set composite for condition. In this case, all
                    // composites will be lines so there is one set per
                    // composite.
                    p->composite.push_back(conditionMap[tag]+1);
                    
                    ++i;
                }
                
                // Finally read surface information.
                mshFile.seekg(sectionMap["SURFACES"]);
                getline(mshFile, line);
                ss.clear(); ss.str(line);
                ss >> word;
                
                tag   = ss.str();
                start = tag.find_first_of('=');
                end   = tag.find_first_of('>');
                nSurf = atoi(tag.substr(start+1,end).c_str());
                
                i = id = 0;
                int elmt, side;
                int periodicTagId = -1;
                
                set<pair<int, int> > visitedPeriodic;
                
                while (i < nSurf)
                {
                    getline(mshFile, line);
                    ss.clear(); ss.str(line);
                    ss >> id >> elmt >> side >> word;
                    elmt--;
                    side--;
                    
                    if (word == "<P>")
                    {
                        // If this is the first periodic boundary condition
                        // encountered, then set up m->condition with two
                        // periodic conditions.
                        if (periodicTagId == -1)
                        {
                            periodicTagId = maxTag;
                            ConditionSharedPtr in  = 
                                ConditionSharedPtr(new Condition());
                            ConditionSharedPtr out = 
                                ConditionSharedPtr(new Condition());
                            for (j = 0; j < m->fields.size(); ++j)
                            {
                                in-> type.push_back(ePeriodic);
                                out->type.push_back(ePeriodic);
                                in-> field.push_back(m->fields[j]);
                                out->field.push_back(m->fields[j]);
                                in-> value.push_back("["+boost::lexical_cast<
                                    string>(periodicTagId+1)+"]");
                                out->value.push_back("["+boost::lexical_cast<
                                    string>(periodicTagId)+"]");
                            }
                            in-> composite.push_back(periodicTagId+1);
                            out->composite.push_back(periodicTagId+2);
                            m->  condition[periodicTagId]   = in;
                            m->  condition[periodicTagId+1] = out;
                        }
                        
                        int elmtB, sideB;
                        
                        ss >> elmtB >> sideB;
                        elmtB--;
                        sideB--;

                        pair<int, int> c1(elmt,  side );
                        pair<int, int> c2(elmtB, sideB);
                        
                        if (visitedPeriodic.count(c1) == 0 && 
                            visitedPeriodic.count(c2) == 0)
                        {
                            visitedPeriodic.insert(make_pair(elmtB, sideB));
                            visitedPeriodic.insert(make_pair(elmt,  side ));
                            insertEdge(elmt,  side,  periodicTagId+1);
                            insertEdge(elmtB, sideB, periodicTagId+2);
                        }
                    }
                    else if (word == "<B>")
                    {
                        ss >> tag;
                        insertEdge(elmt, side, conditionMap[tag]+1);
                    }
                    else
                    {
                        cerr << "Unrecognised or unsupported tag " << word << endl;
                        abort();
                    }
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
        
        void InputSem::insertEdge(int elmt, int side, int tagId)
        {
            EdgeSharedPtr edge = m->element[2][elmt]->GetEdge(side);
            vector<NodeSharedPtr> edgeNodes = edge->edgeNodes;
            edgeNodes.insert(edgeNodes.begin(),edge->n2);
            edgeNodes.insert(edgeNodes.begin(),edge->n1);
            int order = edgeNodes.size()-1;
            
            vector<int> tags;
            tags.push_back(tagId);
            
            ElementType seg = eLine;
            ElmtConfig conf(eLine, order, order > 1, false, true,
                            LibUtilities::eGaussLobattoLegendre);
            ElementSharedPtr E = GetElementFactory().
                CreateInstance(eLine,conf,edgeNodes,tags);
            m->element[1].push_back(E);
        }
    }
}
