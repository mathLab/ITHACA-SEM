////////////////////////////////////////////////////////////////////////////////
//
//  File: OutputGmsh.cpp
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
//  Description: Gmsh file format output.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
using namespace std;

#include "MeshElements.h"
#include "OutputGmsh.h"
#include "InputGmsh.h"

namespace Nektar
{
    namespace Utilities
    {
        ModuleKey OutputGmsh::className = 
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey("msh", eOutputModule), OutputGmsh::create);

        OutputGmsh::OutputGmsh(MeshSharedPtr m) : OutputModule(m)
        {
            map<unsigned int, ElmtConfig>::iterator it;
            
            // Populate #InputGmsh::elmMap and use this to construct an
            // inverse mapping from %ElmtConfig to Gmsh ID.
            for (it = InputGmsh::elmMap.begin(); it != InputGmsh::elmMap.end(); ++it)
            {
                elmMap[it->second] = it->first;
            }
        }

        OutputGmsh::~OutputGmsh()
        {
            
        }
        
        void OutputGmsh::Process()
        {
            // Write MSH header
            mshFile << "$MeshFormat" << endl
                    << "2.2 0 8" << endl
                    << "$EndMeshFormat" << endl;
            
            // Write nodes section.
            mshFile << "$Nodes"            << endl
                    << m->vertexSet.size() << endl;
            
            // Write out nodes.
            NodeSet::iterator it;
            for (it = m->vertexSet.begin(); it != m->vertexSet.end(); ++it)
            {
                mshFile << (*it)->id << " " << (*it)->x << " " 
                        << (*it)->y  << " " << (*it)->z 
                        << endl;
            }
            
            mshFile << "$EndNodes" << endl;
            
            // Write elements section. All other sections are not currently
            // supported (physical names etc).
            mshFile << "$Elements" << endl;
            mshFile << m->GetNumEntities() << endl;
            
            int id = 0;
            
            for (int d = 1; d <= 3; ++d)
            {
                for (int i = 0; i < m->element[d].size(); ++i, ++id)
                {
                    ElementSharedPtr e = m->element[d][i];
                    
                    // First output element ID and type.
                    mshFile << id                   << " " 
                            << elmMap[e->GetConf()] << " ";
                    
                    // Write out number of element tags and then the tags
                    // themselves.
                    vector<int> tags = e->GetTagList();
                    mshFile << tags.size() << " ";
                    
                    for (int j = 0; j < tags.size(); ++j)
                    {
                        mshFile << tags[j] << " ";
                    }
                    
                    // Finally write out node list. First write vertices, then
                    // internal edge nodes, then face nodes.
                    vector<NodeSharedPtr> nodeList = e->GetVertexList();
                    vector<EdgeSharedPtr> edgeList = e->GetEdgeList  ();
                    vector<FaceSharedPtr> faceList = e->GetFaceList  ();
                    
                    tags.clear();
                    
                    for (int j = 0; j < nodeList.size(); ++j)
                    {
                        tags.push_back(nodeList[j]->id);
                    }

                    if (e->GetConf().order > 1)
                    {
                        for (int j = 0; j < edgeList.size(); ++j)
                        {
                            nodeList = edgeList[j]->edgeNodes;
                            for (int k = 0; k < nodeList.size(); ++k)
                            {
                                tags.push_back(nodeList[k]->id);
                            }
                        }
                        
                        for (int j = 0; j < faceList.size(); ++j)
                        {
                            nodeList = faceList[j]->faceNodes;
                            for (int k = 0; k < nodeList.size(); ++k)
                            {
                                tags.push_back(nodeList[k]->id);
                            }
                        }
                    }
                    
                    // Re-order prism vertices.
                    if (e->GetConf().e == ePrism)
                    {
                        // Mirror first in uv plane to swap around
                        // triangular faces
                        swap(tags[0], tags[3]);
                        swap(tags[1], tags[4]);
                        swap(tags[2], tags[5]);
                        // Reorder base points so that face/vertices map
                        // correctly.
                        swap(tags[4], tags[2]);
                        
                        if (e->GetConf().order == 2)
                        {
                            vector<int> nodemap(18);
                            
                            // Vertices remain unchanged.
                            nodemap[ 0] = tags[ 0];
                            nodemap[ 1] = tags[ 1];
                            nodemap[ 2] = tags[ 2];
                            nodemap[ 3] = tags[ 3];
                            nodemap[ 4] = tags[ 4];
                            nodemap[ 5] = tags[ 5];
                            // Reorder edge nodes: first mirror in uv
                            // plane and then place in Nektar++ ordering.
                            nodemap[12] = tags[ 6];
                            nodemap[10] = tags[ 7];
                            nodemap[ 6] = tags[ 8];
                            nodemap[ 8] = tags[ 9];
                            nodemap[13] = tags[10];
                            nodemap[14] = tags[11];
                            nodemap[ 9] = tags[12];
                            nodemap[ 7] = tags[13];
                            nodemap[11] = tags[14];
                            // Face vertices remain unchanged.
                            nodemap[15] = tags[15];
                            nodemap[16] = tags[16];
                            nodemap[17] = tags[17];
                            
                            tags = nodemap;
                        }
                        else if (e->GetConf().order > 2)
                        {
                            cerr << "Error: gmsh prisms only supported up "
                                 << "to second order." << endl;
                            abort();
                        }
                    }
                    
                    // Finally write element nodes.
                    for (int j = 0; j < tags.size(); ++j)
                    {
                        mshFile << tags[j] << " ";
                    }
                    
                    mshFile << endl;
                }
            }
            
            mshFile << "$EndElements" << endl;
        }
    }
}
