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

#include "MeshElements.h"
#include "OutputGmsh.h"
#include "InputGmsh.h"

#include <map>
#include <vector>
#include <string>
using namespace std;

namespace Nektar
{
    namespace Utilities
    {
        ModuleKey OutputGmsh::className = 
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eOutputModule, "msh"), OutputGmsh::create,
                "Writes Gmsh msh file.");

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

        /**
         * @brief Process a mesh to output to Gmsh MSH format.
         * 
         * Gmsh output is fairly straightforward. The file first contains a
         * list of nodes, followed by a list of elements. Since
         * Mesh::vertexSet only contains vertices of the linear elements, we
         * first loop over the elements so that any high-order vertices can be
         * enumerated and then added to the node list. We then print out the
         * list of nodes and finally print the element list.
         */
        void OutputGmsh::Process()
        {
            if (m->verbose)
            {
                cout << "OutputGmsh: Writing file..." << endl;
            }

            // Open the file stream.
            OpenStream();
            
            // Write MSH header
            mshFile << "$MeshFormat" << endl
                    << "2.2 0 8" << endl
                    << "$EndMeshFormat" << endl;
            
            int id = m->vertexSet.size();
            vector<ElementSharedPtr> toComplete;
            
            int maxOrder = -1;
            
            // Do first pass over elements of expansion dimension to determine
            // which elements need completion.
            for (int i = 0; i < m->element[m->expDim].size(); ++i)
            {
                ElementSharedPtr e = m->element[m->expDim][i];
                if (e->GetMaxOrder() > maxOrder)
                {
                    maxOrder = e->GetMaxOrder();
                }
            }
            
            //maxOrder = 2;
            
            for (int d = 1; d <= 3; ++d)
            {
                for (int i = 0; i < m->element[d].size(); ++i)
                {
                    ElementSharedPtr e = m->element[d][i];
                    if ((e->GetConf().order <= 1        && maxOrder > 1) ||
                        (e->GetConf().order == maxOrder && e->GetConf().faceNodes == false))
                    {
                        toComplete.push_back(e);
                    }
                    // Generate geometry information for this element. This will
                    // be stored locally inside each element.
                    SpatialDomains::GeometrySharedPtr geom =
                        m->element[d][i]->GetGeom(m->spaceDim);
                }
            }
            
            // Complete these elements.
            for (int i = 0; i < toComplete.size(); ++i)
            {
                toComplete[i]->Complete(maxOrder);
            }
            
            // Do second pass over elements to enumerate high-order vertices.
            for (int d = 1; d <= 3; ++d)
            {
                //cout << "D = " << d << endl;
                for (int i = 0; i < m->element[d].size(); ++i)
                {
                    // Keep track of faces and edges to ensure that high-order
                    // nodes are only added once on common faces/edges.
                    boost::unordered_set<int> edgesDone;
                    boost::unordered_set<int> facesDone;
                    ElementSharedPtr e = m->element[d][i];
                    
                    //cout << "Element " << i << ": ";
                    
                    if (e->GetConf().order > 1)
                    {
                        vector<NodeSharedPtr> tmp;
                        vector<EdgeSharedPtr> edgeList = e->GetEdgeList();
                        vector<FaceSharedPtr> faceList = e->GetFaceList();
                        vector<NodeSharedPtr> volList  = e->GetVolumeNodes();
                        
                        /*
                        cout << " edge = " << edgeList.size() << " "
                             << " face = " << faceList.size() << " "
                             << " vol  = " << volList.size()  << endl;
                        */
                        
                        for (int j = 0; j < edgeList.size(); ++j)
                        {
                            boost::unordered_set<int>::iterator it = 
                                edgesDone.find(edgeList[j]->id);
                            if (it == edgesDone.end() || d != 3)
                            {
                                tmp.insert(tmp.end(), 
                                           edgeList[j]->edgeNodes.begin(),
                                           edgeList[j]->edgeNodes.end());
                                edgesDone.insert(edgeList[j]->id);
                            }
                        }
                        
                        for (int j = 0; j < faceList.size(); ++j)
                        {
                            boost::unordered_set<int>::iterator it = 
                                facesDone.find(faceList[j]->id);
                            if (it == facesDone.end() || d != 3)
                            {
                                tmp.insert(tmp.end(), 
                                           faceList[j]->faceNodes.begin(),
                                           faceList[j]->faceNodes.end());
                                facesDone.insert(faceList[j]->id);
                            }
                        }
                        
                        tmp.insert(tmp.end(), volList.begin(), volList.end());
                        
                        // Even though faces/edges are at this point unique
                        // across the mesh, still need to test inserts since
                        // high-order nodes may already have been inserted into
                        // the list from an adjoining element or a boundary
                        // element.
                        for (int j = 0; j < tmp.size(); ++j)
                        {
                            pair<NodeSet::iterator, bool> testIns =
                                m->vertexSet.insert(tmp[j]);
                            
                            if (testIns.second)
                            {
                                (*(testIns.first))->id = id++;
                            }
                            else
                            {
                                tmp[j]->id = (*(testIns.first))->id;
                            }
                        }
                    }
                }
            }
            
            // Create ordered set of nodes - not required but looks nicer.
            std::set<NodeSharedPtr>::iterator it;
            std::set<NodeSharedPtr> tmp(m->vertexSet.begin(), m->vertexSet.end());

            // Write out nodes section.
            mshFile << "$Nodes"            << endl
                    << m->vertexSet.size() << endl;
            
            for (it = tmp.begin(); it != tmp.end(); ++it)
            {
                mshFile << (*it)->id << " " << scientific << setprecision(10)
                        << (*it)->x << " " 
                        << (*it)->y  << " " << (*it)->z 
                        << endl;
            }
            
            mshFile << "$EndNodes" << endl;
            
            // Write elements section. All other sections are not currently
            // supported (physical names etc).
            mshFile << "$Elements" << endl;
            mshFile << m->GetNumEntities() << endl;
            
            id = 0;
            
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
                    
                    if (tags.size() == 1)
                    {
                        tags.push_back(tags[0]);
                        tags.push_back(0);
                    }
                    
                    mshFile << tags.size() << " ";
                    
                    for (int j = 0; j < tags.size(); ++j)
                    {
                        mshFile << tags[j] << " ";
                    }
                    
                    // Finally write out node list. First write vertices, then
                    // internal edge nodes, then face nodes.
                    vector<NodeSharedPtr> nodeList = e->GetVertexList ();
                    vector<EdgeSharedPtr> edgeList = e->GetEdgeList   ();
                    vector<FaceSharedPtr> faceList = e->GetFaceList   ();
                    vector<NodeSharedPtr> volList  = e->GetVolumeNodes();
                    
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
                                //cout << "EDGENODE" << endl;
                            }
                        }
                        
                        for (int j = 0; j < faceList.size(); ++j)
                        {
                            nodeList = faceList[j]->faceNodes;
                            for (int k = 0; k < nodeList.size(); ++k)
                            {
                                //cout << "FACENODE" << endl;
                                tags.push_back(nodeList[k]->id);
                            }
                        }
                        
                        for (int j = 0; j < volList.size(); ++j)
                        {
                            //cout << "VOLNODE" << endl;
                            tags.push_back(volList[j]->id);
                        }
                    }

                    // Re-order tetrahedral vertices.
                    if (e->GetConf().e == eTetrahedron)
                    {
                        int order = e->GetConf().order;
                        if (order > 4)
                        {
                            cerr << "Temporary error: Gmsh tets only supported "
                                 << "up to 4th order - will fix soon!" << endl;
                            abort();
                        }
                        int pos = 4;
                        // Swap edge 1->3 nodes with edge 2->3 nodes.
                        pos = 4 + 4*(order-1);
                        for (int j = 0; j < order-1; ++j)
                        {
                            swap(tags[j+pos], tags[j+pos+order-1]);
                        }
                        // Reverse ordering of other vertical edge-interior
                        // nodes.
                        reverse(tags.begin()+4+3*(order-1), tags.begin()+4+4*(order-1));
                        reverse(tags.begin()+4+4*(order-1), tags.begin()+4+5*(order-1));
                        reverse(tags.begin()+4+5*(order-1), tags.begin()+4+6*(order-1));
                        
                        // Swap face 2 nodes with face 3.
                        pos = 4 + 6*(order-1) + 2*(order-2)*(order-1)/2;
                        for (int j = 0; j < (order-2)*(order-1)/2; ++j)
                        {
                            swap(tags[j+pos], tags[j+pos+(order-2)*(order-1)/2]);
                        }
                        
                        // Re-order face points. Gmsh ordering (node->face) is:
                        //
                        // Face 0: 0->2->1
                        // Face 1: 0->1->3
                        // Face 2: 0->3->2
                        // Face 3: 3->1->2
                        //
                        // Therefore need to reorder nodes for faces 0, 2 and
                        // 3 to match nodal ordering.
                        
                        // Re-order face 0: transpose
                        vector<int> tmp((order-2)*(order-1)/2);
                        int a = 0;
                        pos = 4 + 6*(order-1);
                        for (int j = 0; j < order-2; ++j)
                        {
                            for (int k = 0; k < order-2-j; ++k, ++a)
                            {
                                tmp[a] = tags[pos+j+k*(2*(order-2)+1-k)/2];
                            }
                        }
                        for (int j = 0; j < (order-1)*(order-2)/2; ++j)
                        {
                            tags[pos+j] = tmp[j];
                        }
                        
                        // Re-order face 2: transpose
                        pos = 4 + 6*(order-1) + 2*(order-2)*(order-1)/2;
                        a = 0;
                        for (int j = 0; j < order-2; ++j)
                        {
                            for (int k = 0; k < order-2-j; ++k, ++a)
                            {
                                tmp[a] = tags[pos+j+k*(2*(order-2)+1-k)/2];
                            }
                        }
                        for (int j = 0; j < (order-1)*(order-2)/2; ++j)
                        {
                            tags[pos+j] = tmp[j];
                        }
                        
                        // Re-order face 3: reflect in y direction
                        pos = 4 + 6*(order-1)+3*(order-2)*(order-1)/2;
                        a = 0;
                        for (int j = 0; j < order-2; ++j)
                        {
                            for (int k = order-3-j; k >= 0; --k, ++a)
                            {
                                tmp[a] = tags[pos+j+k*(2*(order-2)+1-k)/2];
                            }
                        }
                        
                        for (int j = 0; j < (order-1)*(order-2)/2; ++j)
                        {
                            tags[pos+j] = tmp[j];
                        }
                    }
                    // Re-order prism vertices.
                    else if (e->GetConf().e == ePrism)
                    {
                        int order = e->GetConf().order;
                        if (order > 2)
                        {
                            cerr << "Temporary error: Gmsh prisms only "
                                 << "supported up to 2nd order!" << endl;
                            abort();
                        }
                        
                        // Swap nodes.
                        vector<int> temp (18);
                        temp[0] = tags[0];
                        temp[1] = tags[1];
                        temp[2] = tags[4];
                        temp[3] = tags[3];
                        temp[4] = tags[2];
                        temp[5] = tags[5];
                        temp[6] = tags[6];
                        temp[7] = tags[10];
                        temp[8] = tags[9];
                        temp[9] = tags[11];
                        temp[10] = tags[7];
                        temp[11] = tags[14];
                        temp[12] = tags[8];
                        temp[13] = tags[13];
                        temp[14] = tags[12];
                        temp[15] = tags[15];
                        temp[16] = tags[17];
                        temp[17] = tags[16];
                        for (int k = 0; k < 18; ++k)
                        {
                            tags[k] = temp[k];
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
