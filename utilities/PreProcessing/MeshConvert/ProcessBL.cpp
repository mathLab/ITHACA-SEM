////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessBL.cpp
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
//  Description: Refine boundary layer of elements.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
using namespace std;

#include "MeshElements.h"
#include "ProcessBL.h"

#include <SpatialDomains/MeshGraph.h>
#include <LocalRegions/QuadExp.h>

namespace Nektar
{
    namespace Utilities
    {
        ModuleKey ProcessBL::className = 
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey("bl", eProcessModule), ProcessBL::create);

        ProcessBL::ProcessBL(MeshSharedPtr m) : ProcessModule(m)
        {
            
        }

        ProcessBL::~ProcessBL()
        {
            
        }
        
        void ProcessBL::Process()
        {
            // Create a duplicate of the element list.
            vector<ElementSharedPtr> el = m->element[m->expDim];
            const int nLayers = 6;
            const int layerWidth = 3;
            int nodeId = m->vertexSet.size();
            
            // Erase all elements from the element list.
            m->element[m->expDim].clear();
            
            // Iterate over list of elements of expansion dimension.
            for (int i = 0; i < el.size(); ++i)
            {
                // Hacky hack: this class is currently a test example which
                // takes the bottom row of quadrilaterals and refines them.
                if (el[i]->GetConf().e != eQuadrilateral || i % 2 == 1)
                {
                    m->element[m->expDim].push_back(el[i]);
                    continue;
                }
                
                // Get elemental geometry object.
                SpatialDomains::QuadGeomSharedPtr geom = 
                    boost::dynamic_pointer_cast<SpatialDomains::QuadGeom>(
                        el[i]->GetGeom(m->spaceDim));
                
                // Create basis.
                LibUtilities::BasisKey B0(LibUtilities::eModified_A, 4,
                    LibUtilities::PointsKey(layerWidth, LibUtilities::ePolyEvenlySpaced));
                LibUtilities::BasisKey B1(LibUtilities::eModified_A, 4,
                    LibUtilities::PointsKey(nLayers+1,LibUtilities::ePolyEvenlySpaced));
                
                // Create local region.
                LocalRegions::QuadExpSharedPtr q = 
                    MemoryManager<LocalRegions::QuadExp>::AllocateSharedPtr(B0,B1,geom);
                
                // Grab co-ordinates.
                Array<OneD, NekDouble> x(layerWidth*(nLayers+1));
                Array<OneD, NekDouble> y(layerWidth*(nLayers+1));
                q->GetCoords(x,y);
                
                // Create element layers.
                for (int j = 0; j < nLayers; ++j)
                {
                    // Create corner vertices.
                    vector<NodeSharedPtr> nodeList(4);
                    int offset = j*layerWidth;
                    
                    nodeList[0] = NodeSharedPtr(
                        new Node(nodeId++, 
                                 x[offset], 
                                 y[offset], 0.0));
                    nodeList[1] = NodeSharedPtr(
                        new Node(nodeId++, 
                                 x[offset+layerWidth-1], 
                                 y[offset+layerWidth-1], 0.0));
                    nodeList[2] = NodeSharedPtr(
                        new Node(nodeId++, 
                                 x[offset+2*layerWidth-1], 
                                 y[offset+2*layerWidth-1], 0.0));
                    nodeList[3] = NodeSharedPtr(
                        new Node(nodeId++, 
                                 x[offset+layerWidth], 
                                 y[offset+layerWidth], 0.0));
                    
                    // Create element tags - 0 indicates place the element in
                    // the general domain.
                    vector<int> tags;
                    tags.push_back(0);
                    
                    ElmtConfig conf(eQuadrilateral, 1, false, false);
                    ElementSharedPtr elmt = GetElementFactory().
                        CreateInstance(eQuadrilateral,conf,nodeList,tags);
                    
                    // Add high order information to bottom edge.
                    EdgeSharedPtr e = elmt->GetEdge(0);
                    for (int k = 1; k < layerWidth-1; ++k)
                    {
                        e->edgeNodes.push_back(
                            NodeSharedPtr(
                                new Node(nodeId++,
                                         x[offset+k], 
                                         y[offset+k], 0.0)));
                    }

                    // Add element to mesh.
                    m->element[m->expDim].push_back(elmt);
                }
            }
            
            // Re-process mesh to eliminate duplicate vertices and edges.
            ProcessVertices();
            ProcessEdges();
            ProcessFaces();
            ProcessElements();
            ProcessComposites();
        }
    }
}
