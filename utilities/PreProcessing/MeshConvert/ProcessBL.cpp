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
#include <LocalRegions/PrismExp.h>
#include <LocalRegions/QuadExp.h>
#include <LocalRegions/TriExp.h>

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
            const int nLayers = 10;
            const int layerWidth = 3;
            int nodeId = m->vertexSet.size();
            
            // Erase all elements from the element list.
            m->element[m->expDim].clear();
	     
	   // Iterate over list of elements of expansion dimension.
            for (int i = 0; i < el.size(); ++i)
            {
	      if (el[i]->GetConf().e != ePrism)
		{
		  m->element[m->expDim].push_back(el[i]);		  
		  continue;
		}
	      
              // Find boundary faces if any	
	       int bl0 = el[i]->GetBoundaryLink(0);		 
	       
	      // Get elemental geometry object.
	      SpatialDomains::PrismGeomSharedPtr geom = 
		boost::dynamic_pointer_cast<SpatialDomains::PrismGeom>(
								       el[i]->GetGeom(m->spaceDim));

	      // Default PointsType.
	      LibUtilities::PointsType pt = LibUtilities::eGaussLobattoLegendre;
	      
	      // Create basis.
	      LibUtilities::BasisKey B0(LibUtilities::eModified_A, layerWidth,
					LibUtilities::PointsKey(layerWidth,pt));
	      LibUtilities::BasisKey B1(LibUtilities::eModified_A, 2,
					LibUtilities::PointsKey(nLayers+1,LibUtilities::eBoundaryLayerPoints));
	      LibUtilities::BasisKey B2(LibUtilities::eModified_A, layerWidth,
					LibUtilities::PointsKey(layerWidth,pt));
                
	      // Create local region.
	      LocalRegions::PrismExpSharedPtr q = 
		MemoryManager<LocalRegions::PrismExp>::AllocateSharedPtr(B0,B1,B2,geom);
                
	      // Grab co-ordinates.
	      Array<OneD, NekDouble> x(layerWidth*layerWidth*(nLayers+1));
	      Array<OneD, NekDouble> y(layerWidth*layerWidth*(nLayers+1));
	      Array<OneD, NekDouble> z(layerWidth*layerWidth*(nLayers+1));
	      q->GetCoords(x,y,z);
                
	      // Create element layers.
	      for (int j = 0; j < nLayers; ++j)
		{
		  // Create corner vertices.
		  vector<NodeSharedPtr> nodeList(6);
		  int offset = j*layerWidth;
		  if (j == 0)
		    {
		      nodeList[0] = el[i]->GetVertex(0);
		      
		      nodeList[1] = el[i]->GetVertex(1);

		      nodeList[2] = NodeSharedPtr(
						  new Node(nodeId++, 
							   x[offset+2*layerWidth-1], 
							   y[offset+2*layerWidth-1], 
							   z[offset+2*layerWidth-1]));
		      nodeList[3] = NodeSharedPtr(
						  new Node(nodeId++, 
							   x[offset+layerWidth], 
							   y[offset+layerWidth], z[offset+layerWidth]));
		      nodeList[4] = el[i]->GetVertex(4);

		      nodeList[5] = NodeSharedPtr(
						  new Node(nodeId++, 
							   x[offset+1+layerWidth*(layerWidth-1)*(nLayers+1)+layerWidth], 
							   y[offset+1+layerWidth*(layerWidth-1)*(nLayers+1)+layerWidth],
							   z[offset+1+layerWidth*(layerWidth-1)*(nLayers+1)+layerWidth]));
		    }
		  else if (j == nLayers-1)
		    {
		      nodeList[0] = NodeSharedPtr(
						  new Node(nodeId++, 
							   x[offset], 
							   y[offset], z[offset]));
		      nodeList[1] = NodeSharedPtr(
						  new Node(nodeId++, 
							   x[offset+layerWidth-1], 
							   y[offset+layerWidth-1], 
							   z[offset+layerWidth-1]));
		      nodeList[2] = el[i]->GetVertex(2);

		      nodeList[3] = el[i]->GetVertex(3);

		      nodeList[4] = NodeSharedPtr(
						  new Node(nodeId++, 
							   x[offset+1+layerWidth*(layerWidth-1)*(nLayers+1)], 
							   y[offset+1+layerWidth*(layerWidth-1)*(nLayers+1)],
							   z[offset+1+layerWidth*(layerWidth-1)*(nLayers+1)]));
		      nodeList[5] = el[i]->GetVertex(5);
		    }
		  else
		    {
		      nodeList[0] = NodeSharedPtr(
						  new Node(nodeId++, 
							   x[offset], 
							   y[offset], z[offset]));
		      nodeList[1] = NodeSharedPtr(
						  new Node(nodeId++, 
							   x[offset+layerWidth-1], 
							   y[offset+layerWidth-1], 
							   z[offset+layerWidth-1]));
		      nodeList[2] = NodeSharedPtr(
						  new Node(nodeId++, 
							   x[offset+2*layerWidth-1], 
							   y[offset+2*layerWidth-1], 
							   z[offset+2*layerWidth-1]));
		      nodeList[3] = NodeSharedPtr(
						  new Node(nodeId++, 
							   x[offset+layerWidth], 
							   y[offset+layerWidth], z[offset+layerWidth]));
		      nodeList[4] = NodeSharedPtr(
						  new Node(nodeId++, 
							   x[offset+1+layerWidth*(layerWidth-1)*(nLayers+1)], 
							   y[offset+1+layerWidth*(layerWidth-1)*(nLayers+1)],
							   z[offset+1+layerWidth*(layerWidth-1)*(nLayers+1)]));
		      nodeList[5] = NodeSharedPtr(
						  new Node(nodeId++, 
							   x[offset+1+layerWidth*(layerWidth-1)*(nLayers+1)+layerWidth], 
							   y[offset+1+layerWidth*(layerWidth-1)*(nLayers+1)+layerWidth],
							   z[offset+1+layerWidth*(layerWidth-1)*(nLayers+1)+layerWidth]));
		    }		  
                    
		  // Create element tags - 0 indicates place the element in
		  // the general domain.
		  vector<int> tags;
		  tags.push_back(0);
                    
		  ElmtConfig conf(ePrism, 1, false, false,false);
		  ElementSharedPtr elmt = GetElementFactory().
		    CreateInstance(ePrism,conf,nodeList,tags);
                    
		      // Add high order information to bottom edge.
		      EdgeSharedPtr e0 = elmt->GetEdge(0);
		      EdgeSharedPtr e1 = elmt->GetEdge(2);
                      EdgeSharedPtr e2 = elmt->GetEdge(4);
		      for (int k = 1; k < layerWidth-1; ++k)
		  	{
		  	  e0->edgeNodes.push_back(
		  				 NodeSharedPtr(
		  					       new Node(nodeId++,
		  							x[offset+k], 
		  							y[offset+k], z[offset+k])));
		  	  e1->edgeNodes.push_back(
		  				 NodeSharedPtr(
		  					       new Node(nodeId++,
		  							x[offset+k*(layerWidth)*(nLayers+1)], 
		  							y[offset+k*(layerWidth)*(nLayers+1)],
									z[offset+k*(layerWidth)*(nLayers+1)])));
		  	  e2->edgeNodes.push_back(
		  				 NodeSharedPtr(
		  					       new Node(nodeId++,
		  							x[offset+k*(layerWidth)*(nLayers+1)+layerWidth-1], 
		  							y[offset+k*(layerWidth)*(nLayers+1)+layerWidth-1],
									z[offset+k*(layerWidth)*(nLayers+1)+layerWidth-1])));
		  	}
		      e0->curveType = pt;
		      e1->curveType = pt;
		      e2->curveType = pt;
		      if (bl0 != -1)
			{
			  if (j == 0)
			    {			 
			      ElementSharedPtr e = m->element[m->expDim-1][bl0];
			      e->SetVertex(0, nodeList[0]);
			      e->SetVertex(1, nodeList[1]);
			      e->SetVertex(2, nodeList[2]);
			      e->SetVertex(3, nodeList[3]);
			    }
			  else
			    {
			      vector<NodeSharedPtr> qNodeList(4);
			      for (int k = 0; k < 4; ++k )
				{
				  qNodeList[k] = nodeList[k];
				}
			      vector<int> tagBE;
			      tagBE =  m->element[m->expDim-1][bl0]->GetTagList();
			      ElmtConfig bconf(eQuadrilateral, 1, false, false);
			      ElementSharedPtr boundaryElmt = GetElementFactory().
				CreateInstance(eQuadrilateral,bconf,qNodeList,tagBE);
			      m->element[m->expDim-1].push_back(boundaryElmt);
			    }
			}

		      //Add element to mesh.
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
