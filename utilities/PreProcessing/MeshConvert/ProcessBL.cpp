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
#include <LibUtilities/BasicUtils/SharedArray.hpp>

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

      // Stores nodes which have been created in the refinement process.
      vector<NodeSharedPtr> createdNodes;
      
      // Very expensive function - search for this node to see if it has
      // already been created.
      NodeSharedPtr checkNode(Node *n)
      {
          int i;
          for (i = 0; i < createdNodes.size(); ++i)
          {
              NodeSharedPtr n1 = createdNodes[i];
              if (((n1->x-n->x)*(n1->x-n->x) + (n1->y-n->y)*(n1->y-n->y) + (n1->z-n->z)*(n1->z-n->z)) < 1.0e-14)
              {
                  return n1;
              }
          }
          createdNodes.push_back(NodeSharedPtr(n));
          return createdNodes.back();
      }
        
    void ProcessBL::Process()
    {
      // Create a duplicate of the element list.
      vector<ElementSharedPtr> el = m->element[m->expDim];
      const int nLayers = 10;
      const int layerWidth = 3;
      const int tetsOn = 1; // Set to 0 for prismatic boundary layer and 1 for tetrahedral boundary layer.
      int nodeId = m->vertexSet.size();
      int offsetSurf = 0; // Counter for surface elements to be removed in Case 5 of tetrahedral spliting.
	    	    
      // Erase all elements from the element list.
      m->element[m->expDim].clear();
	     
      // Iterate over list of elements of expansion dimension.
      for (int i = 0; i < el.size(); ++i)
	{
            cout << i << endl;
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
	      if (j == 0) // For the first layer use nodes of the original prism at the bottom.
		{
		  nodeList[0] = el[i]->GetVertex(0);
		      
		  nodeList[1] = el[i]->GetVertex(1);

		  nodeList[2] = checkNode(
					      new Node(nodeId++, 
						       x[offset+2*layerWidth-1], 
						       y[offset+2*layerWidth-1], 
						       z[offset+2*layerWidth-1]));
		  nodeList[3] = checkNode(
					      new Node(nodeId++, 
						       x[offset+layerWidth], 
						       y[offset+layerWidth], z[offset+layerWidth]));
		  nodeList[4] = el[i]->GetVertex(4);

		  nodeList[5] = checkNode(
					      new Node(nodeId++, 
						       x[offset+1+layerWidth*(layerWidth-1)*(nLayers+1)+layerWidth], 
						       y[offset+1+layerWidth*(layerWidth-1)*(nLayers+1)+layerWidth],
						       z[offset+1+layerWidth*(layerWidth-1)*(nLayers+1)+layerWidth]));
		}
	      else if (j == nLayers-1) // For the last layer use nodes of the original prism at the top.
		{
		  nodeList[0] = checkNode(
					      new Node(nodeId++, 
						       x[offset], 
						       y[offset], z[offset]));
		  nodeList[1] = checkNode(
					      new Node(nodeId++, 
						       x[offset+layerWidth-1], 
						       y[offset+layerWidth-1], 
						       z[offset+layerWidth-1]));
		  nodeList[2] = el[i]->GetVertex(2);

		  nodeList[3] = el[i]->GetVertex(3);

		  nodeList[4] = checkNode(
					      new Node(nodeId++, 
						       x[offset+1+layerWidth*(layerWidth-1)*(nLayers+1)], 
						       y[offset+1+layerWidth*(layerWidth-1)*(nLayers+1)],
						       z[offset+1+layerWidth*(layerWidth-1)*(nLayers+1)]));
		  nodeList[5] = el[i]->GetVertex(5);
		}
	      else // All the intermediate layers.
		{
		  nodeList[0] = checkNode(
					      new Node(nodeId++, 
						       x[offset], 
						       y[offset],
						       z[offset]));
		  nodeList[1] = checkNode(
					      new Node(nodeId++, 
						       x[offset+layerWidth-1], 
						       y[offset+layerWidth-1], 
						       z[offset+layerWidth-1]));
		  nodeList[2] = checkNode(
					      new Node(nodeId++, 
						       x[offset+2*layerWidth-1], 
						       y[offset+2*layerWidth-1], 
						       z[offset+2*layerWidth-1]));
		  nodeList[3] = checkNode(
					      new Node(nodeId++, 
						       x[offset+layerWidth], 
						       y[offset+layerWidth],
						       z[offset+layerWidth]));
		  nodeList[4] = checkNode(
					      new Node(nodeId++, 
						       x[offset+1+layerWidth*(layerWidth-1)*(nLayers+1)], 
						       y[offset+1+layerWidth*(layerWidth-1)*(nLayers+1)],
						       z[offset+1+layerWidth*(layerWidth-1)*(nLayers+1)]));
		  nodeList[5] = checkNode(
					      new Node(nodeId++, 
						       x[offset+1+layerWidth*(layerWidth-1)*(nLayers+1)+layerWidth], 
						       y[offset+1+layerWidth*(layerWidth-1)*(nLayers+1)+layerWidth],
						       z[offset+1+layerWidth*(layerWidth-1)*(nLayers+1)+layerWidth]));
		}		  
                    
	      // Create element tags - 0 indicates place the element in
	      // the general domain.
	      vector<int> tags;
	      tags.push_back(0);
                  
	      // Create the element.
	      ElmtConfig conf(ePrism, 1, false, false,false);
	      ElementSharedPtr elmt = GetElementFactory().
		CreateInstance(ePrism,conf,nodeList,tags); 

	      // Add high order information to bottom edge of each layer.
	      EdgeSharedPtr e0 = elmt->GetEdge(0);
	      EdgeSharedPtr e2 = elmt->GetEdge(2);
	      EdgeSharedPtr e4 = elmt->GetEdge(4);
	      EdgeSharedPtr e5 = elmt->GetEdge(5);
	      EdgeSharedPtr e6 = elmt->GetEdge(6);
	      EdgeSharedPtr e7 = elmt->GetEdge(7);
	      for (int k = 1; k < layerWidth-1; ++k)
		{
		  e0->edgeNodes.push_back(
					  checkNode(
							new Node(nodeId++,
								 x[offset+k], 
								 y[offset+k],
								 z[offset+k])));
		  e2->edgeNodes.push_back(
					  checkNode(
							new Node(nodeId++,
								 x[offset+k*(layerWidth)*(nLayers+1)], 
								 y[offset+k*(layerWidth)*(nLayers+1)],
								 z[offset+k*(layerWidth)*(nLayers+1)])));
		  e4->edgeNodes.push_back(
					  checkNode(
							new Node(nodeId++,
								 x[offset+k*(layerWidth)*(nLayers+1)+layerWidth-1], 
								 y[offset+k*(layerWidth)*(nLayers+1)+layerWidth-1],
								 z[offset+k*(layerWidth)*(nLayers+1)+layerWidth-1])));
		  e5->edgeNodes.push_back(
					  checkNode(
							new Node(nodeId++,
								 x[offset+layerWidth+k*(layerWidth)*(nLayers+1)+layerWidth-1], 
								 y[offset+layerWidth+k*(layerWidth)*(nLayers+1)+layerWidth-1],
								 z[offset+layerWidth+k*(layerWidth)*(nLayers+1)+layerWidth-1])));
		  e6->edgeNodes.push_back(
					  checkNode(
							new Node(nodeId++,
								 x[offset+layerWidth+k], 
								 y[offset+layerWidth+k],
								 z[offset+layerWidth+k])));
		  e7->edgeNodes.push_back(
					  checkNode(
							new Node(nodeId++,
								 x[offset+layerWidth+k*(layerWidth)*(nLayers+1)], 
								 y[offset+layerWidth+k*(layerWidth)*(nLayers+1)],
								 z[offset+layerWidth+k*(layerWidth)*(nLayers+1)])));		     
		}
	      e0->curveType = pt;
	      e2->curveType = pt;
	      e4->curveType = pt;
	      e5->curveType = pt;
	      e6->curveType = pt;
	      e7->curveType = pt;

	      if (tetsOn != 1)  // If a mixed prism/tet mesh, add prism element to mesh.
		{
		  // Change the surface elements of the quad face on the symettry plane to match the layers of prisims.
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
		  m->element[m->expDim].push_back(elmt);
		}

	      else if (tetsOn == 1) // If only tet mesh, split each prim into 3 tetrahedrons and then push back tets to the mesh.
		{			  			  			 			  
		  // Get geom of the layered element.
		  SpatialDomains::PrismGeomSharedPtr geomLayer = 
		    boost::dynamic_pointer_cast<SpatialDomains::PrismGeom>(
									   elmt->GetGeom(m->spaceDim));
		  // Create local region.
		  LocalRegions::PrismExpSharedPtr qs = 
		    MemoryManager<LocalRegions::PrismExp>::AllocateSharedPtr(B0,B0,B0,geomLayer);
			   
		  // Get the coordiantes of the high order layered prism element.
		  Array<OneD, NekDouble> xs(layerWidth*layerWidth*layerWidth);
		  Array<OneD, NekDouble> ys(layerWidth*layerWidth*layerWidth);
		  Array<OneD, NekDouble> zs(layerWidth*layerWidth*layerWidth);
		  qs->GetCoords(xs,ys,zs);

		  vector<NodeSharedPtr> nodeList0(4);
		  vector<NodeSharedPtr> nodeList1(4);
		  vector<NodeSharedPtr> nodeList2(4);
                  m->splitMap[i].first = 0;
                  m->splitMap[i].second = 1;
		  if (m->splitMap[i].first == 0)
		    {
		      if (m->splitMap[i].second == -1) // Case 1
			{
			  // Assign prism's vertices to 3 tetrahebrons according to the triangular orientation.
			  nodeList0[0] = elmt->GetVertex(0);
			  nodeList0[1] = elmt->GetVertex(4);
			  nodeList0[2] = elmt->GetVertex(2);
			  nodeList0[3] = elmt->GetVertex(5);
				  
			  nodeList1[0] = elmt->GetVertex(0);
			  nodeList1[1] = elmt->GetVertex(4);
			  nodeList1[2] = elmt->GetVertex(1);
			  nodeList1[3] = elmt->GetVertex(2);
				  
			  nodeList2[0] = elmt->GetVertex(0);
			  nodeList2[1] = elmt->GetVertex(5);
			  nodeList2[2] = elmt->GetVertex(2);
			  nodeList2[3] = elmt->GetVertex(3);
				  
			  // Create elements.
			  ElmtConfig confTet(eTetrahedron, 1, false, false,false);
			  ElementSharedPtr elmtTet0 = GetElementFactory().
			    CreateInstance(eTetrahedron,confTet,nodeList0,tags);	
			  ElementSharedPtr elmtTet1 = GetElementFactory().
			    CreateInstance(eTetrahedron,confTet,nodeList1,tags);
			  ElementSharedPtr elmtTet2 = GetElementFactory().
			    CreateInstance(eTetrahedron,confTet,nodeList2,tags);
				  
			  // Add high order information to element 0 edges 0-4, 0-2, 0-5, 4-2 & 2-5
			  EdgeSharedPtr e00 = elmtTet0->GetEdge(0);
			  EdgeSharedPtr e01 = elmtTet0->GetEdge(1);
			  EdgeSharedPtr e02 = elmtTet0->GetEdge(2);
			  EdgeSharedPtr e03 = elmtTet0->GetEdge(3);
			  EdgeSharedPtr e05 = elmtTet0->GetEdge(5);				  
			  // Add high order information to element 1 edges 0-4, 0-1, 0-2, 4-1 & 4-2
			  EdgeSharedPtr e10 = elmtTet1->GetEdge(0);
			  EdgeSharedPtr e11 = elmtTet1->GetEdge(1);
			  EdgeSharedPtr e12 = elmtTet1->GetEdge(2);
			  EdgeSharedPtr e13 = elmtTet1->GetEdge(3);				  
			  EdgeSharedPtr e14 = elmtTet1->GetEdge(4);				  
			  // Add high order information to element 2 edges 0-5, 0-2, 5-2, 5-3 & 2-3
			  EdgeSharedPtr e20 = elmtTet2->GetEdge(0);
			  EdgeSharedPtr e21 = elmtTet2->GetEdge(1);
			  EdgeSharedPtr e23 = elmtTet2->GetEdge(3);
			  EdgeSharedPtr e24 = elmtTet2->GetEdge(4);
			  EdgeSharedPtr e25 = elmtTet2->GetEdge(5);
			  
			  for (int k = 1; k < layerWidth-1; ++k)
			    {
			      e00->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+k*(layerWidth)*(nLayers+1)], 
									      y[offset+k*(layerWidth)*(nLayers+1)],
									      z[offset+k*(layerWidth)*(nLayers+1)])));
			      e01->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[(layerWidth+1)*k],
									      ys[(layerWidth+1)*k],
									      zs[(layerWidth+1)*k])));
			      e02->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[layerWidth*(layerWidth+1)*k],
									      ys[layerWidth*(layerWidth+1)*k],
									      zs[layerWidth*(layerWidth+1)*k])));
			      e03->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[(layerWidth-1)+layerWidth*(layerWidth-1)*(k+1)],
									      ys[(layerWidth-1)+layerWidth*(layerWidth-1)*(k+1)],
									      zs[(layerWidth-1)+layerWidth*(layerWidth-1)*(k+1)])));
			      e05->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+layerWidth+k*(layerWidth)*(nLayers+1)+layerWidth-1], 
									      y[offset+layerWidth+k*(layerWidth)*(nLayers+1)+layerWidth-1],
									      z[offset+layerWidth+k*(layerWidth)*(nLayers+1)+layerWidth-1])));
			      //
			      e10->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+k*(layerWidth)*(nLayers+1)], 
									      y[offset+k*(layerWidth)*(nLayers+1)],
									      z[offset+k*(layerWidth)*(nLayers+1)])));
			      e11->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+k], 
									      y[offset+k],
									      z[offset+k])));
			      e12->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[(layerWidth+1)*k],
									      ys[(layerWidth+1)*k],
									      zs[(layerWidth+1)*k])));
			      e13->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+k*(layerWidth)*(nLayers+1)+layerWidth-1], 
									      y[offset+k*(layerWidth)*(nLayers+1)+layerWidth-1],
									      z[offset+k*(layerWidth)*(nLayers+1)+layerWidth-1])));				      
			      e14->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[(layerWidth-1)+layerWidth*(layerWidth-1)*(k+1)],
									      ys[(layerWidth-1)+layerWidth*(layerWidth-1)*(k+1)],
									      zs[(layerWidth-1)+layerWidth*(layerWidth-1)*(k+1)])));				      
			      //
			      e20->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[layerWidth*(layerWidth+1)*k],
									      ys[layerWidth*(layerWidth+1)*k],
									      zs[layerWidth*(layerWidth+1)*k])));
			      e21->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[(layerWidth+1)*k],
									      ys[(layerWidth+1)*k],
									      zs[(layerWidth+1)*k])));
			      e23->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+layerWidth+k*(layerWidth)*(nLayers+1)+layerWidth-1], 
									      y[offset+layerWidth+k*(layerWidth)*(nLayers+1)+layerWidth-1],
									      z[offset+layerWidth+k*(layerWidth)*(nLayers+1)+layerWidth-1])));
			      e24->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+layerWidth+k*(layerWidth)*(nLayers+1)], 
									      y[offset+layerWidth+k*(layerWidth)*(nLayers+1)],
									      z[offset+layerWidth+k*(layerWidth)*(nLayers+1)])));
			      e25->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+layerWidth+k], 
									      y[offset+layerWidth+k],
									      z[offset+layerWidth+k])));				      
			    }
			  reverse(e03->edgeNodes.begin(), e03->edgeNodes.end());
			  reverse(e13->edgeNodes.begin(), e13->edgeNodes.end());
			  reverse(e14->edgeNodes.begin(), e14->edgeNodes.end());
			  reverse(e24->edgeNodes.begin(), e24->edgeNodes.end());
			  reverse(e25->edgeNodes.begin(), e25->edgeNodes.end());				  
			  e00->curveType = pt;
			  e01->curveType = pt;
			  e02->curveType = pt;
			  e03->curveType = pt;
			  e05->curveType = pt;
			  e10->curveType = pt;
			  e11->curveType = pt;
			  e12->curveType = pt;
			  e13->curveType = pt;				  
			  e14->curveType = pt;
			  e20->curveType = pt;
			  e21->curveType = pt;
			  e23->curveType = pt;
			  e24->curveType = pt;
			  e25->curveType = pt;

			  // Push back tet elements to the mesh.
			  m->element[m->expDim].push_back(elmtTet0);
			  m->element[m->expDim].push_back(elmtTet1);
			  m->element[m->expDim].push_back(elmtTet2);
			}
		      else if (m->splitMap[i].second == 1) // Case 2
			{
			  // Assign prism's vertices to 3 tetrahebrons according to the triangular orientation.
			  nodeList0[0] = elmt->GetVertex(0);
			  nodeList0[1] = elmt->GetVertex(4);
			  nodeList0[2] = elmt->GetVertex(1);
			  nodeList0[3] = elmt->GetVertex(5);
				  
			  nodeList1[0] = elmt->GetVertex(0);
			  nodeList1[1] = elmt->GetVertex(5);
			  nodeList1[2] = elmt->GetVertex(1);
			  nodeList1[3] = elmt->GetVertex(2);
				  
			  nodeList2[0] = elmt->GetVertex(0);
			  nodeList2[1] = elmt->GetVertex(5);
			  nodeList2[2] = elmt->GetVertex(2);
			  nodeList2[3] = elmt->GetVertex(3);

			  // Create elements.
			  ElmtConfig confTet(eTetrahedron, 1, false, false,false);
			  ElementSharedPtr elmtTet0 = GetElementFactory().
			    CreateInstance(eTetrahedron,confTet,nodeList0,tags);	
			  ElementSharedPtr elmtTet1 = GetElementFactory().
			    CreateInstance(eTetrahedron,confTet,nodeList1,tags);
			  ElementSharedPtr elmtTet2 = GetElementFactory().
			    CreateInstance(eTetrahedron,confTet,nodeList2,tags);				  

			  // Add high order information to element 0 edges 0-4, 0-1, 0-5, 4-1 & 1-5
			  EdgeSharedPtr e00 = elmtTet0->GetEdge(0);
			  EdgeSharedPtr e01 = elmtTet0->GetEdge(1);
			  EdgeSharedPtr e02 = elmtTet0->GetEdge(2);
			  EdgeSharedPtr e03 = elmtTet0->GetEdge(3);
			  EdgeSharedPtr e05 = elmtTet0->GetEdge(5);
			  // Add high order information to element 1 edges 0-5, 0-1, 0-2, 5-1 & 5-2
			  EdgeSharedPtr e10 = elmtTet1->GetEdge(0);
			  EdgeSharedPtr e11 = elmtTet1->GetEdge(1);
			  EdgeSharedPtr e12 = elmtTet1->GetEdge(2);
			  EdgeSharedPtr e13 = elmtTet1->GetEdge(3);
			  EdgeSharedPtr e14 = elmtTet1->GetEdge(4);				  
			  // Add high order information to element 2 edges 0-5, 0-2, 5-2, 5-3 & 2-3
			  EdgeSharedPtr e20 = elmtTet2->GetEdge(0);
			  EdgeSharedPtr e21 = elmtTet2->GetEdge(1);
			  EdgeSharedPtr e23 = elmtTet2->GetEdge(3);
			  EdgeSharedPtr e24 = elmtTet2->GetEdge(4);
			  EdgeSharedPtr e25 = elmtTet2->GetEdge(5);
				  
			  for (int k = 1; k < layerWidth-1; ++k)
			    {
			      e00->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+k*(layerWidth)*(nLayers+1)], 
									      y[offset+k*(layerWidth)*(nLayers+1)],
									      z[offset+k*(layerWidth)*(nLayers+1)])));
			      e01->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+k], 
									      y[offset+k], z[offset+k])));
			      e02->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[layerWidth*(layerWidth+1)*k],
									      ys[layerWidth*(layerWidth+1)*k],
									      zs[layerWidth*(layerWidth+1)*k])));
			      e03->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+k*(layerWidth)*(nLayers+1)+layerWidth-1], 
									      y[offset+k*(layerWidth)*(nLayers+1)+layerWidth-1],
									      z[offset+k*(layerWidth)*(nLayers+1)+layerWidth-1])));
			      e05->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[(layerWidth-1)+layerWidth*(layerWidth+1)*k],
									      ys[(layerWidth-1)+layerWidth*(layerWidth+1)*k],
									      zs[(layerWidth-1)+layerWidth*(layerWidth+1)*k])));
			      //
			      e10->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[layerWidth*(layerWidth+1)*k],
									      ys[layerWidth*(layerWidth+1)*k],
									      zs[layerWidth*(layerWidth+1)*k])));
			      e11->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+k], 
									      y[offset+k], z[offset+k])));
			      e12->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[(layerWidth+1)*k], 
									      ys[(layerWidth+1)*k],
									      zs[(layerWidth+1)*k])));
			      e13->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[(layerWidth-1)+layerWidth*(layerWidth+1)*k],
									      ys[(layerWidth-1)+layerWidth*(layerWidth+1)*k],
									      zs[(layerWidth-1)+layerWidth*(layerWidth+1)*k])));
			      e14->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+layerWidth+k*(layerWidth)*(nLayers+1)+layerWidth-1], 
									      y[offset+layerWidth+k*(layerWidth)*(nLayers+1)+layerWidth-1],
									      z[offset+layerWidth+k*(layerWidth)*(nLayers+1)+layerWidth-1])));
			      //
			      e20->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[layerWidth*(layerWidth+1)*k],
									      ys[layerWidth*(layerWidth+1)*k],
									      zs[layerWidth*(layerWidth+1)*k])));
			      e21->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[(layerWidth+1)*k], 
									      ys[(layerWidth+1)*k],
									      zs[(layerWidth+1)*k])));
			      e23->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+layerWidth+k*(layerWidth)*(nLayers+1)+layerWidth-1], 
									      y[offset+layerWidth+k*(layerWidth)*(nLayers+1)+layerWidth-1],
									      z[offset+layerWidth+k*(layerWidth)*(nLayers+1)+layerWidth-1])));
			      e24->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+layerWidth+k*(layerWidth)*(nLayers+1)], 
									      y[offset+layerWidth+k*(layerWidth)*(nLayers+1)],
									      z[offset+layerWidth+k*(layerWidth)*(nLayers+1)])));
			      e25->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+layerWidth+k], 
									      y[offset+layerWidth+k], z[offset+layerWidth+k])));
			    }
			  reverse(e03->edgeNodes.begin(), e03->edgeNodes.end());
			  reverse(e13->edgeNodes.begin(), e13->edgeNodes.end());
			  reverse(e14->edgeNodes.begin(), e14->edgeNodes.end());
			  reverse(e24->edgeNodes.begin(), e24->edgeNodes.end());
			  reverse(e25->edgeNodes.begin(), e25->edgeNodes.end());
			  e00->curveType = pt;
			  e01->curveType = pt;
			  e02->curveType = pt;
			  e03->curveType = pt;
			  e05->curveType = pt;
			  e10->curveType = pt;
			  e11->curveType = pt;
			  e12->curveType = pt;
			  e13->curveType = pt;
			  e14->curveType = pt;
			  e20->curveType = pt;
			  e21->curveType = pt;
			  e23->curveType = pt;
			  e24->curveType = pt;
			  e25->curveType = pt;
				  
			  // Push back tet elements to the mesh.
			  m->element[m->expDim].push_back(elmtTet0);
			  m->element[m->expDim].push_back(elmtTet1);
			  m->element[m->expDim].push_back(elmtTet2);				  
			}
		    }

		  if (m->splitMap[i].first == 1)
		    {
		      if (m->splitMap[i].second == -1) // Case 3
			{
			  // Assign prism's vertices to 3 tetrahebrons according to the triangular orientation.
			  nodeList0[0] = elmt->GetVertex(1);
			  nodeList0[1] = elmt->GetVertex(0);
			  nodeList0[2] = elmt->GetVertex(5);
			  nodeList0[3] = elmt->GetVertex(3);
				  
			  nodeList1[0] = elmt->GetVertex(1);
			  nodeList1[1] = elmt->GetVertex(0);
			  nodeList1[2] = elmt->GetVertex(4);
			  nodeList1[3] = elmt->GetVertex(5);
				  
			  nodeList2[0] = elmt->GetVertex(1);
			  nodeList2[1] = elmt->GetVertex(3);
			  nodeList2[2] = elmt->GetVertex(5);
			  nodeList2[3] = elmt->GetVertex(2);
				  
			  // Create elements.
			  ElmtConfig confTet(eTetrahedron, 1, false, false,false);
			  ElementSharedPtr elmtTet0 = GetElementFactory().
			    CreateInstance(eTetrahedron,confTet,nodeList0,tags);	
			  ElementSharedPtr elmtTet1 = GetElementFactory().
			    CreateInstance(eTetrahedron,confTet,nodeList1,tags);
			  ElementSharedPtr elmtTet2 = GetElementFactory().
			    CreateInstance(eTetrahedron,confTet,nodeList2,tags);
				  
			  // Add high order information to element 0 edges 1-0, 1-5, 1-3, 0-5 & 5-3
			  EdgeSharedPtr e00 = elmtTet0->GetEdge(0);
			  EdgeSharedPtr e01 = elmtTet0->GetEdge(1);
			  EdgeSharedPtr e02 = elmtTet0->GetEdge(2);
			  EdgeSharedPtr e03 = elmtTet0->GetEdge(3);
			  EdgeSharedPtr e05 = elmtTet0->GetEdge(5);				  
			  // Add high order information to element 1 edges 1-0, 1-4, 1-5, 0-4 & 0-5
			  EdgeSharedPtr e10 = elmtTet1->GetEdge(0);
			  EdgeSharedPtr e11 = elmtTet1->GetEdge(1);
			  EdgeSharedPtr e12 = elmtTet1->GetEdge(2);
			  EdgeSharedPtr e13 = elmtTet1->GetEdge(3);				  
			  EdgeSharedPtr e14 = elmtTet1->GetEdge(4);				  
			  // Add high order information to element 2 edges 1-3, 1-5, 3-5, 3-2 & 5-2
			  EdgeSharedPtr e20 = elmtTet2->GetEdge(0);
			  EdgeSharedPtr e21 = elmtTet2->GetEdge(1);
			  EdgeSharedPtr e23 = elmtTet2->GetEdge(3);
			  EdgeSharedPtr e24 = elmtTet2->GetEdge(4);
			  EdgeSharedPtr e25 = elmtTet2->GetEdge(5);				  				  
				  				 				  
			  for (int k = 1; k < layerWidth-1; ++k)
			    {
			      e00->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+k], 
									      y[offset+k], z[offset+k])));
			      e01->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[(layerWidth-1)+layerWidth*(layerWidth+1)*k],
									      ys[(layerWidth-1)+layerWidth*(layerWidth+1)*k],
									      zs[(layerWidth-1)+layerWidth*(layerWidth+1)*k])));
			      e02->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[(layerWidth-1)*(k+1)],
									      ys[(layerWidth-1)*(k+1)],
									      zs[(layerWidth-1)*(k+1)])));
			      e03->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[layerWidth*(layerWidth+1)*k],
									      ys[layerWidth*(layerWidth+1)*k],
									      zs[layerWidth*(layerWidth+1)*k])));	
			      e05->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+layerWidth+k*(layerWidth)*(nLayers+1)], 
									      y[offset+layerWidth+k*(layerWidth)*(nLayers+1)],
									      z[offset+layerWidth+k*(layerWidth)*(nLayers+1)])));
				      
			      //
			      e10->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+k], 
									      y[offset+k], z[offset+k])));
			      e11->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+k*(layerWidth)*(nLayers+1)+layerWidth-1], 
									      y[offset+k*(layerWidth)*(nLayers+1)+layerWidth-1],
									      z[offset+k*(layerWidth)*(nLayers+1)+layerWidth-1])));
			      e12->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[(layerWidth-1)+layerWidth*(layerWidth+1)*k],
									      ys[(layerWidth-1)+layerWidth*(layerWidth+1)*k],
									      zs[(layerWidth-1)+layerWidth*(layerWidth+1)*k])));
			      e13->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+k*(layerWidth)*(nLayers+1)], 
									      y[offset+k*(layerWidth)*(nLayers+1)],
									      z[offset+k*(layerWidth)*(nLayers+1)])));				      
			      e14->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[layerWidth*(layerWidth+1)*k],
									      ys[layerWidth*(layerWidth+1)*k],
									      zs[layerWidth*(layerWidth+1)*k])));				      
			      //
			      e20->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[(layerWidth-1)+layerWidth*(layerWidth+1)*k],
									      ys[(layerWidth-1)+layerWidth*(layerWidth+1)*k],
									      zs[(layerWidth-1)+layerWidth*(layerWidth+1)*k])));
			      e21->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[(layerWidth-1)+layerWidth*(layerWidth+1)*k],
									      ys[(layerWidth-1)+layerWidth*(layerWidth+1)*k],
									      zs[(layerWidth-1)+layerWidth*(layerWidth+1)*k])));
			      e23->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+layerWidth+k*(layerWidth)*(nLayers+1)], 
									      y[offset+layerWidth+k*(layerWidth)*(nLayers+1)],
									      z[offset+layerWidth+k*(layerWidth)*(nLayers+1)])));
			      e24->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+layerWidth+k], 
									      y[offset+layerWidth+k],
									      z[offset+layerWidth+k])));
			      e25->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+layerWidth+k*(layerWidth)*(nLayers+1)+layerWidth-1], 
									      y[offset+layerWidth+k*(layerWidth)*(nLayers+1)+layerWidth-1],
									      z[offset+layerWidth+k*(layerWidth)*(nLayers+1)+layerWidth-1])));     
			    }
			  reverse(e00->edgeNodes.begin(), e00->edgeNodes.end());
			  reverse(e05->edgeNodes.begin(), e05->edgeNodes.end());
			  reverse(e10->edgeNodes.begin(), e10->edgeNodes.end());
			  reverse(e25->edgeNodes.begin(), e25->edgeNodes.end());				  
			  e00->curveType = pt;
			  e01->curveType = pt;
			  e02->curveType = pt;
			  e03->curveType = pt;
			  e05->curveType = pt;
			  e10->curveType = pt;
			  e11->curveType = pt;
			  e12->curveType = pt;
			  e13->curveType = pt;				  
			  e14->curveType = pt;
			  e20->curveType = pt;
			  e21->curveType = pt;
			  e23->curveType = pt;
			  e24->curveType = pt;
			  e25->curveType = pt;	 
				  
			  // Push back tet elements to the mesh.
			  m->element[m->expDim].push_back(elmtTet0);
			  m->element[m->expDim].push_back(elmtTet1);
			  m->element[m->expDim].push_back(elmtTet2);
			}
		      else if (m->splitMap[i].second == 1) // Case 4
			{
			  // Assign prism's vertices to 3 tetrahebrons according to the triangular orientation.
			  nodeList0[0] = elmt->GetVertex(1);
			  nodeList0[1] = elmt->GetVertex(0);
			  nodeList0[2] = elmt->GetVertex(4);
			  nodeList0[3] = elmt->GetVertex(3);
				  
			  nodeList1[0] = elmt->GetVertex(1);
			  nodeList1[1] = elmt->GetVertex(3);
			  nodeList1[2] = elmt->GetVertex(4);
			  nodeList1[3] = elmt->GetVertex(5);
				  
			  nodeList2[0] = elmt->GetVertex(1);
			  nodeList2[1] = elmt->GetVertex(3);
			  nodeList2[2] = elmt->GetVertex(5);
			  nodeList2[3] = elmt->GetVertex(2);
				  
			  //Add elements to mesh.
			  ElmtConfig confTet(eTetrahedron, 1, false, false,false);
			  ElementSharedPtr elmtTet0 = GetElementFactory().
			    CreateInstance(eTetrahedron,confTet,nodeList0,tags);
			  ElementSharedPtr elmtTet1 = GetElementFactory().
			    CreateInstance(eTetrahedron,confTet,nodeList1,tags);
			  ElementSharedPtr elmtTet2 = GetElementFactory().
			    CreateInstance(eTetrahedron,confTet,nodeList2,tags);	  

			  // Add high order information to element 0 edges 1-0, 1-4, 1-3, 0-4 & 4-3
			  EdgeSharedPtr e00 = elmtTet0->GetEdge(0);
			  EdgeSharedPtr e01 = elmtTet0->GetEdge(1);
			  EdgeSharedPtr e02 = elmtTet0->GetEdge(2);
			  EdgeSharedPtr e03 = elmtTet0->GetEdge(3);
			  EdgeSharedPtr e05 = elmtTet0->GetEdge(5);				  
			  // Add high order information to element 1 edges 1-3, 1-4, 1-5, 3-4 & 3-5
			  EdgeSharedPtr e10 = elmtTet1->GetEdge(0);
			  EdgeSharedPtr e11 = elmtTet1->GetEdge(1);
			  EdgeSharedPtr e12 = elmtTet1->GetEdge(2);
			  EdgeSharedPtr e13 = elmtTet1->GetEdge(3);				  
			  EdgeSharedPtr e14 = elmtTet1->GetEdge(4);				  
			  // Add high order information to element 2 edges 1-3, 1-5, 3-5, 3-2 & 5-2
			  EdgeSharedPtr e20 = elmtTet2->GetEdge(0);
			  EdgeSharedPtr e21 = elmtTet2->GetEdge(1);
			  EdgeSharedPtr e23 = elmtTet2->GetEdge(3);
			  EdgeSharedPtr e24 = elmtTet2->GetEdge(4);
			  EdgeSharedPtr e25 = elmtTet2->GetEdge(5);				  				  
				  				 				  
			  for (int k = 1; k < layerWidth-1; ++k)
			    {
			      e00->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+k], 
									      y[offset+k],
									      z[offset+k])));
			      e01->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+k*(layerWidth)*(nLayers+1)+layerWidth-1], 
									      y[offset+k*(layerWidth)*(nLayers+1)+layerWidth-1],
									      z[offset+k*(layerWidth)*(nLayers+1)+layerWidth-1])));
			      e02->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[(layerWidth-1)*(k+1)],
									      ys[(layerWidth-1)*(k+1)],
									      zs[(layerWidth-1)*(k+1)])));
			      e03->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+k*(layerWidth)*(nLayers+1)], 
									      y[offset+k*(layerWidth)*(nLayers+1)],
									      z[offset+k*(layerWidth)*(nLayers+1)])));	
			      e05->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[layerWidth*(layerWidth-1)*(k+1)],
									      ys[layerWidth*(layerWidth-1)*(k+1)],
									      zs[layerWidth*(layerWidth-1)*(k+1)])));
				      
			      //
			      e10->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[(layerWidth-1)+layerWidth*(layerWidth+1)*k],
									      ys[(layerWidth-1)+layerWidth*(layerWidth+1)*k],
									      zs[(layerWidth-1)+layerWidth*(layerWidth+1)*k])));
			      e11->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+k*(layerWidth)*(nLayers+1)+layerWidth-1], 
									      y[offset+k*(layerWidth)*(nLayers+1)+layerWidth-1],
									      z[offset+k*(layerWidth)*(nLayers+1)+layerWidth-1])));
			      e12->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[(layerWidth-1)+layerWidth*(layerWidth+1)*k],
									      ys[(layerWidth-1)+layerWidth*(layerWidth+1)*k],
									      zs[(layerWidth-1)+layerWidth*(layerWidth+1)*k])));
			      e13->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[layerWidth*(layerWidth-1)*(k+1)],
									      ys[layerWidth*(layerWidth-1)*(k+1)],
									      zs[layerWidth*(layerWidth-1)*(k+1)])));			      
			      e14->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+layerWidth+k*(layerWidth)*(nLayers+1)], 
									      y[offset+layerWidth+k*(layerWidth)*(nLayers+1)],
									      z[offset+layerWidth+k*(layerWidth)*(nLayers+1)])));				      
			      //
			      e20->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[(layerWidth-1)+layerWidth*(layerWidth+1)*k],
									      ys[(layerWidth-1)+layerWidth*(layerWidth+1)*k],
									      zs[(layerWidth-1)+layerWidth*(layerWidth+1)*k])));
			      e21->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[(layerWidth-1)+layerWidth*(layerWidth+1)*k],
									      ys[(layerWidth-1)+layerWidth*(layerWidth+1)*k],
									      zs[(layerWidth-1)+layerWidth*(layerWidth+1)*k])));
			      e23->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+layerWidth+k*(layerWidth)*(nLayers+1)], 
									      y[offset+layerWidth+k*(layerWidth)*(nLayers+1)],
									      z[offset+layerWidth+k*(layerWidth)*(nLayers+1)])));
			      e24->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+layerWidth+k], 
									      y[offset+layerWidth+k],
									      z[offset+layerWidth+k])));
			      e25->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+layerWidth+k*(layerWidth)*(nLayers+1)+layerWidth-1], 
									      y[offset+layerWidth+k*(layerWidth)*(nLayers+1)+layerWidth-1],
									      z[offset+layerWidth+k*(layerWidth)*(nLayers+1)+layerWidth-1])));     
			    }
			  reverse(e00->edgeNodes.begin(), e00->edgeNodes.end());
			  reverse(e05->edgeNodes.begin(), e05->edgeNodes.end());
			  reverse(e25->edgeNodes.begin(), e25->edgeNodes.end());				  
			  e00->curveType = pt;
			  e01->curveType = pt;
			  e02->curveType = pt;
			  e03->curveType = pt;
			  e05->curveType = pt;
			  e10->curveType = pt;
			  e11->curveType = pt;
			  e12->curveType = pt;
			  e13->curveType = pt;				  
			  e14->curveType = pt;
			  e20->curveType = pt;
			  e21->curveType = pt;
			  e23->curveType = pt;
			  e24->curveType = pt;
			  e25->curveType = pt;

			  // Push back tet elements to the mesh.
			  m->element[m->expDim].push_back(elmtTet0);
			  m->element[m->expDim].push_back(elmtTet1);
			  m->element[m->expDim].push_back(elmtTet2);
			}
		    }

		  if (m->splitMap[i].first == 2)
		    {
		      if (m->splitMap[i].second == -1) // Case 5
			{
			  // Assign prism's vertices to 3 tetrahebrons according to the triangular orientation.
			  nodeList0[0] = elmt->GetVertex(4);
			  nodeList0[1] = elmt->GetVertex(1);
			  nodeList0[2] = elmt->GetVertex(3);
			  nodeList0[3] = elmt->GetVertex(2);
				  
			  nodeList1[0] = elmt->GetVertex(4);
			  nodeList1[1] = elmt->GetVertex(1);
			  nodeList1[2] = elmt->GetVertex(0);
			  nodeList1[3] = elmt->GetVertex(3);
				  
			  nodeList2[0] = elmt->GetVertex(4);
			  nodeList2[1] = elmt->GetVertex(2);
			  nodeList2[2] = elmt->GetVertex(3);
			  nodeList2[3] = elmt->GetVertex(5);
				  				  
			  //Add elements to mesh.
			  ElmtConfig confTet(eTetrahedron, 1, false, false,false);
			  ElementSharedPtr elmtTet0 = GetElementFactory().
			    CreateInstance(eTetrahedron,confTet,nodeList0,tags);
			  ElementSharedPtr elmtTet1 = GetElementFactory().
			    CreateInstance(eTetrahedron,confTet,nodeList1,tags);
			  ElementSharedPtr elmtTet2 = GetElementFactory().
			    CreateInstance(eTetrahedron,confTet,nodeList2,tags);	  

			  // Add high order information to element 0 edges 4-1, 4-3, 4-2, 1-3 & 3-2
			  EdgeSharedPtr e00 = elmtTet0->GetEdge(0);
			  EdgeSharedPtr e01 = elmtTet0->GetEdge(1);
			  EdgeSharedPtr e02 = elmtTet0->GetEdge(2);
			  EdgeSharedPtr e03 = elmtTet0->GetEdge(3);
			  EdgeSharedPtr e05 = elmtTet0->GetEdge(5);				  
			  // Add high order information to element 1 edges 4-1, 4-0, 4-3, 1-0 & 1-3
			  EdgeSharedPtr e10 = elmtTet1->GetEdge(0);
			  EdgeSharedPtr e11 = elmtTet1->GetEdge(1);
			  EdgeSharedPtr e12 = elmtTet1->GetEdge(2);
			  EdgeSharedPtr e13 = elmtTet1->GetEdge(3);				  
			  EdgeSharedPtr e14 = elmtTet1->GetEdge(4);			  
			  // Add high order information to element 2 edges 4-2, 4-3, 2-3, 2-5 & 3-5
			  EdgeSharedPtr e20 = elmtTet2->GetEdge(0);
			  EdgeSharedPtr e21 = elmtTet2->GetEdge(1);
			  EdgeSharedPtr e23 = elmtTet2->GetEdge(3);
			  EdgeSharedPtr e24 = elmtTet2->GetEdge(4);
			  EdgeSharedPtr e25 = elmtTet2->GetEdge(5);				  				  
				  				 				  
			  for (int k = 1; k < layerWidth-1; ++k)
			    {
			      e00->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+k*(layerWidth)*(nLayers+1)+layerWidth-1], 
									      y[offset+k*(layerWidth)*(nLayers+1)+layerWidth-1],
									      z[offset+k*(layerWidth)*(nLayers+1)+layerWidth-1])));
			      e01->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[layerWidth*(layerWidth-1)*(k+1)],
									      ys[layerWidth*(layerWidth-1)*(k+1)],
									      zs[layerWidth*(layerWidth-1)*(k+1)]))); 
			      e02->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[(layerWidth-1)+layerWidth*(layerWidth-1)*(k+1)],
									      ys[(layerWidth-1)+layerWidth*(layerWidth-1)*(k+1)],
									      zs[(layerWidth-1)+layerWidth*(layerWidth-1)*(k+1)])));
			      e03->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[(layerWidth-1)*(k+1)],
									      ys[(layerWidth-1)*(k+1)],
									      zs[(layerWidth-1)*(k+1)])));	
			      e05->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+layerWidth+k], 
									      y[offset+layerWidth+k],
									      z[offset+layerWidth+k])));				      
			      //
			      e10->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+k*(layerWidth)*(nLayers+1)+layerWidth-1], 
									      y[offset+k*(layerWidth)*(nLayers+1)+layerWidth-1],
									      z[offset+k*(layerWidth)*(nLayers+1)+layerWidth-1])));
			      e11->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+k*(layerWidth)*(nLayers+1)], 
									      y[offset+k*(layerWidth)*(nLayers+1)],
									      z[offset+k*(layerWidth)*(nLayers+1)])));
			      e12->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[layerWidth*(layerWidth-1)*(k+1)],
									      ys[layerWidth*(layerWidth-1)*(k+1)],
									      zs[layerWidth*(layerWidth-1)*(k+1)])));
			      e13->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+k], 
									      y[offset+k],
									      z[offset+k])));			      
			      e14->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[(layerWidth-1)*(k+1)],
									      ys[(layerWidth-1)*(k+1)],
									      zs[(layerWidth-1)*(k+1)])));				      
			      //
			      e20->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[(layerWidth-1)+layerWidth*(layerWidth-1)*(k+1)],
									      ys[(layerWidth-1)+layerWidth*(layerWidth-1)*(k+1)],
									      zs[(layerWidth-1)+layerWidth*(layerWidth-1)*(k+1)])));
			      e21->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[layerWidth*(layerWidth-1)*(k+1)],
									      ys[layerWidth*(layerWidth-1)*(k+1)],
									      zs[layerWidth*(layerWidth-1)*(k+1)])));
			      e23->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+layerWidth+k], 
									      y[offset+layerWidth+k], 
									      z[offset+layerWidth+k])));
			      e24->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+layerWidth+k*(layerWidth)*(nLayers+1)+layerWidth-1], 
									      y[offset+layerWidth+k*(layerWidth)*(nLayers+1)+layerWidth-1],
									      z[offset+layerWidth+k*(layerWidth)*(nLayers+1)+layerWidth-1])));
			      e25->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+layerWidth+k*(layerWidth)*(nLayers+1)], 
									      y[offset+layerWidth+k*(layerWidth)*(nLayers+1)],
									      z[offset+layerWidth+k*(layerWidth)*(nLayers+1)]))); 
			    }
			  reverse(e00->edgeNodes.begin(), e00->edgeNodes.end());
			  reverse(e01->edgeNodes.begin(), e01->edgeNodes.end());
			  reverse(e02->edgeNodes.begin(), e02->edgeNodes.end());
			  reverse(e10->edgeNodes.begin(), e10->edgeNodes.end());
			  reverse(e11->edgeNodes.begin(), e11->edgeNodes.end());
			  reverse(e12->edgeNodes.begin(), e12->edgeNodes.end());
			  reverse(e13->edgeNodes.begin(), e13->edgeNodes.end());
			  reverse(e20->edgeNodes.begin(), e20->edgeNodes.end());
			  reverse(e21->edgeNodes.begin(), e21->edgeNodes.end());
			  reverse(e23->edgeNodes.begin(), e23->edgeNodes.end());
			  e00->curveType = pt;
			  e01->curveType = pt;
			  e02->curveType = pt;
			  e03->curveType = pt;
			  e05->curveType = pt;
			  e10->curveType = pt;
			  e11->curveType = pt;
			  e12->curveType = pt;
			  e13->curveType = pt;				  
			  e14->curveType = pt;
			  e20->curveType = pt;
			  e21->curveType = pt;
			  e23->curveType = pt;
			  e24->curveType = pt;
			  e25->curveType = pt;			 				   

			  // Add triangular surface elements on the symmetry plane.
			  if (bl0 != -1)
			    {
			      if (j == 0) // Remove the quad element (once per element).
				{			
				  // Remove element at this inde
				  m->element[m->expDim-1].erase(m->element[m->expDim-1].begin()+(bl0-offsetSurf));
				  offsetSurf++;
				}			      		    			    				  	
			      vector<NodeSharedPtr> tNodeList0(3);
			      vector<NodeSharedPtr> tNodeList1(3);
			      tNodeList0[0] = nodeList[0];
			      tNodeList0[1] = nodeList[1];
			      tNodeList0[2] = nodeList[3];
			      tNodeList1[0] = nodeList[2];
			      tNodeList1[1] = nodeList[3];
			      tNodeList1[2] = nodeList[1];
			      vector<int> tagBE;
			      tagBE =  m->element[m->expDim-1][bl0]->GetTagList();
			      ElmtConfig bconf(eTriangle, 1, false, false);
			      ElementSharedPtr boundaryElmt0 = GetElementFactory().
				CreateInstance(eTriangle,bconf,tNodeList0,tagBE);			     
			      ElementSharedPtr boundaryElmt1 = GetElementFactory().
				CreateInstance(eTriangle,bconf,tNodeList1,tagBE);
			      m->element[m->expDim-1].push_back(boundaryElmt0);
			      m->element[m->expDim-1].push_back(boundaryElmt1);			     
			    }	
				  
			  // Push back tet elements to the mesh.
			  m->element[m->expDim].push_back(elmtTet0);
			  m->element[m->expDim].push_back(elmtTet1);
			  m->element[m->expDim].push_back(elmtTet2);
			}
		      else if (m->splitMap[i].second == 1) // Case 6
			{
			  nodeList0[0] = elmt->GetVertex(4);
			  nodeList0[1] = elmt->GetVertex(1);
			  nodeList0[2] = elmt->GetVertex(0);
			  nodeList0[3] = elmt->GetVertex(2);
				  
			  nodeList1[0] = elmt->GetVertex(4);
			  nodeList1[1] = elmt->GetVertex(2);
			  nodeList1[2] = elmt->GetVertex(0);
			  nodeList1[3] = elmt->GetVertex(3);
				  
			  nodeList2[0] = elmt->GetVertex(4);
			  nodeList2[1] = elmt->GetVertex(2);
			  nodeList2[2] = elmt->GetVertex(3);
			  nodeList2[3] = elmt->GetVertex(5);

			  //Add elements to mesh.
			  ElmtConfig confTet(eTetrahedron, 1, false, false,false);
			  ElementSharedPtr elmtTet0 = GetElementFactory().
			    CreateInstance(eTetrahedron,confTet,nodeList0,tags);
			  ElementSharedPtr elmtTet1 = GetElementFactory().
			    CreateInstance(eTetrahedron,confTet,nodeList1,tags);
			  ElementSharedPtr elmtTet2 = GetElementFactory().
			    CreateInstance(eTetrahedron,confTet,nodeList2,tags);	  

			  // Add high order information to element 0 edges 4-1, 4-0, 4-2, 1-0 & 0-2
			  EdgeSharedPtr e00 = elmtTet0->GetEdge(0);
			  EdgeSharedPtr e01 = elmtTet0->GetEdge(1);
			  EdgeSharedPtr e02 = elmtTet0->GetEdge(2);
			  EdgeSharedPtr e03 = elmtTet0->GetEdge(3);
			  EdgeSharedPtr e05 = elmtTet0->GetEdge(5);				  
			  // Add high order information to element 1 edges 4-2, 4-0, 4-3, 2-0 & 2-3
			  EdgeSharedPtr e10 = elmtTet1->GetEdge(0);
			  EdgeSharedPtr e11 = elmtTet1->GetEdge(1);
			  EdgeSharedPtr e12 = elmtTet1->GetEdge(2);
			  EdgeSharedPtr e13 = elmtTet1->GetEdge(3);				  
			  EdgeSharedPtr e14 = elmtTet1->GetEdge(4);				  
			  // Add high order information to element 2 edges 4-2, 4-3, 2-3, 2-5 & 3-5
			  EdgeSharedPtr e20 = elmtTet2->GetEdge(0);
			  EdgeSharedPtr e21 = elmtTet2->GetEdge(1);
			  EdgeSharedPtr e23 = elmtTet2->GetEdge(3);
			  EdgeSharedPtr e24 = elmtTet2->GetEdge(4);
			  EdgeSharedPtr e25 = elmtTet2->GetEdge(5);				  				  
				  				 				  
			  for (int k = 1; k < layerWidth-1; ++k)
			    {
			      e00->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+k*(layerWidth)*(nLayers+1)+layerWidth-1], 
									      y[offset+k*(layerWidth)*(nLayers+1)+layerWidth-1],
									      z[offset+k*(layerWidth)*(nLayers+1)+layerWidth-1])));
			      e01->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+k*(layerWidth)*(nLayers+1)], 
									      y[offset+k*(layerWidth)*(nLayers+1)],
									      z[offset+k*(layerWidth)*(nLayers+1)]))); 
			      e02->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[(layerWidth-1)+layerWidth*(layerWidth-1)*(k+1)],
									      ys[(layerWidth-1)+layerWidth*(layerWidth-1)*(k+1)],
									      zs[(layerWidth-1)+layerWidth*(layerWidth-1)*(k+1)])));
			      e03->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+k], 
									      y[offset+k],
									      z[offset+k])));	
			      e05->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[(layerWidth+1)*k],
									      ys[(layerWidth+1)*k],
									      zs[(layerWidth+1)*k])));				      
			      //
			      e10->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[(layerWidth-1)+layerWidth*(layerWidth-1)*(k+1)],
									      ys[(layerWidth-1)+layerWidth*(layerWidth-1)*(k+1)],
									      zs[(layerWidth-1)+layerWidth*(layerWidth-1)*(k+1)])));
			      e11->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+k*(layerWidth)*(nLayers+1)], 
									      y[offset+k*(layerWidth)*(nLayers+1)],
									      z[offset+k*(layerWidth)*(nLayers+1)])));
			      e12->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[layerWidth*(layerWidth-1)*(k+1)],
									      ys[layerWidth*(layerWidth-1)*(k+1)],
									      zs[layerWidth*(layerWidth-1)*(k+1)])));
			      e13->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[(layerWidth+1)*k],
									      ys[(layerWidth+1)*k],
									      zs[(layerWidth+1)*k])));			      
			      e14->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+layerWidth+k], 
									      y[offset+layerWidth+k],
									      z[offset+layerWidth+k])));				      
			      //
			      e20->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[(layerWidth-1)+layerWidth*(layerWidth-1)*(k+1)],
									      ys[(layerWidth-1)+layerWidth*(layerWidth-1)*(k+1)],
									      zs[(layerWidth-1)+layerWidth*(layerWidth-1)*(k+1)])));
			      e21->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      xs[layerWidth*(layerWidth-1)*(k+1)],
									      ys[layerWidth*(layerWidth-1)*(k+1)],
									      zs[layerWidth*(layerWidth-1)*(k+1)])));
			      e23->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+layerWidth+k], 
									      y[offset+layerWidth+k], 
									      z[offset+layerWidth+k])));
			      e24->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+layerWidth+k*(layerWidth)*(nLayers+1)+layerWidth-1], 
									      y[offset+layerWidth+k*(layerWidth)*(nLayers+1)+layerWidth-1],
									      z[offset+layerWidth+k*(layerWidth)*(nLayers+1)+layerWidth-1])));
			      e25->edgeNodes.push_back(
						       checkNode(
								     new Node(nodeId++,
									      x[offset+layerWidth+k*(layerWidth)*(nLayers+1)], 
									      y[offset+layerWidth+k*(layerWidth)*(nLayers+1)],
									      z[offset+layerWidth+k*(layerWidth)*(nLayers+1)]))); 
			    }
			  reverse(e00->edgeNodes.begin(), e00->edgeNodes.end());
			  reverse(e01->edgeNodes.begin(), e01->edgeNodes.end());
			  reverse(e02->edgeNodes.begin(), e02->edgeNodes.end());
			  reverse(e03->edgeNodes.begin(), e03->edgeNodes.end());
			  reverse(e10->edgeNodes.begin(), e10->edgeNodes.end());
			  reverse(e11->edgeNodes.begin(), e11->edgeNodes.end());
			  reverse(e12->edgeNodes.begin(), e12->edgeNodes.end());
			  reverse(e13->edgeNodes.begin(), e13->edgeNodes.end());
			  reverse(e14->edgeNodes.begin(), e14->edgeNodes.end());
			  reverse(e20->edgeNodes.begin(), e20->edgeNodes.end());
			  reverse(e21->edgeNodes.begin(), e21->edgeNodes.end());
			  reverse(e23->edgeNodes.begin(), e23->edgeNodes.end());
			  e00->curveType = pt;
			  e01->curveType = pt;
			  e02->curveType = pt;
			  e03->curveType = pt;
			  e05->curveType = pt;
			  e10->curveType = pt;
			  e11->curveType = pt;
			  e12->curveType = pt;
			  e13->curveType = pt;				  
			  e14->curveType = pt;
			  e20->curveType = pt;
			  e21->curveType = pt;
			  e23->curveType = pt;
			  e24->curveType = pt;
			  e25->curveType = pt;			 				   	 
			  
			  // Push back tet elements to the mesh.
			  m->element[m->expDim].push_back(elmtTet0);
			  m->element[m->expDim].push_back(elmtTet1);
			  m->element[m->expDim].push_back(elmtTet2);
			}
		    }
		}		
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
