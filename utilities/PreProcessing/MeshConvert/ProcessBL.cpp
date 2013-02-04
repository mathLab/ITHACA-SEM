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
//  Description: Refine prismatic boundary layer elements.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
using namespace std;

#include "MeshElements.h"
#include "ProcessBL.h"

#include <LibUtilities/Foundations/BLPoints.h>
#include <LocalRegions/PrismExp.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

namespace Nektar
{
    namespace Utilities
    {
        ModuleKey ProcessBL::className = 
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eProcessModule, "bl"), ProcessBL::create,
                "Refines a prismatic boundary layer.");
        
        ProcessBL::ProcessBL(MeshSharedPtr m) : ProcessModule(m)
        {
	    // BL mesh configuration.
            config["layers"]     = ConfigOption(false, "0",       
                "Number of layers to refine.");
	    config["nopoints"]   = ConfigOption(false, "3",       
                "Number of points in high order elements.");
	    config["powercoeff"] = ConfigOption(false, "2",       
                "Initial power law coeficiant for layer spacing.");
	    // Physical parameters of the problem.
            config["Re"]         = ConfigOption(false, "11e6",
                "Reynolds number to adapt to.");
	    config["BLlength"]   = ConfigOption(false, "0.8059",       
                "The length of the cord.");
            config["yplus"]      = ConfigOption(false, "5",
                "y^+ value (i.e. height of first element).");
            config["delta"]      = ConfigOption(false, "0.04",
                "Hight of elements to refine (m).");
            config["rho"]        = ConfigOption(false, "1.205",
                "Density (kg/m^3).");
            config["mu"]         = ConfigOption(false, "1.82e-5",
                "Dynamic viscosity (kg/(m*s)).");
	    config["TetsOff"]    = ConfigOption(true,  "0",
                "Use this option to turn off splitting prisms into tetrahedra");
        }
      
        ProcessBL::~ProcessBL()
        {
          
        }
        
        void ProcessBL::Process()
        {
            if (m->verbose)
            {
                cout << "ProcessBL: Refining prismatic boundary layer..."
                     << endl;
            }

	    // Initialisation of parameters and settings.

            // Physical problem parameters.
	    double Reynolds   = config["Re"].        as<double>();
	    double cord       = config["BLlength"].  as<double>();
	    double y_plus     = config["yplus"].     as<double>();
	    double delta_int  = config["delta"].     as<double>();
	    // Physical constants for air (default at temperature of 293K
	    // (20C)).
	    double visc_mu    = config["mu"].        as<double>(); //kg\(m*s)
	    double dens_rho   = config["rho"].       as<double>(); //kg\m^3
	    
	    // Mesh configuration parameters; powercoefficient only works when
	    // iLayers is on, layers works when iLayers is off.
	    double rr         = config["powercoeff"].as<double>(); 
	    int    nLayers    = config["layers"].    as<int>   (); 
	    int    layerWidth = config["nopoints"].  as<int>   ();

	    // Derived parameters.
	    double U          = visc_mu*Reynolds/(cord*dens_rho);
	    double Cf         = pow(2*log10(Reynolds)-0.65,-2.3);
	    double tau_w      = Cf*0.5*dens_rho*pow(U,2);
	    double delta_y    = visc_mu*y_plus/(sqrt(tau_w/dens_rho)*dens_rho);
	    double delta_star = delta_y*2/delta_int;
            
            // Used only when iLayers is on.
	    int npoints = int(ceil(log10(2/delta_star)/log10(rr)))+2;
            // Used for when iLayers is off.
	    double rn = powf(2.0/delta_star,1/double(nLayers+1));

	    // Options - Can be automaticaly changed based on if "powercoeff" or
	    // "layers" is specified.
	    bool tetsOn = !config["TetsOff"].as<bool>();
            // Set to false for user defined no. of layers and true for user
            // defined scaling factor.
	    bool iLayers = true; 

	    if (nLayers != 0)
            {
		iLayers = false;
            }
	    
	    if (iLayers)
            {
                nLayers = npoints-1;
		rn = powf(2.0/delta_star,1/double(npoints-2));
            }
            
	    // Printouts - enable only in -v option??
	    cerr << "Reynolds no. : " << Reynolds << ",  y plus : " 
                 << y_plus << endl;
	    cerr << "U : " << U << "m/s,  Cf : " << Cf << ",  tau_w : " 
                 << tau_w << endl;
	    cerr << "delta : " << delta_y << "m" << endl;
	    cerr << "delta_star : " << delta_star << endl;
	    cerr << "Number of layers :  " << nLayers << endl;
	    cerr << "geometric factor :  " << rn << endl;

            // Create a duplicate of the element list.
            vector<ElementSharedPtr> el = m->element[m->expDim];
            int nodeId = m->vertexSet.size();
            
            // Erase all elements from the element list. Elements will be
            // re-added as they are split.
            m->element[m->expDim].clear();

            // Set up map which identifies edges (as pairs of vertex ids)
            // including their vertices to the offset/stride in the 3d array
            // of collapsed co-ordinates. Note that this map also includes the
            // diagonal edges along quadrilateral faces which will be used to
            // add high-order information to the split tetrahedra.
            map<pair<int,int>, pair<int,int> > edgeMap;
            map<pair<int,int>, pair<int,int> >::iterator it;
            
            // Standard prismatic edges (0->8)
            int nq = layerWidth;
            edgeMap[pair<int,int>(0,1)] = pair<int,int>(0,            1);
            edgeMap[pair<int,int>(3,2)] = pair<int,int>(nq*(nq-1),    1);
            edgeMap[pair<int,int>(0,3)] = pair<int,int>(0,            nq);
            edgeMap[pair<int,int>(1,2)] = pair<int,int>(nq-1,         nq);
	    edgeMap[pair<int,int>(4,5)] = pair<int,int>((nq-1)*nq*nq, nq);
            edgeMap[pair<int,int>(0,4)] = pair<int,int>(0,            nq*nq);
            edgeMap[pair<int,int>(1,4)] = pair<int,int>(nq-1,         nq*nq);
            edgeMap[pair<int,int>(2,5)] = pair<int,int>(nq*nq-1,      nq*nq);
            edgeMap[pair<int,int>(3,5)] = pair<int,int>(nq*(nq-1),    nq*nq);
            // Face 0 diagonals
            edgeMap[pair<int,int>(0,2)] = pair<int,int>(0,            nq+1);
            edgeMap[pair<int,int>(1,3)] = pair<int,int>(nq-1,         nq-1);
            // Face 2 diagonals
            edgeMap[pair<int,int>(1,5)] = pair<int,int>(nq-1,         nq*nq+nq);
            edgeMap[pair<int,int>(2,4)] = pair<int,int>(nq*nq-1,      nq*nq-nq);
            // Face 4 diagonals
            edgeMap[pair<int,int>(0,5)] = pair<int,int>(0,            nq*nq+nq);
            edgeMap[pair<int,int>(3,4)] = pair<int,int>(nq*(nq-1),    nq*nq-nq);

            // Prismatic node -> face map.
            int prismFaceNodes[5][4] = {
                {0,1,2,3},{0,1,4,-1},{1,2,5,4},{3,2,5,-1},{0,3,5,4}};
            
            // Default PointsType.
            LibUtilities::PointsType pt = LibUtilities::eGaussLobattoLegendre;
            
            // Pass delta_star spacing through to BLPoints.
            LibUtilities::BLPoints::delta_star = delta_star;
            
            // Iterate over list of elements of expansion dimension.
            for (int i = 0; i < el.size(); ++i)
            {
                if (el[i]->GetConf().e != ePrism)
                {
                    m->element[m->expDim].push_back(el[i]);		  
                    continue;
                }
	      
                // Find quadrilateral boundary faces if any
                std::map<int, int> bLink;
                for (int j = 0; j < 5; j += 2)
                {
                    int bl = el[i]->GetBoundaryLink(j);
                    if (bl != -1)
                    {
                        bLink[j] = bl;
                    }
                }
	       
                // Get elemental geometry object.
                SpatialDomains::PrismGeomSharedPtr geom = 
                    boost::dynamic_pointer_cast<SpatialDomains::PrismGeom>(
                        el[i]->GetGeom(m->spaceDim));
                
                // Create basis.
                LibUtilities::BasisKey B0(
                    LibUtilities::eModified_A, layerWidth,
                    LibUtilities::PointsKey(layerWidth,pt));
                LibUtilities::BasisKey B1(
                    LibUtilities::eModified_A, 2,
                    LibUtilities::PointsKey(
                        nLayers+1, LibUtilities::eBoundaryLayerPoints));
                LibUtilities::BasisKey B2(
                    LibUtilities::eModified_A, layerWidth,
                    LibUtilities::PointsKey(layerWidth,pt));
                
                // Create local region.
                LocalRegions::PrismExpSharedPtr q = 
                    MemoryManager<LocalRegions::PrismExp>::AllocateSharedPtr(
                        B0,B1,B2,geom);
                
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

                    // For the first layer use nodes of the original prism at
                    // the bottom.
                    if (j == 0)
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
                    // For the last layer use nodes of the original prism at the
                    // top.
                    else if (j == nLayers-1)
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
                    ElmtConfig conf(ePrism, 1, true, true, false);
                    ElementSharedPtr elmt = GetElementFactory().
                        CreateInstance(ePrism,conf,nodeList,tags); 

                    EdgeSharedPtr e0 = elmt->GetEdge(0);
                    EdgeSharedPtr e2 = elmt->GetEdge(2);
                    EdgeSharedPtr e4 = elmt->GetEdge(4);
                    EdgeSharedPtr e5 = elmt->GetEdge(5);
                    EdgeSharedPtr e6 = elmt->GetEdge(6);
                    EdgeSharedPtr e7 = elmt->GetEdge(7);

                    for (int k = 1; k < layerWidth-1; ++k)
                    {
                        e0->edgeNodes.push_back(
                            NodeSharedPtr(
                                new Node(nodeId++,
                                         x[offset+k], 
                                         y[offset+k],
                                         z[offset+k])));
                        e2->edgeNodes.push_back(
                            NodeSharedPtr(
                                new Node(nodeId++,
                                         x[offset+layerWidth+k], 
                                         y[offset+layerWidth+k],
                                         z[offset+layerWidth+k])));
                        e4->edgeNodes.push_back(
                            NodeSharedPtr(
                                new Node(nodeId++,
                                         x[offset+k*(layerWidth)*(nLayers+1)], 
                                         y[offset+k*(layerWidth)*(nLayers+1)],
                                         z[offset+k*(layerWidth)*(nLayers+1)])));
                        e5->edgeNodes.push_back(
                            NodeSharedPtr(
                                new Node(nodeId++,
                                         x[offset+k*(layerWidth)*(nLayers+1)+layerWidth-1], 
                                         y[offset+k*(layerWidth)*(nLayers+1)+layerWidth-1],
                                         z[offset+k*(layerWidth)*(nLayers+1)+layerWidth-1])));
                        e6->edgeNodes.push_back(
                            NodeSharedPtr(
                                new Node(nodeId++,
                                         x[offset+layerWidth+k*(layerWidth)*(nLayers+1)+layerWidth-1], 
                                         y[offset+layerWidth+k*(layerWidth)*(nLayers+1)+layerWidth-1],
                                         z[offset+layerWidth+k*(layerWidth)*(nLayers+1)+layerWidth-1])));
                        e7->edgeNodes.push_back(
                            NodeSharedPtr(
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

                    // Change the surface elements of the quad face on the
                    // symmetry plane to match the layers of prisms.
                    map<int,int>::iterator it;
                    for (it = bLink.begin(); it != bLink.end(); ++it)
                    {
                        int fid = it->first;
                        int bl  = it->second;
                        if (j == 0)
                        {
                            ElementSharedPtr e = m->element[m->expDim-1][bl];
                            e->SetVertex(0, nodeList[prismFaceNodes[fid][0]]);
                            e->SetVertex(1, nodeList[prismFaceNodes[fid][1]]);
                            e->SetVertex(2, nodeList[prismFaceNodes[fid][2]]);
                            e->SetVertex(3, nodeList[prismFaceNodes[fid][3]]);
                            elmt->SetBoundaryLink(fid,bl);
                        }
                        else
                        {
                            vector<NodeSharedPtr> qNodeList(4);
                            for (int k = 0; k < 4; ++k)
                            {
                                qNodeList[k] = nodeList[prismFaceNodes[fid][k]];
                            }
                            vector<int> tagBE;
                            tagBE = m->element[m->expDim-1][bl]->GetTagList();
                            ElmtConfig bconf(eQuadrilateral, 1, true, true, false);
                            ElementSharedPtr boundaryElmt = GetElementFactory().
                                CreateInstance(eQuadrilateral,bconf,qNodeList,tagBE);
                            elmt->SetBoundaryLink(fid,m->element[m->expDim-1].size());
                            m->element[m->expDim-1].push_back(boundaryElmt);
                        }
                    }
                    
                    m->element[m->expDim].push_back(elmt);
                }
            }

            /*
             * Split all element types into tetrahedra. This is based on the
             * algorithm found in:
             * 
             * "How to Subdivide Pyramids, Prisms and Hexahedra into
             * Tetrahedra", J. Dompierre et al.
             */
            if (tetsOn)
            {
                // Denotes a set of indices inside m->element[m->expDim-1]
                // which are to be removed. These are precisely the
                // quadrilateral boundary faces which will be replaced by two
                // triangular faces.
                set<int> toRemove;
                
                // Represents table 2 of paper; each row i represents a
                // rotation of the prism nodes such that vertex i is placed at
                // position 0.
                static int indir[6][6] = {
                    {0,1,2,3,4,5},
                    {1,2,0,4,5,3},
                    {2,0,1,5,3,4},
                    {3,5,4,0,2,1},
                    {4,3,5,1,0,2},
                    {5,4,3,2,1,0}
                };
                
                // Represents table 3 of paper; the first three rows are the
                // three tetrahedra if the first condition is met; the latter
                // three rows are the three tetrahedra if the second condition
                // is met.
                static int prismTet[6][4] = {
                    {0,1,2,5},
                    {0,1,5,4},
                    {0,4,5,3},
                    {0,1,2,4},
                    {0,4,2,5},
                    {0,4,5,3}
                };
                
                // Represents the order of tetrahedral edges (in Nektar++
                // ordering).
                static int tetEdges[6][2] = {
                    {0,1}, {1,2},
                    {0,2}, {0,3},
                    {1,3}, {2,3}};
                
                // A tetrahedron nodes -> faces map.
                static int tetFaceNodes[4][3] = {
                    {0,1,2},{0,1,3},{1,2,3},{0,2,3}};

                // Make a copy of the element list.
                el = m->element[m->expDim];
                m->element[m->expDim].clear();
                
                for (int i = 0; i < el.size(); ++i)
                {
                    if (el[i]->GetConf().e != ePrism)
                    {
                        m->element[m->expDim].push_back(el[i]);
                        continue;
                    }                    
                    
                    vector<NodeSharedPtr> nodeList(6);
                    
                    // Map Nektar++ ordering (vertices 0,1,2,3 are base quad)
                    // to paper ordering (vertices 0,1,2 are first triangular
                    // face).
                    int mapPrism[6] = {0,1,4,3,2,5};
                    for (int j = 0; j < 6; ++j)
                    {
                        nodeList[j] = el[i]->GetVertex(mapPrism[j]);
                    }
                    
                    // Determine minimum ID of the nodes in this prism.
                    int minElId = nodeList[0]->id;
                    int minId   = 0;
                    for (int j = 1; j < 6; ++j)
                    {
                        int curId = nodeList[j]->id;
                        if (curId < minElId)
                        {
                            minElId = curId;
                            minId   = j;
                        }
                    }
                    
                    int offset;
                    
                    // Split prism using paper criterion.
                    if (min(nodeList[indir[minId][1]]->id, nodeList[indir[minId][5]]->id) <
                        min(nodeList[indir[minId][2]]->id, nodeList[indir[minId][4]]->id))
                    {
                        offset = 0;
                    }
                    else if (min(nodeList[indir[minId][1]]->id, nodeList[indir[minId][5]]->id) >
                             min(nodeList[indir[minId][2]]->id, nodeList[indir[minId][4]]->id))
                    {
                        offset = 3;
                    }
                    else
                    {
                        cerr << "Connectivity issue with prism->tet splitting." << endl;
                        abort();
                    }

                    // Create local prismatic region so that co-ordinates of
                    // the mapped element can be read from.
                    SpatialDomains::PrismGeomSharedPtr geomLayer = 
                        boost::dynamic_pointer_cast<SpatialDomains::PrismGeom>(
                            el[i]->GetGeom(m->spaceDim));
                    LibUtilities::BasisKey B0(LibUtilities::eModified_A, layerWidth,
                                              LibUtilities::PointsKey(layerWidth,pt));
                    LocalRegions::PrismExpSharedPtr qs = 
                        MemoryManager<LocalRegions::PrismExp>::AllocateSharedPtr(
                            B0,B0,B0,geomLayer);
                    
                    // Get the coordiantes of the high order prismatic element.
                    Array<OneD, NekDouble> xs(layerWidth*layerWidth*layerWidth);
                    Array<OneD, NekDouble> ys(layerWidth*layerWidth*layerWidth);
                    Array<OneD, NekDouble> zs(layerWidth*layerWidth*layerWidth);
                    qs->GetCoords(xs,ys,zs);
                    
                    for (int j = 0; j < 3; ++j)
                    {
                        vector<NodeSharedPtr> tetNodes(4);
                        
                        for (int k = 0; k < 4; ++k)
                        {
                            tetNodes[k] = nodeList[indir[minId][prismTet[j+offset][k]]];
                        }
                        
                        // Add high order information to tetrahedral edges.
                        for (int k = 0; k < 6; ++k)
                        {
                            // Determine prismatic nodes which correspond with
                            // this edge. Apply prism map to this to get
                            // Nektar++ ordering (this works since as a
                            // permutation, prismMap^{-1} = prismMap).
                            int n1 = mapPrism[indir[minId][prismTet[j+offset][tetEdges[k][0]]]];
                            int n2 = mapPrism[indir[minId][prismTet[j+offset][tetEdges[k][1]]]];
                            
                            // Find offset/stride
                            it = edgeMap.find(pair<int,int>(n1,n2));
                            if (it == edgeMap.end())
                            {
                                it = edgeMap.find(pair<int,int>(n2,n1));
                                if (it == edgeMap.end())
                                {
                                    cerr << "Couldn't find prism edges " << n1 << " " << n2 << endl;
                                    abort();
                                }
                                // Extract vertices -- reverse order.
                                for (int l = layerWidth-2; l >= 1; --l)
                                {
                                    int pos = it->second.first + l*it->second.second;
                                    tetNodes.push_back(
                                        NodeSharedPtr(
                                            new Node(nodeId++, xs[pos], ys[pos], zs[pos])));
                                }
                            }
                            else
                            {
                                // Extract vertices -- forwards order.
                                for (int l = 1; l < layerWidth-1; ++l)
                                {
                                    int pos = it->second.first + l*it->second.second;
                                    tetNodes.push_back(
                                        NodeSharedPtr(
                                            new Node(nodeId++, xs[pos], ys[pos], zs[pos])));
                                }
                            }
                        }
                        
                        vector<int> tags = el[i]->GetTagList();
                        //tags.push_back(0);
                        
                        ElmtConfig conf(eTetrahedron, layerWidth-1, false, false);
                        ElementSharedPtr elmt = GetElementFactory().
                            CreateInstance(eTetrahedron,conf,tetNodes,tags);
                        
                        m->element[m->expDim].push_back(elmt);
                    }

                    // Now check to see if this one of the quadrilateral faces
                    // is associated with a boundary condition. If it is, we
                    // split the face into two triangles and mark the existing
                    // face for removal.
                    //
                    // Note that this algorithm has significant room for
                    // improvement and is likely one of the least optimal
                    // approachs - however implementation is simple.
                    for (int fid = 0; fid < 5; fid += 2)
                    {
                        int bl = el[i]->GetBoundaryLink(fid);
                        
                        if (bl == -1)
                        {
                            continue;
                        }
                        
                        vector<NodeSharedPtr> triNodeList(3);
                        vector<int>           faceNodes  (3);
                        vector<int>           tmp;
                        vector<int>           tagBE;
                        ElmtConfig            bconf(eTriangle, 1, true, true);
                        ElementSharedPtr      elmt;
                        
                        // Mark existing boundary face for removal.
                        toRemove.insert(bl);
                        tagBE =  m->element[m->expDim-1][bl]->GetTagList();
                        
                        // First loop over tets.
                        for (int j = 0; j < 3; ++j)
                        {
                            // Now loop over faces.
                            for (int k = 0; k < 4; ++k)
                            {
                                // Finally determine Nektar++ local node
                                // numbers for this face.
                                for (int l = 0; l < 3; ++l)
                                {
                                    faceNodes[l] = mapPrism[indir[minId][prismTet[j+offset][tetFaceNodes[k][l]]]];
                                }
                                
                                tmp = faceNodes;
                                sort(faceNodes.begin(), faceNodes.end());

                                // If this face matches a triple denoting a
                                // split quad face, add the face to the
                                // expansion list.
                                if ((fid == 0 && (
                                        (faceNodes[0] == 0 && faceNodes[1] == 1 && faceNodes[2] == 2)   ||
                                        (faceNodes[0] == 0 && faceNodes[1] == 2 && faceNodes[2] == 3)   ||
                                        (faceNodes[0] == 0 && faceNodes[1] == 1 && faceNodes[2] == 3)   ||
                                        (faceNodes[0] == 1 && faceNodes[1] == 2 && faceNodes[2] == 3))) ||
                                    (fid == 2 && (
                                        (faceNodes[0] == 1 && faceNodes[1] == 2 && faceNodes[2] == 5)   ||
                                        (faceNodes[0] == 1 && faceNodes[1] == 4 && faceNodes[2] == 5)   ||
                                        (faceNodes[0] == 1 && faceNodes[1] == 2 && faceNodes[2] == 4)   ||
                                        (faceNodes[0] == 2 && faceNodes[1] == 4 && faceNodes[2] == 5))) ||
                                    (fid == 4 && (
                                        (faceNodes[0] == 0 && faceNodes[1] == 3 && faceNodes[2] == 5) ||
                                        (faceNodes[0] == 0 && faceNodes[1] == 4 && faceNodes[2] == 5) ||
                                        (faceNodes[0] == 0 && faceNodes[1] == 3 && faceNodes[2] == 4) ||
                                        (faceNodes[0] == 3 && faceNodes[1] == 4 && faceNodes[2] == 5))))
                                {
                                    triNodeList[0] = nodeList[mapPrism[tmp[0]]];
                                    triNodeList[1] = nodeList[mapPrism[tmp[1]]];
                                    triNodeList[2] = nodeList[mapPrism[tmp[2]]];
                                    elmt           = GetElementFactory().
                                        CreateInstance(eTriangle,bconf,triNodeList,tagBE);
                                    m->element[m->expDim-1].push_back(elmt);
                                }
                            }
                        }
                    }
                }
                
                // Remove 2D elements.
                vector<ElementSharedPtr> tmp;
                for (int i = 0; i < m->element[m->expDim-1].size(); ++i)
                {
                    set<int>::iterator it = toRemove.find(i);
                    if (it == toRemove.end())
                    {
                        tmp.push_back(m->element[m->expDim-1][i]);
                    }
                }
                
                m->element[m->expDim-1] = tmp;
            }
            
            // Re-process mesh to eliminate duplicate vertices and edges.
            ProcessVertices();
            ProcessEdges();
            ProcessFaces();
            ProcessElements();
            ProcessComposites();
        }

        // Very expensive function - search for this node to see if it has
        // already been created. This should be replaced by an appropriate
        // data structure for nearest-neighbour searches (e.g. kd-tree).
        NodeSharedPtr ProcessBL::checkNode(Node *n)
        {
            int i;
            for (i = 0; i < createdNodes.size(); ++i)
            {
                NodeSharedPtr n1 = createdNodes[i];
                if (((n1->x-n->x)*(n1->x-n->x) + 
                     (n1->y-n->y)*(n1->y-n->y) + 
                     (n1->z-n->z)*(n1->z-n->z)) < 1.0e-14)
                {
                    delete n;
                    return n1;
                }
            }
            createdNodes.push_back(NodeSharedPtr(n));
            return createdNodes.back();
        }
    }
}
