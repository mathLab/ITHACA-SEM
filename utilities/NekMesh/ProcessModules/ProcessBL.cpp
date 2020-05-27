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
//  Description: Refine prismatic or quadrilateral boundary layer elements.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>

#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Foundations/BLPoints.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Interpreter/Interpreter.h>
#include <LocalRegions/HexExp.h>
#include <LocalRegions/PrismExp.h>
#include <LocalRegions/QuadExp.h>

#include "ProcessBL.h"
#include <NekMeshUtils/MeshElements/Element.h>

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{
ModuleKey ProcessBL::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "bl"), ProcessBL::create,
    "Refines a prismatic or quadrilateral boundary layer.");

int **helper2d(int lda, int arr[][2])
{
    int **ret = new int *[lda];
    for (int i = 0; i < lda; ++i)
    {
        ret[i]    = new int[2];
        ret[i][0] = arr[i][0];
        ret[i][1] = arr[i][1];
    }
    return ret;
}

int **helper2d(int lda, int arr[][4])
{
    int **ret = new int *[lda];
    for (int i = 0; i < lda; ++i)
    {
        ret[i]    = new int[4];
        ret[i][0] = arr[i][0];
        ret[i][1] = arr[i][1];
        ret[i][2] = arr[i][2];
        ret[i][3] = arr[i][3];
    }
    return ret;
}

int **helper2d(int lda, int arr[][3])
{
    int **ret = new int *[lda];
    for (int i = 0; i < lda; ++i)
    {
        ret[i]    = new int[3];
        ret[i][0] = arr[i][0];
        ret[i][1] = arr[i][1];
        ret[i][2] = arr[i][2];
    }
    return ret;
}

struct SplitMapHelper
{
    int size;
    int dir;
    int oppositeFace;
    int bfacesSize;
    int *bfaces;

    int nEdgeToSplit;
    int *edgesToSplit;
    int **edgeVert;

    int **conn;

    int nEdgeToCurve;
    int *edgesToCurve;
    int **intEdgeFace;
    int **extEdgeFace;
    int blpDir;
    int **gll;
};

ProcessBL::ProcessBL(MeshSharedPtr m) : ProcessModule(m)
{
    // BL mesh configuration.
    m_config["layers"] =
        ConfigOption(false, "2", "Number of layers to refine.");
    m_config["nq"] =
        ConfigOption(false, "5", "Number of points in high order elements.");
    m_config["surf"] =
        ConfigOption(false, "", "Tag identifying surface connected to prism.");
    m_config["r"] =
        ConfigOption(false, "2.0", "Ratio to use in geometry progression.");
}

ProcessBL::~ProcessBL()
{
}

void ProcessBL::Process()
{
    if (m_mesh->m_verbose)
    {
        cout << "ProcessBL: Refining boundary layer..." << endl;
    }
    int dim = m_mesh->m_expDim;
    switch (dim)
    {
        case 2:
            BoundaryLayer2D();
            break;

        case 3:
            BoundaryLayer3D();
            break;

        default:
            ASSERTL0(0, "Dimension not supported")
            break;
    }

    ProcessVertices();
    ProcessEdges();
    ProcessFaces();
    ProcessElements();
    ProcessComposites();
}

void ProcessBL::BoundaryLayer2D()
{
    // This implementation of 2D does not support bidirectional splitting

    int nodeId  = m_mesh->m_vertexSet.size();
    int nl      = m_config["layers"].as<int>();
    int nq      = m_config["nq"].    as<int>();

    // determine if geometric ratio is string or a constant.
    LibUtilities::Interpreter rEval;
    NekDouble r             =  1;
    int       rExprId       = -1;
    bool      ratioIsString = false;

    if (m_config["r"].isType<NekDouble>())
    {
        r = m_config["r"].as<NekDouble>();
    }
    else
    {
        std::string rstr = m_config["r"].as<string>();
        rExprId = rEval.DefineFunction("x y z", rstr);
        ratioIsString = true;
    }

    // Default PointsType.
    LibUtilities::PointsType pt = LibUtilities::eGaussLobattoLegendre;

    // Map which takes element ID to edge on surface. This enables
    // splitting to occur in either y-direction of the prism.
    map<int, int> splitEls;

    // edgeMap associates geometry edge IDs to the (nl+1) vertices which
    // are generated along that edge when a prism is split, and is used
    // to avoid generation of duplicate vertices. It is stored as an
    // unordered map for speed.
    std::unordered_map<int, vector<NodeSharedPtr> > edgeMap;

    string surf = m_config["surf"].as<string>();
    if (surf.size() > 0)
    {
        vector<unsigned int> surfs;
        ParseUtils::GenerateSeqVector(surf, surfs);
        sort(surfs.begin(), surfs.end());

        // If surface is defined, process list of elements to find those
        // that are connected to it.
        for (int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); ++i)
        {
            ElementSharedPtr el = m_mesh->m_element[m_mesh->m_expDim][i];
            int nSurf = el->GetEdgeCount();

            for (int j = 0; j < nSurf; ++j)
            {
                int bl = el->GetBoundaryLink(j);
                if (bl == -1)
                {
                    continue;
                }

                ElementSharedPtr bEl  = m_mesh->m_element[m_mesh->m_expDim-1][bl];
                vector<int>      tags = bEl->GetTagList();
                vector<int>      inter;

                sort(tags.begin(), tags.end());
                set_intersection(surfs.begin(), surfs.end(),
                                 tags .begin(), tags .end(),
                                 back_inserter(inter));
                ASSERTL0(inter.size() <= 1,
                         "Intersection of surfaces wrong");

                if (inter.size() == 1)
                {
                    if (el->GetConf().m_e != LibUtilities::eQuadrilateral)
                    {
                        cerr << "WARNING: Found non-quad element "
                             << "to split in surface " << surf
                             << "; ignoring" << endl;
                        continue;
                    }

                    if (splitEls.count(el->GetId()) > 0)
                    {
                        cerr << "WARNING: quad already found; "
                             << "ignoring" << endl;
                        continue;
                    }

                    splitEls[el->GetId()] = j;
                }
            }
        }
    }
    else
    {
        ASSERTL0(false, "Surface must be specified.");
    }

    if (splitEls.size() == 0)
    {
        cerr << "WARNING: No elements detected to split." << endl;
        return;
    }

    // Erase all elements from the element list. Elements will be
    // re-added as they are split.
    vector<ElementSharedPtr> el = m_mesh->m_element[m_mesh->m_expDim];
    m_mesh->m_element[m_mesh->m_expDim].clear();

    // Iterate over list of elements of expansion dimension.
    for (int i = 0; i < el.size(); ++i)
    {


        if (splitEls.count(el[i]->GetId()) == 0)
        {
            m_mesh->m_element[m_mesh->m_expDim].push_back(el[i]);
            continue;
        }

        // Find other boundary faces if any
        std::map<int, int> bLink;
        for (int j = 0; j < 4; j++)
        {
            int bl = el[i]->GetBoundaryLink(j);
            if ( (bl != -1) && ( j != splitEls[el[i]->GetId()]) )
            {
                bLink[j] = bl;
            }
        }

        // Get elemental geometry object.
        SpatialDomains::QuadGeomSharedPtr geom =
            std::dynamic_pointer_cast<SpatialDomains::QuadGeom>(
                el[i]->GetGeom(m_mesh->m_spaceDim));

        // Determine whether to use reverse points.
        // (if edges 1 or 2 are on the surface)
        LibUtilities::PointsType t = ( (splitEls[el[i]->GetId()]+1) %4) < 2 ?
            LibUtilities::eBoundaryLayerPoints :
            LibUtilities::eBoundaryLayerPointsRev;


        if(ratioIsString) // determine value of r base on geom
        {
            NekDouble x,y,z;
            NekDouble x1,y1,z1;
            int nverts = geom->GetNumVerts();

            x = y = z = 0.0;

            for(int i = 0; i < nverts; ++i)
            {
                geom->GetVertex(i)->GetCoords(x1,y1,z1);
                x += x1; y += y1; z += z1;
            }
            x /= (NekDouble) nverts;
            y /= (NekDouble) nverts;
            z /= (NekDouble) nverts;
            r = rEval.Evaluate(rExprId,x,y,z,0.0);
        }

        // Create basis.
        LibUtilities::BasisKey B0(
            LibUtilities::eModified_A, nq,
            LibUtilities::PointsKey(nq,pt));
        LibUtilities::BasisKey B1(
            LibUtilities::eModified_A, 2,
            LibUtilities::PointsKey(nl+1, t, r));

        // Create local region.
        LocalRegions::QuadExpSharedPtr q;
        if (splitEls[el[i]->GetId()] % 2)
        {
            q = MemoryManager<LocalRegions::QuadExp>::AllocateSharedPtr(
                B1,B0,geom);
        }
        else
        {
            q = MemoryManager<LocalRegions::QuadExp>::AllocateSharedPtr(
                B0,B1,geom);
        }

        // Grab co-ordinates.
        Array<OneD, NekDouble> x(nq*(nl+1), 0.0);
        Array<OneD, NekDouble> y(nq*(nl+1), 0.0);
        Array<OneD, NekDouble> z(nq*(nl+1), 0.0);
        q->GetCoords(x,y,z);

        vector<vector<NodeSharedPtr> > edgeNodes(2);

        // Loop over edges to be split.
        for (int j = 0; j < 2; ++j)
        {
            int locEdge = (splitEls[el[i]->GetId()]+1+2*j)%4;
            int edgeId  = el[i]->GetEdge(locEdge)->m_id;

            // Determine whether we have already generated vertices
            // along this edge.
            auto eIt = edgeMap.find(edgeId);

            if (eIt == edgeMap.end())
            {
                // If not then resize storage to hold new points.
                edgeNodes[j].resize(nl+1);

                // Re-use existing vertices at endpoints of edge to
                // avoid duplicating the existing vertices.
                edgeNodes[j][0]  = el[i]->GetVertex(locEdge);
                edgeNodes[j][nl] = el[i]->GetVertex((locEdge+1)%4);

                 // Variable geometric ratio
                if(ratioIsString)
                {
                    NekDouble x0,y0;
                    NekDouble x1,y1;
                    NekDouble xm,ym,zm=0.0;

                    // -> Find edge end and mid points
                    x0 = edgeNodes[j][0]->m_x;
                    y0 = edgeNodes[j][0]->m_y;

                    x1 = edgeNodes[j][nl]->m_x;
                    y1 = edgeNodes[j][nl]->m_y;

                    xm = 0.5*(x0+x1);
                    ym = 0.5*(y0+y1);

                    // evaluate r factor based on mid point value
                    NekDouble rnew;
                    rnew = rEval.Evaluate(rExprId,xm,ym,zm,0.0);

                    // Get basis with new r;
                    t =  (j==0) ? LibUtilities::eBoundaryLayerPoints :
                                 LibUtilities::eBoundaryLayerPointsRev;
                    LibUtilities::PointsKey Pkey(nl+1, t, rnew);
                    LibUtilities::PointsSharedPtr newP
                        = LibUtilities::PointsManager()[Pkey];

                    const Array<OneD, const NekDouble> z = newP->GetZ();

                    // Create new interior nodes based on this new blend
                    for (int k = 1; k < nl; ++k)
                    {
                        xm = 0.5*(1+z[k])*(x1-x0) + x0;
                        ym = 0.5*(1+z[k])*(y1-y0) + y0;
                        zm = 0.0;
                        edgeNodes[j][k] = NodeSharedPtr(
                                new Node(nodeId++, xm,ym,zm));
                    }
                }
                else
                {
                    // Create new interior nodes.
                    int pos = 0;
                    for (int k = 1; k < nl; ++k)
                    {
                        switch (locEdge)
                        {
                            case 0:
                                pos = k;
                                break;
                            case 1:
                                pos = nq -1 + k*nq;
                                break;
                            case 2:
                                pos = nq*(nl+1) -1 - k;
                                break;
                            case 3:
                                pos = nq*nl - k*nq;
                                break;
                            default:
                                NEKERROR(ErrorUtil::efatal,
                                         "Quad edge should be < 4.");
                            break;
                        }
                        edgeNodes[j][k] = NodeSharedPtr(
                            new Node(nodeId++, x[pos], y[pos], z[pos]));
                    }
                }

                // Store these edges in edgeMap.
                edgeMap[edgeId] = edgeNodes[j];
            }
            else
            {
                // Check orientation
                if (eIt->second[0] == el[i]->GetVertex(locEdge))
                {
                    // Same orientation: copy nodes
                    edgeNodes[j] = eIt->second;
                }
                else
                {
                    // Reversed orientation: copy in reversed order
                    edgeNodes[j].resize(nl+1);
                    for (int k = 0; k < nl+1; ++k)
                    {
                        edgeNodes[j][k] = eIt->second[nl-k];
                    }
                }

            }
        }

        // Create element layers.
        for (int j = 0; j < nl; ++j)
        {
            // Get corner vertices.
            vector<NodeSharedPtr> nodeList(4);
            switch (splitEls[el[i]->GetId()])
            {
                case 0:
                {
                    nodeList[0] = edgeNodes[1][nl-j  ];
                    nodeList[1] = edgeNodes[0][j  ];
                    nodeList[2] = edgeNodes[0][j+1];
                    nodeList[3] = edgeNodes[1][nl-j-1];
                    break;
                }
                case 1:
                {
                    nodeList[0] = edgeNodes[1][j];
                    nodeList[1] = edgeNodes[1][j+1];
                    nodeList[2] = edgeNodes[0][nl-j-1];
                    nodeList[3] = edgeNodes[0][nl-j];
                    break;
                }
                case 2:
                {
                    nodeList[0] = edgeNodes[0][nl-j  ];
                    nodeList[1] = edgeNodes[1][j  ];
                    nodeList[2] = edgeNodes[1][j+1];
                    nodeList[3] = edgeNodes[0][nl-j-1];
                    break;
                }
                case 3:
                {
                    nodeList[0] = edgeNodes[0][j];
                    nodeList[1] = edgeNodes[0][j+1];
                    nodeList[2] = edgeNodes[1][nl-j-1];
                    nodeList[3] = edgeNodes[1][nl-j];
                    break;
                }
            }
            // Create the element.
            ElmtConfig conf(LibUtilities::eQuadrilateral, 1, true, false, true);
            ElementSharedPtr elmt = GetElementFactory().
                CreateInstance(
              LibUtilities::eQuadrilateral,conf,nodeList,el[i]->GetTagList());

            // Add high order nodes to split edges.
            for (int l = 0; l < 2; ++l)
            {
                int locEdge = (splitEls[el[i]->GetId()]+2*l)%4;
                EdgeSharedPtr HOedge = elmt->GetEdge(
                    locEdge);
                int pos = 0;
                for (int k = 1; k < nq-1; ++k)
                {
                    switch (locEdge)
                    {
                        case 0:
                            pos = j*nq+ k;
                            break;
                        case 1:
                            pos = j+1 + k*(nl+1);
                            break;
                        case 2:
                            pos = (j+1)*nq + (nq-1) - k;
                            break;
                        case 3:
                            pos = (nl+1)*(nq-1) + j - k*(nl+1);
                            break;
                        default:
                            NEKERROR(ErrorUtil::efatal,
                                     "Quad edge should be < 4.");
                            break;
                    }
                    HOedge->m_edgeNodes.push_back(
                        NodeSharedPtr(
                            new Node(nodeId++,x[pos],y[pos],0.0)));
                }
                HOedge->m_curveType = pt;
            }

            // Change the elements on the boundary
            // to match the layers
            for (auto &it : bLink)
            {
                int eid = it.first;
                int bl  = it.second;

                if (j == 0)
                {
                    // For first layer reuse existing 2D element.
                    ElementSharedPtr e = m_mesh->m_element[m_mesh->m_expDim-1][bl];
                    for (int k = 0; k < 2; ++k)
                    {
                        e->SetVertex(
                            k, nodeList[(eid+k)%4]);
                    }
                }
                else
                {
                    // For all other layers create new element.
                    vector<NodeSharedPtr> qNodeList(2);
                    for (int k = 0; k < 2; ++k)
                    {
                        qNodeList[k] = nodeList[(eid+k)%4];
                    }
                    vector<int> tagBE;
                    tagBE = m_mesh->m_element[m_mesh->m_expDim-1][bl]->GetTagList();
                    ElmtConfig bconf(LibUtilities::eSegment,1,true,true,false);
                    ElementSharedPtr boundaryElmt = GetElementFactory().
                        CreateInstance(LibUtilities::eSegment,bconf,
                                       qNodeList,tagBE);
                    m_mesh->m_element[m_mesh->m_expDim-1].push_back(boundaryElmt);
                }
            }

            m_mesh->m_element[m_mesh->m_expDim].push_back(elmt);
        }
    }
}

void ProcessBL::BoundaryLayer3D()
{
    // A set containing all element types which are valid.
    set<LibUtilities::ShapeType> validElTypes;
    validElTypes.insert(LibUtilities::ePrism);
    validElTypes.insert(LibUtilities::eHexahedron);

    int nodeId = m_mesh->m_vertexSet.size();
    int nl     = m_config["layers"].as<int>();
    int nq     = m_config["nq"].as<int>();

    // determine if geometric ratio is string or a constant.
    LibUtilities::Interpreter rEval;
    NekDouble r        = 1;
    int rExprId        = -1;
    bool ratioIsString = false;

    if (m_config["r"].isType<NekDouble>())
    {
        r = m_config["r"].as<NekDouble>();
    }
    else
    {
        std::string rstr = m_config["r"].as<string>();
        rExprId = rEval.DefineFunction("x y z", rstr);
        ratioIsString = true;
    }

    // Prismatic node -> face map.
    int prismFaceNodes[5][4] = {
        {0, 1, 2, 3}, {0, 1, 4, -1}, {1, 2, 5, 4}, {3, 2, 5, -1}, {0, 3, 5, 4}};
    int hexFaceNodes[6][4] = {{0, 1, 2, 3}, {0, 1, 5, 4}, {1, 2, 6, 5},
                              {3, 2, 6, 7}, {0, 3, 7, 4}, {4, 5, 6, 7}};
    map<LibUtilities::ShapeType, int **> faceNodeMap;
    faceNodeMap[LibUtilities::ePrism]      = helper2d(5, prismFaceNodes);
    faceNodeMap[LibUtilities::eHexahedron] = helper2d(6, hexFaceNodes);

    // Default PointsType.
    LibUtilities::PointsKey ekey(nq, LibUtilities::eGaussLobattoLegendre);
    Array<OneD, NekDouble> gll;
    LibUtilities::PointsManager()[ekey]->GetPoints(gll);

    // Map which takes element ID to face on surface. This enables
    // splitting to occur in either y-direction of the prism.
    unordered_map<int, int> splitEls;
    unordered_map<int, int>::iterator sIt;

    // Set up maps which takes an edge (in nektar++ ordering) and return
    // their offset and stride in the 3d array of collapsed quadrature
    // points. Note that this map includes only the edges that are on
    // the triangular faces as the edges in the normal direction are
    // linear.
    map<LibUtilities::ShapeType, map<int, SplitMapHelper> > splitMap;

    ////////////////////////////////////////
    // HEX DIR X
    ////////////////////////////////////////

    SplitMapHelper splitHex0;
    int splitMapBFacesHex0[4]   = {1, 2, 3, 4};
    int splitedgehex0[4]        = {4, 5, 6, 7};
    int splitHex0EdgeVert[4][2] = {{0, 4}, {1, 5}, {2, 6}, {3, 7}};
    int splitMapConnHex0[8][2]  = {{0, 0}, {1, 0}, {2, 0}, {3, 0},
                                  {0, 1}, {1, 1}, {2, 1}, {3, 1}};
    int splitedgestocurvehex0[8] = {0, 1, 2, 3, 8, 9, 10, 11};
    int splithex0gll[8][3]       = {{-1, -1, -1}, {1, -1, -1}, {1, 1, -1},
                              {-1, 1, -1},  {-1, -1, 1}, {1, -1, 1},
                              {1, 1, 1},    {-1, 1, 1}};
    int splitHex0IntEdgeFace[4][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 4}};
    int splitHex0ExtEdgeFace[4][2] = {{8, 1}, {9, 2}, {10, 3}, {11, 4}};

    splitHex0.size            = 8;
    splitHex0.dir             = 0;
    splitHex0.oppositeFace    = 5;
    splitHex0.nEdgeToSplit    = 4;
    splitHex0.edgesToSplit    = splitedgehex0;
    splitHex0.edgeVert        = helper2d(4, splitHex0EdgeVert);
    splitHex0.conn            = helper2d(8, splitMapConnHex0);
    splitHex0.bfacesSize      = 4;
    splitHex0.bfaces          = splitMapBFacesHex0;
    splitHex0.nEdgeToCurve    = 8;
    splitHex0.edgesToCurve    = splitedgestocurvehex0;
    splitHex0.intEdgeFace     = helper2d(4, splitHex0IntEdgeFace);
    splitHex0.extEdgeFace     = helper2d(4, splitHex0ExtEdgeFace);
    splitHex0.blpDir          = 2;
    splitHex0.gll             = helper2d(8, splithex0gll);

    splitMap[LibUtilities::eHexahedron][0] = splitHex0;
    int splitMapConnHex0rev[8][2]  = {{0, 1}, {1, 1}, {2, 1}, {3, 1}, {0, 0}, {1, 0}, {2, 0}, {3, 0}};

    SplitMapHelper splitHex5;
    splitHex5.size                         = 8;
    splitHex5.dir                          = 1;
    splitHex5.oppositeFace                 = 0;
    splitHex5.nEdgeToSplit                 = 4;
    splitHex5.edgesToSplit                 = splitedgehex0;
    splitHex5.edgeVert                     = helper2d(4, splitHex0EdgeVert);
    splitHex5.conn                         = helper2d(8, splitMapConnHex0rev);
    splitHex5.bfacesSize                   = 4;
    splitHex5.bfaces                       = splitMapBFacesHex0;
    splitHex5.nEdgeToCurve                 = 8;
    splitHex5.edgesToCurve                 = splitedgestocurvehex0;
    splitHex5.intEdgeFace                  = helper2d(4, splitHex0ExtEdgeFace);
    splitHex5.extEdgeFace                  = helper2d(4, splitHex0IntEdgeFace);
    splitHex5.blpDir                       = 2;
    splitHex5.gll                          = helper2d(8, splithex0gll);
    splitMap[LibUtilities::eHexahedron][5] = splitHex5;

    ////////////////////////////////////////
    // HEX DIR Y
    ////////////////////////////////////////

    SplitMapHelper splitHex1;
    int splitMapBFacesHex1[4]   = {0, 2, 5, 4};
    int splitedgehex1[4]        = {11, 9, 1, 3};
    int splitHex1EdgeVert[4][2] = {{4, 7}, {5, 6}, {1, 2}, {0, 3}};
    int splitMapConnHex1[8][2]  = {{3, 0}, {2, 0}, {2, 1}, {3, 1},
                                  {0, 0}, {1, 0}, {1, 1}, {0, 1}};
    int splitedgestocurvehex1[8] = {4, 8, 5, 0, 7, 10, 6, 2};
    int splitHex1IntEdgeFace[4][2] = {{0, 0}, {4, 4}, {5, 2}, {8, 5}};
    int splitHex1ExtEdgeFace[4][2] = {{2, 0}, {6, 2}, {7, 4}, {10, 5}};

    splitHex1.size                         = 8;
    splitHex1.dir                          = 0;
    splitHex1.oppositeFace                 = 3;
    splitHex1.nEdgeToSplit                 = 4;
    splitHex1.edgesToSplit                 = splitedgehex1;
    splitHex1.edgeVert                     = helper2d(4, splitHex1EdgeVert);
    splitHex1.conn                         = helper2d(8, splitMapConnHex1);
    splitHex1.bfacesSize                   = 4;
    splitHex1.bfaces                       = splitMapBFacesHex1;
    splitHex1.nEdgeToCurve                 = 8;
    splitHex1.edgesToCurve                 = splitedgestocurvehex1;
    splitHex1.intEdgeFace                  = helper2d(4, splitHex1IntEdgeFace);
    splitHex1.extEdgeFace                  = helper2d(4, splitHex1ExtEdgeFace);
    splitHex1.blpDir                       = 1;
    splitHex1.gll                          = helper2d(8, splithex0gll);
    splitMap[LibUtilities::eHexahedron][1] = splitHex1;

    SplitMapHelper splitHex3;
    int splitMapConnHex1rev[8][2]  = {{3, 1}, {2, 1}, {2, 0}, {3, 0}, {0, 1},
                                      {1, 1}, {1, 0}, {0, 0}};

    splitHex3.size                         = 8;
    splitHex3.dir                          = 1;
    splitHex3.oppositeFace                 = 1;
    splitHex3.nEdgeToSplit                 = 4;
    splitHex3.edgesToSplit                 = splitedgehex1;
    splitHex3.edgeVert                     = helper2d(4, splitHex1EdgeVert);
    splitHex3.conn                         = helper2d(8, splitMapConnHex1rev);
    splitHex3.bfacesSize                   = 4;
    splitHex3.bfaces                       = splitMapBFacesHex1;
    splitHex3.nEdgeToCurve                 = 8;
    splitHex3.edgesToCurve                 = splitedgestocurvehex1;
    splitHex3.intEdgeFace                  = helper2d(4, splitHex1ExtEdgeFace);
    splitHex3.extEdgeFace                  = helper2d(4, splitHex1IntEdgeFace);
    splitHex3.blpDir                       = 1;
    splitHex3.gll                          = helper2d(8, splithex0gll);
    splitMap[LibUtilities::eHexahedron][3] = splitHex3;

    ////////////////////////////////////////
    // HEX DIR Z
    ////////////////////////////////////////

    SplitMapHelper splitHex4;
    int splitMapBFacesHex4[4]   = {0, 1, 5, 3};
    int splitedgehex4[4]        = {8, 0, 2, 10};
    int splitHex4EdgeVert[4][2] = {{4, 5}, {0, 1}, {3, 2}, {7, 6}};
    int splitMapConnHex4[8][2]  = {{1, 0}, {1, 1}, {2, 1}, {2, 0},
                                  {0, 0}, {0, 1}, {3, 1}, {3, 0}};
    int splitedgestocurvehex4[8] = {4, 11, 7, 3, 5, 9, 6, 1};
    int splitHex4IntEdgeFace[4][2] = {{3, 0}, {4, 1}, {7, 3}, {11, 5}};
    int splitHex4ExtEdgeFace[4][2] = {{1, 0}, {5, 1}, {6, 3}, {9, 5}};

    splitHex4.size                         = 8;
    splitHex4.dir                          = 0;
    splitHex4.oppositeFace                 = 2;
    splitHex4.nEdgeToSplit                 = 4;
    splitHex4.edgesToSplit                 = splitedgehex4;
    splitHex4.edgeVert                     = helper2d(4, splitHex4EdgeVert);
    splitHex4.conn                         = helper2d(8, splitMapConnHex4);
    splitHex4.bfacesSize                   = 4;
    splitHex4.bfaces                       = splitMapBFacesHex4;
    splitHex4.nEdgeToCurve                 = 8;
    splitHex4.edgesToCurve                 = splitedgestocurvehex4;
    splitHex4.intEdgeFace                  = helper2d(4, splitHex4IntEdgeFace);
    splitHex4.extEdgeFace                  = helper2d(4, splitHex4ExtEdgeFace);
    splitHex4.blpDir                       = 0;
    splitHex4.gll                          = helper2d(8, splithex0gll);
    splitMap[LibUtilities::eHexahedron][4] = splitHex4;

    SplitMapHelper splitHex2;
    int splitMapConnHex4rev[8][2]  = {{1, 1}, {1, 0}, {2, 0}, {2, 1},
                                      {0, 1}, {0, 0}, {3, 0}, {3, 1}};

    splitHex2.size                         = 8;
    splitHex2.dir                          = 1;
    splitHex2.oppositeFace                 = 4;
    splitHex2.nEdgeToSplit                 = 4;
    splitHex2.edgesToSplit                 = splitedgehex4;
    splitHex2.edgeVert                     = helper2d(4, splitHex4EdgeVert);
    splitHex2.conn                         = helper2d(8, splitMapConnHex4rev);
    splitHex2.bfacesSize                   = 4;
    splitHex2.bfaces                       = splitMapBFacesHex4;
    splitHex2.nEdgeToCurve                 = 8;
    splitHex2.edgesToCurve                 = splitedgestocurvehex4;
    splitHex2.intEdgeFace                  = helper2d(4, splitHex4ExtEdgeFace);
    splitHex2.extEdgeFace                  = helper2d(4, splitHex4IntEdgeFace);
    splitHex2.blpDir                       = 0;
    splitHex2.gll                          = helper2d(8, splithex0gll);
    splitMap[LibUtilities::eHexahedron][2] = splitHex2;

    ////////////////////////////////////////
    // PRISM DIR Y
    ////////////////////////////////////////

    SplitMapHelper splitprism1;
    int splitMapBFacesPrism1[3]   = {0, 2, 4};
    int splitedgeprism1[3]        = {3, 1, 8};
    int splitPrism1EdgeVert[3][2] = {{0, 3}, {1, 2}, {4, 5}};
    int splitMapConnPrism1[6][2]  = {{0, 0}, {1, 0}, {1, 1}, {0, 1},
                                     {2, 0}, {2, 1}};
    int splitedgestocurveprism1[6] = {0, 4, 5, 2, 6, 7};
    int splitprism1gll[6][3]       = {{-1, -1, -1}, {1, -1, -1}, {1, 1, -1},
                              {-1, 1, -1},  {-1, -1, 1}, {-1, 1, 1}};
    int splitPrism1IntEdgeFace[3][2] = {{0, 0}, {4, 4}, {5, 2}};
    int splitPrism1ExtEdgeFace[3][2] = {{2, 0}, {6, 2}, {7, 4}};

    splitprism1.size                  = 6;
    splitprism1.dir                   = 0;
    splitprism1.oppositeFace          = 3;
    splitprism1.nEdgeToSplit          = 3;
    splitprism1.edgesToSplit          = splitedgeprism1;
    splitprism1.edgeVert              = helper2d(3, splitPrism1EdgeVert);
    splitprism1.conn                  = helper2d(6, splitMapConnPrism1);
    splitprism1.bfacesSize            = 3;
    splitprism1.bfaces                = splitMapBFacesPrism1;
    splitprism1.nEdgeToCurve          = 6;
    splitprism1.edgesToCurve          = splitedgestocurveprism1;
    splitprism1.intEdgeFace           = helper2d(3, splitPrism1IntEdgeFace);
    splitprism1.extEdgeFace           = helper2d(3, splitPrism1ExtEdgeFace);
    splitprism1.blpDir                = 1;
    splitprism1.gll                   = helper2d(8, splitprism1gll);
    splitMap[LibUtilities::ePrism][1] = splitprism1;

    SplitMapHelper splitprism3;
    int splitMapConnPrism1rev[6][2]  = {{0, 1}, {1, 1}, {1, 0}, {0, 0},
                                        {2, 1}, {2, 0}};
    splitprism3.size                  = 6;
    splitprism3.dir                   = 1;
    splitprism3.oppositeFace          = 1;
    splitprism3.nEdgeToSplit          = 3;
    splitprism3.edgesToSplit          = splitedgeprism1;
    splitprism3.edgeVert              = helper2d(3, splitPrism1EdgeVert);
    splitprism3.conn                  = helper2d(6, splitMapConnPrism1rev);
    splitprism3.bfacesSize            = 3;
    splitprism3.bfaces                = splitMapBFacesPrism1;
    splitprism3.nEdgeToCurve          = 6;
    splitprism3.edgesToCurve          = splitedgestocurveprism1;
    splitprism3.intEdgeFace           = helper2d(3, splitPrism1ExtEdgeFace);
    splitprism3.extEdgeFace           = helper2d(3, splitPrism1IntEdgeFace);
    splitprism3.blpDir                = 1;
    splitprism3.gll                   = helper2d(8, splitprism1gll);
    splitMap[LibUtilities::ePrism][3] = splitprism3;

    // edgeMap associates geometry edge IDs to the (nl+1) vertices which are
    // generated along that edge when a prism is split, and is used to avoid
    // generation of duplicate vertices. It is stored as an unordered map for
    // speed.
    unordered_map<int, vector<NodeSharedPtr> > edgeMap;
    unordered_map<int, vector<NodeSharedPtr> >::iterator eIt;

    string surf = m_config["surf"].as<string>();
    if (surf.size() == 0)
    {
        cout << "no surfaces to split" << endl;
        return;
    }
    vector<unsigned int> surfs;
    ParseUtils::GenerateSeqVector(surf.c_str(), surfs);
    sort(surfs.begin(), surfs.end());

    // If surface is defined, process list of elements to find those
    // that are connected to it.
    for (int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); ++i)
    {
        ElementSharedPtr el = m_mesh->m_element[m_mesh->m_expDim][i];
        int nSurf           = el->GetFaceCount();

        for (int j = 0; j < nSurf; ++j)
        {
            int bl = el->GetBoundaryLink(j);
            if (bl == -1)
            {
                continue;
            }

            ElementSharedPtr bEl = m_mesh->m_element[m_mesh->m_expDim - 1][bl];
            vector<int> tags     = bEl->GetTagList();
            vector<int> inter;

            sort(tags.begin(), tags.end());
            set_intersection(surfs.begin(), surfs.end(), tags.begin(),
                             tags.end(), back_inserter(inter));
            ASSERTL0(inter.size() <= 1, "Intersection of surfaces wrong");

            if (inter.size() == 1)
            {
                if (el->GetConf().m_e == LibUtilities::eHexahedron)
                {
                    map<int, SplitMapHelper>::iterator f =
                        splitMap[LibUtilities::eHexahedron].find(j);
                    if (f == splitMap[LibUtilities::eHexahedron].end())
                    {
                        cout << "WARNING: hex split on face " << j
                             << " unsupported"  << endl;
                        continue;
                    }

                    if (splitEls.count(el->GetId()) > 0)
                    {
                        cerr << "WARNING: hex already found; "
                             << "ignoring" << endl;
                    }

                    splitEls[el->GetId()] = j;
                }
                else if (el->GetConf().m_e == LibUtilities::ePrism)
                {
                    map<int, SplitMapHelper>::iterator f =
                        splitMap[LibUtilities::ePrism].find(j);
                    if (f == splitMap[LibUtilities::ePrism].end())
                    {
                        cout << "WARNING: Prism split on face " << j
                             << " unsupported" << endl;
                        continue;
                    }

                    if (splitEls.count(el->GetId()) > 0)
                    {
                        cerr << "WARNING: prism already found; "
                             << "ignoring" << endl;
                    }

                    splitEls[el->GetId()] = j;
                }
                else if (validElTypes.count(el->GetConf().m_e) == 0)
                {
                    cerr << "WARNING: Unsupported element type "
                         << "found in surface " << j << "; "
                         << "ignoring" << endl;
                    continue;
                }
            }
        }
    }

    if (splitEls.size() == 0)
    {
        cerr << "WARNING: No elements detected to split." << endl;
        return;
    }

    // Erase all elements from the element list. Elements will be
    // re-added as they are split.
    vector<ElementSharedPtr> el = m_mesh->m_element[m_mesh->m_expDim];
    m_mesh->m_element[m_mesh->m_expDim].clear();

    map<int, SpatialDomains::Geometry3DSharedPtr> geomMap;
    map<int, SpatialDomains::GeometrySharedPtr> edgeGeomMap;
    for (int i = 0; i < el.size(); ++i)
    {
        const int elId = el[i]->GetId();
        sIt            = splitEls.find(elId);
        if (sIt == splitEls.end())
        {
            continue;
        }

        // Get elemental geometry object and put into map.
        geomMap[elId] = dynamic_pointer_cast<SpatialDomains::Geometry3D>(
            el[i]->GetGeom(m_mesh->m_spaceDim));

        // Get all edge geometry too for evaluations.
        for(int j = 0; j < el[i]->GetEdgeCount(); j++)
        {
            EdgeSharedPtr e = el[i]->GetEdge(j);
            auto f = edgeGeomMap.find(e->m_id);
            if(f == edgeGeomMap.end())
            {
                edgeGeomMap[e->m_id] = e->GetGeom(m_mesh->m_spaceDim);
            }
        }
    }

    // node id to element id to para
    map<int, map<int, NekDouble> > paraElm; //parametric WRT to element
    map<int, map<int, NekDouble> > paraEdg; //parametric WRT edge

    map<int, vector<ElementSharedPtr> > elToStack;

    // Iterate over list of elements of expansion dimension.
    for (int i = 0; i < el.size(); ++i)
    {
        const int elId = el[i]->GetId();
        sIt            = splitEls.find(elId);

        // Any elements we don't process are ignored.
        if (sIt == splitEls.end())
        {
            m_mesh->m_element[m_mesh->m_expDim].push_back(el[i]);
            continue;
        }

        const int faceNum              = sIt->second;
        LibUtilities::ShapeType elType = el[i]->GetConf().m_e;

        SplitMapHelper &sMap = splitMap[elType][faceNum];

        // Find boundary faces if any
        std::map<int, int> bLink;
        for (int j = 0; j < sMap.bfacesSize; ++j)
        {
            int bl = el[i]->GetBoundaryLink(sMap.bfaces[j]);
            if (bl != -1)
            {
                bLink[sMap.bfaces[j]] = bl;
            }
        }

        // Determine whether to use reverse points.
        int nSplitEdge = sMap.nEdgeToSplit;
        vector<vector<NodeSharedPtr> > edgeNodes(nSplitEdge);

        // Loop over edges to be split.
        for (int j = 0; j < nSplitEdge; ++j)
        {
            int locEdge       = sMap.edgesToSplit[j];
            EdgeSharedPtr edg = el[i]->GetEdge(locEdge);
            int edgeId        = edg->m_id;
            SpatialDomains::GeometrySharedPtr geom = edgeGeomMap[edg->m_id];
            geom->FillGeom();

            // Determine value of r based on geometry.
            if (ratioIsString)
            {
                NekDouble x = 0.0, y = 0.0, z = 0.0, x1, y1, z1;
                int nverts = geom->GetNumVerts();

                for (int i = 0; i < nverts; ++i)
                {
                    geom->GetVertex(i)->GetCoords(x1, y1, z1);
                    x += x1;
                    y += y1;
                    z += z1;
                }
                x /= (NekDouble)nverts;
                y /= (NekDouble)nverts;
                z /= (NekDouble)nverts;
                r = rEval.Evaluate(rExprId, x, y, z, 0.0);
            }

            // Grab the boundary layer points distributions.
            LibUtilities::PointsKey bkey(
                nl + 1, LibUtilities::eBoundaryLayerPoints, r);
            Array<OneD, NekDouble> blp;
            LibUtilities::PointsManager()[bkey]->GetPoints(blp);

            LibUtilities::PointsKey brkey(
                nl + 1, LibUtilities::eBoundaryLayerPointsRev, r);
            Array<OneD, NekDouble> blpr;
            LibUtilities::PointsManager()[brkey]->GetPoints(blpr);

            // Determine whether we have already generated vertices along this
            // edge.
            eIt = edgeMap.find(edgeId);

            if (eIt == edgeMap.end())
            {
                // If not then resize storage to hold new points.
                edgeNodes[j].resize(nl + 1);

                // Re-use existing vertices at endpoints of edge to
                // avoid duplicating the existing vertices.
                edgeNodes[j][0]  = el[i]->GetVertex(sMap.edgeVert[j][0]);
                edgeNodes[j][nl] = el[i]->GetVertex(sMap.edgeVert[j][1]);

                // Evaluate nodal positions using the xmap.
                StdRegions::StdExpansionSharedPtr xmap = geom->GetXmap();
                Array<OneD, NekDouble> coeffs0 = geom->GetCoeffs(0);
                Array<OneD, NekDouble> coeffs1 = geom->GetCoeffs(1);
                Array<OneD, NekDouble> coeffs2 = geom->GetCoeffs(2);
                Array<OneD, NekDouble> xc(xmap->GetTotPoints());
                Array<OneD, NekDouble> yc(xmap->GetTotPoints());
                Array<OneD, NekDouble> zc(xmap->GetTotPoints());
                xmap->BwdTrans(coeffs0, xc);
                xmap->BwdTrans(coeffs1, yc);
                xmap->BwdTrans(coeffs2, zc);

                Array<OneD,NekDouble> pts;

                NekDouble tb = -1.0;
                NekDouble te = 1.0;
                if (sMap.dir != 0)
                {
                    swap(edgeNodes[j][0], edgeNodes[j][nl]);
                    swap(tb, te);
                }

                // Store parametric positions for the nodes -- i.e. \xi in
                // [-1,1].
                paraElm[edgeNodes[j][0]->m_id][el[i]->GetId()]  = tb;
                paraElm[edgeNodes[j][nl]->m_id][el[i]->GetId()] = te;

                // edgeNodes[j][0] is guaranteed to be the base of the spltting
                // direction; determine whether it is the end of the node
                bool forward = edgeNodes[j][0] == edg->m_n1;
                NekDouble tb2 = -1.0;
                NekDouble te2 = 1.0;
                if(!forward)
                {
                    swap(tb2,te2);
                }
                paraEdg[edgeNodes[j][0]->m_id][el[i]->GetId()]  = tb2;
                paraEdg[edgeNodes[j][nl]->m_id][el[i]->GetId()] = te2;

                // Create new interior nodes.
                for (int k = 1; k < nl; ++k)
                {
                    NekDouble tElm = tb * (1.0 - blp[k]) / 2.0 +
                        te * (1.0 + blp[k]) / 2.0;
                    paraElm[nodeId][el[i]->GetId()] = tElm;

                    NekDouble tEdge = tb2 * (1.0 - blp[k]) / 2.0 +
                        te2 * (1.0 + blp[k]) / 2.0;
                    paraEdg[nodeId][el[i]->GetId()] = tEdge;

                    Array<OneD, NekDouble> xp(1);

                    xp[0] = tEdge;

                    Array<OneD, NekDouble> loc(3);
                    loc[0]          = xmap->PhysEvaluate(xp, xc);
                    loc[1]          = xmap->PhysEvaluate(xp, yc);
                    loc[2]          = xmap->PhysEvaluate(xp, zc);
                    edgeNodes[j][k] = NodeSharedPtr(
                        new Node(nodeId++, loc[0], loc[1], loc[2]));

                    // If this edge is attached to some CAD, then perform a
                    // reverse projection to determine parametrisation on the
                    // CAD curve/surface.
                    if (edg->m_parentCAD)
                    {
                        if (edg->m_parentCAD->GetType() == CADType::eCurve)
                        {
                            CADCurveSharedPtr c = std::dynamic_pointer_cast<
                                CADCurve>(edg->m_parentCAD);
                            NekDouble t;
                            c->loct(loc, t);
                            edgeNodes[j][k]->SetCADCurve(c, t);
                        }
                        else if (edg->m_parentCAD->GetType() == CADType::eSurf)
                        {
                            CADSurfSharedPtr s = std::dynamic_pointer_cast<
                                CADSurf>(edg->m_parentCAD);
                            Array<OneD, NekDouble> uv = s->locuv(loc);
                            edgeNodes[j][k]->SetCADSurf(s, uv);
                        }
                    }
                }

                // Store these edges in edgeMap.
                edgeMap[edgeId] = edgeNodes[j];
            }
            else
            {
                edgeNodes[j] = eIt->second;
                NekDouble tb = -1.0;
                NekDouble te = 1.0;
                if (sMap.dir == 1)
                {
                    swap(tb, te);
                }
                for (int k = 0; k < eIt->second.size(); ++k)
                {
                    NekDouble tElm = tb * (1.0 - blp[k]) / 2.0 +
                        te * (1.0 + blp[k]) / 2.0;
                    paraElm[eIt->second[k]->m_id][el[i]->GetId()] = tElm;
                }
            }
        }

        // Create element layers.
        for (int j = 0; j < nl; ++j)
        {
            // Get corner vertices.
            vector<NodeSharedPtr> nodeList(sMap.size);
            for (int k = 0; k < sMap.size; ++k)
            {
                nodeList[k] = edgeNodes[sMap.conn[k][0]][j + sMap.conn[k][1]];
            }

            // Create the element.
            ElmtConfig conf(elType, 1, false, false, false);
            ElementSharedPtr elmt = GetElementFactory().CreateInstance(
                elType, conf, nodeList, el[i]->GetTagList());
            elmt->m_parentCAD = el[i]->m_parentCAD;

            // Always copy CAD parency from any side faces that are split.
            for (int k = 0; k < sMap.bfacesSize; ++k)
            {
                elmt->GetFace(sMap.bfaces[k])->m_parentCAD =
                    el[i]->GetFace(sMap.bfaces[k])->m_parentCAD;
            }

            // Copy face CAD from interior surface.
            if (j == 0)
            {
                elmt->GetFace(faceNum)->m_parentCAD =
                    el[i]->GetFace(faceNum)->m_parentCAD;
            }

            // Copy face CAD from exterior surface.
            if (j == nl - 1)
            {
                elmt->GetFace(sMap.oppositeFace)->m_parentCAD =
                    el[i]->GetFace(sMap.oppositeFace)->m_parentCAD;
            }

            // Always copy CAD parency from any side edges that are split.
            for (int k = 0; k < sMap.nEdgeToSplit; ++k)
            {
                elmt->GetEdge(sMap.edgesToSplit[k])->m_parentCAD =
                    el[i]->GetEdge(sMap.edgesToSplit[k])->m_parentCAD;
            }

            // For the first layer, copy CAD from interior curves.
            if (j == 0)
            {
                for (int k = 0; k < sMap.nEdgeToCurve / 2; ++k)
                {
                    elmt->GetEdge(sMap.intEdgeFace[k][0])->m_parentCAD =
                        el[i]->GetEdge(sMap.intEdgeFace[k][0])->m_parentCAD;
                }
            }
            // For the other layers, copy CAD of side surfaces.
            else
            {
                for (int k = 0; k < sMap.nEdgeToCurve / 2; ++k)
                {
                    elmt->GetEdge(sMap.intEdgeFace[k][0])->m_parentCAD =
                        el[i]->GetFace(sMap.intEdgeFace[k][1])->m_parentCAD;
                }
            }

            // For the last layer, copy CAD from exterior curves.
            if (j == nl - 1)
            {
                for (int k = 0; k < sMap.nEdgeToCurve / 2; ++k)
                {
                    elmt->GetEdge(sMap.extEdgeFace[k][0])->m_parentCAD =
                        el[i]->GetEdge(sMap.extEdgeFace[k][0])->m_parentCAD;
                }
            }
            // For the other layers, copy CAD of side surfaces.
            else
            {
                for (int k = 0; k < sMap.nEdgeToCurve / 2; ++k)
                {
                    elmt->GetEdge(sMap.extEdgeFace[k][0])->m_parentCAD =
                        el[i]->GetFace(sMap.extEdgeFace[k][1])->m_parentCAD;
                }
            }

            //if edge already exists (i.e top of layer) grab it
            for (int k = 0; k < elmt->GetEdgeCount(); ++k)
            {
                EdgeSharedPtr ed = elmt->GetEdge(k);
                auto fi  = m_mesh->m_edgeSet.find(ed);

                if (fi != m_mesh->m_edgeSet.end())
                {
                    elmt->SetEdge(k, *fi);
                }
            }

            elToStack[elId].push_back(elmt);

            // Change the surface elements to match the layers of
            // elements on the boundary of the domain.
            map<int, int>::iterator it;
            for (it = bLink.begin(); it != bLink.end(); ++it)
            {
                int fid = it->first;
                int bl  = it->second;

                vector<NodeSharedPtr> qNodeList(4);
                for (int k = 0; k < 4; ++k)
                {
                    qNodeList[k] = nodeList[faceNodeMap[elType][fid][k]];
                }
                vector<int> tagBE;
                tagBE =
                    m_mesh->m_element[m_mesh->m_expDim - 1][bl]->GetTagList();
                ElmtConfig bconf(LibUtilities::eQuadrilateral, 1, false, false,
                                 false);
                ElementSharedPtr boundaryElmt =
                    GetElementFactory().CreateInstance(
                        LibUtilities::eQuadrilateral, bconf, qNodeList, tagBE);

                // Overwrite first layer boundary element with new
                // boundary element, otherwise push this back to end of
                // the boundary list
                if (j == 0)
                {
                    m_mesh->m_element[m_mesh->m_expDim - 1][bl] = boundaryElmt;
                }
                else
                {
                    m_mesh->m_element[m_mesh->m_expDim - 1].push_back(
                        boundaryElmt);
                }
            }

            m_mesh->m_element[m_mesh->m_expDim].push_back(elmt);
        }
    }

    // At this point the split stack has been split linearly and the stack is is
    // saved in elToStack.  The vertices are unique but edges are not, but we
    // need the existing edgeset.

    EdgeSet oldEdgeSet = m_mesh->m_edgeSet;
    ProcessEdges();

    // Edges are now unique.

    for (int i = 0; i < el.size(); ++i)
    {
        const int elId = el[i]->GetId();
        sIt            = splitEls.find(elId);

        if (sIt == splitEls.end())
        {
            continue;
        }

        SpatialDomains::Geometry3DSharedPtr gm = geomMap[elId];
        gm->FillGeom();
        StdRegions::StdExpansionSharedPtr xmape = gm->GetXmap();
        Array<OneD, NekDouble> coeffs0e = gm->GetCoeffs(0);
        Array<OneD, NekDouble> coeffs1e = gm->GetCoeffs(1);
        Array<OneD, NekDouble> coeffs2e = gm->GetCoeffs(2);
        Array<OneD, NekDouble> xce(xmape->GetTotPoints());
        Array<OneD, NekDouble> yce(xmape->GetTotPoints());
        Array<OneD, NekDouble> zce(xmape->GetTotPoints());
        xmape->BwdTrans(coeffs0e, xce);
        xmape->BwdTrans(coeffs1e, yce);
        xmape->BwdTrans(coeffs2e, zce);

        const int faceNum              = sIt->second;
        LibUtilities::ShapeType elType = el[i]->GetConf().m_e;

        SplitMapHelper &sMap = splitMap[elType][faceNum];

        int nSplitEdge = sMap.nEdgeToSplit;
        vector<vector<NodeSharedPtr> > edgeNodes(nSplitEdge);
        vector<ElementSharedPtr> stack = elToStack[elId];

        // make layers high-order
        for(int j = 0; j < sMap.nEdgeToSplit; j++)
        {
            int locEdge         = sMap.edgesToSplit[j];
            EdgeSharedPtr edg   = el[i]->GetEdge(locEdge);

            for(int k = 0; k < nl; k++)
            {
                EdgeSharedPtr nwEdg = stack[k]->GetEdge(locEdge);

                // If the edge was in the old mesh we dont want to touch it;
                // likewise if it wasn't but has already been curved, we don't
                // want it so we will add done edge to the set.
                auto f = oldEdgeSet.find(nwEdg);
                if (f != oldEdgeSet.end())
                {
                    continue;
                }

                if (edgeGeomMap.find(edg->m_id) == edgeGeomMap.end())
                {
                    continue;
                }

                SpatialDomains::GeometrySharedPtr geom = edgeGeomMap[edg->m_id];
                geom->FillGeom();
                StdRegions::StdExpansionSharedPtr xmap = geom->GetXmap();
                Array<OneD, NekDouble> coeffs0 = geom->GetCoeffs(0);
                Array<OneD, NekDouble> coeffs1 = geom->GetCoeffs(1);
                Array<OneD, NekDouble> coeffs2 = geom->GetCoeffs(2);
                Array<OneD, NekDouble> xc(xmap->GetTotPoints());
                Array<OneD, NekDouble> yc(xmap->GetTotPoints());
                Array<OneD, NekDouble> zc(xmap->GetTotPoints());
                xmap->BwdTrans(coeffs0, xc);
                xmap->BwdTrans(coeffs1, yc);
                xmap->BwdTrans(coeffs2, zc);

                NekDouble tb = paraEdg[nwEdg->m_n1->m_id][el[i]->GetId()];
                NekDouble te = paraEdg[nwEdg->m_n2->m_id][el[i]->GetId()];

                for (int l = 1; l < nq - 1; l++)
                {
                    Array<OneD, NekDouble> xp(1);
                    xp[0] =
                        tb * (1.0 - gll[l]) / 2.0 + te * (1.0 + gll[l]) / 2.0;
                    Array<OneD, NekDouble> loc(3);
                    loc[0] = xmap->PhysEvaluate(xp, xc);
                    loc[1] = xmap->PhysEvaluate(xp, yc);
                    loc[2] = xmap->PhysEvaluate(xp, zc);
                    nwEdg->m_edgeNodes.push_back(
                        NodeSharedPtr(new Node(0, loc[0], loc[1], loc[2])));

                    if (nwEdg->m_parentCAD)
                    {
                        if (nwEdg->m_parentCAD->GetType() == CADType::eCurve)
                        {
                            CADCurveSharedPtr c = std::dynamic_pointer_cast<
                                CADCurve>(nwEdg->m_parentCAD);
                            NekDouble t;
                            c->loct(loc, t);
                            nwEdg->m_edgeNodes.back()->SetCADCurve(c, t);
                        }
                        else if (nwEdg->m_parentCAD->GetType() == CADType::eSurf)
                        {
                            CADSurfSharedPtr s = std::dynamic_pointer_cast<
                                CADSurf>(nwEdg->m_parentCAD);
                            Array<OneD, NekDouble> uv = s->locuv(loc);
                            nwEdg->m_edgeNodes.back()->SetCADSurf(s, uv);
                        }
                    }
                }

                nwEdg->m_curveType = LibUtilities::eGaussLobattoLegendre;
                oldEdgeSet.insert(nwEdg);
            }
        }

        for (int j = 0; j < sMap.nEdgeToCurve; j++)
        {
            for (int k = 0; k < nl; k++)
            {
                map<int, int *> gllMap;
                for (int l = 0; l < sMap.size; l++)
                {
                    NodeSharedPtr n = stack[k]->GetVertex(l);
                    gllMap[n->m_id] = sMap.gll[l];
                }

                int locEdge          = sMap.edgesToCurve[j];
                EdgeSharedPtr nwEdg  = stack[k]->GetEdge(locEdge);
                auto f = oldEdgeSet.find(nwEdg);
                if (f != oldEdgeSet.end())
                {
                    continue;
                }

                nwEdg->m_edgeNodes.clear();

                Array<OneD, NekDouble> tb(3), te(3);
                tb[0] = gllMap[nwEdg->m_n1->m_id][0];
                tb[1] = gllMap[nwEdg->m_n1->m_id][1];
                tb[2] = gllMap[nwEdg->m_n1->m_id][2];
                te[0] = gllMap[nwEdg->m_n2->m_id][0];
                te[1] = gllMap[nwEdg->m_n2->m_id][1];
                te[2] = gllMap[nwEdg->m_n2->m_id][2];

                tb[sMap.blpDir] = paraElm[nwEdg->m_n1->m_id][el[i]->GetId()];
                te[sMap.blpDir] = paraElm[nwEdg->m_n2->m_id][el[i]->GetId()];

                for (int l = 1; l < nq - 1; l++)
                {
                    Array<OneD, NekDouble> xp(3);
                    for (int m = 0; m < 3; m++)
                    {
                        xp[m] = tb[m] * (1.0 - gll[l]) / 2.0 +
                                        te[m] * (1.0 + gll[l]) / 2.0;

                        if(xp[m] < -1.0 || xp[m] > 1.0)
                        {
                            cout << "WARNING: lies outside parameter range: "
                                 << xp[m] << endl;
                        }
                    }


                    Array<OneD, NekDouble> loc(3);
                    loc[0] = xmape->PhysEvaluate(xp, xce);
                    loc[1] = xmape->PhysEvaluate(xp, yce);
                    loc[2] = xmape->PhysEvaluate(xp, zce);
                    nwEdg->m_edgeNodes.push_back(
                        NodeSharedPtr(new Node(0, loc[0], loc[1], loc[2])));

                    if (nwEdg->m_parentCAD)
                    {
                        if (nwEdg->m_parentCAD->GetType() == CADType::eCurve)
                        {
                            CADCurveSharedPtr c =
                                std::dynamic_pointer_cast<CADCurve>(
                                    nwEdg->m_parentCAD);
                            NekDouble t;
                            c->loct(loc, t);
                            nwEdg->m_edgeNodes.back()->SetCADCurve(c, t);
                        }
                        else if (nwEdg->m_parentCAD->GetType() ==
                                 CADType::eSurf)
                        {
                            CADSurfSharedPtr s =
                                std::dynamic_pointer_cast<CADSurf>(
                                    nwEdg->m_parentCAD);
                            Array<OneD, NekDouble> uv = s->locuv(loc);
                            nwEdg->m_edgeNodes.back()->SetCADSurf(s, uv);
                        }
                    }
                }

                nwEdg->m_curveType = LibUtilities::eGaussLobattoLegendre;

                oldEdgeSet.insert(nwEdg);
            }
        }
    }

    ProcessVertices();
    ProcessEdges();
    ProcessFaces();
    ProcessElements();
    ProcessComposites();
}
}
}
