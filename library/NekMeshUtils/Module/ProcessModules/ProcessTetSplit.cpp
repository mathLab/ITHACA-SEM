////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessTetSplit.cpp
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
//  Description: Split prisms -> tets
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <StdRegions/StdNodalPrismExp.h>
#include <LocalRegions/PrismExp.h>

#include <NekMeshUtils/MeshElements/Element.h>

#include "ProcessTetSplit.h"

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{

typedef std::pair<int, int> ipair;

ModuleKey ProcessTetSplit::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "tetsplit"),
        ProcessTetSplit::create,
        "Split prismatic elements to tetrahedra");

ProcessTetSplit::ProcessTetSplit(MeshSharedPtr m) : ProcessModule(m)
{
    m_config["nq"] =
        ConfigOption(false, "3", "Number of points in high order elements.");
}

ProcessTetSplit::~ProcessTetSplit()
{
}

void ProcessTetSplit::Process()
{
    int nodeId = m_mesh->m_vertexSet.size();

    // Set up map which identifies edges (as pairs of vertex ids)
    // including their vertices to the offset/stride in the 3d array of
    // collapsed co-ordinates. Note that this map also includes the
    // diagonal edges along quadrilateral faces which will be used to
    // add high-order information to the split tetrahedra.
    map<ipair, ipair> edgeMap;

    int nq  = m_config["nq"].as<int>();
    int ne  = nq - 2;            // Edge interior
    int nft = (ne - 1) * ne / 2; // Face interior (triangle)
    int nfq = ne * ne;           // Face interior (quad)
    int o   = 6;

    // Standard prismatic edges (0->8)
    edgeMap[ipair(0, 1)] = ipair(o, 1);
    edgeMap[ipair(1, 2)] = ipair(o + 1 * ne, 1);
    edgeMap[ipair(3, 2)] = ipair(o + 2 * ne, 1);
    edgeMap[ipair(0, 3)] = ipair(o + 3 * ne, 1);
    edgeMap[ipair(0, 4)] = ipair(o + 4 * ne, 1);
    edgeMap[ipair(1, 4)] = ipair(o + 5 * ne, 1);
    edgeMap[ipair(2, 5)] = ipair(o + 6 * ne, 1);
    edgeMap[ipair(3, 5)] = ipair(o + 7 * ne, 1);
    edgeMap[ipair(4, 5)] = ipair(o + 8 * ne, 1);

    // Face 0 diagonals
    o += 9 * ne;
    edgeMap[ipair(0, 2)] = ipair(o, ne + 1);
    edgeMap[ipair(1, 3)] = ipair(o + ne - 1, ne - 1);

    // Face 2 diagonals
    o += nfq + nft;
    edgeMap[ipair(1, 5)] = ipair(o, ne + 1);
    edgeMap[ipair(2, 4)] = ipair(o + ne - 1, ne - 1);

    // Face 4 diagonals
    o += nfq + nft;
    edgeMap[ipair(0, 5)] = ipair(o, ne + 1);
    edgeMap[ipair(3, 4)] = ipair(o + ne - 1, ne - 1);

    // Default PointsType.
    LibUtilities::PointsType pt = LibUtilities::ePolyEvenlySpaced;

    /*
     * Split all element types into tetrahedra. This is based on the
     * algorithm found in:
     *
     * "How to Subdivide Pyramids, Prisms and Hexahedra into
     * Tetrahedra", J. Dompierre et al.
     */

    // Denotes a set of indices inside m_mesh->m_element[m_mesh->m_expDim-1]
    // which are
    // to be removed. These are precisely the quadrilateral boundary
    // faces which will be replaced by two triangular faces.
    set<int> toRemove;

    // Represents table 2 of paper; each row i represents a rotation of
    // the prism nodes such that vertex i is placed at position 0.
    static int indir[6][6] = {{0, 1, 2, 3, 4, 5},
                              {1, 2, 0, 4, 5, 3},
                              {2, 0, 1, 5, 3, 4},
                              {3, 5, 4, 0, 2, 1},
                              {4, 3, 5, 1, 0, 2},
                              {5, 4, 3, 2, 1, 0}};

    // Represents table 3 of paper; the first three rows are the three
    // tetrahedra if the first condition is met; the latter three rows
    // are the three tetrahedra if the second condition is met.
    static int prismTet[6][4] = {{0, 1, 2, 5},
                                 {0, 1, 5, 4},
                                 {0, 4, 5, 3},
                                 {0, 1, 2, 4},
                                 {0, 4, 2, 5},
                                 {0, 4, 5, 3}};

    // Represents the order of tetrahedral edges (in Nektar++ ordering).
    static int tetEdges[6][2] = {
        {0, 1}, {1, 2}, {0, 2}, {0, 3}, {1, 3}, {2, 3}};

    // A tetrahedron nodes -> faces map.
    static int tetFaceNodes[4][3] = {
        {0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {0, 2, 3}};

    static NodeSharedPtr stdPrismNodes[6] = {
        NodeSharedPtr(new Node(0, -1, -1, -1)),
        NodeSharedPtr(new Node(1, 1, -1, -1)),
        NodeSharedPtr(new Node(2, 1, 1, -1)),
        NodeSharedPtr(new Node(3, -1, 1, -1)),
        NodeSharedPtr(new Node(4, -1, -1, 1)),
        NodeSharedPtr(new Node(5, -1, 1, 1))};

    // Generate equally spaced nodal points on triangle.
    Array<OneD, NekDouble> rp(nq * (nq + 1) / 2), sp(nq * (nq + 1) / 2);
    LibUtilities::PointsKey elec(nq, LibUtilities::eNodalTriEvenlySpaced);
    LibUtilities::PointsManager()[elec]->GetPoints(rp, sp);

    // Make a copy of the element list.
    vector<ElementSharedPtr> el = m_mesh->m_element[m_mesh->m_expDim];
    m_mesh->m_element[m_mesh->m_expDim].clear();

    for (int i = 0; i < el.size(); ++i)
    {
        if (el[i]->GetConf().m_e != LibUtilities::ePrism)
        {
            m_mesh->m_element[m_mesh->m_expDim].push_back(el[i]);
            continue;
        }

        vector<NodeSharedPtr> nodeList(6);

        // Map Nektar++ ordering (vertices 0,1,2,3 are base quad) to
        // paper ordering (vertices 0,1,2 are first triangular face).
        int mapPrism[6] = {0, 1, 4, 3, 2, 5};
        for (int j = 0; j < 6; ++j)
        {
            nodeList[j] = el[i]->GetVertex(mapPrism[j]);
        }

        // Determine minimum ID of the nodes in this prism.
        int minElId = nodeList[0]->m_id;
        int minId = 0;
        for (int j = 1; j < 6; ++j)
        {
            int curId = nodeList[j]->m_id;
            if (curId < minElId)
            {
                minElId = curId;
                minId   = j;
            }
        }

        int offset;

        // Split prism using paper criterion.
        int id1 = min(nodeList[indir[minId][1]]->m_id,
                      nodeList[indir[minId][5]]->m_id);
        int id2 = min(nodeList[indir[minId][2]]->m_id,
                      nodeList[indir[minId][4]]->m_id);

        if (id1 < id2)
        {
            offset = 0;
        }
        else if (id1 > id2)
        {
            offset = 3;
        }
        else
        {
            cerr << "Connectivity issue with prism->tet splitting." << endl;
            abort();
        }

        // Create local prismatic region so that co-ordinates of the
        // mapped element can be read from.
        SpatialDomains::PrismGeomSharedPtr geomLayer =
            std::dynamic_pointer_cast<SpatialDomains::PrismGeom>(
                el[i]->GetGeom(m_mesh->m_spaceDim));
        LibUtilities::BasisKey B0(
            LibUtilities::eOrtho_A,
            nq,
            LibUtilities::PointsKey(nq, LibUtilities::eGaussLobattoLegendre));
        LibUtilities::BasisKey B1(
            LibUtilities::eOrtho_B,
            nq,
            LibUtilities::PointsKey(nq, LibUtilities::eGaussRadauMAlpha1Beta0));
        LocalRegions::PrismExpSharedPtr qs =
            MemoryManager<LocalRegions::PrismExp>::AllocateSharedPtr(
                B0, B0, B1, geomLayer);

        // Get the coordinates of the high order prismatic element.
        Array<OneD, NekDouble> wsp(3 * nq * nq * nq);
        Array<OneD, Array<OneD, NekDouble> > x(3);
        x[0] = wsp;
        x[1] = wsp + 1 * nq * nq * nq;
        x[2] = wsp + 2 * nq * nq * nq;
        qs->GetCoords(x[0], x[1], x[2]);

        LibUtilities::BasisKey NB0(
            LibUtilities::eModified_A, nq, LibUtilities::PointsKey(nq, pt));
        LibUtilities::BasisKey NB1(
            LibUtilities::eModified_B, nq, LibUtilities::PointsKey(nq, pt));

        // Process face data. Initially put coordinates into equally
        // spaced nodal distribution.
        StdRegions::StdNodalPrismExpSharedPtr nodalPrism =
            MemoryManager<StdRegions::StdNodalPrismExp>::AllocateSharedPtr(
                B0, B0, B1, LibUtilities::eNodalPrismEvenlySpaced);

        int nCoeffs = nodalPrism->GetNcoeffs();
        Array<OneD, NekDouble> wsp2(3 * nCoeffs);
        Array<OneD, Array<OneD, NekDouble> > xn(3);
        xn[0] = wsp2;
        xn[1] = wsp2 + 1 * nCoeffs;
        xn[2] = wsp2 + 2 * nCoeffs;

        for (int j = 0; j < 3; ++j)
        {
            Array<OneD, NekDouble> tmp(nodalPrism->GetNcoeffs());
            qs->FwdTrans(x[j], tmp);
            nodalPrism->ModalToNodal(tmp, xn[j]);
        }

        // Store map of nodes.
        map<int, int> prismVerts;
        for (int j = 0; j < 6; ++j)
        {
            prismVerts[el[i]->GetVertex(j)->m_id] = j;
        }

        for (int j = 0; j < 3; ++j)
        {
            vector<NodeSharedPtr> tetNodes(4);

            // Extract vertices for tetrahedron.
            for (int k = 0; k < 4; ++k)
            {
                tetNodes[k] = nodeList[indir[minId][prismTet[j + offset][k]]];
            }

            // Add high order information to tetrahedral edges.
            for (int k = 0; k < 6; ++k)
            {
                // Determine prismatic nodes which correspond with this
                // edge. Apply prism map to this to get Nektar++
                // ordering (this works since as a permutation,
                // prismMap^2 = id).
                int n1 = mapPrism
                    [indir[minId][prismTet[j + offset][tetEdges[k][0]]]];
                int n2 = mapPrism
                    [indir[minId][prismTet[j + offset][tetEdges[k][1]]]];

                // Find offset/stride
                auto it = edgeMap.find(pair<int, int>(n1, n2));
                if (it == edgeMap.end())
                {
                    it = edgeMap.find(pair<int, int>(n2, n1));
                    if (it == edgeMap.end())
                    {
                        cerr << "Couldn't find prism edges " << n1 << " " << n2
                             << endl;
                        abort();
                    }
                    // Extract vertices -- reverse order.
                    for (int l = ne - 1; l >= 0; --l)
                    {
                        int pos = it->second.first + l * it->second.second;
                        tetNodes.push_back(NodeSharedPtr(new Node(
                            nodeId++, xn[0][pos], xn[1][pos], xn[2][pos])));
                    }
                }
                else
                {
                    // Extract vertices -- forwards order.
                    for (int l = 0; l < ne; ++l)
                    {
                        int pos = it->second.first + l * it->second.second;
                        tetNodes.push_back(NodeSharedPtr(new Node(
                            nodeId++, xn[0][pos], xn[1][pos], xn[2][pos])));
                    }
                }
            }

            // Create new tetrahedron with edge curvature.
            vector<int> tags = el[i]->GetTagList();
            ElmtConfig conf(LibUtilities::eTetrahedron, nq - 1, false, false);
            ElementSharedPtr elmt = GetElementFactory().CreateInstance(
                LibUtilities::eTetrahedron, conf, tetNodes, tags);

            // Extract interior face data.
            for (int k = 0; k < 4; ++k)
            {
                // First determine which nodes of prism are being used
                // for this face.
                FaceSharedPtr face = elmt->GetFace(k);
                vector<NodeSharedPtr> triNodes(3);

                for (int l = 0; l < 3; ++l)
                {
                    NodeSharedPtr v = face->m_vertexList[l];
                    triNodes[l]     = stdPrismNodes[prismVerts[v->m_id]];
                }

                // Create a triangle with the standard nodes of the
                // prism.
                vector<int> tags;
                ElmtConfig conf(LibUtilities::eTriangle, 1, false, false);
                ElementSharedPtr elmt = GetElementFactory().CreateInstance(
                    LibUtilities::eTriangle, conf, triNodes, tags);
                SpatialDomains::GeometrySharedPtr triGeom = elmt->GetGeom(3);
                triGeom->FillGeom();
                o                 = 3 + 3 * ne;
                face->m_curveType = LibUtilities::eNodalTriEvenlySpaced;

                for (int l = 0; l < nft; ++l)
                {
                    Array<OneD, NekDouble> tmp1(2), tmp2(3);
                    tmp1[0] = rp[o + l];
                    tmp1[1] = sp[o + l];
                    tmp2[0] = triGeom->GetCoord(0, tmp1);
                    tmp2[1] = triGeom->GetCoord(1, tmp1);
                    tmp2[2] = triGeom->GetCoord(2, tmp1);

                    // PhysEvaluate on prism geometry.
                    NekDouble xc = geomLayer->GetCoord(0, tmp2);
                    NekDouble yc = geomLayer->GetCoord(1, tmp2);
                    NekDouble zc = geomLayer->GetCoord(2, tmp2);
                    face->m_faceNodes.push_back(
                        NodeSharedPtr(new Node(nodeId++, xc, yc, zc)));
                }
            }

            m_mesh->m_element[m_mesh->m_expDim].push_back(elmt);
        }

        // Now check to see if this one of the quadrilateral faces is
        // associated with a boundary condition. If it is, we split the
        // face into two triangles and mark the existing face for
        // removal.
        //
        // Note that this algorithm has significant room for improvement
        // and is likely one of the least optimal approachs - however
        // implementation is simple.
        for (int fid = 0; fid < 5; fid += 2)
        {
            int bl = el[i]->GetBoundaryLink(fid);

            if (bl == -1)
            {
                continue;
            }

            vector<NodeSharedPtr> triNodeList(3);
            vector<int> faceNodes(3);
            vector<int> tmp;
            vector<int> tagBE;
            ElmtConfig bconf(LibUtilities::eTriangle, 1, true, true);
            ElementSharedPtr elmt;

            // Mark existing boundary face for removal.
            toRemove.insert(bl);
            tagBE = m_mesh->m_element[m_mesh->m_expDim - 1][bl]->GetTagList();

            // First loop over tets.
            for (int j = 0; j < 3; ++j)
            {
                // Now loop over faces.
                for (int k = 0; k < 4; ++k)
                {
                    // Finally determine Nektar++ local node numbers for
                    // this face.
                    for (int l = 0; l < 3; ++l)
                    {
                        faceNodes[l] = mapPrism
                            [indir[minId]
                                  [prismTet[j + offset][tetFaceNodes[k][l]]]];
                    }

                    tmp = faceNodes;
                    sort(faceNodes.begin(), faceNodes.end());

                    // If this face matches a triple denoting a split
                    // quad face, add the face to the expansion list.
                    if ((fid == 0 && ((faceNodes[0] == 0 && faceNodes[1] == 1 &&
                                       faceNodes[2] == 2) ||
                                      (faceNodes[0] == 0 && faceNodes[1] == 2 &&
                                       faceNodes[2] == 3) ||
                                      (faceNodes[0] == 0 && faceNodes[1] == 1 &&
                                       faceNodes[2] == 3) ||
                                      (faceNodes[0] == 1 && faceNodes[1] == 2 &&
                                       faceNodes[2] == 3))) ||
                        (fid == 2 && ((faceNodes[0] == 1 && faceNodes[1] == 2 &&
                                       faceNodes[2] == 5) ||
                                      (faceNodes[0] == 1 && faceNodes[1] == 4 &&
                                       faceNodes[2] == 5) ||
                                      (faceNodes[0] == 1 && faceNodes[1] == 2 &&
                                       faceNodes[2] == 4) ||
                                      (faceNodes[0] == 2 && faceNodes[1] == 4 &&
                                       faceNodes[2] == 5))) ||
                        (fid == 4 && ((faceNodes[0] == 0 && faceNodes[1] == 3 &&
                                       faceNodes[2] == 5) ||
                                      (faceNodes[0] == 0 && faceNodes[1] == 4 &&
                                       faceNodes[2] == 5) ||
                                      (faceNodes[0] == 0 && faceNodes[1] == 3 &&
                                       faceNodes[2] == 4) ||
                                      (faceNodes[0] == 3 && faceNodes[1] == 4 &&
                                       faceNodes[2] == 5))))
                    {
                        triNodeList[0] = nodeList[mapPrism[tmp[0]]];
                        triNodeList[1] = nodeList[mapPrism[tmp[1]]];
                        triNodeList[2] = nodeList[mapPrism[tmp[2]]];
                        elmt = GetElementFactory().CreateInstance(
                            LibUtilities::eTriangle, bconf, triNodeList, tagBE);
                        m_mesh->m_element[m_mesh->m_expDim - 1].push_back(elmt);
                    }
                }
            }
        }
    }

    // Remove 2D elements.
    vector<ElementSharedPtr> tmp;
    for (int i = 0; i < m_mesh->m_element[m_mesh->m_expDim - 1].size(); ++i)
    {
        if (toRemove.find(i) == toRemove.end())
        {
            tmp.push_back(m_mesh->m_element[m_mesh->m_expDim - 1][i]);
        }
    }

    m_mesh->m_element[m_mesh->m_expDim - 1] = tmp;

    // Re-process mesh to eliminate duplicate vertices and edges.
    ProcessVertices();
    ProcessEdges();
    ProcessFaces();
    ProcessElements();
    ProcessComposites();
}
}
}
