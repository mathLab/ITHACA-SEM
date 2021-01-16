////////////////////////////////////////////////////////////////////////////////
//
//  File: Mesh.cpp
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
//  Description: Mesh object.
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/ManagerAccess.h>
#include <NekMesh/MeshElements/Mesh.h>

#include <NekMesh/CADSystem/CADCurve.h>
#include <NekMesh/CADSystem/CADSurf.h>

using namespace std;

namespace Nektar
{
namespace NekMesh
{

/**
 * @brief Return the number of elements of the expansion dimension.
 */
unsigned int Mesh::GetNumElements()
{
    return m_element[m_expDim].size();
}

/**
 * @brief Return the number of boundary elements (i.e. one below the
 * expansion dimension).
 */
unsigned int Mesh::GetNumBndryElements()
{
    unsigned int i, nElmt = 0;

    for (i = 0; i < m_expDim; ++i)
        nElmt += m_element[i].size();

    return nElmt;
}

/**
 * @brief Return the total number of entities in the mesh (i.e. all
 * elements, regardless of dimension).
 */
unsigned int Mesh::GetNumEntities()
{
    unsigned int nEnt = 0;

    for (unsigned int d = 0; d <= m_expDim; ++d)
    {
        nEnt += m_element[d].size();
    }

    return nEnt;
}

/**
 * @brief Convert this mesh into a mesh of uniform polynomial order @p order
 * with a curve point distribution @p distType.
 *
 * This routine adds curvature points into a mesh so that the resulting elements
 * are all of a uniform order @p order and all high-order vertices are
 * consistently ordered. It proceeds in a bottom-up fashion:
 *
 * - First construct all edge, face and elemental geometry mappings.
 * - Then call the local MakeOrder functions on each edge, face and element of
 *   dimension Mesh::m_expDim.
 * - Finally, any boundary elements are updated so that they have the same
 *   interior degrees of freedom as their corresponding edge or face links.
 */
void Mesh::MakeOrder(int order, LibUtilities::PointsType distType, Logger &log)
{
    // Going to make a copy of the curavture information, since this is cheaper
    // than using Nektar's Geometry objects. Currently, the geometry objects
    // which make up a 3D element dont use the volume nodes, they are just
    // stored, so we can get away without copying them.

    int id = m_vertexSet.size();

    EdgeSet::iterator eit;
    FaceSet::iterator fit;

    std::unordered_map<int, EdgeSharedPtr> edgeCopies;
    std::unordered_map<int, FaceSharedPtr> faceCopies;

    // Decide on distribution of points to use for each shape type based on the
    // input we've been supplied.
    std::map<LibUtilities::ShapeType, LibUtilities::PointsType> pTypes;
    if (distType == LibUtilities::ePolyEvenlySpaced)
    {
        pTypes[LibUtilities::eSegment]  = LibUtilities::ePolyEvenlySpaced;
        pTypes[LibUtilities::eTriangle] = LibUtilities::eNodalTriEvenlySpaced;
        pTypes[LibUtilities::eQuadrilateral] = LibUtilities::ePolyEvenlySpaced;
        pTypes[LibUtilities::eTetrahedron] =
            LibUtilities::eNodalTetEvenlySpaced;
        pTypes[LibUtilities::ePrism] = LibUtilities::eNodalPrismEvenlySpaced;
        pTypes[LibUtilities::eHexahedron] = LibUtilities::ePolyEvenlySpaced;
    }
    else if (distType == LibUtilities::eGaussLobattoLegendre)
    {
        pTypes[LibUtilities::eSegment]  = LibUtilities::eGaussLobattoLegendre;
        pTypes[LibUtilities::eTriangle] = LibUtilities::eNodalTriElec;
        pTypes[LibUtilities::eQuadrilateral] =
            LibUtilities::eGaussLobattoLegendre;
        pTypes[LibUtilities::ePrism]       = LibUtilities::eNodalPrismElec;
        pTypes[LibUtilities::eTetrahedron] = LibUtilities::eNodalTetElec;
        pTypes[LibUtilities::eHexahedron] = LibUtilities::eGaussLobattoLegendre;
    }
    else
    {
        ASSERTL1(false, "Mesh::MakeOrder does not support this points type.");
    }

    // Begin by copying mesh objects for edges and faces so that we don't affect
    // any neighbouring elements in the mesh as we process each element. At the
    // same time we delete the curvature from the original edge and face, which
    // will be re-added with the MakeOrder routine.

    // First, we fill in the volume-interior nodes. This preserves the original
    // curvature of the mesh.
    const int nElmt = m_element[m_expDim].size();
    int tmpId       = 0;
    for (int i = 0; i < nElmt; ++i)
    {
        log(VERBOSE).Progress(i, nElmt, "MakeOrder: Elements");
        ElementSharedPtr el                    = m_element[m_expDim][i];
        SpatialDomains::GeometrySharedPtr geom = el->GetGeom(m_spaceDim);
        geom->FillGeom();
        el->MakeOrder(order, geom, pTypes[el->GetConf().m_e], m_spaceDim,
                      tmpId);
    }

    log(VERBOSE).Newline();

    // Now make copies of each of the edges.
    for (eit = m_edgeSet.begin(); eit != m_edgeSet.end(); eit++)
    {
        edgeCopies[(*eit)->m_id] = EdgeSharedPtr(new Edge(*(*eit)));
        (*eit)->m_edgeNodes.clear();
    }

    // Now copy faces. Make sure that this is a "deep copy", so that the face's
    // edge list corresponds to the copied edges, otherwise we end up in a
    // non-consistent state.
    for (fit = m_faceSet.begin(); fit != m_faceSet.end(); fit++)
    {
        FaceSharedPtr tmpFace = FaceSharedPtr(new Face(*(*fit)));

        for (int i = 0; i < tmpFace->m_edgeList.size(); ++i)
        {
            tmpFace->m_edgeList[i] = edgeCopies[tmpFace->m_edgeList[i]->m_id];
        }

        faceCopies[(*fit)->m_id] = tmpFace;
        (*fit)->m_faceNodes.clear();
    }

    std::unordered_set<int> processedEdges, processedFaces, processedVolumes;

    // note if CAD previously existed on the face or edge, the new points need
    // to be projected onto the CAD entity.

    // Call MakeOrder with our generated geometries on each edge to fill in edge
    // interior nodes.
    int ct = 0;
    for (eit = m_edgeSet.begin(); eit != m_edgeSet.end(); eit++, ct++)
    {
        log(VERBOSE).Progress(ct, m_edgeSet.size(), "MakeOrder: Edges");

        int edgeId = (*eit)->m_id;

        if (processedEdges.find(edgeId) != processedEdges.end())
        {
            continue;
        }

        EdgeSharedPtr cpEdge                   = edgeCopies[edgeId];
        SpatialDomains::GeometrySharedPtr geom = cpEdge->GetGeom(m_spaceDim);
        geom->FillGeom();

        (*eit)->MakeOrder(order, geom, pTypes[LibUtilities::eSegment],
                          m_spaceDim, id);
        processedEdges.insert(edgeId);
    }

    if (m_edgeSet.size() > 0)
    {
        log(VERBOSE).Newline();
    }

    // Call MakeOrder with our generated geometries on each face to fill in face
    // interior nodes.
    ct = 0;
    for (fit = m_faceSet.begin(); fit != m_faceSet.end(); fit++, ct++)
    {
        log(VERBOSE).Progress(ct, m_faceSet.size(), "MakeOrder: Faces");

        int faceId = (*fit)->m_id;

        if (processedFaces.find(faceId) != processedFaces.end())
        {
            continue;
        }

        FaceSharedPtr cpFace                   = faceCopies[faceId];
        SpatialDomains::GeometrySharedPtr geom = cpFace->GetGeom(m_spaceDim);
        geom->FillGeom();

        LibUtilities::ShapeType type = (*fit)->m_vertexList.size() == 3
                                           ? LibUtilities::eTriangle
                                           : LibUtilities::eQuadrilateral;
        (*fit)->MakeOrder(order, geom, pTypes[type], m_spaceDim, id);
        processedFaces.insert(faceId);
    }

    if (m_faceSet.size() > 0)
    {
        log(VERBOSE).Newline();
    }

    // Copy curvature into boundary conditions
    for (int i = 0; i < m_element[1].size(); ++i)
    {
        ElementSharedPtr el = m_element[1][i];
        EdgeSharedPtr edge  = el->GetEdgeLink();

        if (!edge)
        {
            continue;
        }

        // Copy face curvature
        el->MakeOrder(order, SpatialDomains::GeometrySharedPtr(),
                      pTypes[el->GetConf().m_e], m_spaceDim, id, true);
        el->SetVolumeNodes(edge->m_edgeNodes);
    }

    for (int i = 0; i < m_element[2].size(); ++i)
    {
        ElementSharedPtr el = m_element[2][i];
        FaceSharedPtr face  = el->GetFaceLink();

        if (!face)
        {
            continue;
        }

        // Copy face curvature
        el->MakeOrder(order, SpatialDomains::GeometrySharedPtr(),
                      pTypes[el->GetConf().m_e], m_spaceDim, id, true);
        el->SetVolumeNodes(face->m_faceNodes);
    }

    for (int i = 0; i < nElmt; ++i)
    {
        vector<NodeSharedPtr> tmp = m_element[m_expDim][i]->GetVolumeNodes();
        for (int j = 0; j < tmp.size(); ++j)
        {
            tmp[j]->m_id = id++;
        }
    }
}

/**
 * @brief Print out basic statistics of this mesh.
 */
void Mesh::PrintStats(Logger &log)
{
    log << "Mesh statistics:" << std::endl;

    log << "  - Mesh dimension       : " << m_spaceDim << std::endl
        << "  - Element dimension    : " << m_expDim << std::endl
        << "  - Has CAD attached     : " << (m_cad ? "yes" : "no") << std::endl
        << "  - Node count           : " << m_vertexSet.size() << std::endl;

    if (m_edgeSet.size() > 0)
    {
        log << "  - Edge count           : " << m_edgeSet.size() << std::endl;
    }

    if (m_faceSet.size() > 0)
    {
        log << "  - Face count           : " << m_faceSet.size() << std::endl;
    }

    log << "  - Elements             : " << m_element[m_expDim].size()
        << std::endl;
    log << "  - Bnd elements         : " << m_element[m_expDim - 1].size()
        << std::endl;

    // Print out number of composites
    log << "  - Number of composites : " << m_composite.size() << std::endl;

    // Calculate domain extent
    auto extent = m_element[m_expDim][0]->GetBoundingBox();
    for (int i = 1; i < m_element[m_expDim].size(); ++i)
    {
        auto el           = m_element[m_expDim][i]->GetBoundingBox();
        extent.first.m_x  = std::min(extent.first.m_x, el.first.m_x);
        extent.first.m_y  = std::min(extent.first.m_y, el.first.m_y);
        extent.first.m_z  = std::min(extent.first.m_z, el.first.m_z);
        extent.second.m_x = std::max(extent.second.m_x, el.second.m_x);
        extent.second.m_y = std::max(extent.second.m_y, el.second.m_y);
        extent.second.m_z = std::max(extent.second.m_z, el.second.m_z);
    }

    log << "  - Lower mesh extent    : " << extent.first.m_x << " "
        << extent.first.m_y << " " << extent.first.m_z << std::endl
        << "  - Upper mesh extent    : " << extent.second.m_x << " "
        << extent.second.m_y << " " << extent.second.m_z << std::endl;

    std::map<LibUtilities::ShapeType, std::pair<int, int>> elmtCounts;

    for (int i = 1; i < LibUtilities::SIZE_ShapeType; ++i)
    {
        elmtCounts[(LibUtilities::ShapeType)i] = std::make_pair(0, 0);
    }

    for (int dim = 0; dim <= 3; ++dim)
    {
        for (auto &elmt : m_element[dim])
        {
            auto &counts = elmtCounts[elmt->GetShapeType()];

            if (elmt->IsDeformed())
            {
                counts.second++;
            }
            else
            {
                counts.first++;
            }
        }
    }

    log << "Element counts (regular/deformed/total):" << std::endl;
    for (int i = 1; i < LibUtilities::SIZE_ShapeType; ++i)
    {
        auto shapeType = (LibUtilities::ShapeType)i;
        auto counts    = elmtCounts[shapeType];

        if (counts.first + counts.second == 0)
        {
            continue;
        }

        log << "  - " << std::setw(12) << std::left
            << LibUtilities::ShapeTypeMap[(LibUtilities::ShapeType)i] << ": "
            << setw(12) << counts.first << "  " << setw(12) << counts.second
            << "  " << setw(12) << counts.first + counts.second << std::endl;
    }
}

} // namespace NekMesh
} // namespace Nektar
