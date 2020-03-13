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
//  Description: Mesh object.
//
////////////////////////////////////////////////////////////////////////////////

#include <NekMeshUtils/MeshElements/Mesh.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/BasicUtils/Progressbar.hpp>

using namespace std;

namespace Nektar
{
namespace NekMeshUtils
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
void Mesh::MakeOrder(int                      order,
                     LibUtilities::PointsType distType)
{
    int id = m_vertexSet.size();

    EdgeSet::iterator eit;
    FaceSet::iterator fit;

    boost::unordered_map<int, SpatialDomains::Geometry1DSharedPtr> edgeGeoms;
    boost::unordered_map<int, SpatialDomains::Geometry2DSharedPtr> faceGeoms;
    boost::unordered_map<int, SpatialDomains::GeometrySharedPtr> volGeoms;

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
        pTypes[LibUtilities::ePrism] = LibUtilities::eNodalPrismElec;
        pTypes[LibUtilities::eTetrahedron] = LibUtilities::eNodalTetElec;
        pTypes[LibUtilities::eHexahedron] = LibUtilities::eGaussLobattoLegendre;
    }
    else
    {
        ASSERTL1(false, "Mesh::MakeOrder does not support this points type.");
    }

    // Begin by generating Nektar++ geometry objects for edges, faces and
    // elements so that we don't affect any neighbouring elements in the mesh as
    // we process each element.
    for(eit = m_edgeSet.begin(); eit != m_edgeSet.end(); eit++)
    {
        SpatialDomains::Geometry1DSharedPtr geom =
            (*eit)->GetGeom(m_spaceDim);
        geom->FillGeom();
        edgeGeoms[(*eit)->m_id] = geom;
    }

    for(fit = m_faceSet.begin(); fit != m_faceSet.end(); fit++)
    {
        SpatialDomains::Geometry2DSharedPtr geom =
            (*fit)->GetGeom(m_spaceDim);
        geom->FillGeom();
        faceGeoms[(*fit)->m_id] = geom;
    }

    for(int i = 0; i < m_element[m_expDim].size(); i++)
    {
        ElementSharedPtr el = m_element[m_expDim][i];
        SpatialDomains::GeometrySharedPtr geom =
            el->GetGeom(m_spaceDim);
        geom->FillGeom();
        volGeoms[el->GetId()] = geom;
    }

    boost::unordered_set<int> processedEdges, processedFaces, processedVolumes;

    // note if CAD previously existed on the face or edge, the new points need
    // to be projected onto the CAD entity.

    // Call MakeOrder with our generated geometries on each edge to fill in edge
    // interior nodes.
    int ct = 0;
    for (eit = m_edgeSet.begin(); eit != m_edgeSet.end(); eit++, ct++)
    {
        if (m_verbose)
        {
            LibUtilities::PrintProgressbar(
                ct, m_edgeSet.size(), "MakeOrder: Edges: ");
        }
        int edgeId = (*eit)->m_id;

        if (processedEdges.find(edgeId) != processedEdges.end())
        {
            continue;
        }

        (*eit)->MakeOrder(order, edgeGeoms[edgeId],
                          pTypes[LibUtilities::eSegment], m_spaceDim, id);
        processedEdges.insert(edgeId);
    }

    // Call MakeOrder with our generated geometries on each face to fill in face
    // interior nodes.
    ct = 0;
    for (fit = m_faceSet.begin(); fit != m_faceSet.end(); fit++, ct++)
    {
        if (m_verbose)
        {
            LibUtilities::PrintProgressbar(
                ct, m_faceSet.size(), "MakeOrder: Faces: ");
        }
        int faceId = (*fit)->m_id;

        if (processedFaces.find(faceId) != processedFaces.end())
        {
            continue;
        }

        LibUtilities::ShapeType type = (*fit)->m_vertexList.size() == 3 ?
            LibUtilities::eTriangle : LibUtilities::eQuadrilateral;
        (*fit)->MakeOrder(order, faceGeoms[faceId], pTypes[type], m_spaceDim,
                          id);
        processedFaces.insert(faceId);
    }

    // Copy curvature into boundary conditions
    for (int i = 0; i < m_element[1].size(); ++i)
    {
        ElementSharedPtr el = m_element[1][i];
        EdgeSharedPtr edge = el->GetEdgeLink();

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
        FaceSharedPtr face = el->GetFaceLink();

        if (!face)
        {
            continue;
        }

        // Copy face curvature
        el->MakeOrder(order, SpatialDomains::GeometrySharedPtr(),
                      pTypes[el->GetConf().m_e], m_spaceDim, id, true);
        el->SetVolumeNodes(face->m_faceNodes);
    }

    // Finally, fill in volumes.
    const int nElmt = m_element[m_expDim].size();
    for (int i = 0; i < nElmt; ++i)
    {
        if (m_verbose)
        {
            LibUtilities::PrintProgressbar(i, nElmt, "MakeOrder: Elements: ");
        }
        ElementSharedPtr el = m_element[m_expDim][i];
        el->MakeOrder(order, volGeoms[el->GetId()], pTypes[el->GetConf().m_e],
                      m_spaceDim, id);
    }

    if (m_verbose)
    {
        cout << endl;
    }
}

}
}
