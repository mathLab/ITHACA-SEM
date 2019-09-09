////////////////////////////////////////////////////////////////////////////////
//
//  File: Mesh.h
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

#ifndef NEKMESHUTILS_MESHELEMENTS_MESH
#define NEKMESHUTILS_MESHELEMENTS_MESH

#include <set>

#include <LibUtilities/BasicUtils/FieldIO.h>

#include <NekMeshUtils/NekMeshUtilsDeclspec.h>
#include <NekMeshUtils/MeshElements/Element.h>
#include <NekMeshUtils/MeshElements/Composite.h>

namespace Nektar
{
namespace NekMeshUtils
{

class Octree;
typedef std::shared_ptr<Octree> OctreeSharedPtr;

/**
 * Enumeration of condition types (Dirichlet, Neumann, etc).
 */
enum ConditionType
{
    eDirichlet,
    eNeumann,
    eRobin,
    ePeriodic,
    eHOPCondition,
    SIZE_ConditionType
};

/**
 * @brief Defines a boundary condition.
 *
 * A boundary condition is defined by its type (e.g. Dirichlet), the
 * field it applies to, the value imposed on this field and the
 * composite which the boundary condition is applied to.
 */
struct Condition
{
    Condition() : type(), field(), value(), m_composite()
    {
    }
    std::vector<ConditionType> type;
    std::vector<std::string> field;
    std::vector<std::string> value;
    std::vector<int> m_composite;
};

typedef std::shared_ptr<Condition> ConditionSharedPtr;
typedef std::map<int, ConditionSharedPtr> ConditionMap;

NEKMESHUTILS_EXPORT bool operator==(ConditionSharedPtr const &c1,
                                    ConditionSharedPtr const &c2);

class Mesh
{
public:
    NEKMESHUTILS_EXPORT Mesh() : m_verbose(false), m_nummode(0)
    {
    }

    /// Verbose flag
    bool                            m_verbose;
    /// Dimension of the expansion.
    unsigned int                    m_expDim;
    /// Dimension of the space in which the mesh is defined.
    unsigned int                    m_spaceDim;
    /// a order tag to aid output, a bit of a hack
    unsigned int                    m_nummode;
    /// List of mesh nodes.
    std::vector<NodeSharedPtr>      m_node;
    /// Set of element vertices.
    NodeSet                         m_vertexSet;
    /// used for meshing purposes to keep trac of ids
    int                             m_numNodes;
    /// Set of element edges.
    EdgeSet                         m_edgeSet;
    /// Set of element faces.
    FaceSet                         m_faceSet;
    /// Map for elements.
    ElementMap                      m_element;
    /// Map for composites.
    CompositeMap                    m_composite;
    /// Boundary conditions maps tag to condition.
    ConditionMap                    m_condition;
    /// List of fields names.
    std::vector<std::string>        m_fields;
    /// Map of vertex normals.
    std::unordered_map<int, Node>   m_vertexNormals;
    /// Set of all pairs of element ID and edge/face number on which to
    /// apply spherigon surface smoothing.
    std::set<std::pair<int, int> >  m_spherigonSurfs;
    /// List of face labels for composite annotation
    std::map<int, std::string>      m_faceLabels;
    /// CAD system pointer, if there is no cad its empty
    CADSystemSharedPtr              m_cad;
    /// Octree system pointer, if there is no octree its empty
    OctreeSharedPtr                 m_octree;
    /// Metadata map for storing any mesh generation parameters
    LibUtilities::FieldMetaDataMap  m_metadata;
    /// MPI communicator in case we end up using MPI multiple times from
    /// Nektar++ SessionReader object.
    LibUtilities::CommSharedPtr     m_comm;

    /// Returns the total number of elements in the mesh with
    /// dimension expDim.
    NEKMESHUTILS_EXPORT unsigned int GetNumElements();
    /// Returns the total number of elements in the mesh with
    /// dimension < expDim.
    NEKMESHUTILS_EXPORT unsigned int GetNumBndryElements();
    /// Returns the total number of entities in the mesh.
    NEKMESHUTILS_EXPORT unsigned int GetNumEntities();

    NEKMESHUTILS_EXPORT void MakeOrder(int                      order,
                                       LibUtilities::PointsType distType);

    NEKMESHUTILS_EXPORT void PrintStats(std::ostream &out);
};
/// Shared pointer to a mesh.
typedef std::shared_ptr<Mesh> MeshSharedPtr;
}
}

#endif
