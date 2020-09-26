///////////////////////////////////////////////////////////////////////////////
//
// File Octree.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Implementation of an octree sorting algorithm.
//
///////////////////////////////////////////////////////////////////////////////

#include "Octree.h"
#include <stdexcept>
#include <math.h>

namespace Nektar
{
namespace FieldUtils
{

/**
 * @brief Construct a new octree object
 *
 */
Octree::Octree() : m_maxPts(0), m_nMshPts(0), m_nNodes(0), m_nLeaves(0),
                   m_maxDepth(0), m_root(nullptr)
{
}

/**
 * @brief Construct a new octree object
 *
 * @param pts
 * @param maxPts
 * @param bounds
 */
Octree::Octree(const Array<OneD, Array<OneD, NekDouble> > &pts, int maxPts,
               const Array<OneD, NekDouble> &bounds)
{
    // Set some values
    m_maxPts  = maxPts;
    m_nMshPts = pts.size();

    // Create first (root) node
    std::vector<int> indices(m_nMshPts);
    for (int i = 0; i < m_nMshPts; ++i)
    {
        indices[i] = i;
    }
    m_root = std::make_shared<Octant>(0, 1, 0, bounds);
    m_root->SetIndices(indices);
    m_nodes.push_back(m_root);

    // Create the tree
    m_root->Subdivide(m_maxPts, pts, m_nodes);

    // Get some data after the tree is computed
    m_maxDepth = 0;
    AdvanceToStats(m_root->GetID());
    m_nNodes = m_nodes.size();

    // Set the pointers to the neighbouring nodes
    SetNeighbours(m_root->GetID());
}

/**
 * @brief Construct a new octree object
 *
 * @param pts
 * @param maxPts
 */
Octree::Octree(const Array<OneD, Array<OneD, NekDouble> > &pts, int maxPts)
{
    // Small margin to avoid rounding errors
    NekDouble margin = 1e-10;

    // Find coordinates of the bounding box
    Array<OneD, NekDouble> bounds(6);
    bounds[0] = pts[0][0];
    bounds[1] = pts[0][0];
    bounds[2] = pts[0][1];
    bounds[3] = pts[0][1];
    bounds[4] = pts[0][2];
    bounds[5] = pts[0][2];
    for (int i = 1; i < m_nMshPts; ++i)
    {
        bounds[0] = (bounds[0] < pts[i][0]) ? bounds[0] : pts[i][0];
        bounds[1] = (bounds[1] > pts[i][0]) ? bounds[1] : pts[i][0];
        bounds[2] = (bounds[2] < pts[i][1]) ? bounds[2] : pts[i][1];
        bounds[3] = (bounds[3] > pts[i][1]) ? bounds[3] : pts[i][1];
        bounds[4] = (bounds[4] < pts[i][2]) ? bounds[4] : pts[i][2];
        bounds[5] = (bounds[5] > pts[i][2]) ? bounds[5] : pts[i][2];
    }

    // Add the margin
    for (int i = 0; i < 6; ++i)
    {
        bounds[i] -= pow(-1,i) * margin;
    }

    // Call the octree constructor
    Octree(pts, maxPts, bounds);
}

/**
 * @brief Given the coordinates 'coords' of a point, returns the leaf octant
 * that contains it. If 'depth' is specified, returs the node that contains the
 * point such that its depth is lower or equal to the one indicated. If the
 * point lies outside the tree, it returns -1
 *
 * @param coords
 * @param depth
 * @return int
 */
int Octree::QueryNode(const Array<OneD, NekDouble> &coords, int depth)
{
    int nodeID = -1;

    if (coords.size())
    {
        // First, check if inside the octree
        Array<OneD, NekDouble> bounds = m_root->GetBounds();
        if ((coords[0] >= bounds[0]) && (coords[0] <= bounds[1]) &&
            (coords[1] >= bounds[2]) && (coords[1] <= bounds[3]) &&
            (coords[2] >= bounds[4]) && (coords[2] <= bounds[5]))
        {
            // Initialise 'node'
            OctantSharedPtr node = m_root;

            // Keep advancing towards the end of the branch
            while (!node->IsLeaf() && node->GetDepth() < depth)
            {
                int loc = node->GetLocInNode(coords);
                node = node->GetChildren()[loc-1];
            }

            nodeID = node->GetID();
        }
    }

    return nodeID;
}

/**
 * @brief Finds the ID of the closest point in 'pts' to the one specified by
 * 'coords'. It also returns the distance between both points in 'distance'
 *
 * @param pts
 * @param coords
 * @param distance
 * @param pointInd
 * @return int
 */
int Octree::QueryClosest(const Array<OneD, Array<OneD, NekDouble> > &pts,
                         const Array<OneD, NekDouble> &coords,
                         NekDouble &distance, int pointInd)
{
    int index = -1;
    distance  = std::numeric_limits<NekDouble>::max();

    if (coords.size())
    {
        // Find the corresponding node to 'coords'
        int nodeInd;
        if (pointInd > 0)
        {
            nodeInd = pointInd;
        }
        else
        {
            nodeInd = QueryNode(coords);
        }

        // List the indices of all the candidate points
        std::vector<int> indices(m_nodes[nodeInd]->GetIndices());
        for (OctantWeakPtr neigh : m_nodes[nodeInd]->GetNeighbours())
        {
            for (int i : neigh.lock()->GetIndices())
            {
                indices.push_back(i);
            }
        }

        // Check the distances with all the nodes
        for (int i : indices)
        {
            NekDouble sub = pts[i][0]-coords[0];
            NekDouble tmpDistance = sub * sub;
            for (int j = 1; j < 3; ++j)
            {
                sub = pts[i][j]-coords[j];
                tmpDistance += sub * sub;
            }
            tmpDistance = std::sqrt(tmpDistance);

            if (distance > tmpDistance)
            {
                distance = tmpDistance;
                index = i;
            }
        }
    }

    return index;
}

/**
 * @brief Returns the indices of the points of the mesh contained in the tree
 *
 * @param nodeID
 * @return std::vector<int>
 */
std::vector<int> Octree::QueryPoints(int nodeID)
{
    return m_nodes[nodeID]->GetIndices();
}

/**
 * @brief Returns the IDs of the octants that surround the queried node. First,
 * it finds the neighbouring nodes with the same or lower depth and, if they
 * are not leaf nodes, return all the leaf octants contained in them. This
 * means that, for octants of depth 2, this function returns all the leaf nodes
 * in the tree except those lying inside the queried octant
 *
 * @param nodeID
 * @return std::vector<int>
 */
std::vector<int> Octree::QueryNeighbours(int nodeID)
{
    std::vector<int> indices;
    for (const OctantWeakPtr &node : m_nodes[nodeID]->GetNeighbours())
    {
        indices.push_back(node.lock()->GetID());
    }
    return indices;
}

/**
 * @brief Returns some characteristic values of the tree.
 *
 * @param maxPts
 * @param nPts
 * @param nNodes
 * @param nLeaves
 * @param depth
 */
void Octree::GetStats(int &maxPts, int &nPts, int &nNodes,
                      int &nLeaves, int &depth)
{
    maxPts  = m_maxPts;
    nPts    = m_nMshPts;
    nNodes  = m_nNodes;
    nLeaves = m_nLeaves;
    depth   = m_maxDepth;
}

/**
 * @brief Goes through all the nodes of the octree counting the number of
 * octants and the maximum depth reached
 *
 * @param nodeID
 */
void Octree::AdvanceToStats(int nodeID)
{
    // Update stats if we reached the end of the branch
    if (m_nodes[nodeID]->IsLeaf())
    {
        m_nLeaves++;
        m_maxDepth = (m_maxDepth > m_nodes[nodeID]->GetDepth()) ? m_maxDepth :
                                                m_nodes[nodeID]->GetDepth();
    }
    // In any other case, dig into the tree
    else
    {
        Array<OneD, OctantSharedPtr> children = m_nodes[nodeID]->GetChildren();
        for (OctantSharedPtr child : children)
        {
            AdvanceToStats(child->GetID());
        }
    }
}

/**
 * @brief Once the nodes of the octree are created, sets their neighbours as
 * explained in 'Octree::QueryNeighbours'
 *
 * @param nodeID
 */
void Octree::SetNeighbours(int nodeID)
{
    // Array with the different steps
    static int steps[26][3] = {{-1,-1,-1}, {1,0,0}, {1,0,0}, {0,1,0}, {0,1,0},
                               {-1,0,0}, {-1,0,0}, {0,-1,0}, {1,0,0},
                               {-1,-1,1}, {1,0,0}, {1,0,0}, {0,1,0}, {0,1,0},
                               {-1,0,0}, {-1,0,0}, {0,-1,0}, {0,-1,1}, {1,0,0},
                               {1,0,0}, {0,1,0}, {0,1,0}, {-1,0,0}, {-1,0,0},
                               {0,-1,0}, {1,0,0}};

    // Advance to the leaves of the octree
    if (!m_nodes[nodeID]->IsLeaf())
    {
        for (OctantSharedPtr child : m_nodes[nodeID]->GetChildren())
        {
            SetNeighbours(child->GetID());
        }
    }
    else
    {
        // delta * steps
        Array<OneD, NekDouble> probeCoords(m_nodes[nodeID]->GetCentre());
        for (int step = 0; step < 26; ++step)
        {
            for (int i = 0; i < 3; ++i)
            {
                probeCoords[i] += m_nodes[nodeID]->GetDelta()*steps[step][i];
            }

            // For each neighbour, find the leaves and add them
            int neighInd = QueryNode(probeCoords, m_nodes[nodeID]->GetDepth());
            if (neighInd > -1)
            {
                std::vector<OctantSharedPtr> leaves;
                m_nodes[neighInd]->GetLeaves(leaves);
                m_nodes[nodeID]->AddNeighbours(leaves);
            }
        }
    }
}

/**
 * @brief Construct a new Octree::Octant object
 *
 */
Octree::Octant::Octant() : m_nPts(-1), m_loc(-1), m_depth(-1), m_id(-1),
                        m_delta(-1), m_centre(3), m_bounds(6), m_isLeaf(true)
{
}

/**
 * @brief Construct a new Octree::Octant object
 *
 * @param loc
 * @param depth
 * @param id
 * @param bounds
 */
Octree::Octant::Octant(int loc, int depth, int id,
                       const Array<OneD, NekDouble> &bounds) :
                            m_nPts(0), m_loc(loc), m_depth(depth),
                            m_id(id), m_isLeaf(true)
{
    // Check the size of 'bounds'
    if (bounds.size() != 6)
    {
        throw std::out_of_range("Size of bounds must be 6.");
    }

    // If all deltas are not equal, use the largest ones
    NekDouble deltaX = bounds[1] - bounds[0];
    NekDouble deltaY = bounds[3] - bounds[2];
    NekDouble deltaZ = bounds[5] - bounds[4];
    if (deltaX != deltaY || deltaY != deltaZ)
    {
        m_delta = (deltaX > deltaY) ? deltaX : deltaY;
        m_delta = (m_delta > deltaZ) ? m_delta : deltaZ;
    }
    else
    {
        m_delta = deltaX;
    }

    // Fill in the rest of the data
    m_centre = Array<OneD, NekDouble>(3);
    m_bounds = Array<OneD, NekDouble>(6);
    for (int i = 0; i < 3; ++i)
    {
        m_centre[i]     = (bounds[2*i+1] + bounds[2*i])/2.0;
        m_bounds[2*i]   = m_centre[i] - m_delta/2.0;
        m_bounds[2*i+1] = m_centre[i] + m_delta/2.0;
    }
}

/**
 * @brief Construct a new Octree::Octant object
 *
 * @param loc
 * @param parent
 */
Octree::Octant::Octant(int loc, Octant &parent) : m_nPts(0), m_loc(loc),
                                                  m_id(-1), m_isLeaf(true)
{
    // Set depth
    m_depth = parent.GetDepth() + 1;

    // Set delta
    m_delta = parent.GetDelta()/2.0;

    // Set centre
    NekDouble centreDX;
    NekDouble centreDY;
    NekDouble centreDZ;
    switch (loc)
    {
        case 1:  // x-, y-, z-
            centreDX = -m_delta/2.0;
            centreDY = -m_delta/2.0;
            centreDZ = -m_delta/2.0;
            break;
        case 2:  // x+, y-, z-
            centreDX =  m_delta/2.0;
            centreDY = -m_delta/2.0;
            centreDZ = -m_delta/2.0;
            break;
        case 3:  // x+, y+, z-
            centreDX =  m_delta/2.0;
            centreDY =  m_delta/2.0;
            centreDZ = -m_delta/2.0;
            break;
        case 4:  // x-, y+, z-
            centreDX = -m_delta/2.0;
            centreDY =  m_delta/2.0;
            centreDZ = -m_delta/2.0;
            break;
        case 5:  // x-, y-, z+
            centreDX = -m_delta/2.0;
            centreDY = -m_delta/2.0;
            centreDZ =  m_delta/2.0;
            break;
        case 6:  // x+, y-, z+
            centreDX =  m_delta/2.0;
            centreDY = -m_delta/2.0;
            centreDZ =  m_delta/2.0;
            break;
        case 7:  // x+, y+, z+
            centreDX =  m_delta/2.0;
            centreDY =  m_delta/2.0;
            centreDZ =  m_delta/2.0;
            break;
        case 8:  // x-, y+, z+
            centreDX = -m_delta/2.0;
            centreDY =  m_delta/2.0;
            centreDZ =  m_delta/2.0;
            break;
        default:
            throw std::out_of_range("Loc must be in the range (1,8).");
    }
    Array<OneD, NekDouble> pCentre = parent.GetCentre();
    m_centre = Array<OneD, NekDouble>(3);
    m_centre[0] = pCentre[0] + centreDX;
    m_centre[1] = pCentre[1] + centreDY;
    m_centre[2] = pCentre[2] + centreDZ;

    // Set bounds
    m_bounds = Array<OneD, NekDouble>(6);
    for (int i = 0; i < 3; ++i)
    {
        m_bounds[2*i]   = m_centre[i] - m_delta/2.0;
        m_bounds[2*i+1] = m_centre[i] + m_delta/2.0;
    }
}

/**
 * @brief Updates 'leaves' so that it contains all the leaf nodes belonging to
 * the octant
 *
 * @param leaves
 */
void Octree::Octant::GetLeaves(std::vector<OctantSharedPtr>& leaves)
{
    if (m_isLeaf)
    {
        leaves.push_back(shared_from_this());
    }
    else
    {
        for (OctantSharedPtr child : m_children)
        {
            child->GetLeaves(leaves);
        }
    }
}

/**
 * @brief Sets the values of 'm_pointInd' to those in 'indices'
 *
 * @param indices
 */
void Octree::Octant::SetIndices(const std::vector<int> &indices)
{
    for (int i : indices)
    {
        m_pointInd.push_back(i);
    }
    m_nPts = indices.size();
}

/**
 * @brief Adds to 'm_neighbours' the octants that are not already in the list
 *
 * @param neighbours
 */
void Octree::Octant::AddNeighbours(
        const std::vector<OctantSharedPtr> &neighbours)
{
    for (const OctantSharedPtr &neighbour : neighbours)
    {
        bool equal = false;
        for (const OctantWeakPtr &neigh: m_neighbours)
        {
            if (neigh.lock()->GetID() == neighbour->GetID())
            {
                equal = true;
                break;
            }
        }
        if (!equal)
        {
            m_neighbours.push_back(neighbour);
        }
    }
}

/**
 * @brief Adds to 'm_pointInd' the IDs of the points in 'pts' that fall inside
 * the octant
 *
 * @param pts
 * @param indices
 */
void Octree::Octant::AddPoints(const Array<OneD, Array<OneD, NekDouble> > &pts,
                               const std::vector<int> &indices)
{
    for (int i : indices)
    {
        // Check if the point is inside the node
        Array<OneD, NekDouble> pt = pts[i];
        if ((pt[0] < m_bounds[0]) || (pt[0] > m_bounds[1]))
        {
            continue;
        }
        if ((pt[1] < m_bounds[2]) || (pt[1] > m_bounds[3]))
        {
            continue;
        }
        if ((pt[2] < m_bounds[4]) || (pt[2] > m_bounds[5]))
        {
            continue;
        }

        // If so, add it to the list
        m_nPts++;
        m_pointInd.push_back(i);

        // Flag it as a leaf node
        m_isLeaf = true;
    }
}

/**
 * @brief Recursively divides the octant into 8 children and fills the leaf
 * nodes with their corresponding points. Does NOT add neighbours
 *
 * @param maxPts
 * @param pts
 * @param nodes
 */
void Octree::Octant::Subdivide(int maxPts,
                               const Array<OneD, Array<OneD, NekDouble> > &pts,
                               std::vector<OctantSharedPtr> &nodes)
{
    // For a non-leaf node
    if (m_nPts > maxPts)
    {
        // Create and fill children RECURSIVELY
        m_children = Array<OneD, OctantSharedPtr>(8);
        for (int i = 0; i < 8; ++i)
        {
            OctantSharedPtr newChild =
                std::make_shared<Octant>(i+1, *shared_from_this());
            newChild->AddPoints(pts, m_pointInd);
            newChild->SetID(nodes.size());  // ID's start from 0

            // Add it to the list
            m_children[i] = newChild;
            nodes.push_back(newChild);

            // Keep dividing
            newChild->Subdivide(maxPts, pts, nodes);  // Recursion
        }

        // Not a leaf node anymore
        m_pointInd.clear();
        m_nPts   = 0;
        m_isLeaf = false;
    }
}

/**
 * @brief Returns the position inside an octant in the range (1-8). The name
 * convention is as follows: let \f$\Delta\f$ be a value between 0 and the
 * length of the side of the father octant, and \f$x_c,y_c,z_c\f$ be the
 * coordinates of the centre of the father node. Then, position 1 corresponds
 * to \f$x=x_c-\Delta\f$, \f$y=y_c-\Delta\f$ and \f$y=y_c-\Delta\f$. The next
 * positions are obtained by rotating counter-clockwise around the Z axis and
 * then, making the same rotation for \f$z=z_c+\Delta\f$
 *
 * @param coords
 * @return int
 */
int Octree::Octant::GetLocInNode(const Array<OneD, NekDouble> &coords)
{
    // Different positions as bits in 'posByte'
    // MSB <==> LSB
    unsigned char posByte;

    if (coords[0] <= m_centre[0])  // x-
    {
        posByte = 153;   //0b10011001;
    }
    else                           // x+
    {
        posByte = 102;   //0b01100110;
    }
    if (coords[1] <= m_centre[1])  // y-
    {
        posByte &= 51;   //0b00110011;
    }
    else                           // y+
    {
        posByte &= 204;  //0b11001100;
    }
    if (coords[2] <= m_centre[2])  // z-
    {
        posByte &= 15;   //0b00001111;
    }
    else                           // z+
    {
        posByte &= 240;  //0b11110000;
    }

    // Transform into a position in the range (1,8)
    int position = 1;
    while (posByte > 1)
    {
        posByte = posByte >> 1;
        position++;
    }

    return position;
}
}
}
