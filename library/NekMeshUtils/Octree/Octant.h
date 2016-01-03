////////////////////////////////////////////////////////////////////////////////
//
//  File: Octant.h
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
//  Description: class of individal octree octants
//
////////////////////////////////////////////////////////////////////////////////


#ifndef NEKTAR_MESHUTILS_OCTREE_OCTANT_H
#define NEKTAR_MESHUTILS_OCTREE_OCTANT_H

#include <NekMeshUtils/Octree/CurvaturePoint.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

#include <boost/unordered_set.hpp>

namespace Nektar
{
namespace NekMeshUtils
{

enum OctantFace
{
    eUp,
    eDown,
    eForward,
    eBack,
    eLeft,
    eRight
};

enum OctantLocation
{
    eInside,
    eOutside,
    eOnBoundary,
    eUnknown
};

class Octant; //have to forward declare the class for the sharedptr
typedef boost::shared_ptr<Octant> OctantSharedPtr;
typedef std::set<OctantSharedPtr> OctantSet;

/**
 * @brief this class contains the infomration and methods for individal octants
 * in the Octree
 */
class Octant
{
    public:
        friend class MemoryManager<Octant>;

        /**
         * @brief Defualt constructor
         */
        Octant(int i, OctantSharedPtr p, Array<OneD, OctantFace> dir);

        //master constructor
        Octant(int i, NekDouble x, NekDouble y, NekDouble z, NekDouble dx,
                       const std::vector<CurvaturePointSharedPtr> &cplist);

        void Subdivide(OctantSharedPtr p, int &numoct);

        void CompileLeaves(std::vector<OctantSharedPtr> &Octants)
        {
            for(int i = 0; i < 8; i++)
            {
                if(m_children[i]->IsLeaf())
                {
                    Octants.push_back(m_children[i]);
                }
                else
                {
                    m_children[i]->CompileLeaves(Octants);
                }
            }
        }

        int GetId()
        {
            return m_id;
        }

        Array<OneD, NekDouble> GetLoc()
        {
            return m_loc;
        }

        NekDouble DX()
        {
            return m_hd;
        }

        std::vector<CurvaturePointSharedPtr> GetCPList()
        {
            return m_localCPIDList;
        }

        int NumCurvePoint()
        {
            return m_localCPIDList.size();
        }

        int NumValidCurvePoint()
        {
            return m_numValidPoints;
        }

        void SetDelta(NekDouble d)
        {
            m_delta.first = true;
            m_delta.second = d;
        }

        NekDouble GetDelta()
        {
            ASSERTL0(m_delta.first, "Tried to acsess delta of octant"
                                    "which has not been set");
            return m_delta.second;
        }

        void SetChildren(Array<OneD, OctantSharedPtr> c)
        {
            m_children = c;
        }

        bool IsLeaf()
        {
            return m_leaf;
        }

        NekDouble FX(OctantFace f)
        {
            switch (f)
            {
                case eUp:
                    return m_loc[1] + m_hd;
                    break;
                case eDown:
                    return m_loc[1] - m_hd;
                    break;
                case eForward:
                    return m_loc[0] + m_hd;
                    break;
                case eBack:
                    return m_loc[0] - m_hd;
                    break;
                case eLeft:
                    return m_loc[2] + m_hd;
                    break;
                case eRight:
                    return m_loc[2] - m_hd;
                    break;
            }
        }

        void RemoveNeigbour(int id, OctantFace f);

        void SetNeigbour(OctantSharedPtr o, OctantFace f)
        {
            m_neigbours[f].push_back(o);
        }

        std::map<OctantFace, std::vector<OctantSharedPtr> > GetNeigbours()
        {
            return m_neigbours;
        }

        bool NeedDivide()
        {
            return m_needToDivide;
        }

        NekDouble Distance(OctantSharedPtr o)
        {
            Array<OneD, NekDouble> loc = o->GetLoc();
            return sqrt((loc[0] - m_loc[0])*(loc[0] - m_loc[0]) +
                        (loc[1] - m_loc[1])*(loc[1] - m_loc[1]) +
                        (loc[2] - m_loc[2])*(loc[2] - m_loc[2]));
        }

        bool IsDeltaKnown()
        {
            return m_delta.first;
        }

        void SetLocation(OctantLocation l)
        {
            m_location = l;
        }

        CurvaturePointSharedPtr GetCPPoint()
        {
            ASSERTL0(m_localCPIDList.size() > 0, "tried to get cp point where there is none");
            return m_localCPIDList[0];
        }

        OctantLocation GetLocation()
        {
            return m_location;
        }

    private:

        ///id
        int m_id;
        /// leaf identifer
        bool m_leaf;
        /// parent id
        OctantSharedPtr m_parent;
        /// list of child ids
        Array<OneD, OctantSharedPtr> m_children;
        /// x,y,z location of the octant
        Array<OneD, NekDouble> m_loc;
        /// half dimension of the octant
        NekDouble m_hd;
        /// curvature sampling point list
        std::vector<CurvaturePointSharedPtr> m_localCPIDList;
        int m_numValidPoints;
        /// mesh sizing parameter
        std::pair<bool, NekDouble> m_delta;
        /// idenify if division is needed
        bool m_needToDivide; //asume no need to divide
        /// idenify if delta has ben set
        OctantLocation m_location;
        /// list of neighbours
        std::map<OctantFace, std::vector<OctantSharedPtr> > m_neigbours;
};

bool operator==(OctantSharedPtr const &p1, OctantSharedPtr const &p2);

}
}


#endif
