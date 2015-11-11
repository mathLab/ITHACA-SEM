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

#include <MeshUtils/Octree/CurvaturePoint.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

#include <boost/unordered_set.hpp>

namespace Nektar
{
namespace MeshUtils
{

class Octant; //have to forward declare the class for the sharedptr
typedef boost::shared_ptr<Octant> OctantSharedPtr;
typedef boost::unordered_set<OctantSharedPtr> OctantSet;

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
        Octant(OctantSharedPtr p, Array<OneD, NekDouble> dir);

        //master constructor
        Octant(NekDouble x, NekDouble y, NekDouble z, NekDouble dx,
                       const std::vector<CurvaturePointSharedPtr> &cplist);

        /**
         * @brief scans over all octants in the octantlist and finds neighouring
         * octants
         */
        void CreateNeighbourList(OctantSet OctantList);

        /**
         * @brief get boolean on whether the octant needs to divide based on
         * geometry
         */
        bool GetDivide(){return m_needToDivide;}

        /**
         * @brief returns the id of the ith child of 8 for a non leaf octant
         */
        OctantSharedPtr GetChild(int i){return m_children[i];}

        /**
         * @brief get boolean on whether the octant is a leaf or not
         */
        bool IsLeaf(){return m_leaf;}

        /**
         * @brief alter leaf status
         */
        void SetLeaf(bool l){m_leaf = l;}

        /**
         * @brief get boolean on whether the octant has curvature sample points
         * within its volume
         */
        bool HasPoints()
        {
            if(m_localCPIDList.size()>0)
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        /**
         * @brief return mesh sizing specifcaion of the octant
         */
        NekDouble GetDelta(){return m_delta;}

        /**
         * @brief alter mesh sizing specification
         */
        void SetDelta(NekDouble d)
        {
            m_delta = d;
            m_deltaSet=true;
            ASSERTL0(m_delta > 0.0, "delta assignement less than 0");
        }

        /**
         * @brief get boolean on whether the octant has had a mesh spacing
         * assigned
         */
        bool IsDeltaKnown(){return m_deltaSet;}

        /**
         * @brief get the neigbour list of the octant
         */
        std::vector<OctantSharedPtr> GetNeighbourList(){return m_neighbourList;}

        /**
         * @brief get the location of the octant
         */
        Array<OneD, NekDouble> GetLoc(){return m_loc;}

        /**
         * @brief get the half dimension of the octant
         */
        NekDouble DX(){return m_hd;}

        /**
         * @brief get the far x coordiate of the octant volume
         * backwards or forards depending on dir (should be -1 or 1)
         */
        NekDouble FX(NekDouble dir){return m_loc[0]+dir*m_hd;}

        /**
         * @brief get the far y coordiate of the octant volume
         * backwards or forards depending on dir (should be -1 or 1)
         */
        NekDouble FY(NekDouble dir){return m_loc[1]+dir*m_hd;}

        /**
         * @brief get the far z coordiate of the octant volume
         * backwards or forards depending on dir (should be -1 or 1)
         */
        NekDouble FZ(NekDouble dir){return m_loc[2]+dir*m_hd;}

        /**
         * @brief get the number of curvature sampling points within the octant
         */
        int NumCurvePoint(){return m_localCPIDList.size();}

        /**
         * @brief get the list of curvature sampling points within the octant
         */
        std::vector<CurvaturePointSharedPtr> GetCPList(){return m_localCPIDList;}

        /**
         * @brief get the number of valid curvature sampling points
         */
        int NumValidCurvePoint(){return m_numValidPoints;}

        /**
         * @brief level of the octant within the octree
         */
        int GetLevel(){return m_level;}

        /**
         * @brief set the list of child octants
         */
        void SetChildren(Array<OneD, OctantSharedPtr> i){m_children = i;}

        /**
         * @brief clear neigbour list
         */
        void DeleteNeighbourList(){m_neighbourList.clear();}

        /**
         * @brief get the parent of this octant
         */
        OctantSharedPtr GetParent(){return m_parent;}

        /**
         * @brief get the distance between this octant and another
         */
        NekDouble Distance(const OctantSharedPtr &oct)
        {
            Array<OneD, NekDouble> octloc = oct->GetLoc();
            NekDouble r = sqrt((m_loc[0]-octloc[0])*(m_loc[0]-octloc[0])+
                               (m_loc[1]-octloc[1])*(m_loc[1]-octloc[1])+
                               (m_loc[2]-octloc[2])*(m_loc[2]-octloc[2]));
            return r;
        }

        /**
         * @brief get the distance between this octant and a curavture sampling
         * point
         */
        NekDouble CPDistance(const CurvaturePointSharedPtr &cu)
        {
            Array<OneD, NekDouble> cploc = cu->GetLoc();
            NekDouble r = sqrt((m_loc[0]-cploc[0])*(m_loc[0]-cploc[0])+
                               (m_loc[1]-cploc[0])*(m_loc[1]-cploc[0])+
                               (m_loc[2]-cploc[0])*(m_loc[2]-cploc[0]));
            return r;
        }

        /**
         * @brief get the diagnal length of the octant (from corner to center)
         */
        NekDouble DiagonalDim()
        {
            return sqrt(3.0*m_hd*m_hd);
        }

        int GetLocation(){return m_location;}
        bool KnowsLocation(){return m_locationKnown;}

        void SetLocation(int l)
        {
            m_locationKnown = true;
            m_location = l;
        }

    private:

        /// leaf identifer
        bool m_leaf;
        /// parent id
        OctantSharedPtr m_parent;
        /// list of child ids
        Array<OneD, OctantSharedPtr> m_children;
        /// level of the octant
        int m_level;
        /// x,y,z location of the octant
        Array<OneD, NekDouble> m_loc;
        /// half dimension of the octant
        NekDouble m_hd;
        /// curvature sampling point list
        std::vector<CurvaturePointSharedPtr> m_localCPIDList;
        /// mesh sizing parameter
        NekDouble m_delta;
        /// list of ids of neigbours
        std::vector<OctantSharedPtr> m_neighbourList;
        /// idenify if division is needed
        bool m_needToDivide; //asume no need to divide
        /// idenify if delta has ben set
        bool m_deltaSet; //will not know delta
        /// idenify if orientation has been set
        int m_numValidPoints;
        /// location with respect to the domain
        bool m_locationKnown;
        int m_location; //1 is interior 2 is boundary 3 is outside

};

bool operator==(OctantSharedPtr const &p1, OctantSharedPtr const &p2);

}
}


#endif
