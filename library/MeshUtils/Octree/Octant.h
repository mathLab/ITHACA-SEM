////////////////////////////////////////////////////////////////////////////////
//
//  File: Octree.h
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
//  Description: cad object methods.
//
////////////////////////////////////////////////////////////////////////////////


#ifndef NEKTAR_MESHUTILS_OCTREE_OCTANT_H
#define NEKTAR_MESHUTILS_OCTREE_OCTANT_H

#include <MeshUtils/Octree/CurvaturePoint.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>


namespace Nektar {
namespace MeshUtils {

class Octant; //have to forward declare the class for the sharedptr
typedef boost::shared_ptr<Octant> OctantSharedPtr;

class Octant
{
    public:
        friend class MemoryManager<Octant>;

        Octant(NekDouble x, NekDouble y, NekDouble z, NekDouble dx,
               int p, int l,
               const std::vector<CurvaturePointSharedPtr> &CurvaturePointList);

        void CreateNeighbourList(const std::vector<OctantSharedPtr> &OctantList);

        bool Divide(){return m_needToDivide;}
        int GetChild(int i){return m_children[i];}
        bool GetLeaf(){return m_leaf;}
        void SetLeaf(bool l){m_leaf = l;}
        bool hasPoints()
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
        NekDouble GetDelta(){return m_delta;}
        void SetDelta(NekDouble d)
        {
            m_delta = d;
            m_deltaSet=true;
        }
        bool GetDeltaKnown(){return m_deltaSet;}
        int GetOrient(){return m_orientation;}
        void SetOrient(int i)
        {
            m_orientation = i;
            m_orientSet = true;
        }
        bool GetOrientKnown(){return m_orientSet;}
        std::vector<int> GetNeighbourList(){return m_neighbourList;}
        Array<OneD, NekDouble> GetLoc(){return m_loc;}
        NekDouble DX(){return m_hd;}
        NekDouble FX(NekDouble dir){return m_loc[0]+dir*m_hd;}
        NekDouble FY(NekDouble dir){return m_loc[1]+dir*m_hd;}
        NekDouble FZ(NekDouble dir){return m_loc[2]+dir*m_hd;}
        int NumCurvePoint(){return m_localCPIDList.size();}
        std::vector<CurvaturePointSharedPtr> GetCPList(){return m_localCPIDList;}
        int NumValidCurvePoint(){return m_numValidPoints;}
        int GetLevel(){return m_level;}
        void SetChildren(Array<OneD, int> i){m_children = i;}
        void DeleteNeighbourList(){m_neighbourList.clear();}
        int GetParent(){return m_parent;}

        NekDouble Distance(const OctantSharedPtr &oct)
        {
            Array<OneD, NekDouble> octloc = oct->GetLoc();
            NekDouble r = sqrt((m_loc[0]-octloc[0])*(m_loc[0]-octloc[0])+
                               (m_loc[1]-octloc[1])*(m_loc[1]-octloc[1])+
                               (m_loc[2]-octloc[2])*(m_loc[2]-octloc[2]));
            return r;
        }

        NekDouble CPDistance(const CurvaturePointSharedPtr &cu)
        {
            NekDouble r = sqrt((m_loc[0]-cu->X())*(m_loc[0]-cu->X())+
                               (m_loc[1]-cu->Y())*(m_loc[1]-cu->Y())+
                               (m_loc[2]-cu->Z())*(m_loc[2]-cu->Z()));
            return r;
        }

        NekDouble DiagonalDim()
        {
            return sqrt(3.0*m_hd*m_hd);
        }



    private:

        bool m_leaf; //assume leaf
        int m_parent;
        Array<OneD, int> m_children;
        int m_level;
        Array<OneD, NekDouble> m_loc;
        NekDouble m_hd;

        std::vector<CurvaturePointSharedPtr> m_localCPIDList;
        NekDouble m_delta;
        std::vector<int> m_neighbourList;
        bool m_needToDivide; //asume no need to divide
        bool m_deltaSet; //will not know delta
        bool m_orientSet; //does not know orient
        int m_orientation; //1 is in 2 is partial (haspoints) 3 is out
        int m_numValidPoints;

};

}
}


#endif
