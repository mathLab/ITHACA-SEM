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


#ifndef NEKTAR_LIB_UTILITIES_MESHUTILS_OCTREE_OCTANT_H
#define NEKTAR_LIB_UTILITIES_MESHUTILS_OCTREE_OCTANT_H

#include <LibUtilities/MeshUtils/CurvaturePoint.h>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>


namespace Nektar {
    namespace LibUtilities {
        namespace MeshUtils {
            
            class Octant; //have to forward declare the class for the sharedptr
            typedef boost::shared_ptr<Octant> OctantSharedPtr;
            
            class Octant
            {
            public:
                friend class MemoryManager<Octant>;
                
                Octant(NekDouble x, NekDouble y, NekDouble z,
                       NekDouble dx, NekDouble dy, NekDouble dz,
                       int p, int l,
                       const std::vector<CurvaturePointSharedPtr> &CurvaturePointList,
                       const std::vector<int> &CPList);
                
                void AddCurvaturePoint(int i, bool valid, NekDouble delta,
                                       NekDouble &maxDif, NekDouble &minDif,
                                       NekDouble &av);
                
                
                bool Divide(){return m_needToDivide;}
                int GetChild(int i){return m_children[i];}
                bool isLeaf(){return m_leaf;}
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
                bool isDeltaKnown(){return m_deltaSet;}
                int GetOrient(){return m_orientation;}
                void SetOrient(int i)
                {
                    m_orientation = i;
                    m_orientSet = true;
                }
                bool isOrientKnown(){return m_orientSet;}
                std::vector<int> GetNeighbourList(){return m_neighbourList;}
                NekDouble X(){return m_x;}
                NekDouble Y(){return m_y;}
                NekDouble Z(){return m_z;}
                NekDouble DX(){return m_dx;}
                NekDouble DY(){return m_dy;}
                NekDouble DZ(){return m_dz;}
                NekDouble FX(NekDouble dir){return m_x+dir*m_dx;}
                NekDouble FY(NekDouble dir){return m_y+dir*m_dy;}
                NekDouble FZ(NekDouble dir){return m_z+dir*m_dz;}
                int NumCurvePoint(){return m_localCPIDList.size();}
                std::vector<int> GetCPList(){return m_localCPIDList;}
                int NumValidCurvePoint(){return m_numValidPoints;}
                int GetLevel(){return m_level;}
                void SetChildren(Array<OneD, int> i){m_children = i;}
                void LeafFalse(){m_leaf = false;}
                
                
                
                
                
            private:
                
                bool m_leaf; //assume leaf
                int m_parent;
                Array<OneD, int> m_children;
                int m_level;
                NekDouble m_x;
                NekDouble m_y;
                NekDouble m_z;
                NekDouble m_dx;
                NekDouble m_dy;
                NekDouble m_dz;
                std::vector<int> m_localCPIDList;
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
}


#endif