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

#include <LibUtilities/MeshUtils/Octant.h>

using namespace std;
namespace Nektar {
namespace LibUtilities {
namespace MeshUtils {
    
    Octant::Octant(NekDouble x, NekDouble y, NekDouble z,
                   NekDouble dx, NekDouble dy, NekDouble dz,
                   int p, int l,
                   vector<CurvaturePointSharedPtr> &CurvaturePointList,
                   vector<OctantSharedPtr> &OctantList):
                   m_parent(p),m_level(l),
                   m_x(x),m_y(y),m_z(z),m_dx(dx),m_dy(dy),m_dz(dz)
    {
        m_leaf = true;
        m_needToDivide = false;
        m_deltaSet = false;
        m_orientSet = false;
        m_orientation = -1;
        int m_numValidPoints = 0;
        NekDouble av=0;
        NekDouble maxDif=0;
        NekDouble minDif=10000.0;
        
        if(m_parent!=-1)
        {
            vector<int> CPList = OctantList[m_parent]->GetCPList();
            for(int i = 0; i<CPList.size(); i++)
            {
                if(CurvaturePointList[CPList[i]]->X()>=FX(-1) &&
                   CurvaturePointList[CPList[i]]->X()<=FX(+1) &&
                   CurvaturePointList[CPList[i]]->Y()>=FY(-1) &&
                   CurvaturePointList[CPList[i]]->Y()<=FY(+1) &&
                   CurvaturePointList[CPList[i]]->Z()>=FZ(-1) &&
                   CurvaturePointList[CPList[i]]->Z()<=FX(+1))
                {
                    AddCurvaturePoint(CPList[i],
                                      CurvaturePointList[CPList[i]]->IsValid(),
                                      CurvaturePointList[CPList[i]]->GetDelta(),
                                      maxDif, minDif, av);
                }
            }
        }else
        {
            for(int i = 0; i<CurvaturePointList.size(); i++)
            {
                if(CurvaturePointList[i]->X() >= FX(-1) &&
                   CurvaturePointList[i]->X() <= FX(+1) &&
                   CurvaturePointList[i]->Y() >= FY(-1) &&
                   CurvaturePointList[i]->Y() <= FY(+1) &&
                   CurvaturePointList[i]->Z() >= FZ(-1) &&
                   CurvaturePointList[i]->Z() <= FZ(+1))
                {
                    AddCurvaturePoint(i, CurvaturePointList[i]->IsValid(),
                                      CurvaturePointList[i]->GetDelta(),
                                      maxDif, minDif, av);
                }
            }
        }
        
        if(NumCurvePoint()>0)
        {
            SetOrient(2);
        }
        if(NumValidCurvePoint()>0)
        {
            if(maxDif/minDif>1.1)
            {
                m_needToDivide=true;
            }
            
            SetDelta(av/NumValidCurvePoint());
            
            ASSERTL0(GetDelta()>0, "negative delta assignment");
            
        }
    }
    
    void Octant::AddCurvaturePoint(int i, bool valid, NekDouble delta,
                                   NekDouble &maxDif, NekDouble &minDif,
                                   NekDouble &av)
    {
        if(valid)
        {
            if(delta>maxDif)
            {
                maxDif = delta;
            }
            
            if(delta<minDif)
            {
                minDif = delta;
            }
            av+=delta;
            m_localCPIDList.push_back(i);
            m_numValidPoints++;
        }
        else
        {
            m_localCPIDList.push_back(i);
        }
    }
    
    
    
}
}
}