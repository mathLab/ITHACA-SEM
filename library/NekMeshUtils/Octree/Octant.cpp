////////////////////////////////////////////////////////////////////////////////
//
//  File: Octant.cpp
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
//  Description: octant object methods.
//
////////////////////////////////////////////////////////////////////////////////

#include <limits>

#include <NekMeshUtils/Octree/Octant.h>

using namespace std;
namespace Nektar
{
namespace NekMeshUtils
{

Octant::Octant(int i, OctantSharedPtr p, Array<OneD, NekDouble> dir) : m_id(i), m_parent(p)
{
    cout << "made: " << i << endl;

    //initialise variables to defualt states
    m_leaf = true;
    m_needToDivide = false;
    m_numValidPoints = 0;
    m_delta = pair<bool, NekDouble>(false, 0.0);
    NekDouble av = 0;
    NekDouble maxDif = 0;
    NekDouble minDif=numeric_limits<double>::max();
    m_location = eUnknown;

    //pull information from parent
    Array<OneD, NekDouble> parentloc = m_parent->GetLoc();
    m_loc = Array<OneD, NekDouble>(3);
    m_loc[0] = parentloc[0] + dir[0] * m_parent->DX() / 2.0;
    m_loc[1] = parentloc[1] + dir[1] * m_parent->DX() / 2.0;
    m_loc[2] = parentloc[2] + dir[2] * m_parent->DX() / 2.0;

    m_hd = m_parent->DX() / 2.0;
    vector<CurvaturePointSharedPtr> CurvaturePointList = m_parent->GetCPList();

    //setup complete

    //look over the curvature point list provided by the parent,
    //firstly look to see if it is in the new octant and if so
    //add it to the conserdation of the delta specification
    for(int i = 0; i<CurvaturePointList.size(); i++)
    {
        Array<OneD, NekDouble> cploc = CurvaturePointList[i]->GetLoc();
        if(!(cploc[0] > m_loc[0] + m_hd ||
             cploc[0] < m_loc[0] - m_hd ||
             cploc[1] > m_loc[1] + m_hd ||
             cploc[1] < m_loc[1] - m_hd ||
             cploc[2] > m_loc[2] + m_hd ||
             cploc[2] < m_loc[2] - m_hd ))
        {
            m_localCPIDList.push_back(CurvaturePointList[i]);
            if(CurvaturePointList[i]->IsValid())
            {
                if(CurvaturePointList[i]->GetDelta()>maxDif)
                {
                    maxDif = CurvaturePointList[i]->GetDelta();
                }

                if(CurvaturePointList[i]->GetDelta()<minDif)
                {
                    minDif = CurvaturePointList[i]->GetDelta();
                }
                av+=CurvaturePointList[i]->GetDelta();

                m_numValidPoints++;
            }
        }
    }

    //if it has valid points delta can be calculated
    if(NumValidCurvePoint()>0)
    {
        //geometrically octant should subdivide
        if(maxDif/minDif>1.1)
        {
            m_needToDivide=true;
        }

        SetDelta(minDif);

        if(GetDelta() < m_hd*10)
            m_needToDivide = true;

        ASSERTL0(GetDelta()>0, "negative delta assignment");

    }

    if(NumCurvePoint()>0)
    {
        m_location = eOnBoundary;
    }
}

//constructor for the master octant
Octant::Octant(int i, NekDouble x, NekDouble y, NekDouble z, NekDouble dx,
               const vector<CurvaturePointSharedPtr> &cplist)
               : m_id(i), m_hd(dx)
{
    cout << "made: " << i << endl;

    //initialise variables to defualt states
    m_leaf = false;
    m_needToDivide = true;
    m_numValidPoints = 0;
    m_delta = pair<bool, NekDouble>(false, 0.0);

    m_loc = Array<OneD, NekDouble>(3);
    m_loc[0] = x;
    m_loc[1] = y;
    m_loc[2] = z;

    m_localCPIDList = cplist;

    for(int i = 0; i < m_localCPIDList.size(); i++)
    {
        if(m_localCPIDList[i]->IsValid())
        {
            m_numValidPoints++;
        }
    }

    m_location = eOnBoundary;
}

void Octant::Subdivide(OctantSharedPtr p, NekDouble minDelta, int &numoct)
{
    if(!p->m_needToDivide) return;

    m_leaf = false; //set as not leaf and make children

    Array<OneD, OctantSharedPtr> children(8);

    for(int i = 0; i < 8; i++)
    {
        //set up x,y,z ordering of the 8 octants
        Array<OneD, NekDouble> dir(3);
        if(i<4)
        {
            dir[2] = +1.0;
            if(i<2)
            {
                dir[0] = +1.0;
            }
            else
            {
                dir[0] = -1.0;
            }
            if(i==0||i==3)
            {
                dir[1] = +1.0;
            }
            else
            {
                dir[1] = -1.0;
            }
        }
        else
        {
            dir[2] = -1.0;
            if(i<6)
            {
                dir[0] = +1.0;
            }
            else
            {
                dir[0] = -1.0;
            }
            if(i==4||i==7)
            {
                dir[1] = +1.0;
            }
            else
            {
                dir[1] = -1.0;
            }
        }

        children[i] = boost::shared_ptr<Octant>(new Octant(numoct++, p, dir));
    }

    this->SetChildren(children);

    //need to figure out neigbours here
    for(int i = 0; i < 8; i++)
    {

    }

    for(int i = 0; i < 8; i++)
    {
        if(children[i]->DX() / 4.0 > minDelta)
        {
            children[i]->Subdivide(children[i], minDelta, numoct);
        }
    }
}

bool operator==(OctantSharedPtr const &p1, OctantSharedPtr const &p2)
{
    if(p1->GetId() == p2->GetId())
        return true;
    else
        return false;
}


}
}
