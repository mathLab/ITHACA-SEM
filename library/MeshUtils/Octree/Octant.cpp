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

#include <MeshUtils/Octree/Octant.h>

using namespace std;
namespace Nektar
{
namespace MeshUtils
{

Octant::Octant(OctantSharedPtr p, Array<OneD, NekDouble> dir) : m_parent(p)
{
    //initialise variables to defualt states
    m_leaf = true;
    m_needToDivide = false;
    m_deltaSet = false;
    m_numValidPoints = 0;
    m_delta = -1;
    NekDouble av=0;
    NekDouble maxDif=0;
    NekDouble minDif=numeric_limits<double>::max();

    //pull information from parent
    Array<OneD, NekDouble> parentloc = m_parent->GetLoc();
    m_loc = Array<OneD, NekDouble>(3);
    m_loc[0] = parentloc[0] + dir[0] * m_parent->DX() / 2.0;
    m_loc[1] = parentloc[1] + dir[1] * m_parent->DX() / 2.0;
    m_loc[2] = parentloc[2] + dir[2] * m_parent->DX() / 2.0;

    m_hd = m_parent->DX() / 2.0;
    m_level = m_parent->GetLevel() + 1;
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

        SetDelta(av/NumValidCurvePoint());

        ASSERTL0(GetDelta()>0, "negative delta assignment");

    }

    if(NumCurvePoint()>0)
    {
        //location should be on the boundary
        m_locationKnown = true;
        m_location = 2;
    }
    else
    {
        m_location = false;
    }
}

//constructor for the master octant
Octant::Octant(NekDouble x, NekDouble y, NekDouble z, NekDouble dx,
               const vector<CurvaturePointSharedPtr> &cplist)
               : m_hd(dx)
{
    //initialise variables to defualt states
    m_leaf = false;
    m_needToDivide = true;
    m_deltaSet = true;
    m_numValidPoints = 0;
    m_delta = -1;

    m_loc = Array<OneD, NekDouble>(3);
    m_loc[0] = x;
    m_loc[1] = y;
    m_loc[2] = z;

    m_level = 0;
    m_localCPIDList = cplist;

    for(int i = 0; i < m_localCPIDList.size(); i++)
    {
        if(m_localCPIDList[i]->IsValid())
        {
            m_numValidPoints++;
        }
    }

    m_location = 2;
}

void Octant::CreateNeighbourList(OctantSet Octants)
{
    //clear old list
    DeleteNeighbourList();

    OctantSet::iterator it;

    //look over all octants and consider if they are neigbours
    for(it = Octants.begin(); it != Octants.end(); it++)
    {
        OctantSharedPtr oct = *it;

        if(oct == boost::shared_ptr<Octant>(this))
            continue;

        if(oct->IsLeaf())
        {
            //work out the max distance between the two octants if they were
            //joined corner to corner
            NekDouble rmax = DiagonalDim() + oct->DiagonalDim();

            NekDouble ractual = Distance(oct);

            //if the actucal distance between them is greater that this
            //this octant should not be considered, massive speed up
            if(ractual > 1.1*rmax)
            {
                continue;
            }

            //check overlapping in 2-D for all three dimensions
            if(fabs(FX(-1) - oct->FX(+1)) < 1E-6 ||
               fabs(FX(+1) - oct->FX(-1)) < 1E-6 )
            {
                //check yz rects
                bool Cond1 = false;
                bool Cond2 = false;
                bool Cond3 = false;
                bool Cond4 = false;
                if(FY(+1) < oct->FY(-1))
                    Cond1 = true;
                if(FY(-1) > oct->FY(+1))
                    Cond1 = true;
                if(FZ(+1) < oct->FZ(-1))
                    Cond1 = true;
                if(FZ(-1) > oct->FZ(+1))
                    Cond1 = true;
                if(!Cond1 && !Cond2 && !Cond3 && !Cond4)
                {
                    m_neighbourList.push_back(oct);
                }
            }
            else if(fabs(FY(-1) - oct->FY(+1)) <1E-6 ||
                    fabs(FY(+1) - oct->FY(-1)) <1E-6)
            {
                //check xz rects
                bool Cond1 = false;
                bool Cond2 = false;
                bool Cond3 = false;
                bool Cond4 = false;
                if(FX(+1) < oct->FX(-1))
                    Cond1 = true;
                if(FX(-1) > oct->FX(+1))
                    Cond1 = true;
                if(FZ(+1) < oct->FZ(-1))
                    Cond1 = true;
                if(FZ(-1) > oct->FZ(+1))
                    Cond1 = true;
                if(!Cond1 && !Cond2 && !Cond3 && !Cond4)
                {
                    m_neighbourList.push_back(oct);
                }
            }
            else if(fabs(FZ(-1) - oct->FZ(+1)) <1E-6 ||
                    fabs(FZ(+1) - oct->FZ(-1)) <1E-6)
            {
                //check xy rects
                bool Cond1 = false;
                bool Cond2 = false;
                bool Cond3 = false;
                bool Cond4 = false;
                if(FX(+1) < oct->FX(-1))
                    Cond1 = true;
                if(FX(-1) > oct->FX(+1))
                    Cond1 = true;
                if(FY(+1) < oct->FY(-1))
                    Cond1 = true;
                if(FY(-1) > oct->FY(+1))
                    Cond1 = true;
                if(!Cond1 && !Cond2 && !Cond3 && !Cond4)
                {
                    m_neighbourList.push_back(oct);
                }
            }
        }
    }
}

bool operator==(OctantSharedPtr const &p1, OctantSharedPtr const &p2)
{
    Array<OneD, NekDouble> loc1 = p1->GetLoc();
    Array<OneD, NekDouble> loc2 = p2->GetLoc();
    if(loc1[0] == loc2[0] && loc1[1] == loc2[1] && loc1[2] == loc2[2])
        return true;
    else
        return false;
}

}
}
