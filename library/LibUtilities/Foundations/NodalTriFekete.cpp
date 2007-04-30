///////////////////////////////////////////////////////////////////////////////
//
// File NodalTriFekete.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
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
// Description: 2D Nodal Triangle Fekete Point Definitions
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <algorithm>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/Foundations/Foundations.hpp>

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/Polylib/Polylib.h>
#include <LibUtilities/Foundations/NodalTriFekete.h>
#include <LibUtilities/Foundations/NodalTriFeketeData.h>

namespace Nektar
{
    namespace LibUtilities 
    {
        void NodalTriFekete::CalculatePoints()
        {
            // Allocate the storage for points
            PointsBaseType::CalculatePoints();

            int index=0,isum=0;
            const int offset = 3; //offset to match Datafile
            NekDouble a,b,c;
            unsigned int numPoints = GetNumPoints();

            // initialize values
            for(unsigned int i=0; i < numPoints-2; ++i)
            {
                index += NodalTriFeketeNPTS[i];
            }

            for(unsigned int i=0; i < NodalTriFeketeNPTS[numPoints-2]; ++i, ++index)
            {
                if(int(NodalTriFeketeData[index][0]))
                {
                    a = NodalTriFeketeData[index][3]; 
                    b = NodalTriFeketeData[index][4]; 
                    c = NodalTriFeketeData[index][5]; 

                    (*m_points[0])[isum] = 2.0*b - 1.0;
                    (*m_points[1])[isum] = 2.0*c - 1.0;
                    isum++;
                    continue;
                }//end symmetry1


                if(int(NodalTriFeketeData[index][1]) == 1)
                {
                    for(unsigned int j=0; j < 3; ++j)
                    {
                        a = NodalTriFeketeData[index][offset+perm3A_2d[j][0]];
                        b = NodalTriFeketeData[index][offset+perm3A_2d[j][1]];
                        c = NodalTriFeketeData[index][offset+perm3A_2d[j][2]];
                        (*m_points[0])[isum] = 2.0*b - 1.0;
                        (*m_points[1])[isum] = 2.0*c - 1.0;
                        isum++;
                    }//end j
                    continue;
                }//end symmetry3a

                if(int(NodalTriFeketeData[index][1]) == 2)
                {
                    for(unsigned int j=0; j < 3; ++j)
                    {
                        a = NodalTriFeketeData[index][offset+perm3B_2d[j][0]];
                        b = NodalTriFeketeData[index][offset+perm3B_2d[j][1]];
                        c = NodalTriFeketeData[index][offset+perm3B_2d[j][2]];
                        (*m_points[0])[isum] = 2.0*b - 1.0;
                        (*m_points[1])[isum] = 2.0*c - 1.0;
                        isum++;
                    }//end j
                    continue;   
                }//end symmetry3b


                if(int(NodalTriFeketeData[index][2]))
                {
                    for(unsigned int j=0; j < 6; ++j)
                    {
                        a = NodalTriFeketeData[index][offset+perm6_2d[j][0]];
                        b = NodalTriFeketeData[index][offset+perm6_2d[j][1]];
                        c = NodalTriFeketeData[index][offset+perm6_2d[j][2]];
                        (*m_points[0])[isum] = 2.0*b - 1.0;
                        (*m_points[1])[isum] = 2.0*c - 1.0;
                        isum++;
                    }//end j
                    continue;   
                }//end symmetry6
            }//end npts

            NodalPointReorder2d();

            ASSERTL1((isum==m_pointsKey.GetTotNumPoints()),"sum not equal to npts");
        }

        void NodalTriFekete::CalculateWeights()
        {
            // No weights computed
        }

        void NodalTriFekete::CalculateDerivMatrix()
        {
            // No derivative matrix computed
        }

        boost::shared_ptr<NodalTriFekete::PointsBaseType> NodalTriFekete::Create(const PointsKey &key)
        {
            boost::shared_ptr<PointsBaseType> returnval(MemoryManager::AllocateSharedPtr<NodalTriFekete>(key));
            returnval->Initialize();
            return returnval;
        }

        void NodalTriFekete::NodalPointReorder2d()
        {
            int istart,iend,isum=0;
            const int numvert = 3;
            const int numepoints = GetNumPoints()-2;

            if(numepoints==0)
            {
                return;
            }

            // bubble sort for first edge
            istart = numvert + isum;
            iend = istart + numepoints;
            for(int i=istart; i<iend; ++i)
            {
                for(int j=istart; j<iend-1; ++j)
                {
                    if((*m_points[0])[j+1] < (*m_points[0])[j])
                    {
                        std::swap((*m_points[0])[j+1], (*m_points[0])[j]);
                        std::swap((*m_points[1])[j+1], (*m_points[1])[j]);
                    }
                }
            }
            isum += numepoints;

            // bubble sort for second edge
            istart = numvert + isum;
            iend = istart + numepoints;
            for(int i=istart; i<iend; ++i)
            {
                for(int j=istart;j<iend-1; ++j)
                {
                    if((*m_points[0])[j+1] > (*m_points[0])[j])
                    {
                        std::swap((*m_points[0])[j+1], (*m_points[0])[j]);
                        std::swap((*m_points[1])[j+1], (*m_points[1])[j]);
                    }
                }
            }
            isum += numepoints;

            // bubble sort for third edge
            istart = numvert + isum;
            iend = istart + numepoints;
            for(int i=istart; i<iend; ++i)
            {
                for(int j=istart; j<iend-1; ++j)
                {
                    if((*m_points[1])[j+1] > (*m_points[1])[j])
                    {
                        std::swap((*m_points[0])[j+1], (*m_points[0])[j]);
                        std::swap((*m_points[1])[j+1], (*m_points[1])[j]);
                    }
                }
            }

            return;
        }     
    } // end of namespace stdregion
} // end of namespace stdregion

/**
* $Log: NodalTriFekete.cpp,v $
* Revision 1.10  2007/04/29 00:31:57  jfrazier
* Updated to use multi_arrays.
*
*/
