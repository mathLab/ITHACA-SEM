///////////////////////////////////////////////////////////////////////////////
//
// File NodalTetElec.cpp
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
// Description: 3D Nodal Tet Electrostatic Point Definitions
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <algorithm>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/Foundations/Foundations.hpp>

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/Polylib/Polylib.h>
#include <LibUtilities/Foundations/NodalTetElec.h>
#include <LibUtilities/Foundations/NodalTetElecData.h>

namespace Nektar
{
    namespace LibUtilities 
    {
        void NodalTetElec::CalculatePoints()
        {
            // Allocate the storage for points
            Points<double>::CalculatePoints();

            int index=0,isum=0;
            const int offset = 5; //offset to match Datafile
            NekDouble a,b,c,d;
            unsigned int numPoints = GetNumPoints();

            // initialize values
            for(unsigned int i=0; i < numPoints-2; ++i)
            {
                index += NodalTetElecNPTS[i];
            }

            for(unsigned int i=0; i < NodalTetElecNPTS[numPoints-2]; ++i, ++index)
            {
                // 1 Point Symmetry: aaaa
                if(int(NodalTetElecData[index][0]))
                {
                    a = NodalTetElecData[index][5];
                    b = NodalTetElecData[index][6];
                    c = NodalTetElecData[index][7];
                    d = NodalTetElecData[index][8];

                    m_points[0][isum] = 2.0*b - 1.0;
                    m_points[1][isum] = 2.0*c - 1.0;
                    m_points[2][isum] = 2.0*d - 1.0;
                    isum++;
                    continue;
                }//end symmetry 1

                
                // 4 Point symmetry: aaab or abbb
                if(int(NodalTetElecData[index][1]))
                {
                    for(unsigned int j=0; j < 4; ++j)
                    {
                        a = NodalTetElecData[index][offset + perm4_3d[j][0]];
                        b = NodalTetElecData[index][offset + perm4_3d[j][1]];
                        c = NodalTetElecData[index][offset + perm4_3d[j][2]];
                        d = NodalTetElecData[index][offset + perm4_3d[j][3]];
                        
                        m_points[0][isum] = 2.0*b - 1.0;
                        m_points[1][isum] = 2.0*c - 1.0;
                        m_points[2][isum] = 2.0*d - 1.0;
                        isum++;
                    }//end j
                    continue;
                }//end symmetry 4

                
                // 6 Point symmetry: aabb
                if(int(NodalTetElecData[index][2]))
                {
                    for(unsigned int j=0; j < 6; ++j)
                    {
                        a = NodalTetElecData[index][offset + perm6_3d[j][0]];
                        b = NodalTetElecData[index][offset + perm6_3d[j][1]];
                        c = NodalTetElecData[index][offset + perm6_3d[j][2]];
                        d = NodalTetElecData[index][offset + perm6_3d[j][3]];
                        
                        m_points[0][isum] = 2.0*b - 1.0;
                        m_points[1][isum] = 2.0*c - 1.0;
                        m_points[2][isum] = 2.0*d - 1.0;
                        isum++;
                    }//end j
                    continue;   
                }//end symmetry6                
                

                // 12 Point symmetry: case aabc
                if(int(NodalTetElecData[index][3]) == 1)
                {
                    for(unsigned int j=0; j < 12; ++j)
                    {
                        a = NodalTetElecData[index][offset + perm12A_3d[j][0]];
                        b = NodalTetElecData[index][offset + perm12A_3d[j][1]];
                        c = NodalTetElecData[index][offset + perm12A_3d[j][2]];
                        d = NodalTetElecData[index][offset + perm12A_3d[j][3]];
                        
                        m_points[0][isum] = 2.0*b - 1.0;
                        m_points[1][isum] = 2.0*c - 1.0;
                        m_points[2][isum] = 2.0*d - 1.0;
                        isum++;
                    }//end j
                    continue;
                }//end symmetry 12 aabc

                
                // 12 Point symmetry: case abcc
                if(int(NodalTetElecData[index][3]) == 2)
                {
                    for(unsigned int j=0; j < 12; ++j)
                    {
                        a = NodalTetElecData[index][offset + perm12B_3d[j][0]];
                        b = NodalTetElecData[index][offset + perm12B_3d[j][1]];
                        c = NodalTetElecData[index][offset + perm12B_3d[j][2]];
                        d = NodalTetElecData[index][offset + perm12B_3d[j][3]];
                        
                        m_points[0][isum] = 2.0*b - 1.0;
                        m_points[1][isum] = 2.0*c - 1.0;
                        m_points[2][isum] = 2.0*d - 1.0;
                        isum++;
                    }//end j
                    continue;
                }//end symmetry 12 abcc


                // 12 Point symmetry: case abbc
                if(int(NodalTetElecData[index][3]) == 3)
                {
                    for(unsigned int j=0; j < 12; ++j)
                    {
                        a = NodalTetElecData[index][offset + perm12C_3d[j][0]];
                        b = NodalTetElecData[index][offset + perm12C_3d[j][1]];
                        c = NodalTetElecData[index][offset + perm12C_3d[j][2]];
                        d = NodalTetElecData[index][offset + perm12C_3d[j][3]];
                        
                        m_points[0][isum] = 2.0*b - 1.0;
                        m_points[1][isum] = 2.0*c - 1.0;
                        m_points[2][isum] = 2.0*d - 1.0;
                        isum++;
                    }//end j
                    continue;
                }//end symmetry 12 abbc

                
                // 24 Point symmetry: case abcd
                if(int(NodalTetElecData[index][4]))
                {
                    for(unsigned int j=0; j < 24; ++j)
                    {
                        a = NodalTetElecData[index][offset + perm24_3d[j][0]];
                        b = NodalTetElecData[index][offset + perm24_3d[j][1]];
                        c = NodalTetElecData[index][offset + perm24_3d[j][2]];
                        d = NodalTetElecData[index][offset + perm24_3d[j][3]];
                        
                        m_points[0][isum] = 2.0*b - 1.0;
                        m_points[1][isum] = 2.0*c - 1.0;
                        m_points[2][isum] = 2.0*d - 1.0;
                        isum++;
                    }//end j
                    continue;
                }//end symmetry24abcd
                                

            }//end npts

            NodalPointReorder3d();

            ASSERTL1((isum==m_pointsKey.GetTotNumPoints()),"sum not equal to npts");            
        }

        void NodalTetElec::CalculateWeights()
        {
            // No weights computed
        }

        void NodalTetElec::CalculateDerivMatrix()
        {
            // No derivative matrix computed
        }

        boost::shared_ptr<NodalTetElec::PointsBaseType> NodalTetElec::Create(const PointsKey &key)
        {
            boost::shared_ptr<PointsBaseType> returnval(MemoryManager<NodalTetElec>::AllocateSharedPtr(key));
            returnval->Initialize();
            return returnval;
        }

        void NodalTetElec::NodalPointReorder3d()
        {
        }     

    } // end of namespace stdregion
} // end of namespace stdregion


/**
* $Log: NodalTetElec.cpp,v $
* Revision 1.9  2007/07/22 23:03:26  bnelson
* Backed out Nektar::ptr.
*
* Revision 1.8  2007/07/20 00:28:26  bnelson
* Replaced boost::shared_ptr with Nektar::ptr
*
* Revision 1.7  2007/05/15 03:37:24  bnelson
* Updated to use the new Array object.
*
* Revision 1.6  2007/04/30 23:29:09  jfrazier
* More conversion to multi_array.
*
* Revision 1.5  2007/04/29 00:31:57  jfrazier
* Updated to use multi_arrays.
*
*/
