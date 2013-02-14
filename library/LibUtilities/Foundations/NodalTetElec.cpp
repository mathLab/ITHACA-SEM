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

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/Foundations/NodalTetElec.h>
#include <LibUtilities/Foundations/NodalTetElecData.h>
#include <LibUtilities/Foundations/NodalUtil.h>
#include <LibUtilities/LinearAlgebra/NekVector.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>

namespace Nektar
{
    namespace LibUtilities 
    {
    
        // ////////////////////////////////////////////////////////
        //  Coordinate the nodal tetrahedron electrostatic points
        
        void NodalTetElec::CalculatePoints()
        {
            // Allocate the storage for points
            Points<NekDouble>::CalculatePoints();

            int index=0,isum=0;
            const int offset = 5; //offset to match Datafile
            NekDouble b,c,d;
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

            ASSERTL1((static_cast<unsigned int>(isum)==m_pointsKey.GetTotNumPoints()),"sum not equal to npts");
        }

        void NodalTetElec::CalculateWeights()
        {
            // Allocate the storage for points
            PointsBaseType::CalculateWeights();

            typedef DataType T;
            
            // Solve the Vandermonde system of integrals for the weight vector
            NekVector<T> w = MakeTetWeights(NekVector<T>(m_points[0]), NekVector<T>(m_points[1]), NekVector<T>(m_points[2]));
            
            m_weights = Array<OneD,T>( w.GetRows(), w.GetPtr() );
        }

          // ////////////////////////////////////////
        //        CalculateInterpMatrix()
        void NodalTetElec::CalculateInterpMatrix(const Array<OneD, const NekDouble>& xia, const Array<OneD, const NekDouble>& yia,
                                                 const Array<OneD, const NekDouble>& zia, Array<OneD, NekDouble>& interp)
                                                   
        {
             NekVector<NekDouble>  x( m_points[0] );
             NekVector<NekDouble>  y( m_points[1] );
             NekVector<NekDouble>  z( m_points[2] );
             NekVector<NekDouble> xi( xia );
             NekVector<NekDouble> yi( yia );
             NekVector<NekDouble> zi( zia );
             NekMatrix<NekDouble> interMat = GetTetInterpolationMatrix(x, y, z, xi, yi, zi);

             int rows = xi.GetRows(), cols = GetTotNumPoints();
             for( int i = 0; i < rows; ++i ) {
                for( int j = 0; j < cols; ++j ) {
                    interp[j + i*cols] = interMat(i,j);
                }
             }             
         }

        void NodalTetElec::CalculateDerivMatrix()
        {
             // Allocate the derivative matrix.
            PointsBaseType::CalculateDerivMatrix();

            NekVector<NekDouble> x( m_points[0] );
            NekVector<NekDouble> y( m_points[1] );
            NekVector<NekDouble> z( m_points[2] );
            NekVector<NekDouble> xi = x;
            NekVector<NekDouble> yi = y;
            NekVector<NekDouble> zi = z;

            *m_derivmatrix[0] = *GetTetXDerivativeMatrix(x,y,z,xi,yi,zi);

            *m_derivmatrix[1] = *GetTetYDerivativeMatrix(x,y,z,xi,yi,zi);

            *m_derivmatrix[2] = *GetTetZDerivativeMatrix(x,y,z,xi,yi,zi);
        }

        boost::shared_ptr<PointsBaseType> NodalTetElec::Create(const PointsKey &key)
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


