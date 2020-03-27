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
        bool NodalTetElec::initPointsManager[] = {
            PointsManager().RegisterCreator(PointsKey(0, eNodalTetElec),         NodalTetElec::Create)
        };

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

            m_util = MemoryManager<NodalUtilTetrahedron>::AllocateSharedPtr(
                numPoints - 1, m_points[0], m_points[1], m_points[2]);
        }

        void NodalTetElec::CalculateWeights()
        {
            // Allocate the storage for points
            PointsBaseType::CalculateWeights();

            typedef DataType T;

            // Solve the Vandermonde system of integrals for the weight vector
            NekVector<T> w = m_util->GetWeights();
            m_weights = Array<OneD,T>( w.GetRows(), w.GetPtr() );
        }

          // ////////////////////////////////////////
        //        CalculateInterpMatrix()
        void NodalTetElec::CalculateInterpMatrix(const Array<OneD, const NekDouble>& xia, const Array<OneD, const NekDouble>& yia,
                                                 const Array<OneD, const NekDouble>& zia, Array<OneD, NekDouble>& interp)

        {
             Array<OneD, Array<OneD, NekDouble> > xi(3);
             xi[0] = xia;
             xi[1] = yia;
             xi[2] = zia;

             std::shared_ptr<NekMatrix<NekDouble> > mat =
                 m_util->GetInterpolationMatrix(xi);
             Vmath::Vcopy(mat->GetRows() * mat->GetColumns(), mat->GetRawPtr(),
                          1, &interp[0], 1);
         }

        void NodalTetElec::CalculateDerivMatrix()
        {
            // Allocate the derivative matrix.
            PointsBaseType::CalculateDerivMatrix();

            m_derivmatrix[0] = m_util->GetDerivMatrix(0);
            m_derivmatrix[1] = m_util->GetDerivMatrix(1);
            m_derivmatrix[2] = m_util->GetDerivMatrix(2);
        }

        std::shared_ptr<PointsBaseType> NodalTetElec::Create(const PointsKey &key)
        {
            std::shared_ptr<PointsBaseType> returnval(MemoryManager<NodalTetElec>::AllocateSharedPtr(key));
            returnval->Initialize();
            return returnval;
        }

        void NodalTetElec::NodalPointReorder3d()
        {
            int cnt;
            int istart,iend;

            const int nVerts = 4;
            const int nEdgeInteriorPoints = GetNumPoints()-2;
            const int nFaceInteriorPoints = (GetNumPoints()-3)*(GetNumPoints()-2)/2;
            //const int nBoundaryPoints = 4 + 6*nEdgeInteriorPoints + 4*nFaceInteriorPoints;
            const int nAllPoints = GetNumPoints()*(GetNumPoints()+1)*(GetNumPoints()+2)/6;
            if(nEdgeInteriorPoints==0)
            {
                return;
            }

            //group all edge 1 points
            istart = nVerts;
            for(int i = cnt = istart; i < nAllPoints; i++)
            {
                if( fabs(m_points[1][i] + 1.0) < NekConstants::kNekZeroTol &&
                    fabs(m_points[2][i] + 1.0) < NekConstants::kNekZeroTol)
                {
                    std::swap(m_points[0][cnt], m_points[0][i]);
                    std::swap(m_points[1][cnt], m_points[1][i]);
                    std::swap(m_points[2][cnt], m_points[2][i]);
                    cnt++;
                }
            }

            // bubble sort edge 1 (counterclockwise numbering)
            iend = istart + nEdgeInteriorPoints;
            for(int i = istart; i < iend; i++)
            {
                for(int j = istart+1; j < iend; j++)
                {
                    if(m_points[0][j] < m_points[0][j-1])
                    {
                        std::swap(m_points[0][j], m_points[0][j-1]);
                        std::swap(m_points[1][j], m_points[1][j-1]);
                        std::swap(m_points[2][j], m_points[2][j-1]);
                    }
                }
            }

            // group the points of edge 2 together;
            istart = iend;
            for(int i = cnt = istart; i < nAllPoints; i++)
            {
                if( fabs(m_points[1][i]+m_points[0][i]) < NekConstants::kNekZeroTol &&
                    fabs(m_points[2][i] + 1.0) < NekConstants::kNekZeroTol)
                {
                    std::swap(m_points[0][cnt], m_points[0][i]);
                    std::swap(m_points[1][cnt], m_points[1][i]);
                    std::swap(m_points[2][cnt], m_points[2][i]);
                    cnt++;
                }
            }

            // bubble sort edge 2 (counterclockwise numbering)
            iend = istart + nEdgeInteriorPoints;
            for(int i = istart; i < iend; i++)
            {
                for(int j = istart+1; j < iend; j++)
                {
                    if(m_points[1][j] < m_points[1][j-1])
                    {
                        std::swap(m_points[0][j], m_points[0][j-1]);
                        std::swap(m_points[1][j], m_points[1][j-1]);
                        std::swap(m_points[2][j], m_points[2][j-1]);
                    }
                }
            }

            // group the points of edge 3 together;
            istart = iend;
            for(int i = cnt = istart; i < nAllPoints; i++)
            {
                if( fabs(m_points[0][i] + 1.0) < NekConstants::kNekZeroTol &&
                    fabs(m_points[2][i] + 1.0) < NekConstants::kNekZeroTol)
                {
                    std::swap(m_points[0][cnt], m_points[0][i]);
                    std::swap(m_points[1][cnt], m_points[1][i]);
                    std::swap(m_points[2][cnt], m_points[2][i]);
                    cnt++;
                }
            }

            // bubble sort edge 3 (counterclockwise numbering)
            iend = istart + nEdgeInteriorPoints;
            for(int i = istart; i < iend; i++)
            {
                for(int j = istart+1; j < iend; j++)
                {
                    if(m_points[1][j] > m_points[1][j-1])
                    {
                        std::swap(m_points[0][j], m_points[0][j-1]);
                        std::swap(m_points[1][j], m_points[1][j-1]);
                        std::swap(m_points[2][j], m_points[2][j-1]);
                    }
                }
            }

            // group the points of edge 4 together;
            istart = iend;
            for(int i = cnt = istart; i < nAllPoints; i++)
            {
                if( fabs(m_points[0][i] + 1.0) < NekConstants::kNekZeroTol &&
                    fabs(m_points[1][i] + 1.0) < NekConstants::kNekZeroTol)
                {
                    std::swap(m_points[0][cnt], m_points[0][i]);
                    std::swap(m_points[1][cnt], m_points[1][i]);
                    std::swap(m_points[2][cnt], m_points[2][i]);
                    cnt++;
                }
            }

            // bubble sort edge 3 (counterclockwise numbering)
            iend = istart + nEdgeInteriorPoints;
            for(int i = istart; i < iend; i++)
            {
                for(int j = istart+1; j < iend; j++)
                {
                    if(m_points[2][j] < m_points[2][j-1])
                    {
                        std::swap(m_points[0][j], m_points[0][j-1]);
                        std::swap(m_points[1][j], m_points[1][j-1]);
                        std::swap(m_points[2][j], m_points[2][j-1]);
                    }
                }
            }

            // group the points of edge 5 together;
            istart = iend;
            for(int i = cnt = istart; i < nAllPoints; i++)
            {
                if( fabs(m_points[0][i]+m_points[2][i]) < NekConstants::kNekZeroTol &&
                    fabs(m_points[1][i] + 1.0) < NekConstants::kNekZeroTol)
                {
                    std::swap(m_points[0][cnt], m_points[0][i]);
                    std::swap(m_points[1][cnt], m_points[1][i]);
                    std::swap(m_points[2][cnt], m_points[2][i]);
                    cnt++;
                }
            }

            // bubble sort edge 5 (counterclockwise numbering)
            iend = istart + nEdgeInteriorPoints;
            for(int i = istart; i < iend; i++)
            {
                for(int j = istart+1; j < iend; j++)
                {
                    if(m_points[2][j] < m_points[2][j-1])
                    {
                        std::swap(m_points[0][j], m_points[0][j-1]);
                        std::swap(m_points[1][j], m_points[1][j-1]);
                        std::swap(m_points[2][j], m_points[2][j-1]);
                    }
                }
            }

            // group the points of edge 6 together;
            istart = iend;
            for(int i = cnt = istart; i < nAllPoints; i++)
            {
                if( fabs(m_points[1][i]+m_points[2][i]) < NekConstants::kNekZeroTol &&
                    fabs(m_points[0][i] + 1.0) < NekConstants::kNekZeroTol)
                {
                    std::swap(m_points[0][cnt], m_points[0][i]);
                    std::swap(m_points[1][cnt], m_points[1][i]);
                    std::swap(m_points[2][cnt], m_points[2][i]);
                    cnt++;
                }
            }

            // bubble sort edge 6 (counterclockwise numbering)
            iend = istart + nEdgeInteriorPoints;
            for(int i = istart; i < iend; i++)
            {
                for(int j = istart+1; j < iend; j++)
                {
                    if(m_points[2][j] < m_points[2][j-1])
                    {
                        std::swap(m_points[0][j], m_points[0][j-1]);
                        std::swap(m_points[1][j], m_points[1][j-1]);
                        std::swap(m_points[2][j], m_points[2][j-1]);
                    }
                }
            }

            if(GetNumPoints() < 4)
            {
                //no face points
                return;
            }

            // group the points of face 1 together;
            istart = iend;
            for(int i = cnt = istart; i < nAllPoints; i++)
            {
                if(fabs(m_points[2][i] + 1.0) < NekConstants::kNekZeroTol)
                {
                    std::swap(m_points[0][cnt], m_points[0][i]);
                    std::swap(m_points[1][cnt], m_points[1][i]);
                    std::swap(m_points[2][cnt], m_points[2][i]);
                    cnt++;
                }
            }

            // bubble sort face1 (tensor numbering)
            iend = istart + nFaceInteriorPoints;
            bool repeat = true;
            while(repeat)
            {
                repeat = false;
                for(int i = istart; i < iend - 1; i++)
                {
                    if(m_points[1][i] > m_points[1][i+1])
                    {
                        std::swap(m_points[0][i+1], m_points[0][i]);
                        std::swap(m_points[1][i+1], m_points[1][i]);
                        std::swap(m_points[2][i+1], m_points[2][i]);
                        repeat = true;
                    }
                }
            }
            int offset = 0;
            int npl = GetNumPoints() - 3;
            while(npl > 1)
            {
                repeat = true;
                while(repeat)
                {
                    repeat = false;
                    for(int i = offset+istart; i < offset+istart + npl - 1; i++)
                    {
                        if(m_points[0][i] > m_points[0][i+1])
                        {
                            std::swap(m_points[0][i+1], m_points[0][i]);
                            std::swap(m_points[1][i+1], m_points[1][i]);
                            std::swap(m_points[2][i+1], m_points[2][i]);
                            repeat = true;
                        }
                    }
                }
                offset += npl;
                npl--;
            }

            // group the points of face 2 together;
            istart = iend;
            for(int i = cnt = istart; i < nAllPoints; i++)
            {
                if(fabs(m_points[1][i] + 1.0) < NekConstants::kNekZeroTol)
                {
                    std::swap(m_points[0][cnt], m_points[0][i]);
                    std::swap(m_points[1][cnt], m_points[1][i]);
                    std::swap(m_points[2][cnt], m_points[2][i]);
                    cnt++;
                }
            }

            // bubble sort face2 (tensor numbering)
            iend = istart + nFaceInteriorPoints;
            repeat = true;
            while(repeat)
            {
                repeat = false;
                for(int i = istart; i < iend - 1; i++)
                {
                    if(m_points[2][i] > m_points[2][i+1])
                    {
                        std::swap(m_points[0][i+1], m_points[0][i]);
                        std::swap(m_points[1][i+1], m_points[1][i]);
                        std::swap(m_points[2][i+1], m_points[2][i]);
                        repeat = true;
                    }
                }
            }
            offset = 0;
            npl = GetNumPoints() - 3;
            while(npl > 1)
            {
                repeat = true;
                while(repeat)
                {
                    repeat = false;
                    for(int i = offset+istart; i < offset+istart + npl - 1; i++)
                    {
                        if(m_points[0][i] > m_points[0][i+1])
                        {
                            std::swap(m_points[0][i+1], m_points[0][i]);
                            std::swap(m_points[1][i+1], m_points[1][i]);
                            std::swap(m_points[2][i+1], m_points[2][i]);
                            repeat = true;
                        }
                    }
                }
                offset += npl;
                npl--;
            }

            // group the points of face 3 together;
            istart = iend;
            for(int i = cnt = istart; i < nAllPoints; i++)
            {
                if(fabs(m_points[1][i]+m_points[0][i]+m_points[2][i]+1.0) < 1E-9) //nek zero tol too small
                {
                    std::swap(m_points[0][cnt], m_points[0][i]);
                    std::swap(m_points[1][cnt], m_points[1][i]);
                    std::swap(m_points[2][cnt], m_points[2][i]);
                    cnt++;
                }
            }

            // bubble sort face3 (tensor numbering)
            iend = istart + nFaceInteriorPoints;
            repeat = true;
            while(repeat)
            {
                repeat = false;
                for(int i = istart; i < iend - 1; i++)
                {
                    if(m_points[2][i] > m_points[2][i+1])
                    {
                        std::swap(m_points[0][i+1], m_points[0][i]);
                        std::swap(m_points[1][i+1], m_points[1][i]);
                        std::swap(m_points[2][i+1], m_points[2][i]);
                        repeat = true;
                    }
                }
            }
            offset = 0;
            npl = GetNumPoints() - 3;
            while(npl > 1)
            {
                repeat = true;
                while(repeat)
                {
                    repeat = false;
                    for(int i = offset+istart; i < offset+istart + npl - 1; i++)
                    {
                        if(m_points[1][i] > m_points[1][i+1])
                        {
                            std::swap(m_points[0][i+1], m_points[0][i]);
                            std::swap(m_points[1][i+1], m_points[1][i]);
                            std::swap(m_points[2][i+1], m_points[2][i]);
                            repeat = true;
                        }
                    }
                }
                offset += npl;
                npl--;
            }

            // group the points of face 4 together;
            istart = iend;
            for(int i = cnt = istart; i < nAllPoints; i++)
            {
                if(fabs(m_points[0][i] + 1.0) < NekConstants::kNekZeroTol)
                {
                    std::swap(m_points[0][cnt], m_points[0][i]);
                    std::swap(m_points[1][cnt], m_points[1][i]);
                    std::swap(m_points[2][cnt], m_points[2][i]);
                    cnt++;
                }
            }

            // bubble sort face4 (tensor numbering)
            iend = istart + nFaceInteriorPoints;
            repeat = true;
            while(repeat)
            {
                repeat = false;
                for(int i = istart; i < iend - 1; i++)
                {
                    if(m_points[2][i] > m_points[2][i+1])
                    {
                        std::swap(m_points[0][i+1], m_points[0][i]);
                        std::swap(m_points[1][i+1], m_points[1][i]);
                        std::swap(m_points[2][i+1], m_points[2][i]);
                        repeat = true;
                    }
                }
            }
            offset = 0;
            npl = GetNumPoints() - 3;
            while(npl > 1)
            {
                repeat = true;
                while(repeat)
                {
                    repeat = false;
                    for(int i = offset+istart; i < offset+istart + npl - 1; i++)
                    {
                        if(m_points[1][i] > m_points[1][i+1])
                        {
                            std::swap(m_points[0][i+1], m_points[0][i]);
                            std::swap(m_points[1][i+1], m_points[1][i]);
                            std::swap(m_points[2][i+1], m_points[2][i]);
                            repeat = true;
                        }
                    }
                }
                offset += npl;
                npl--;
            }

        }

    } // end of namespace stdregion
} // end of namespace stdregion


