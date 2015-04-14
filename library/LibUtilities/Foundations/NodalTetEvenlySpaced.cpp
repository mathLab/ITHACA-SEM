///////////////////////////////////////////////////////////////////////////////
//
// File NodalTetEvenlySpaced.cpp
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
// Description: 3D Nodal Tetrahedron Evenly Spaced Point Definitions
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Foundations/NodalTetEvenlySpaced.h>
#include <LibUtilities/Foundations/NodalUtil.h>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <LibUtilities/LinearAlgebra/NekVector.hpp>
#include <vector>



namespace Nektar
{
    namespace LibUtilities
    {
        namespace
        {
           // construct the geometory and set the coordinate of tetrahedron
           // edges and vertices are ordered as anticlockwise
           
            bool isVertex(int x, int y, int z, int npts){
                return (x==0 && y==0 && z==0) || (x==(npts-1) && y==0 && z==0) || (x==0 && y==(npts-1) && z==0) || (x==0 && y==0 && z==(npts-1));
            }

            bool isEdge_01(int x, int y, int z, int npts){  // edge 0
                return y==0 && z==0;
            }

            bool isEdge_12(int x, int y, int z, int npts){  // edge 1
                return z==0 && x + y == npts -1;
            }

            bool isEdge_20(int x, int y, int z, int npts){  // edge 2
                return x==0 && z==0;
            }

            bool isEdge_03(int x, int y, int z, int npts){  // edge 3
                return x==0 && y==0;
            }

            bool isEdge_13(int x, int y, int z, int npts){  // edge 4
                return y==0 && x + z == npts -1;
            }

            bool isEdge_23(int x, int y, int z, int npts){  // edge 5
                return x==0 && y + z == npts -1;
            }

            bool isEdge(int x, int y, int z, int npts){
               return isEdge_01(x,y,z,npts)||isEdge_12(x,y,z,npts)||isEdge_20(x,y,z,npts)
                      ||isEdge_03(x,y,z,npts)||isEdge_13(x,y,z,npts)||isEdge_23(x,y,z,npts);
            }

            bool isFace_012(int x, int y, int z, int npts){  // bottom face (face 0)
                return z==0;
            }

            bool isFace_013(int x, int y, int z, int npts){  // face 1
                return y==0;
            }

            bool isFace_123(int x, int y, int z, int npts){  // face 2
                return x+y+z==npts-1;
            }

            bool isFace_203(int x, int y, int z, int npts){  // face 3
                return x==0;
            }

            bool isFace(int x, int y, int z, int npts){
                return isFace_012(x, y, z, npts) || isFace_013(x, y, z, npts)
                       || isFace_123(x, y, z, npts) || isFace_203(x, y, z, npts);
            }
        }

        // Calculate evenly spaced number of points
        void NodalTetEvenlySpaced::CalculatePoints()
        {
            // Allocate the storage for points
            PointsBaseType::CalculatePoints();
            
            // Populate m_points
            unsigned int npts = GetNumPoints();
            NekDouble delta = 2.0/(npts - 1.0);
            for(unsigned int z=0, index=0; z<npts; ++z){
                for(int y=0; y<npts-z; ++y){
                    for(int x=0; x<npts-z-y; ++x, ++index){
                        NekDouble xi = -1.0 + x*delta;
                        NekDouble yi = -1.0 + y*delta;
                        NekDouble zi = -1.0 + z*delta;

                        m_points[0][index] = xi;
                        m_points[1][index] = yi;
                        m_points[2][index] = zi;                        
                    }
                }                
            }

            NodalPointReorder3d();            
        }

        void NodalTetEvenlySpaced::NodalPointReorder3d()
        {
            unsigned int npts = GetNumPoints();
            using std::vector;
            vector<int> vertex;
            vector<int> iEdge_01;  // interior edge 0
            vector<int> iEdge_12;  // interior edge 1
            vector<int> iEdge_20;  // interior edge 2
            vector<int> iEdge_03;  // interior edge 3
            vector<int> iEdge_13;  // interior edge 4
            vector<int> iEdge_23;  // interior edge 5
            vector<int> iFace_012; // interior face 0
            vector<int> iFace_013; // interior face 1 
            vector<int> iFace_123; // interior face 2
            vector<int> iFace_203; // interior face 3
            vector<int> interiorVolumePoints; // interior volume points
            vector<int> map;

            // Build the lattice tetrahedron left to right - bottom to top
            for(int z=0, index=0; z<npts; ++z){
                for(int y=0; y<npts-z; ++y){
                    for(int x=0; x<npts-z-y; ++x, ++index){

                        if( isVertex(x,y,z,npts) ){ // vertex
                        
                            vertex.push_back(index);
                            
                        } else if( isEdge(x,y,z,npts) ){ // interior edge
                        
                            if( isEdge_01(x,y,z,npts) ){  // interior edge 0
                            
                                iEdge_01.push_back(index);
                                
                            } else if( isEdge_12(x,y,z,npts) ){  // interior edge 1

                                iEdge_12.push_back(index);
                                
                            } else if( isEdge_20(x,y,z,npts) ){  // interior edge 2

                                iEdge_20.insert(iEdge_20.begin(), index);
                                
                            } else if( isEdge_03(x,y,z,npts) ){ // interior edge 3

                                    iEdge_03.push_back(index);
                                
                            } else if( isEdge_13(x,y,z,npts) ){ // interior edge 4

                                iEdge_13.push_back(index);
                                
                            } else if( isEdge_23(x,y,z,npts) ){ // interior edge 5

                                  iEdge_23.push_back(index);
                                
                            }
                            
                        } else if( isFace(x,y,z,npts) ) {  // interior face

                            if( isFace_012(x,y,z,npts) ){  // interior face 0

                                iFace_012.push_back(index);
                            
                            } else if( isFace_013(x,y,z,npts) ){  // interior face 1

                                iFace_013.push_back(index);
                            
                            } else if( isFace_123(x,y,z,npts) ){ // interior face 2

                                iFace_123.push_back(index);

                            } else if( isFace_203(x,y,z,npts) ){  // interior face 3

                                iFace_203.push_back(index);

                            }
                        } else {  // interior volume points

                            interiorVolumePoints.push_back(index);                            
                        }
                    }
                }
            }
           
            // Mapping the vertex, edges, faces, interior volume points using the permutation matrix,
            // so the points are ordered anticlockwise.
            for(unsigned int n=0; n<vertex.size(); ++n){

                map.push_back(vertex[n]);
            }

            for(unsigned int n=0; n<iEdge_01.size(); ++n){

                map.push_back(iEdge_01[n]);
            }

            for(unsigned int n=0; n<iEdge_12.size(); ++n){

                map.push_back(iEdge_12[n]);
            }

            for(unsigned int n=0; n<iEdge_20.size(); ++n){

                map.push_back(iEdge_20[n]);
            }

            for(unsigned int n=0; n<iEdge_03.size(); ++n){

                map.push_back(iEdge_03[n]);
            }

            for(unsigned int n=0; n<iEdge_13.size(); ++n){

                map.push_back(iEdge_13[n]);
            }

            for(unsigned int n=0; n<iEdge_23.size(); ++n){

                map.push_back(iEdge_23[n]);
            }

            for(unsigned int n=0; n<iFace_012.size(); ++n){

                map.push_back(iFace_012[n]);
            }

            for(unsigned int n=0; n<iFace_013.size(); ++n){

                map.push_back(iFace_013[n]);
            }

            for(unsigned int n=0; n<iFace_123.size(); ++n){

                map.push_back(iFace_123[n]);
            }

            for(unsigned int n=0; n<iFace_203.size(); ++n){

                map.push_back(iFace_203[n]);
            }

            for(unsigned int n=0; n<interiorVolumePoints.size(); ++n){

                map.push_back(interiorVolumePoints[n]);
            }


            Array<OneD, NekDouble> points[3];
            points[0] = Array<OneD, NekDouble>(GetTotNumPoints());
            points[1] = Array<OneD, NekDouble>(GetTotNumPoints());
            points[2] = Array<OneD, NekDouble>(GetTotNumPoints());
            for(unsigned int index=0; index<map.size(); ++index){

                points[0][index] = m_points[0][index];
                points[1][index] = m_points[1][index];
                points[2][index] = m_points[2][index];

            }

            for(unsigned int index=0; index<map.size(); ++index){

                m_points[0][index] = points[0][map[index]];
                m_points[1][index] = points[1][map[index]];
                m_points[2][index] = points[2][map[index]];

            }
                                   
        }

        

        void NodalTetEvenlySpaced::CalculateWeights()
        {            
            // Allocate storage for points
            PointsBaseType::CalculateWeights();

            typedef DataType T;

            // Solve the Vandermonde system of integrals for the weight vector
            NekVector<T> w = MakeTetWeights(NekVector<T>(m_points[0]), NekVector<T>(m_points[1]), NekVector<T>(m_points[2]));
            
            m_weights = Array<OneD,T>( w.GetRows(), w.GetPtr() );

        }


        // ////////////////////////////////////////
        //        CalculateInterpMatrix()
        void NodalTetEvenlySpaced::CalculateInterpMatrix(const Array<OneD, const NekDouble>& xia,
                                                         const Array<OneD, const NekDouble>& yia,
                                                         const Array<OneD, const NekDouble>& zia,
                                                         Array<OneD, NekDouble>& interp)
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

        // ////////////////////////////////////////
        //        CalculateDerivMatrix()
        void NodalTetEvenlySpaced::CalculateDerivMatrix()
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

        boost::shared_ptr<PointsBaseType> NodalTetEvenlySpaced::Create(const PointsKey &key)
        {
            boost::shared_ptr<PointsBaseType> returnval(MemoryManager<NodalTetEvenlySpaced>::AllocateSharedPtr(key));

            returnval->Initialize();

            return returnval;
        }
        


    } // end of namespace 
} // end of namespace 

