///////////////////////////////////////////////////////////////////////////////
//
// File NodalTriEvenlySpaced.cpp
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
// Description: 2D Nodal Triangle Evenly Spaced Point Definitions
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/Foundations/NodalTriEvenlySpaced.h>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/Foundations/NodalUtil.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <LibUtilities/LinearAlgebra/NekVector.hpp>
#include <vector>



namespace Nektar
{
    namespace LibUtilities
    {
        bool NodalTriEvenlySpaced::initPointsManager[] = {
            PointsManager().RegisterCreator(PointsKey(0, eNodalTriEvenlySpaced), NodalTriEvenlySpaced::Create)
        };

        namespace
        {
           // construct the geometory and set the coordinate of triangle
           // edges and vertices are ordered as anticlockwise
            bool isVertex(int i, int j, int npts){
                return (i==0 && j==0) || (i==(npts-1) && j==0) || (i==0 && j==(npts-1));
            }

            bool isEdge(int i, int j, int npts){
                return i==0 || j==0 || i+j==npts-1; //i+j=tot num of steps
            }
            
            bool isEdge_1(int i, int j, int npts){
                boost::ignore_unused(j, npts);
                return i==0;
            }
            
            bool isEdge_2(int i, int j, int npts){
                return i+j==npts-1;
            }
        }

        
        void NodalTriEvenlySpaced::CalculatePoints()
        {
            // Allocate the storage for points
            PointsBaseType::CalculatePoints();

            // Populate m_points
            unsigned int npts = GetNumPoints();
            NekDouble delta = 2.0/(npts - 1.0);
            for(int i=0, index=0; i<npts; ++i){ // y-direction
                for(int j=0; j<npts-i; ++j,++index){ // x-direction
                    NekDouble    x = -1.0 + j*delta;
                    NekDouble    y = -1.0 + i*delta;
                    m_points[0][index] = x;
                    m_points[1][index] = y;
                }
            }

            NodalPointReorder2d();

            m_util = MemoryManager<NodalUtilTriangle>::AllocateSharedPtr(
                npts - 1, m_points[0], m_points[1]);
        }


        void NodalTriEvenlySpaced::CalculateWeights()
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
        void NodalTriEvenlySpaced::CalculateInterpMatrix(const Array<OneD, const NekDouble>& xia,
                                                         const Array<OneD, const NekDouble>& yia,
                                                         Array<OneD, NekDouble>& interp)
        {
            Array<OneD, Array<OneD, NekDouble> > xi(2);
            xi[0] = xia;
            xi[1] = yia;

            std::shared_ptr<NekMatrix<NekDouble> > mat =
                m_util->GetInterpolationMatrix(xi);
            Vmath::Vcopy(mat->GetRows() * mat->GetColumns(), mat->GetRawPtr(),
                         1, &interp[0], 1);
        }

        // ////////////////////////////////////////
        //        CalculateDerivMatrix()
        void NodalTriEvenlySpaced::CalculateDerivMatrix()
        {
            
            // Allocate the derivative matrix.
            PointsBaseType::CalculateDerivMatrix();

            m_derivmatrix[0] = m_util->GetDerivMatrix(0);
            m_derivmatrix[1] = m_util->GetDerivMatrix(1);
        } 

        std::shared_ptr<PointsBaseType> NodalTriEvenlySpaced::Create(const PointsKey &key)
        {
            std::shared_ptr<PointsBaseType> returnval(MemoryManager<NodalTriEvenlySpaced>::AllocateSharedPtr(key));

            returnval->Initialize();

            return returnval;
        }
        
        void NodalTriEvenlySpaced::NodalPointReorder2d()
        {
            unsigned int npts = GetNumPoints();
            using std::vector;
            vector<int> vertex;
            vector<int> iEdge_1; // interior edge points on the bottom triangle edge
            vector<int> iEdge_2; // interior edge points on the right triangle edge
            vector<int> iEdge_3; // interior edge points on the left triangle edge
            vector<int> interiorPoints;
            vector<int> map;

               // Build the lattice triangle left to right - bottom to top
            for(int i=0, index=0; i<npts; ++i){ // y-direction
                for(int j=0; j<npts-i; ++j,++index){ // x-direction

                    if( isVertex(i,j,npts) ) {

                        vertex.push_back(index);

                    } else if( isEdge(i,j,npts) ) { // interior edge

                        if(isEdge_1(i,j,npts)){  // bottom edge

                            iEdge_1.push_back(index);

                        }else if(isEdge_2(i,j,npts)){  // right edge

                            iEdge_2.push_back(index);

                        }else   // left edge
                        {
                            // Add backwards.  This is because of counter clockwise.
                            iEdge_3.insert(iEdge_3.begin(), index);
                        }

                    } else { // Interior points

                        interiorPoints.push_back(index);
                    }
                }
            }

            // Mapping the vertex, edges, and interior points using the permutation matrix,
            // so the points are ordered anticlockwise.
            for(unsigned int k=0; k<vertex.size(); ++k){

                map.push_back(vertex[k]);
            }

            for(unsigned int k=0; k<iEdge_1.size(); ++k){

                map.push_back(iEdge_1[k]);
            }

            for(unsigned int k=0; k<iEdge_2.size(); ++k){

                map.push_back(iEdge_2[k]);
            }

            for(unsigned int k=0; k<iEdge_3.size(); ++k){

                map.push_back(iEdge_3[k]);
            }
            
            for(unsigned int k=0; k<interiorPoints.size(); ++k){

                map.push_back(interiorPoints[k]);
            }


            Array<OneD,NekDouble> points[2];
            points[0] = Array<OneD,NekDouble>(GetTotNumPoints());
            points[1] = Array<OneD,NekDouble>(GetTotNumPoints());
            for(unsigned int index=0; index<map.size(); ++index){
                points[0][index] = m_points[0][index];
                points[1][index] = m_points[1][index];
            }
            
            for(unsigned int index=0; index<map.size(); ++index){
                m_points[0][index] = points[0][map[index]];
                m_points[1][index] = points[1][map[index]];
            }

        }

    } // end of namespace 
} // end of namespace 

