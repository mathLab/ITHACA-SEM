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
// Description: 2D Nodal Triangle Evenly Spaced Point Definitions
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/NodalQuadElec.h>
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
        namespace
        {
            bool isVertex(int i, int j, int npts){
                return (i==0 && j==0) ||
                       (i==(npts-1) && j==0) ||
                       (i==0 && j==(npts-1)) ||
                       (i==(npts-1) && j==(npts-1));
            }

            bool isEdge(int i, int j, int npts){
                return i==0 || j==0 || i==(npts-1) || j==(npts-1);
            }

            bool isEdge_01(int i, int j, int npts){
                return i==0;
            }

            bool isEdge_12(int i, int j, int npts){
                return j==(npts-1);
            }

            bool isEdge_23(int i, int j, int npts){
                return i==(npts-1);
            }

            bool isEdge_30(int i, int j, int npts){
                return j==0;
            }
        }


        void NodalQuadElec::CalculatePoints()
        {
            // Allocate the storage for points
            PointsBaseType::CalculatePoints();

            // Populate m_points
            unsigned int npts = GetNumPoints();
            LibUtilities::PointsKey pkey(npts, LibUtilities::eGaussLobattoLegendre);
            Array<OneD, NekDouble> u;
            LibUtilities::PointsManager()[pkey]->GetPoints(u);

            for(int i=0, index=0; i<npts; ++i){
                for(int j=0; j<npts; ++j,++index){
                    NekDouble    x = u[i];
                    NekDouble    y = u[j];
                    m_points[0][index] = x;
                    m_points[1][index] = y;
                }
            }

            NodalPointReorder2d();

            m_util = MemoryManager<NodalUtilQuad>::AllocateSharedPtr(
                npts - 1, m_points[0], m_points[1]);
        }


        void NodalQuadElec::CalculateWeights()
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
        void NodalQuadElec::CalculateInterpMatrix(const Array<OneD, const NekDouble>& xia,
                                                         const Array<OneD, const NekDouble>& yia,
                                                         Array<OneD, NekDouble>& interp)
        {
            Array<OneD, Array<OneD, NekDouble> > xi(2);
            xi[0] = xia;
            xi[1] = yia;

            boost::shared_ptr<NekMatrix<NekDouble> > mat =
                m_util->GetInterpolationMatrix(xi);
            Vmath::Vcopy(mat->GetRows() * mat->GetColumns(), mat->GetRawPtr(),
                         1, &interp[0], 1);
        }

        // ////////////////////////////////////////
        //        CalculateDerivMatrix()
        void NodalQuadElec::CalculateDerivMatrix()
        {

            // Allocate the derivative matrix.
            PointsBaseType::CalculateDerivMatrix();

            m_derivmatrix[0] = m_util->GetDerivMatrix(0);
            m_derivmatrix[1] = m_util->GetDerivMatrix(1);
        }

        boost::shared_ptr<PointsBaseType> NodalQuadElec::Create(const PointsKey &key)
        {
            boost::shared_ptr<PointsBaseType> returnval(MemoryManager<NodalQuadElec>::AllocateSharedPtr(key));

            returnval->Initialize();

            return returnval;
        }

        void NodalQuadElec::NodalPointReorder2d()
        {
            unsigned int npts = GetNumPoints();
            using std::vector;
            vector<int> vertex;
            vector<int> iEdge_01;
            vector<int> iEdge_12;
            vector<int> iEdge_23;
            vector<int> iEdge_30;
            vector<int> interiorPoints;
            vector<int> map;

               // Build the lattice triangle left to right - bottom to top
            for(int i=0, index=0; i<npts; ++i){ // y-direction
                for(int j=0; j<npts; ++j,++index){ // x-direction

                    if( isVertex(i,j,npts) ) {

                        vertex.push_back(index);

                    } else if( isEdge(i,j,npts) ) { // interior edge

                        if(isEdge_01(i,j,npts)){  // bottom edge

                            iEdge_01.push_back(index);

                        }else if(isEdge_12(i,j,npts)){  // right edge

                            iEdge_12.push_back(index);

                        }else if(isEdge_23(i,j,npts)){  // right edge

                            iEdge_23.push_back(index);

                        }else if(isEdge_30(i,j,npts)){  // right edge

                            iEdge_30.push_back(index);
                        }

                    } else { // Interior points

                        interiorPoints.push_back(index);
                    }
                }
            }

            swap(vertex[2],vertex[3]);
            reverse(iEdge_23.begin(),iEdge_23.end());
            reverse(iEdge_30.begin(),iEdge_30.end());

            // Mapping the vertex, edges, and interior points using the permutation matrix,
            // so the points are ordered anticlockwise.
            for(unsigned int k=0; k<vertex.size(); ++k){
                map.push_back(vertex[k]);
            }

            for(unsigned int k=0; k<iEdge_01.size(); ++k){

                map.push_back(iEdge_01[k]);
            }

            for(unsigned int k=0; k<iEdge_12.size(); ++k){

                map.push_back(iEdge_12[k]);
            }

            for(unsigned int k=0; k<iEdge_23.size(); ++k){

                map.push_back(iEdge_23[k]);
            }

            for(unsigned int k=0; k<iEdge_30.size(); ++k){

                map.push_back(iEdge_30[k]);
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

