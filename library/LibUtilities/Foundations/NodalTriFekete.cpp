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

 //#include <LibUtilities/LinearAlgebra/NekMatrixFwd.hpp>
 #include <LibUtilities/LinearAlgebra/Lapack.hpp>
 //#include <LibUtilities/ExpressionTemplates/ExpressionTemplates.hpp>
 //#include <LibUtilities/Memory/NekMemoryManager.hpp>
//  
//  #include <boost/call_traits.hpp>
//  
//  #include <algorithm>

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

                    m_points[0][isum] = 2.0*b - 1.0;
                    m_points[1][isum] = 2.0*c - 1.0;
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
                        m_points[0][isum] = 2.0*b - 1.0;
                        m_points[1][isum] = 2.0*c - 1.0;
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
                        m_points[0][isum] = 2.0*b - 1.0;
                        m_points[1][isum] = 2.0*c - 1.0;
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
                        m_points[0][isum] = 2.0*b - 1.0;
                        m_points[1][isum] = 2.0*c - 1.0;
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
            // Allocate the storage for points
            PointsBaseType::CalculateWeights();
            unsigned int npts = m_pointsKey.GetNumPoints();
            
            PointsKey nodalTriFeketeKey(npts, eNodalTriFekete);
            ptr<PointsBaseType> ptr = PointsManager()[nodalTriFeketeKey];
                ConstArray<TwoD, NekDouble> z;
                ConstArray<OneD, NekDouble> w;

              //  ptr->GetZW(z,w);
                
            
        }


        Array<OneD, NekDouble> LegendrePoly(int degree, Array<OneD, NekDouble>& x);
        Array<OneD, NekDouble> JacobiPoly(int degree, Array<OneD, NekDouble>& x, int alpha, int beta);
        Array<OneD, NekDouble> DubinerPoly(int p, int q, Array<OneD, NekDouble>& x, Array<OneD, NekDouble>& y);
        Array<TwoD, NekDouble> getVandermonde(Array<OneD, NekDouble>& x, Array<OneD, NekDouble>& y, int degree);
        Array<TwoD, NekDouble> getVandermonde(Array<OneD, NekDouble>& x, Array<OneD, NekDouble>& y);
        //static void Invert(unsigned int rows, unsigned int columns, Array<OneD, NekDouble>& data);
        Array<OneD, NekDouble> vectorizeMatrix(const Array<TwoD, NekDouble> & A, int M, int N);
        int getDegree(int nBasisFunctions);
        

        //Array<TwoD, NekDouble> getInterpolationMatrix(Array<OneD, NekDouble>& x, Array<OneD, NekDouble>& y,Array<OneD, NekDouble>& xi, Array<OneD, NekDouble>& yi){
        //    int nNodes = x.num_elements();
        //    int degree = getDegree(nNodes);

        //    int M, N; // TODO: set to row/col
        //                
        //    
        //    Array<TwoD, NekDouble> S = getVandermonde(x, y); // Square 'short' matrix
        //    Array<TwoD, NekDouble> T = getVandermonde(xi, yi, degree); // Coefficient interpolation matrix (tall)

        //    Array<OneD, NekDouble> invMatrix = vectorizeMatrix(S, M, N);
        //    //Invert(M,N,invMatrix);
        //    
        //    for(int i=0; i<M; ++i){
        //      for(int j=0; j<N; ++j){
        //      //  v[CalculateIndex(i,j,M,N)] = T[i][j];
        //     }
        //    }
        //    // Get the interpolation matrix
        ////    return T * invMatrix;
        //    return S; // TODO fix
        //}

        static unsigned int CalculateIndex(unsigned int row, unsigned int column, unsigned int matrixRows, unsigned int matrixColumns) {
                 return row*matrixColumns + column;
             }
        Array<OneD, NekDouble> vectorizeMatrix(const Array<TwoD, NekDouble> & A, int M, int N) {
            Array<OneD, NekDouble> v(M*N);
            for( int i=0; i<M; ++i ) {
                for( int j=0; j<N; ++j ) {
                    v[CalculateIndex(i,j,M,N)] = A[i][j];
                }
            }
            return v;
        }


        static void Transpose(unsigned int& rows, unsigned int& columns, Array<OneD, NekDouble>& data)  {
                 Array<OneD, NekDouble> temp(data.num_elements());

                 for(unsigned int row = 0; row < rows; ++row)
                 {
                     for(unsigned int column = 0; column < columns; ++column)
                     {
                         unsigned int firstIndex = CalculateIndex(row, column, rows, columns);
                         unsigned int secondIndex = CalculateIndex(column, row, columns, rows);

                         temp[secondIndex] = data[firstIndex];
                     }
                 }

                 std::swap(rows, columns);
                 std::swap(data, temp);
             }
        
      
//      static void Invert(unsigned int rows, unsigned int columns, Array<OneD, NekDouble>& data) {
////#ifdef NEKTAR_USING_LAPACK
//                 ASSERTL0(rows == columns, "Matrix Inversion only works for square arrays.");
//
//                 /// Incoming data is row major, make it column major for lapack calls.
//                 Transpose(rows, columns, data);
//
//                 int m = rows;
//                 int n = columns;
//                 int pivotSize = std::max(1, std::min(m, n));
//
//                 Array<OneD, int> ipivot(pivotSize);
//                 int info = 0;
//                 Lapack::Dgetrf(m, n, data.get(), m, ipivot.get(), info);
//
//                 if( info < 0 )
//                 {
//                     std::string message = "ERROR: The " + boost::lexical_cast<std::string>(-info) +
//                             "th parameter had an illegal parameter for dgetrf";
//                     ASSERTL0(false, message.c_str());
//                 }
//                 else if( info > 0 )
//                 {
//                     std::string message = "ERROR: Element u_" + boost::lexical_cast<std::string>(info) +
//                             boost::lexical_cast<std::string>(info) + " is 0 from dgetrf";
//                     ASSERTL0(false, message.c_str());
//                 }
//
//                 unsigned int workSize = 64*n;
//
//                 Array<OneD, NekDouble> work(workSize);
//                 Lapack::Dgetri(n, data.get(), n, ipivot.get(), work.get(), workSize, info);
//
//                 if( info < 0 )
//                 {
//                     std::string message = "ERROR: The " + boost::lexical_cast<std::string>(-info) +
//                             "th parameter had an illegal parameter for dgetri";
//                     ASSERTL0(false, message.c_str());
//                 }
//                 else if( info > 0 )
//                 {
//                     std::string message = "ERROR: Element u_" + boost::lexical_cast<std::string>(info) +
//                             boost::lexical_cast<std::string>(info) + " is 0 from dgetri";
//                     ASSERTL0(false, message.c_str());
//                 }
//
//                 // Put it back to row major form.
//                 Transpose(rows, columns, data);
//
// //#else
//                 // TODO
//            //     BOOST_STATIC_ASSERT(sizeof(DataType) == 0);
//// #endif //NEKTAR_USING_LAPACK
//        }

        
        Array<TwoD, NekDouble> getVandermonde(Array<OneD, NekDouble>& x, Array<OneD, NekDouble>& y, int degree){
            int M = x.num_elements();
            int N = (degree + 1) * (degree + 2) / 2;
            
            Array<TwoD, NekDouble> V(M, N);            
            for(int i=0; i<M; ++i){
                for(int j=0; j<N; ++j){
                    V[i][j] = 0.0;
                }
            }
            int n = 0;
            for(int d=0; d<=degree; ++d){
                for(int p=0; p<=d; ++p){
                    int q = d - p;
                    Array<OneD, NekDouble> columnVector = DubinerPoly(p, q, x, y);

                    // Set n-th column of V to the DubinerPoly vector
                    for(int i=0; i<M; ++i){
                        V[i][n] = columnVector[n];
                    }
                    n = n+1;
                }
            }

            return V;
        }

        Array<TwoD, NekDouble> getVandermonde(Array<OneD, NekDouble>& x, Array<OneD, NekDouble>& y){
            int M = x.num_elements();
            int degree = getDegree(M);
            return getVandermonde( x, y, degree );
        }
        
        Array<OneD, NekDouble> DubinerPoly(int p, int q, Array<OneD, NekDouble>& x, Array<OneD, NekDouble>& y){
            // Get the coordinate transform
            int size = y.num_elements();
            Array<OneD, NekDouble> eta_1(size);

            // Initialize the horizontal coordinate of the triangle (beta in Barycentric coordinates)
            for(int el=0; el<size; ++el){
                if( y[el] < 1.0 - 1e-16 ){
                    eta_1[el] = 2.0*(1.0 + x[el]) / (1.0 - y[el]) - 1.0;
                } else {
                    eta_1[el] = -1.0; // When y is close to 1, then we have a removeable singularity
                }
            }

            // Initialize the vertical coordinate of the triangle (gamma in Barycentric coordinates)
            Array<OneD, NekDouble> eta_2 = y;

            // Orthogonal Dubiner polynomial
            int alpha = 2 * p + 1; int beta = 0;
            Array<OneD, NekDouble> psi(size);
            Array<OneD, NekDouble> psi_a = LegendrePoly(p, eta_1);
            Array<OneD, NekDouble> upsilon = JacobiPoly(q, eta_2, alpha, beta);
            NekDouble w = sqrt((2.0 * p + 1) * (p + q + 1) / 2); // Normalizing Orthonormal weight

            for(int el=0; el<size; ++el){
               NekDouble zeta   = pow((1.0 - eta_2[el])/2.0, p);
               NekDouble psi_b  = zeta * upsilon[el];
               psi[el]          = w * psi_a[el] * psi_b;
            }
            return psi;
        }

        Array<OneD, NekDouble> JacobiPoly(int degree, Array<OneD, NekDouble>& x, int alpha, int beta){
            int size = x.num_elements();
            Array<OneD, NekDouble> y(size);
            
            if (degree == 0){
                // Set y to ones
                for(int el=0; el<size; ++el){
                    y[el] = 1.0;                
                }
                
            } else if(degree == 1){
            
                for(int el=0; el<size; ++el){
                    y[el] = 0.5*(alpha - beta + (alpha + beta + 2.0) * x[el]);   
                }
                
            } else if(degree > 1){
            
                NekDouble degm1 = degree - 1.0;
                NekDouble tmp = 2.0 * degm1 + alpha + beta;
                NekDouble a1 = 2.0 * (degm1 + 1.0) * (degm1 + alpha + beta + 1.0) * tmp;
                NekDouble a2 = (tmp + 1.0) * (alpha * alpha - beta * beta);
                NekDouble a3 = tmp * (tmp + 1.0) * (tmp + 2.0);
                NekDouble a4 = 2.0 * (degm1 + alpha) * (degm1 + beta) * (tmp + 2.0);

                Array<OneD, NekDouble> z1 = JacobiPoly(degree-1, x, alpha, beta);
                Array<OneD, NekDouble> z2 = JacobiPoly(degree-2, x, alpha, beta);
                for(int el=0; el<size; ++el){
                   y[el] = ((a2 + a3 * x[el]) * z1[el] - a4 * z2[el])/a1;
                }
                
            } else {
                cerr << "Bad degree" << endl;
            }

            return y;
        }


        Array<OneD, NekDouble> LegendrePoly(int degree, Array<OneD, NekDouble>& x){
            int size = x.num_elements();
            Array<OneD, NekDouble> y(size);
            
            if(degree > 1){              
                Array<OneD, NekDouble> a0(size);
                Array<OneD, NekDouble> a1 = x;
                Array<OneD, NekDouble> a2(size);

                // Set a0 to ones
                for(int el=0; el<size; ++el){
                    a0[el] = 1.0;
                }
                for(int i=2; i<=degree; ++i){
                    NekDouble b = NekDouble(2.0*i-1.0)/i;
                    NekDouble c = NekDouble(i-1.0)/i;

                    // multiply each elements in matrix
                    for(int el=0; el<size; ++el){
                        a2[el] = b*x[el]*a1[el] - c*a0[el];
                    }
                    a0 = a1;
                    a1 = a2;
                }
                
                y = a2;

            } else if( degree == 1 ) {
                y = x;
            } else {
                // set y to ones
                for(int el=0; el<size; ++el){
                    y[el] = 1.0;
                }
            }
            return y;
        }

        int getDegree(int nBasisFunctions){
           return (-3 + int(sqrt(1.0 + 8*nBasisFunctions)))/2;
        }
        
        NekDouble  LagrangePoly(NekDouble x, int pt, int npts, const ConstArray<OneD, NekDouble>& xpts) {
            NekDouble h=1.0;

            for(int i=0;i<pt; ++i)   {
                h = h * (x - xpts[i])/(xpts[pt]-xpts[i]);
            }

            for(int i=pt+1;i<npts;++i) {
                h = h * (x - xpts[i])/(xpts[pt]-xpts[i]);
            }

            return h;
        }

//         void NodalTriFekete::CalculateWeights()
//         {
//             // Allocate the storage for points
//             CalculateWeights();
// 
//             // Compute the intitial weights for integration
//             int nPts = GetTotNumPoints();
//             SharedArray< DataType > w(nPts);
//             for( int k = 0; k < nPts; ++k ) {
//                 w[k] = 0;
//                 for( int i = 0; i <= k; ++i ) {
//                     int s0 = 1 - 2*((i+1)&1);
//                     int s1 = 1 - 2*(k&1);
//                     int s2 = 1 - 2*((k-i)&1);
//                     
//                     w[k] += double(s0) / (i+1) * ( (1 - s1)/(k + 2) - (1 + s2)/(k - i + 1));
//                 }
//             }
// 
//             // Form the Vandermonde matrix from the basis functions and point locations
//             ptr<NekMatrix< DataType, eFull > > vandermondePtr( new NekMatrix< DataType, eFull >(nPts,nPts) );
//             NekMatrix< DataType, eFull > & vandermonde = *vandermondePtr;
// 
//             
//             for( int k = 0; k < nPts; ++k ) {
//                 for( int i = 0; i < nPts; ++i ) {
//                     double x = (*m_points[0])[i], y = (*m_points[1])[i];
//                     vandermonde(i, k) = 0;
//                     for( int p = 0; p <= k; ++p ) {
//                         vandermonde(i, k) += powf(x, p) * powf(y, k-p);
//                     }
//                 }
//             }
//             // Get the inverse and transpose of the Vandermonde matrix
//             // NOTE: this is where the code fails to link
// //             vandermonde.Invert();
// //             vandermonde.Transpose();
// //             LinearSystem<NekMatrix<DataType, eFull> > linsys(vandermondePtr);
// //             //DataType *x = linsys.SolveTranspose(b).GetPtr();
// //             //NekVector<DataType> wHat( nPts, x );
// //             //ptr<NekVector<DataType> > wHat( NekVector<DataType>(nPts, w ) );
// //             //linsys.SolveTranspose(wHat, b);
// //             ptr<NekVector<DataType> > b( new NekVector<DataType>(nPts, w) );
// //             //NekVector<DataType> wHat = linsys.SolveTranspose(b);
// //             ptr<NekVector<DataType> > wHat( new NekVector<DataType>(nPts, (DataType*)0) );
// //             //NekVector<DataType> wHat(nPts, (DataType*)0);
// //             linsys.SolveTranspose(wHat, b);
// 
//             double matrix_buf[] = {81, -5, -28, 4};
//             double b_buf[] = {-941, 348};
//             ptr<NekMatrix<double, eFull> > A(new NekMatrix<double, eFull>(2, 2, matrix_buf));
//             ptr<NekVector<double> > b(new NekVector<double>(2, b_buf));
//             LinearSystem<NekMatrix<double, eFull> > linsys(A);
// //             NekVector<double> result = linsys.Solve(b);
//             
//             
//             for( int k = 0; k < nPts; ++k ) {
//                 //(*m_weights)[k] = (*wHat)[k];
//             }
// 
//             // Compute the final integration weights
// //             for( int k = 0; k < nPts; ++k ) {
// //                 m_weights[k] = 0;
// //                 for( int i = 0; i <= nPts; ++i ) {
// //                     m_weights[k] += vandermonde(k,i) * w[k];
// //                 }
// //             }
//         }

        void NodalTriFekete::CalculateDerivMatrix()
        {
            // No derivative matrix computed
        }

        ptr<NodalTriFekete::PointsBaseType> NodalTriFekete::Create(const PointsKey &key)
        {
            ptr<PointsBaseType> returnval(MemoryManager<NodalTriFekete>::AllocateSharedPtr(key));
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
                    if(m_points[0][j+1] < m_points[0][j])
                    {
                        std::swap(m_points[0][j+1], m_points[0][j]);
                        std::swap(m_points[1][j+1], m_points[1][j]);
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
                    if(m_points[0][j+1] > m_points[0][j])
                    {
                        std::swap(m_points[0][j+1], m_points[0][j]);
                        std::swap(m_points[1][j+1], m_points[1][j]);
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
                    if(m_points[1][j+1] > m_points[1][j])
                    {
                        std::swap(m_points[0][j+1], m_points[0][j]);
                        std::swap(m_points[1][j+1], m_points[1][j]);
                    }
                }
            }

            return;
        }     
    } // end of namespace stdregion
} // end of namespace stdregion

/**
* $Log: NodalTriFekete.cpp,v $
* Revision 1.15  2007/07/11 16:33:13  ehan
* Fixed bugs in Visual Studio 8
*
* Revision 1.14  2007/07/11 09:12:13  ehan
* Fixed bug
*
* Revision 1.12  2007/05/15 03:37:24  bnelson
* Updated to use the new Array object.
*
* Revision 1.11  2007/04/30 23:29:10  jfrazier
* More conversion to multi_array.
*
* Revision 1.10  2007/04/29 00:31:57  jfrazier
* Updated to use multi_arrays.
*
*/

