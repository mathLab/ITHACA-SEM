///////////////////////////////////////////////////////////////////////////////
//
// File StdHexExp.cpp
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
// Description: Heaxhedral methods
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdHexExp.h>

#ifdef max
#undef max
#endif

namespace Nektar
{
    namespace StdRegions
    {

        StdHexExp::StdHexExp() // Deafult construct of standard expansion directly called. 
        {
        }

        StdHexExp::StdHexExp(const LibUtilities::BasisKey &Ba, const LibUtilities::BasisKey &Bb, const LibUtilities::BasisKey &Bc):
        StdExpansion3D(Ba.GetNumModes()*Bb.GetNumModes()*Bc.GetNumModes(), Ba, Bb, Bc)
        {    
        }

        StdHexExp::StdHexExp(const StdHexExp &T):
        StdExpansion3D(T)
        {
        }


        // Destructor
        StdHexExp::~StdHexExp()
        {   
        }  

        //////////////////////////////
        /// Integration Methods
        //////////////////////////////

        namespace 
        {
            void TripleTensorProduct(   const ConstArray<OneD, NekDouble>& fx, 
                                        const ConstArray<OneD, NekDouble>& gy, 
                                        const ConstArray<OneD, NekDouble>& hz, 
                                        const ConstArray<OneD, NekDouble>& inarray, 
                                        Array<OneD, NekDouble> & outarray )
            {
            
            // Using matrix operation, not sum-factorization.
            // Regarding the 3D array, inarray[k][j][i], x is changing the fastest and z the slowest.
            // Thus, the first x-vector of points refers to the first row of the first stack. The first y-vector
            // refers to the first column of the first stack. The first z-vector refers to the vector of stacks
            // intersecting the first row and first column. So in C++, i refers to column, j to row, and k to stack.
            // Contrasting this with the usual C++ matrix convention, note that i does not refer to a C++ row, nor j to C++ column.

                int nx = fx.num_elements();
                int ny = gy.num_elements();
                int nz = hz.num_elements();


                // Multiply by integration constants...
                // Hadamard multiplication refers to elementwise multiplication of two vectors.

                // Hadamard each row with the first vector (x-vector); the index i is changing the fastest.
                for (int jk = 0; jk < ny*nz; ++jk)  // For each j and k, iterate over each row in all of the stacks at once
                {
                    Vmath::Vmul(
                        nx,                         // Size of first weight vector
                        &inarray[0] + jk*nx, 1,     // Offset and stride of each row-vector (x is changing fastest)
                        fx.get(), 1,                // First weight vector (with stride of 1)
                        &outarray[0] + jk*nx, 1     // Output has same offset and stride as input
                    );
                }

                // Hadamard each column with the second vector (y-vector)
                for (int k = 0; k < nz; ++k)                    // For each stack in the 3D-array, do the following...
                {
                    for (int i = 0; i < nx; ++i)                // Iterate over each column in the current stack
                    {
                        Vmath::Vmul(
                            ny,                                 // Size of second weight vector
                            &outarray[0] + i + nx*ny*k, nx,     // Offset and stride of each column-vector
                            gy.get(), 1,                        // second weight vector (with stride of 1)
                            &outarray[0] + i + nx*ny*k, nx      // Output has same offset and stride as input
                        );
                    }
                }

                // Hadamard each stack-vector with the third vector (z-vector)
                for (int ij = 0; ij < nx*ny; ++ij)              // Iterate over each element in the topmost stack
                {
                    Vmath::Vmul(
                        nz,                                     // Size of third weight vector
                        &outarray[0] + ij, nx*ny,               // Offset and stride of each stack-vector
                        hz.get(), 1,                            // Third weight vector (with stride of 1)
                        &outarray[0] + ij, nx*ny                // Output has same offset and stride as input
                    );
                }

            }
            

            // Inner-Product with respect to the weights: i.e., this is the triple sum of the product 
            // of the four inputs over the Hexahedron
            // x-dimension is the row, it is the index that changes the fastest
            // y-dimension is the column
            // z-dimension is the stack, it is the index that changes the slowest
            NekDouble TripleInnerProduct( 
                                        const ConstArray<OneD, NekDouble>& fxyz, 
                                        const ConstArray<OneD, NekDouble>& wx, 
                                        const ConstArray<OneD, NekDouble>& wy, 
                                        const ConstArray<OneD, NekDouble>& wz
                                        )
            {
                int Qx = wx.num_elements();
                int Qy = wy.num_elements();
                int Qz = wz.num_elements();

                if( fxyz.num_elements() != Qx*Qy*Qz ) {
                    cerr << "TripleTetrahedralInnerProduct expected " << fxyz.num_elements() << 
                        " quadrature points from the discretized input function but got " << 
                        Qx*Qy*Qz << " instead." << endl;
                }

                // Sum-factorizing over the stacks
                Array<OneD, NekDouble> A(Qx*Qy, 0.0);
                for( int i = 0; i < Qx; ++i ) {
                    for( int j = 0; j < Qy; ++j ) {
                        for( int k = 0; k < Qz; ++k ) {
                            A[i + Qx*j] +=  fxyz[i + Qx*(j + Qy*k)] * wz[k];
                        }
                    }
                }

                // Sum-factorizing over the columns
                Array<OneD, NekDouble> b(Qx, 0.0);
                for( int i = 0; i < Qx; ++i ) {
                    for( int j = 0; j < Qy; ++j ) {
                        b[i] +=  A[i + Qx*j] * wy[j];
                    }
                }

                // Sum-factorizing over the rows
                NekDouble c = 0;
                for( int i = 0; i < Qx; ++i ) {
                    c +=  b[i] * wx[i];
                }

                return c;
            }
        }
        
        NekDouble StdHexExp::Integral3D(const ConstArray<OneD, NekDouble>& inarray, 
                                        const ConstArray<OneD, NekDouble>& wx,
                                        const ConstArray<OneD, NekDouble>& wy, 
                                        const ConstArray<OneD, NekDouble>& wz)
        {
            return TripleInnerProduct( inarray, wx, wy, wz );

        }

	/** \brief Integrate the physical point list \a inarray over hexahedral region and return the value
            
	Inputs:\n
	
	- \a inarray: definition of function to be returned at quadrature point of expansion. 

	Outputs:\n

	- returns \f$\int^1_{-1}\int^1_{-1}\int^1_{-1} u(\xi_1, \xi_2, \xi_3) J[i,j,k] d  \xi_1 d \xi_2 d \xi_3 \f$ \n
	          \f$ = \sum_{i=0}^{Q_1 - 1} \sum_{j=0}^{Q_2 - 1} \sum_{k=0}^{Q_3 - 1} u(\xi_{1i}, \xi_{2j},\xi_{3k})w_{i} w_{j}  w_{k}   \f$ \n
	          where \f$inarray[i,j, k] = u(\xi_{1i},\xi_{2j}, \xi_{3k}) \f$ \n
	          and \f$ J[i,j,k] \f$ is the  Jacobian evaluated at the quadrature point.

        */
        NekDouble StdHexExp::Integral(const ConstArray<OneD, NekDouble>& inarray)
        {
            ConstArray<OneD, NekDouble> w0, w1, w2;

            w0 = ExpPointsProperties(0)->GetW();
            w1 = ExpPointsProperties(1)->GetW();
            w2 = ExpPointsProperties(2)->GetW();

            return Integral3D(inarray, w0, w1, w2);
        }


	/** \brief  Inner product of \a inarray over region with respect to the 
		expansion basis m_base[0]->GetBdata(),m_base[1]->GetBdata(), m_base[2]->GetBdata() and return in \a outarray 
	
		Wrapper call to StdHexExp::IProductWRTBase
	
		Input:\n
	
		- \a inarray: array of function evaluated at the physical collocation points
	
		Output:\n
	
		- \a outarray: array of inner product with respect to each basis over region

        */
        void StdHexExp::IProductWRTBase(const ConstArray<OneD, NekDouble>& inarray, Array<OneD, NekDouble> &outarray)
        {
            IProductWRTBase(m_base[0]->GetBdata(),m_base[1]->GetBdata(), m_base[2]->GetBdata(),inarray,outarray);
        }


	 /** 
		\brief Calculate the inner product of inarray with respect to
		the basis B=base0*base1*base2 and put into outarray:
		
		\f$ \begin{array}{rcl} I_{pqr} = (\phi_{pqr}, u)_{\delta} & = &
		\sum_{i=0}^{nq_0} \sum_{j=0}^{nq_1} \sum_{k=0}^{nq_2}
		\psi_{p}^{a} (\xi_{1i}) \psi_{q}^{a} (\xi_{2j}) \psi_{r}^{a} (\xi_{3k})
		w_i w_j w_k u(\xi_{1,i} \xi_{2,j} \xi_{3,k})	     
		J_{i,j,k}\\ & = & \sum_{i=0}^{nq_0} \psi_p^a(\xi_{1,i})
		\sum_{j=0}^{nq_1} \psi_{q}^a(\xi_{2,j}) \sum_{k=0}^{nq_2} \psi_{r}^a u(\xi_{1i},\xi_{2j},\xi_{3k})
		J_{i,j,k} \end{array} \f$ \n
		
		where
		
		\f$ \phi_{pqr} (\xi_1 , \xi_2 , \xi_3) = \psi_p^a ( \xi_1) \psi_{q}^a (\xi_2) \psi_{r}^a (\xi_3) \f$ \n
		
		which can be implemented as \n
		\f$f_{r} (\xi_{3k}) = \sum_{k=0}^{nq_3} \psi_{r}^a u(\xi_{1i},\xi_{2j},\xi_{3k})
		J_{i,j,k} = {\bf B_3 U}   \f$ \n
		\f$ g_{q} (\xi_{3k}) = \sum_{j=0}^{nq_1} \psi_{q}^a (\xi_{2j}) f_{r} (\xi_{3k})  = {\bf B_2 F}  \f$ \n
		\f$ (\phi_{pqr}, u)_{\delta} = \sum_{k=0}^{nq_0} \psi_{p}^a (\xi_{3k}) g_{q} (\xi_{3k})  = {\bf B_1 G} \f$

        **/
        void StdHexExp::IProductWRTBase(    const ConstArray<OneD, NekDouble>& bx, 
                                            const ConstArray<OneD, NekDouble>& by, 
                                            const ConstArray<OneD, NekDouble>& bz, 
                                            const ConstArray<OneD, NekDouble>& inarray, 
                                            Array<OneD, NekDouble> & outarray )
        {
            int     Qx = m_base[0]->GetNumPoints();
            int     Qy = m_base[1]->GetNumPoints();
            int     Qz = m_base[2]->GetNumPoints();

            int     P = m_base[0]->GetNumModes() - 1;
            int     Q = m_base[1]->GetNumModes() - 1;
            int     R = m_base[2]->GetNumModes() - 1;

            // Build an index map from the rectangle to the triangle -- This is not necessary for hexahedron
            Array<OneD, int> pq  = Array<OneD, int>( (P+1)*(Q+1), -1 );
            for( int p = 0, mode = 0; p <= P; ++p ) {
                for( int q = 0; q <= Q ; ++q, ++mode ) {
                    pq[q + (Q+1)*p] = mode;
                }
            }

            // Create an index map from the hexahedron to the tetrahedron. -- This is not necessary for hexahedron
            // The actual index is too difficult to compute explicitly.
            Array<OneD, int> pqr = Array<OneD, int>( (P+1)*(Q+1)*(R+1), -1 );
            for( int p = 0, mode = 0; p <= P; ++p ) {
                for( int q = 0; q <= Q; ++q ) {
                    for( int r = 0; r <= R; ++r, ++mode ) {
                        pqr[r + (R+1)*(q + (Q+1)*p)] = mode;
                    }
                }
            }

            // Compute innerproduct over each mode in the Hexahedron domain
            for( int p = 0; p <= P; ++p ) {
                for( int q = 0; q <= Q; ++q ) {
                    for( int r = 0; r <= R; ++r ) {

                        // Determine the index for specifying which mode to use in the basis
                        int mode_p      = p;
                        int mode_pq     = pq[q + (Q+1)*p];
                        int mode_pqr    = pqr[r + (R+1)*(q + (Q+1)*p)];

                        // Compute tensor product of inarray with the 3 basis functions
                        Array<OneD, NekDouble> g_pqr = Array<OneD, NekDouble>( Qx*Qy*Qz, 0.0 );
                        for( int k = 0; k < Qz; ++k ) {
                            for( int j = 0; j < Qy; ++j ) {
                                for( int i = 0; i < Qx; ++i ) {
                                    int s = i + Qx*(j + Qy*k);
                                     g_pqr[s] += inarray[s] * 
                                            bx[i + Qx*p] * 
                                            by[j + Qy*q] * 
                                            bz[k + Qz*r];
                                }
                            }
                        }

                        outarray[r + (R+1)*(q + (Q+1)*p)] = Integral( g_pqr );
                    }
                }
            }
        }
        
        
        ///////////////////////////////
        /// Differentiation Methods
        ///////////////////////////////
        
        void StdHexExp::PhysDeriv( Array<OneD, NekDouble> &out_d0,
                                   Array<OneD, NekDouble> &out_d1,
                                   Array<OneD, NekDouble> &out_d2)
        {
            PhysTensorDeriv(this->m_phys, out_d0, out_d1, out_d2);
        }

	 /** 
            \brief Calculate the derivative of the physical points 
            
            For Hexahedral region can use the PhysTensorDeriv function
            defined under StdExpansion.
	    Following tenserproduct:

            
        **/
        void StdHexExp::PhysDeriv(const ConstArray<OneD, NekDouble>& inarray,
                                   Array<OneD, NekDouble> &out_d0,
                                   Array<OneD, NekDouble> &out_d1,
                                   Array<OneD, NekDouble> &out_d2)
        {
            PhysTensorDeriv(inarray, out_d0, out_d1, out_d2);
        }



        //------------------------------
        /// Evaluation Methods
        //-----------------------------

	/** 
            \brief Backward transformation is evaluated at the quadrature points 
		
	    \f$ u^{\delta} (\xi_{1i}, \xi_{2j}, \xi_{3k}) = \sum_{m(pqr)} \hat u_{pqr} \phi_{pqr} (\xi_{1i}, \xi_{2j}, \xi_{3k})\f$
	    
	     Backward transformation is three dimensional tensorial expansion
		
	    \f$ u (\xi_{1i}, \xi_{2j}, \xi_{3k}) = \sum_{p=0}^{Q_x} \psi_p^a (\xi_{1i}) \lbrace { \sum_{q=0}^{Q_y} \psi_{q}^a (\xi_{2j})
	       \lbrace { \sum_{r=0}^{Q_z} \hat u_{pqr} \psi_{r}^a (\xi_{3k}) \rbrace}
	       \rbrace}. \f$
	       
	       And sumfactorizing step of the form is as:\\
	       \f$ f_{r} (\xi_{3k}) = \sum_{r=0}^{Q_z} \hat u_{pqr} \psi_{r}^a (\xi_{3k}),\\ 
	         g_{p} (\xi_{2j}, \xi_{3k}) = \sum_{r=0}^{Q_y} \psi_{p}^a (\xi_{2j}) f_{r} (\xi_{3k}),\\
		 u(\xi_{1i}, \xi_{2j}, \xi_{3k}) = \sum_{p=0}^{Q_x} \psi_{p}^a (\xi_{1i}) g_{p} (\xi_{2j}, \xi_{3k}).
	       \f$		
        **/
       void StdHexExp::BwdTrans(const ConstArray<OneD, NekDouble>& inarray, 
                                Array<OneD, NekDouble> &outarray)
        {

            ASSERTL1( (m_base[1]->GetBasisType() != LibUtilities::eOrtho_B)  ||
                      (m_base[1]->GetBasisType() != LibUtilities::eModified_B),
                "Basis[1] is not a general tensor type");

            ASSERTL1( (m_base[2]->GetBasisType() != LibUtilities::eOrtho_C) ||
                      (m_base[2]->GetBasisType() != LibUtilities::eModified_C),
                "Basis[2] is not a general tensor type");

            int     Qx = m_base[0]->GetNumPoints();
            int     Qy = m_base[1]->GetNumPoints();
            int     Qz = m_base[2]->GetNumPoints();

            int     P = m_base[0]->GetNumModes() - 1;
            int     Q = m_base[1]->GetNumModes() - 1;
            int     R = m_base[2]->GetNumModes() - 1;

            ConstArray<OneD, NekDouble> xBasis  = m_base[0]->GetBdata();
            ConstArray<OneD, NekDouble> yBasis  = m_base[1]->GetBdata();
            ConstArray<OneD, NekDouble> zBasis  = m_base[2]->GetBdata();


            // Build an index map from the rectangle to the triangle -- This is not necessary for hexahedron
            Array<OneD, int> pq  = Array<OneD, int>( (P+1)*(Q+1), -1 );
            for( int p = 0, mode = 0; p <= P; ++p ) {
                for( int q = 0; q <= Q; ++q, ++mode ) {
                    pq[q + (Q+1)*p] = mode;
                }
            }

            // Create an index map from the hexahedron to the tetrahedron.-- This is not necessary for hexahedron
            // Explicit computation of the actual index is error-prone and too difficult.
            Array<OneD, int> pqr = Array<OneD, int>( (P+1)*(Q+1)*(R+1), -1 );
            for( int p = 0, mode = 0; p <= P; ++p ) {
                for( int q = 0; q <= Q; ++q ) {
                    for( int r = 0; r <= R; ++r, ++mode ) {
                        pqr[r + (R+1)*(q + (Q+1)*p)] = mode;
                    }
                }
            }

            // Sum-factorize the triple summation starting with the z-dimension
            for( int k = 0; k < Qz; ++k ) {

                // Create the matrix of coefficients summed over the z-modes
                Array<OneD, NekDouble> Ak((P+1)*(Q+1), 0.0);
                for( int p = 0; p <= P; ++p ) {
                    for( int q = 0; q <= Q; ++q ) {
                        for( int r = 0; r <= R; ++r ) {
                             int mode = r + (R+1)*(q + (Q+1)*p);
                            Ak[q + (Q+1)*p]   +=   inarray[mode]  *  zBasis[k + Qz*r];
                        }
                    }
                }

                // Factorize the y-dimension
                for( int j = 0; j < Qy; ++j ) {

                    // Create the vector of coefficients summed over the y and z-modes
                    Array<OneD, NekDouble> bjk(P+1, 0.0);
                    for( int p = 0; p <= P; ++p ) {
                        for( int q = 0; q <= Q; ++q ) {
                            int mode = q + (Q+1)*p;
                            bjk[p]   +=   Ak[q + (Q+1)*p]  *  yBasis[j + Qy*q];
                        }
                    }

                    // Factorize the x-dimension
                    for( int i = 0; i < Qx; ++i ) {

                        NekDouble cijk = 0.0;
                        for( int p = 0; p <= P; ++p ) {
                            int mode = p;
                            cijk   +=   bjk[p]  *  xBasis[i + Qx*p];
                        }
                        outarray[i + Qx*(j + Qy*k)] = cijk;
                    }
                }
            }
        }


 	/** \brief Forward transform from physical quadrature space
            stored in \a inarray and evaluate the expansion coefficients and
            store in \a (this)->m_coeffs  
            
            Inputs:\n
            
            - \a inarray: array of physical quadrature points to be transformed
            
            Outputs:\n
            
            - (this)->_coeffs: updated array of expansion coefficients. 
            
        */    
        void StdHexExp::FwdTrans(const ConstArray<OneD, NekDouble>& inarray,
                                 Array<OneD, NekDouble> &outarray)
        {
            if((m_base[0]->Collocation())&&(m_base[1]->Collocation())&&(m_base[2]->Collocation()))
            {
                Vmath::Vcopy(GetNcoeffs(), &inarray[0], 1, &outarray[0], 1);
            }
            else
            {
                IProductWRTBase(inarray,outarray);
                
                // get Mass matrix inverse
                StdMatrixKey      masskey(eInvMass,DetShapeType(),*this);
                DNekMatSharedPtr matsys = GetStdMatrix(masskey);
                
                // copy inarray in case inarray == outarray
                DNekVec in (m_ncoeffs,outarray);
                DNekVec out(m_ncoeffs,outarray,eWrapper);

                out = (*matsys)*in;

            }
        }

        /// Single Point Evaluation
        NekDouble StdHexExp::PhysEvaluate(ConstArray<OneD, NekDouble>& coords)
        {
             return  StdExpansion3D::PhysEvaluate3D(coords);  
        }
        
       void StdHexExp::GetCoords( Array<OneD, NekDouble> & xi_x, Array<OneD, NekDouble> & xi_y, Array<OneD, NekDouble> & xi_z)
        {
            ConstArray<OneD, NekDouble> eta_x = ExpPointsProperties(0)->GetZ();
            ConstArray<OneD, NekDouble> eta_y = ExpPointsProperties(1)->GetZ();
            ConstArray<OneD, NekDouble> eta_z = ExpPointsProperties(2)->GetZ();
            int Qx = GetNumPoints(0);
            int Qy = GetNumPoints(1);
            int Qz = GetNumPoints(2);

            // Convert collapsed coordinates into cartesian coordinates: eta --> xi
            for( int k = 0; k < Qz; ++k ) {
                for( int j = 0; j < Qy; ++j ) {
                    for( int i = 0; i < Qx; ++i ) {
                            int s = i + Qx*(j + Qy*k);
                            xi_x[s] = eta_x[i];
                            xi_y[s] = eta_y[j];
                            xi_z[s] = eta_z[k];

                    }
                }
            }
        }

        void StdHexExp::FillMode(const int mode, Array<OneD, NekDouble> &outarray)
        {
            int    i,j;
            int   nquad0 = m_base[0]->GetNumPoints();
            int   nquad1 = m_base[1]->GetNumPoints();
            int   nquad2 = m_base[2]->GetNumPoints();
            
             ConstArray<OneD, NekDouble> base0  = m_base[0]->GetBdata();
             ConstArray<OneD, NekDouble> base1  = m_base[1]->GetBdata();
             ConstArray<OneD, NekDouble> base2  = m_base[2]->GetBdata();
            
            int   btmp0 = m_base[0]->GetNumModes();
            int   btmp1 = m_base[1]->GetNumModes();
            int   mode2 = mode/(btmp0*btmp1);
            int   mode1 = (mode-mode2*btmp0*btmp1)/btmp0;
            int   mode0 = (mode-mode2*btmp0*btmp1)%btmp0;

            ASSERTL2(mode2 == (int)floor((1.0*mode)/(btmp0*btmp1)),
                "Integer Truncation not Equiv to Floor");
            ASSERTL2(mode1 == (int)floor((1.0*mode-mode2*btmp0*btmp1)/(btmp0*btmp1)),
                "Integer Truncation not Equiv to Floor");
            ASSERTL2(m_ncoeffs <= mode,
                "calling argument mode is larger than total expansion order");

            for(i = 0; i < nquad1*nquad2; ++i)
            {
                Vmath::Vcopy(nquad0,(double *)(base0.get() + mode0*nquad0),1,
                             &outarray[0]+i*nquad0, 1);
            }

            for(j = 0; j < nquad2; ++j)
            {
                for(i = 0; i < nquad0; ++i)
                {
                    Vmath::Vmul(nquad1,(double *)(base1.get() + mode1*nquad1),1,
                        &outarray[0]+i+j*nquad0*nquad1, nquad0,
                        &outarray[0]+i+j*nquad0*nquad1, nquad0);
                }
            }

            for(i = 0; i < nquad2; i++)
            {
                Blas::Dscal(nquad0*nquad1,base2[mode2*nquad2+i],
                            &outarray[0]+i*nquad0*nquad1,1);
            }
        }
   
        DNekMatSharedPtr StdHexExp::GenMatrixHex(const StdMatrixKey &mkey)
        {
            int      i,j;
         
            int      order0    = GetBasisNumModes(0);
            int      order1    = GetBasisNumModes(1);
            int      order2    = GetBasisNumModes(2);
            int      tot_order = GetNcoeffs();

            MatrixType  mtype = mkey.GetMatrixType();

             //StdExpansion::GenerateMassMatrix(outarray);
            DNekMatSharedPtr Mat = StdExpansion::CreateGeneralMatrix(mkey);


        switch(mtype)
        {
        case eMass:
            // For Fourier basis set the imaginary component of mean mode
            // to have a unit diagonal component in mass matrix 
            if(m_base[0]->GetBasisType() == LibUtilities::eFourier)
            {
                for(i = 0; i < order1*order2; ++i)
                {
//                     outarray[(order0*i+1)*tot_order+i*order0+1] = 1.0;
                         //(*Mat)((order0*i+1)*tot_order+i*order0+1) = 1.0;
                }
            }

            if(m_base[1]->GetBasisType() == LibUtilities::eFourier)
            {
                for(j = 0; j < order2; ++j)
                {
                    for(i = 0; i < order0; ++i)
                    {
                        //(*Mat)((order0+i)*tot_order+order0+i+j*(order0*order1)*(tot_order+1)) = 1.0;
                    }
                }
            }

            if(m_base[2]->GetBasisType() == LibUtilities::eFourier)
            {
                for(i = 0; i < order0*order1; ++i)
                {
                    //(*Mat)((order0*order1)*(tot_order+1)+i*tot_order +i) = 1.0;
                }
            }
            break;
            
          }
          
          return Mat;
        }

        void StdHexExp::GenLapMatrix(double * outarray)
        {
            ASSERTL0(false, "Not implemented");
        }




  
    }//end namespace
}//end namespace

/** 
* $Log: StdHexExp.cpp,v $
* Revision 1.12  2008/01/08 22:30:31  ehan
* Clean up the codes.
*
* Revision 1.11  2008/01/03 15:44:38  ehan
* Fixed bug.
*
* Revision 1.10  2007/12/17 13:03:51  sherwin
* Modified StdMatrixKey to contain a list of constants and GenMatrix to take a StdMatrixKey
*
* Revision 1.9  2007/12/01 00:52:12  ehan
* Completed implementing and testing following functions:
* Integral, IProductWRTBase, PhysDeriv. BwdTrans, FwdTrans, and PhysEvaluate.
*
* Revision 1.8  2007/10/15 20:38:41  ehan
* Make changes of column major matrix
*
* Revision 1.7  2007/07/20 02:16:54  bnelson
* Replaced boost::shared_ptr with Nektar::ptr
*
* Revision 1.6  2007/01/18 18:44:45  bnelson
* Updates to compile on Visual Studio 2005.
*
* Revision 1.5  2007/01/17 16:36:57  pvos
* updating doxygen documentation
*
* Revision 1.4  2007/01/17 16:05:40  pvos
* updated doxygen documentation
*
* Revision 1.3  2006/12/10 19:00:54  sherwin
* Modifications to handle nodal expansions
*
* Revision 1.2  2006/06/01 14:13:36  kirby
* *** empty log message ***
*
* Revision 1.1  2006/05/04 18:58:31  kirby
* *** empty log message ***
*
* Revision 1.23  2006/04/25 20:23:33  jfrazier
* Various fixes to correct bugs, calls to ASSERT, etc.
*
* Revision 1.22  2006/04/01 21:59:27  sherwin
* Sorted new definition of ASSERT
*
* Revision 1.21  2006/03/21 09:21:32  sherwin
* Introduced NekMemoryManager
*
* Revision 1.20  2006/03/05 22:11:02  sherwin
*
* Sorted out Project1D, Project2D and Project_Diff2D as well as some test scripts
*
* Revision 1.19  2006/03/01 08:25:03  sherwin
*
* First compiling version of StdRegions
*
* Revision 1.18  2006/02/27 23:47:23  sherwin
*
* Standard coding update upto compilation of StdHexExp.cpp
*
*
**/ 


