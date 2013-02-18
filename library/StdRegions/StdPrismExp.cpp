///////////////////////////////////////////////////////////////////////////////
//
// File StdPrismExp.cpp
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
// Description: Prismatic routines built upon StdExpansion3D
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdPrismExp.h>

namespace Nektar
{
    namespace StdRegions
    {

        StdPrismExp::StdPrismExp() // Deafult construct of standard expansion directly called. 
        {
        }

        StdPrismExp::StdPrismExp(const LibUtilities::BasisKey &Ba, 
                                 const LibUtilities::BasisKey &Bb, 
                                 const LibUtilities::BasisKey &Bc) 
            : StdExpansion  (StdPrismData::getNumberOfCoefficients(
                                 Ba.GetNumModes(), 
                                 Bb.GetNumModes(), 
                                 Bc.GetNumModes()),
                             3,Ba,Bb,Bc),
              StdExpansion3D(StdPrismData::getNumberOfCoefficients(
                                 Ba.GetNumModes(), 
                                 Bb.GetNumModes(), 
                                 Bc.GetNumModes()), 
                             Ba,Bb,Bc)
        {
            ASSERTL0(Ba.GetNumModes() <= Bc.GetNumModes(), 
                     "order in 'a' direction is higher than order in 'c' direction");
        }

        StdPrismExp::StdPrismExp(const StdPrismExp &T)
            : StdExpansion(T),
              StdExpansion3D(T)
        {
        }


        // Destructor
        StdPrismExp::~StdPrismExp()
        {   
        } 


        //---------------------------------------
        // Miscellaneous public 3D function
        //---------------------------------------
        
        void StdPrismExp::WriteCoeffsToFile(std::ofstream &outfile)
        {
            int order0 = m_base[0]->GetNumModes();
            int order1 = m_base[1]->GetNumModes();
            int order2 = m_base[2]->GetNumModes();

            Array<OneD, NekDouble> wsp(order0*order1*order2, 0.0);

            NekDouble *mat = wsp.get(); 

            // put coeffs into matrix and reverse order so that r index is
            // fastest for Prism
            Vmath::Zero(order0*order1*order2, mat, 1);

            for(int i = 0, cnt=0; i < order0; ++i)
            {
                for(int j = 0; j < order1-i; ++j)
                {
                    for(int k = 0; k < order2-i-j; ++k, cnt++)
                    {
                        // mat[i+j*order1] = m_coeffs[cnt];
                        mat[i + order1*(j + order2*k)] = m_coeffs[cnt];
                    }
                }
            }

            outfile <<"Coeffs = [" << " "; 

            for(int k = 0; k < order2; ++k)
            {            
                for(int j = 0; j < order1; ++j)
                {
                    for(int i = 0; i < order0; ++i)
                    {
                        outfile << mat[i + order0*(j + order1*k)] << " ";
                    }
                    outfile << std::endl; 
                }
            }
            outfile << "]"; 
        }
        
        
        //---------------------------------------
        // Integration Methods
        //---------------------------------------
        
        void StdPrismExp::TripleTensorProduct(const Array<OneD, const NekDouble>& fx, 
                                              const Array<OneD, const NekDouble>& gy, 
                                              const Array<OneD, const NekDouble>& hz, 
                                              const Array<OneD, const NekDouble>& inarray, 
                                                    Array<OneD,       NekDouble>& outarray)
        {
            // Using matrix operation, not sum-factorization.  Regarding the
            // 3D array, inarray[k][j][i], x is changing the fastest and z the
            // slowest.  Thus, the first x-vector of points refers to the
            // first row of the first stack. The first y-vector refers to the
            // first column of the first stack. The first z-vector refers to
            // the vector of stacks intersecting the first row and first
            // column. So in C++, i refers to column, j to row, and k to
            // stack.  Contrasting this with the usual C++ matrix convention,
            // note that i does not refer to a C++ row, nor j to C++ column.

            int nx = fx.num_elements();
            int ny = gy.num_elements();
            int nz = hz.num_elements();
            
            // Multiply by integration constants...  Hadamard multiplication
            // refers to elementwise multiplication of two vectors.  Hadamard
            // each row with the first vector (x-vector); the index i is
            // changing the fastest.
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
        
        // Inner-Product with respect to the weights: i.e., this is the triple
        // sum of the product of the four inputs over the prism. x-dimension
        // is the row, it is the index that changes the fastest; y-dimension
        // is the column; z-dimension is the stack, it is the index that
        // changes the slowest.
        NekDouble StdPrismExp::TripleInnerProduct(const Array<OneD, const NekDouble>& fxyz, 
                                                  const Array<OneD, const NekDouble>& wx, 
                                                  const Array<OneD, const NekDouble>& wy, 
                                                  const Array<OneD, const NekDouble>& wz)
        {
            int Qx = wx.num_elements();
            int Qy = wy.num_elements();
            int Qz = wz.num_elements();

            if (fxyz.num_elements() != Qx*Qy*Qz) 
            {
                cerr << "TripleInnerProduct expected " << fxyz.num_elements() 
                     << " quadrature points from the discretized input function but got " 
                     << Qx*Qy*Qz << " instead." << endl;
            }
            
            // Sum-factorizing over the stacks
            Array<OneD, NekDouble> A(Qx*Qy, 0.0);
            for (int i = 0; i < Qx; ++i)
            {
                for (int j = 0; j < Qy; ++j)
                {
                    for (int k = 0; k < Qz; ++k)
                    {
                        A[i + Qx*j] += fxyz[i + Qx*(j + Qy*k)] * wz[k];
                    }
                }
            }
            
            // Sum-factorizing over the columns
            Array<OneD, NekDouble> b(Qx, 0.0);
            for (int i = 0; i < Qx; ++i)
            {
                for (int j = 0; j < Qy; ++j)
                {
                    b[i] += A[i + Qx*j] * wy[j];
                }
            }
            
            // Sum-factorizing over the rows
            NekDouble c = 0;
            for (int i = 0; i < Qx; ++i)
            {
                c += b[i] * wx[i];
            }
            
            return c;
        }

        NekDouble StdPrismExp::Integral3D(
            const Array<OneD, const NekDouble>& inarray, 
            const Array<OneD, const NekDouble>& wx,
            const Array<OneD, const NekDouble>& wy, 
            const Array<OneD, const NekDouble>& wz)
        {
            return TripleInnerProduct(inarray, wx, wy, wz);
        }

        /** 
         * \brief Integrate the physical point list \a inarray over prismatic
         * region and return the value.
         *
         * Inputs:\n
         *
         * - \a inarray: definition of function to be returned at quadrature
         *    point of expansion.
         *
         * Outputs:\n
         *
         * - returns \f$\int^1_{-1}\int^1_{-1}\int^1_{-1} u(\bar \eta_1,
         * \xi_2, \xi_3) J[i,j,k] d \bar \eta_1 d \xi_2 d \xi_3 \f$ \n \f$ =
         * \sum_{i=0}^{Q_1 - 1} \sum_{j=0}^{Q_2 - 1} \sum_{k=0}^{Q_3 - 1}
         * u(\bar \eta_{1i}^{0,0}, \xi_{2j}^{0,0},\xi_{3k}^{1,0})w_{i}^{0,0}
         * w_{j}^{0,0} \hat w_{k}^{1,0} \f$ \n where \f$ inarray[i,j, k] =
         * u(\bar \eta_{1i}^{0,0}, \xi_{2j}^{0,0},\xi_{3k}^{1,0}) \f$, \n
         * \f$\hat w_{i}^{1,0} = \frac {w_{j}^{1,0}} {2} \f$ \n and \f$
         * J[i,j,k] \f$ is the Jacobian evaluated at the quadrature point.
         */
        NekDouble StdPrismExp::v_Integral(
            const Array<OneD, const NekDouble>& inarray)
        {
            // Using implementation from page 146 of Spencer Sherwin's book.
            int Qz = m_base[2]->GetNumPoints();

            // Get the point distributions:
            // * x is assumed to be Gauss-Lobatto-Legendre (incl. -1 and 1)
            // * y is assumed to be Gauss-Lobatto-Legendre (incl. -1 and 1)
            Array<OneD, const NekDouble> z, wx, wy, wz;
            wx = m_base[0]->GetW();
            wy = m_base[1]->GetW();
            m_base[2]->GetZW(z,wz);

            Array<OneD, NekDouble> wz_hat = Array<OneD, NekDouble>(Qz, 0.0);

            // Convert wz into wz_hat, which includes the 1/2 scale factor.
            // Nothing else need be done if the point distribution is Jacobi
            // (1,0) since (1 - xi_z) is already factored into the weights.
            // Note by coincidence, xi_y = eta_y, xi_z = eta_z (xi_z = z
            // according to our notation)
            switch(m_base[2]->GetPointsType())
            {
                // Common case
                case LibUtilities::eGaussRadauMAlpha1Beta0: // (1,0) Jacobi Inner product
                    Vmath::Smul(Qz, 0.5, (NekDouble *)wz.get(), 1, wz_hat.get(), 1);
                    break;
                
                // Corner cases
                case LibUtilities::eGaussLobattoLegendre:
                case LibUtilities::eGaussRadauMLegendre:
                    for (int k = 0; k < Qz; ++k)
                    {
                        wz_hat[k] = 0.5*(1.0 - z[k]) * wz[k];
                    }
                    break;
                    
                default:
                    ASSERTL0(false, "Unsupported quadrature points type.");
                    break;
            }


            return Integral3D(inarray, wx, wy, wz_hat);
        }


        //---------------------------------------
        // Differentiation Methods
        //---------------------------------------

        /**
         * \brief Calculate the derivative of the physical points 
         *  
         * The derivative is evaluated at the nodal physical points.
         * Derivatives with respect to the local Cartesian coordinates.
         *  
         * \f$\begin{Bmatrix} \frac {\partial} {\partial \xi_1} \\ \frac
         * {\partial} {\partial \xi_2} \\ \frac {\partial} {\partial \xi_3}
         * \end{Bmatrix} = \begin{Bmatrix} \frac 2 {(1-\eta_3)} \frac \partial
         * {\partial \bar \eta_1} \\ \frac {\partial} {\partial \xi_2} \ \
         * \frac {(1 + \bar \eta_1)} {(1 - \eta_3)} \frac \partial {\partial
         * \bar \eta_1} + \frac {\partial} {\partial \eta_3} \end{Bmatrix}\f$
         */

        void StdPrismExp::v_PhysDeriv(const Array<OneD, const NekDouble>& u_physical, 
                                      Array<OneD, NekDouble> &out_dxi1, 
                                      Array<OneD, NekDouble> &out_dxi2,
                                      Array<OneD, NekDouble> &out_dxi3 )
        {
            int    Qx = m_base[0]->GetNumPoints();
            int    Qy = m_base[1]->GetNumPoints();
            int    Qz = m_base[2]->GetNumPoints();
            int    Qtot = Qx*Qy*Qz;

            Array<OneD, NekDouble> dEta_bar1(Qtot,0.0);

            Array<OneD, const NekDouble> eta_x, eta_z;
            eta_x = m_base[0]->GetZ();
            eta_z = m_base[2]->GetZ();

            int i, k;

            bool Do_1 = (out_dxi1.num_elements() > 0)? true:false;
            bool Do_3 = (out_dxi3.num_elements() > 0)? true:false;

            // out_dXi2 is just a tensor derivative so is just passed through
            if(Do_3) 
            {
                PhysTensorDeriv(u_physical, dEta_bar1, out_dxi2, out_dxi3);
            }
            else if(Do_1)
            {
                PhysTensorDeriv(u_physical, dEta_bar1, out_dxi2, NullNekDouble1DArray);
            }
            else // case if just require 2nd direction 
            {
                PhysTensorDeriv(u_physical, NullNekDouble1DArray, 
                                out_dxi2, NullNekDouble1DArray);
            }
            
            if(Do_1)
            {
                for (k = 0; k < Qz; ++k)
                {
                    Vmath::Smul(Qx*Qy,2.0/(1.0-eta_z[k]),&dEta_bar1[0] + k*Qx*Qy,1,
                                &out_dxi1[0] + k*Qx*Qy,1);
                }
            }
            
            if(Do_3)
            {
                // divide dEta_Bar1 by (1-eta_z) 
                for (k = 0; k < Qz; ++k)
                {
                    Vmath::Smul(Qx*Qy, 1.0/(1.0-eta_z[k]),&dEta_bar1[0]+k*Qx*Qy,1,
                                &dEta_bar1[0]+k*Qx*Qy,1);
                }

                // Multiply dEta_Bar1 by (1+eta_x) and add ot out_dxi3
                for (i = 0; i < Qx; ++i)
                {
                    Vmath::Svtvp(Qz*Qy,1.0+eta_x[i],&dEta_bar1[0]+i,Qx,
                                 &out_dxi3[0]+i,Qx,&out_dxi3[0]+i,Qx);
                }                

            }
        }

        void StdPrismExp::v_PhysDeriv(const int dir,
                                      const Array<OneD, const NekDouble>& inarray,
                                            Array<OneD,       NekDouble>& outarray)
        {
            switch(dir)
            {
                case 0:
                {
                    v_PhysDeriv(inarray, outarray, NullNekDouble1DArray,
                                NullNekDouble1DArray);
                    break;
                }
                
                case 1:
                {
                    v_PhysDeriv(inarray, NullNekDouble1DArray, outarray,
                                NullNekDouble1DArray);
                    break;
                }
                
                case 2:
                {
                    v_PhysDeriv(inarray, NullNekDouble1DArray,
                                NullNekDouble1DArray, outarray);
                    break;
                }
                
                default:
                {
                    ASSERTL1(false,"input dir is out of range");
                }
                break;
            }
        }

        void StdPrismExp::v_StdPhysDeriv(const Array<OneD, const NekDouble>& inarray, 
                                               Array<OneD,       NekDouble>& out_d0,
                                               Array<OneD,       NekDouble>& out_d1,
                                               Array<OneD,       NekDouble>& out_d2)
        {
            StdPrismExp::v_PhysDeriv(inarray, out_d0, out_d1, out_d2);
        }
        
        
        //---------------------------------------
        // Transforms
        //---------------------------------------
        
	/**
         * @note 'r' (base[2]) runs fastest in this element.
         *
         * Perform backwards transformation at the quadrature points:
         *
	 * \f$ u^{\delta} (\xi_{1i}, \xi_{2j}, \xi_{3k}) = \sum_{m(pqr)} \hat
	 *  u_{pqr} \phi_{pqr} (\xi_{1i}, \xi_{2j}, \xi_{3k})\f$
         *
         * In the prism this expansion becomes:
         *
	 * \f$ u (\xi_{1i}, \xi_{2j}, \xi_{3k}) = \sum_{p=0}^{Q_x} \psi_p^a
         *  (\xi_{1i}) \lbrace { \sum_{q=0}^{Q_y} \psi_{q}^a (\xi_{2j})
         *  \lbrace { \sum_{r=0}^{Q_z} \hat u_{pqr} \psi_{pr}^b (\xi_{3k})
         *  \rbrace} \rbrace}. \f$
         *
         * And sumfactorizing step of the form is as:\\
         *
         * \f$ f_{pr} (\xi_{3k}) = \sum_{r=0}^{Q_z} \hat u_{pqr} \psi_{pr}^b
         * (\xi_{3k}),\\
         *
         * g_{p} (\xi_{2j}, \xi_{3k}) = \sum_{r=0}^{Q_y} \psi_{p}^a (\xi_{2j})
         * f_{pr} (\xi_{3k}),\ \
         *
         * u(\xi_{1i}, \xi_{2j}, \xi_{3k}) = \sum_{p=0}^{Q_x} \psi_{p}^a
         *  (\xi_{1i}) g_{p} (\xi_{2j}, \xi_{3k}).  \f$
         */
        void StdPrismExp::v_BwdTrans(const Array<OneD, const NekDouble>& inarray, 
                                           Array<OneD,       NekDouble>& outarray)
        {
            ASSERTL1((m_base[1]->GetBasisType() != LibUtilities::eOrtho_B)  ||
                     (m_base[1]->GetBasisType() != LibUtilities::eModified_B),
                     "Basis[1] is not a general tensor type");

            ASSERTL1((m_base[2]->GetBasisType() != LibUtilities::eOrtho_C) ||
                     (m_base[2]->GetBasisType() != LibUtilities::eModified_C),
                     "Basis[2] is not a general tensor type");

            if(m_base[0]->Collocation() && 
               m_base[1]->Collocation() &&
               m_base[2]->Collocation())
            {
                Vmath::Vcopy(m_base[0]->GetNumPoints()*
                             m_base[1]->GetNumPoints()*
                             m_base[2]->GetNumPoints(),
                             inarray, 1, outarray, 1);
            }
            else
            {
                StdPrismExp::v_BwdTrans_SumFac(inarray,outarray);
            }
        }

        void StdPrismExp::v_BwdTrans_SumFac(const Array<OneD, const NekDouble>& inarray, 
                                                  Array<OneD,       NekDouble>& outarray)
        {
            int  nquad1 = m_base[1]->GetNumPoints();
            int  nquad2 = m_base[2]->GetNumPoints();
            int  order0 = m_base[0]->GetNumModes();
            int  order1 = m_base[1]->GetNumModes();
            
            Array<OneD, NekDouble> wsp(nquad2*order1*order0 +
                                       nquad1*nquad2*order0);
            
            BwdTrans_SumFacKernel(m_base[0]->GetBdata(),
                                  m_base[1]->GetBdata(),
                                  m_base[2]->GetBdata(),
                                  inarray,outarray,wsp,true,true,true);
        }


        void StdPrismExp::BwdTrans_SumFacKernel(
            const Array<OneD, const NekDouble> &base0,
            const Array<OneD, const NekDouble> &base1,
            const Array<OneD, const NekDouble> &base2,
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
                  Array<OneD,       NekDouble> &wsp,
            bool                                doCheckCollDir0,
            bool                                doCheckCollDir1,
            bool                                doCheckCollDir2)
        { 
            int i, mode;
            int nquad0    = m_base[0]->GetNumPoints();
            int nquad1    = m_base[1]->GetNumPoints();
            int nquad2    = m_base[2]->GetNumPoints();
            int nummodes0 = m_base[0]->GetNumModes();
            int nummodes1 = m_base[1]->GetNumModes();
            int nummodes2 = m_base[2]->GetNumModes();
            Array<OneD, NekDouble> tmp0 = wsp;
            Array<OneD, NekDouble> tmp1 = tmp0 + nquad2*nummodes1*nummodes0;

            for (i = mode = 0; i < nummodes0; ++i)
            {
                Blas::Dgemm('N', 'N', nquad2, nummodes1, nummodes2-i, 
                            1.0, base2.get()   + mode*nquad2,        nquad2, 
                                 inarray.get() + mode*nummodes1,     nummodes2-i, 
                            0.0, tmp0.get()    + i*nquad2*nummodes1, nquad2);
                mode += nummodes2-i;
            }
            
            if (m_base[0]->GetBasisType() == LibUtilities::eModified_A)
            {
                for(i = 0; i < nummodes1; i++)
                {
                    Blas::Daxpy(nquad2,inarray[1+i*nummodes2],base2.get()+nquad2,1,
                                tmp0.get()+nquad2*(nummodes1+i),1);
                }
            }

            for (i = 0; i < nummodes0; i++)
            {
                Blas::Dgemm('N', 'T', nquad1, nquad2, nummodes1,
                            1.0, base1.get(),                     nquad1,
                                 tmp0.get() + i*nquad2*nummodes1, nquad2, 
                            0.0, tmp1.get() + i*nquad2*nquad1,    nquad1);
            }

            Blas::Dgemm('N', 'T', nquad0, nquad2*nquad1, nummodes0, 
                        1.0, base0.get(),    nquad0, 
                             tmp1.get(),     nquad2*nquad1, 
                        0.0, outarray.get(), nquad0);
        }

	/** 
         * \brief Forward transform from physical quadrature space stored in
         * \a inarray and evaluate the expansion coefficients and store in \a
         * (this)->m_coeffs
         *  
         *  Inputs:\n
         *  - \a inarray: array of physical quadrature points to be transformed
         * 
         * Outputs:\n
         *  - (this)->_coeffs: updated array of expansion coefficients. 
         */
        void StdPrismExp::v_FwdTrans(const Array<OneD, const NekDouble>& inarray,
                                           Array<OneD,       NekDouble>& outarray)
        {
            v_IProductWRTBase(inarray, outarray);

            // Get Mass matrix inverse
            StdMatrixKey      masskey(eInvMass,DetExpansionType(),*this);
            DNekMatSharedPtr  matsys = GetStdMatrix(masskey);
            
            // copy inarray in case inarray == outarray
            DNekVec in (m_ncoeffs, outarray);
            DNekVec out(m_ncoeffs, outarray, eWrapper);

            out = (*matsys)*in;
        }
        
        
        //---------------------------------------
        // Inner product functions
        //---------------------------------------
        
        /** 
         * \brief Calculate the inner product of inarray with respect to the
         * basis B=base0*base1*base2 and put into outarray:
         * 
         * \f$ \begin{array}{rcl} I_{pqr} = (\phi_{pqr}, u)_{\delta} & = &
         * \sum_{i=0}^{nq_0} \sum_{j=0}^{nq_1} \sum_{k=0}^{nq_2} \psi_{p}^{a}
         * (\bar \eta_{1i}) \psi_{q}^{a} (\xi_{2j}) \psi_{pr}^{b} (\xi_{3k})
         * w_i w_j w_k u(\bar \eta_{1,i} \xi_{2,j} \xi_{3,k}) J_{i,j,k}\\ & =
         * & \sum_{i=0}^{nq_0} \psi_p^a(\bar \eta_{1,i}) \sum_{j=0}^{nq_1}
         * \psi_{q}^a(\xi_{2,j}) \sum_{k=0}^{nq_2} \psi_{pr}^b u(\bar
         * \eta_{1i},\xi_{2j},\xi_{3k}) J_{i,j,k} \end{array} \f$ \n
         *  
         * where
         *
         * \f$ \phi_{pqr} (\xi_1 , \xi_2 , \xi_3) = \psi_p^a (\bar \eta_1)
         * \psi_{q}^a (\xi_2) \psi_{pr}^b (\xi_3) \f$ \n
         * 
         * which can be implemented as \n 
         *  
         * \f$f_{pr} (\xi_{3k}) = \sum_{k=0}^{nq_3} \psi_{pr}^b u(\bar
         * \eta_{1i},\xi_{2j},\xi_{3k}) J_{i,j,k} = {\bf B_3 U} \f$ \n \f$
         * g_{q} (\xi_{3k}) = \sum_{j=0}^{nq_1} \psi_{q}^a (\xi_{2j}) f_{pr}
         * (\xi_{3k}) = {\bf B_2 F} \f$ \n \f$ (\phi_{pqr}, u)_{\delta} =
         * \sum_{k=0}^{nq_0} \psi_{p}^a (\xi_{3k}) g_{q} (\xi_{3k}) = {\bf B_1
         * G} \f$
         */
        void StdPrismExp::v_IProductWRTBase(
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            ASSERTL1( (m_base[1]->GetBasisType() != LibUtilities::eOrtho_B)  ||
                      (m_base[1]->GetBasisType() != LibUtilities::eModified_B),
                      "Basis[1] is not a general tensor type");

            ASSERTL1( (m_base[2]->GetBasisType() != LibUtilities::eOrtho_C) ||
                      (m_base[2]->GetBasisType() != LibUtilities::eModified_C),
                      "Basis[2] is not a general tensor type");

            if(m_base[0]->Collocation() && m_base[1]->Collocation())
            {
                MultiplyByQuadratureMetric(inarray,outarray);
            }
            else
            {
                StdPrismExp::v_IProductWRTBase_SumFac(inarray,outarray);
            }
        }
        
        /**
         * Implementation of the local matrix inner product operation.
         */
        void StdPrismExp::v_IProductWRTBase_MatOp(
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            int nq = GetTotPoints();
            StdMatrixKey      iprodmatkey(eIProductWRTBase,DetExpansionType(),*this);
            DNekMatSharedPtr  iprodmat = GetStdMatrix(iprodmatkey);

            Blas::Dgemv('N',m_ncoeffs,nq,1.0,iprodmat->GetPtr().get(),
                        m_ncoeffs, inarray.get(), 1, 0.0, outarray.get(), 1);
        }

        void StdPrismExp::v_IProductWRTBase_SumFac(
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            int nquad1 = m_base[1]->GetNumPoints();
            int nquad2 = m_base[2]->GetNumPoints();
            int order0 = m_base[0]->GetNumModes();
            int order1 = m_base[1]->GetNumModes();

            Array<OneD, NekDouble> tmp(inarray.num_elements());
            Array<OneD, NekDouble> wsp(order0*nquad2*(nquad1+order1));
            
            MultiplyByQuadratureMetric(inarray,tmp);
            
            IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                         m_base[1]->GetBdata(),
                                         m_base[2]->GetBdata(),
                                         tmp,outarray,wsp,
                                         true,true,true);
        }
        
        void StdPrismExp::IProductWRTBase_SumFacKernel(
            const Array<OneD, const NekDouble>& base0,
            const Array<OneD, const NekDouble>& base1,
            const Array<OneD, const NekDouble>& base2,
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble> &outarray,
                  Array<OneD,       NekDouble> &wsp,
            bool                                doCheckCollDir0,
            bool                                doCheckCollDir1,
            bool                                doCheckCollDir2)
        {
            // Interior prism implementation based on Spen's book page
            // 119. and 608.
            const int nquad0 = m_base[0]->GetNumPoints();
            const int nquad1 = m_base[1]->GetNumPoints();
            const int nquad2 = m_base[2]->GetNumPoints();
            const int order0 = m_base[0]->GetNumModes ();
            const int order1 = m_base[1]->GetNumModes ();
            const int order2 = m_base[2]->GetNumModes ();
            
            int i, mode;
            
            ASSERTL1(wsp.num_elements() >= nquad1*nquad2*order0 + 
                                           nquad2*order0*order1,
                     "Insufficient workspace size");
            
            Array<OneD, NekDouble> tmp0 = wsp;
            Array<OneD, NekDouble> tmp1 = wsp + nquad1*nquad2*order0;
            
            // Inner product with respect to the '0' direction
            Blas::Dgemm('T', 'N', nquad1*nquad2, order0, nquad0, 
                        1.0, inarray.get(), nquad0,
                             base0.get(),   nquad0,
                        0.0, tmp0.get(),    nquad1*nquad2);
            
            // Inner product with respect to the '1' direction
            Blas::Dgemm('T', 'N', nquad2*order0, order1, nquad1,
                        1.0, tmp0.get(),  nquad1,
                             base1.get(), nquad1,
                        0.0, tmp1.get(),  nquad2*order0);

            // Inner product with respect to the '2' direction
            for (mode=i=0; i < order0; ++i)
            {
                Blas::Dgemm('T', 'N', order2-i, order1, nquad2,
                            1.0, base2.get() + mode*nquad2,  nquad2,
                                 tmp1.get() + i*nquad2,      nquad2*order0,
                            0.0, outarray.get()+mode*order1, order2-i);
                mode  += order2-i;
            }
            
            // Fix top singular vertices; performs phi_{0,q,1} +=
            // phi_1(xi_1)*phi_q(xi_2)*phi_{01}*phi_r(xi_2).
            if (m_base[0]->GetBasisType() == LibUtilities::eModified_A)
            {
                for (i = 0; i < order1; ++i)
                {
                    mode = GetMode(0,i,1);
                    outarray[mode] += Blas::Ddot(
                        nquad2, base2.get()+nquad2, 1, 
                        tmp1.get()+i*order0*nquad2+nquad2, 1);
                }
            }
        }

        /** 
         * \brief Inner product of \a inarray over region with respect to the
         * object's default expansion basis; output in \a outarray.
         */
        void StdPrismExp::v_IProductWRTDerivBase(
            const int                           dir,
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            v_IProductWRTDerivBase_SumFac(dir,inarray,outarray);
        }

        void StdPrismExp::v_IProductWRTDerivBase_MatOp(
            const int                           dir,
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            ASSERTL0(dir >= 0 && dir <= 2, "input dir is out of range");
            
            int nq = GetTotPoints();
            MatrixType mtype;

            switch (dir)
            {
                case 0:
                    mtype = eIProductWRTDerivBase0;
                    break;
                case 1:
                    mtype = eIProductWRTDerivBase1;
                    break;
                case 2:
                    mtype = eIProductWRTDerivBase2;
                    break;
            }

            StdMatrixKey      iprodmatkey(mtype,DetExpansionType(),*this);
            DNekMatSharedPtr  iprodmat = GetStdMatrix(iprodmatkey);

            Blas::Dgemv('N',m_ncoeffs,nq,1.0,iprodmat->GetPtr().get(),
                        m_ncoeffs, inarray.get(), 1, 0.0, outarray.get(), 1);
        }

        void StdPrismExp::v_IProductWRTDerivBase_SumFac(
            const int                           dir,
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            ASSERTL0(dir >= 0 && dir <= 2, "input dir is out of range");
            
            int i;
            int order0  = m_base[0]->GetNumModes ();
            int order1  = m_base[1]->GetNumModes ();
            int nquad0  = m_base[0]->GetNumPoints();
            int nquad1  = m_base[1]->GetNumPoints();
            int nquad2  = m_base[2]->GetNumPoints();

            const Array<OneD, const NekDouble> &z0 = m_base[0]->GetZ();
            const Array<OneD, const NekDouble> &z2 = m_base[2]->GetZ();
            Array<OneD, NekDouble> gfac0(nquad0);
            Array<OneD, NekDouble> gfac2(nquad2);
            Array<OneD, NekDouble> tmp0 (nquad0*nquad1*nquad2);
            Array<OneD, NekDouble> wsp  (order0*nquad2*(nquad1+order1));

            // set up geometric factor: (1+z0)/2
            for (i = 0; i < nquad0; ++i)
            {
                gfac0[i] = 0.5*(1+z0[i]);
            }
            
            // Set up geometric factor: 2/(1-z2)
            for (i = 0; i < nquad2; ++i)
            {
            	gfac2[i] = 2.0/(1-z2[i]);
            }

            // Scale first derivative term by gfac2.
            if (dir != 1)
            {
                for (i = 0; i < nquad2; ++i)
                {
                    Vmath::Smul(nquad0*nquad1,gfac2[i],
                                &inarray[0]+i*nquad0*nquad1,1,
                                &tmp0   [0]+i*nquad0*nquad1,1);
                }
                MultiplyByQuadratureMetric(tmp0,tmp0);
            }

            switch (dir)
            {
                case 0:
                {
                    IProductWRTBase_SumFacKernel(m_base[0]->GetDbdata(),
                                                 m_base[1]->GetBdata (),
                                                 m_base[2]->GetBdata (),
                                                 tmp0,outarray,wsp,
                                                 true,true,true);
                    break;
                }

                case 1:
                {
                    MultiplyByQuadratureMetric(inarray,tmp0);
                    IProductWRTBase_SumFacKernel(m_base[0]->GetBdata (),
                                                 m_base[1]->GetDbdata(),
                                                 m_base[2]->GetBdata (),
                                                 tmp0,outarray,wsp,
                                                 true,true,true);
                    break;
                }

                case 2:
                {
                    Array<OneD, NekDouble> tmp1(m_ncoeffs);
                    
                    // Scale eta_1 derivative with gfac0.
                    for(i = 0; i < nquad1*nquad2; ++i)
                    {
                        Vmath::Vmul(nquad0,&gfac0[0],1,&tmp0[0]+i*nquad0,1,&tmp0[0]+i*nquad0,1);
                    }
                    
                    IProductWRTBase_SumFacKernel(m_base[0]->GetDbdata(),
                                                 m_base[1]->GetBdata(),
                                                 m_base[2]->GetBdata(),
                                                 tmp0,tmp1,wsp,
                                                 true,true,true);
                    
                    MultiplyByQuadratureMetric(inarray, tmp0);
                    IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                                 m_base[1]->GetBdata(),
                                                 m_base[2]->GetDbdata(),
                                                 tmp0,outarray,wsp,
                                                 true,true,true);

                    Vmath::Vadd(m_ncoeffs,&tmp1[0],1,&outarray[0],1,&outarray[0],1);
                    break;
                }
            }
        }
        
        
        //---------------------------------------
        // Evaluation functions
        //---------------------------------------
        
        NekDouble StdPrismExp::v_PhysEvaluate(
            const Array<OneD, const NekDouble>& xi)
        {
            return StdPrismExp::v_PhysEvaluate(xi,m_phys);
        }

        NekDouble StdPrismExp::v_PhysEvaluate(
            const Array<OneD, const NekDouble>& xi,
            const Array<OneD, const NekDouble>& physvals)
        {
            Array<OneD, NekDouble> eta = Array<OneD, NekDouble>(3);

            if( fabs(xi[2]-1.0) < NekConstants::kNekZeroTol)
            {
                // Very top point of the prism
                eta[0] = -1.0;
                eta[1] = xi[1];
                eta[2] = 1.0;
            }
            else
            {
                // Third basis function collapsed to "pr" direction instead of
                // "qr" direction
                eta[2] = xi[2]; // eta_z = xi_z
                eta[1] = xi[1]; //eta_y = xi_y
                eta[0] = 2.0*(1.0 + xi[0])/(1.0 - xi[2]) - 1.0;
            } 

            return StdExpansion3D::v_PhysEvaluate(eta,physvals);
        }
        
 
        void StdPrismExp::v_GetCoords(Array<OneD, NekDouble>& xi_x,
                                      Array<OneD, NekDouble>& xi_y,
                                      Array<OneD, NekDouble>& xi_z)
        {
            Array<OneD, const NekDouble> etaBar_x = m_base[0]->GetZ();
            Array<OneD, const NekDouble> eta_y    = m_base[1]->GetZ();
            Array<OneD, const NekDouble> eta_z    = m_base[2]->GetZ();
            int Qx = GetNumPoints(0);
            int Qy = GetNumPoints(1);
            int Qz = GetNumPoints(2);

            // Convert collapsed coordinates into cartesian coordinates: eta --> xi
            for (int k = 0; k < Qz; ++k) {
                for (int j = 0; j < Qy; ++j) {
                    for (int i = 0; i < Qx; ++i) {
                        int s = i + Qx*(j + Qy*k);
                        xi_x[s] = (1.0 - eta_z[k])*(1.0 + etaBar_x[i]) / 2.0  -  1.0;
                        xi_y[s] = eta_y[j];
                        xi_z[s] = eta_z[k];
                    }
                }
            }
        }

        void StdPrismExp::v_FillMode(const int mode, Array<OneD, NekDouble> &outarray)
        {
            Array<OneD, NekDouble> tmp(m_ncoeffs,0.0);
            tmp[mode] = 1.0;
            StdPrismExp::v_BwdTrans(tmp, outarray);
        }


        //---------------------------------------
        // Helper functions
        //---------------------------------------

        int StdPrismExp::v_GetNverts() const
        {
            return 6;
        }
        
        int StdPrismExp::v_GetNedges() const
        {
            return 9;
        }
        
        int StdPrismExp::v_GetNfaces() const
        {
            return 5;
        }
        
        /** 
         * \brief Return Shape of region, using ShapeType enum list;
         * i.e. prism.
         */
        ExpansionType StdPrismExp::v_DetExpansionType() const
        {
            return ePrism;
        }
        
        int StdPrismExp::v_NumBndryCoeffs() const
        {
            ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                     GetBasisType(0) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                     GetBasisType(1) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(2) == LibUtilities::eModified_B ||
                     GetBasisType(2) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            
            int P = m_base[0]->GetNumModes() - 1;
            int Q = m_base[1]->GetNumModes() - 1;
            int R = m_base[2]->GetNumModes() - 1;
            
            return (P+1)*(Q+1) + 2*(Q+1)*(R+1) // 3 rect. faces in p-q and q-r planes
                + 2*(R+1) + P*(1 + 2*R - P)    // 2 tri.  faces in p-r plane
                - 3*(Q+1) - 2*(P+1) - 4*(R+1)  // subtract double counted points on edge
                + 6;                           // add vertices
        }

        int StdPrismExp::v_NumDGBndryCoeffs() const
        {
            ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                     GetBasisType(0) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                     GetBasisType(1) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(2) == LibUtilities::eModified_B ||
                     GetBasisType(2) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");

            int P = m_base[0]->GetNumModes()-1;
            int Q = m_base[1]->GetNumModes()-1;
            int R = m_base[2]->GetNumModes()-1;
            
            return (P+1)*(Q+1)               // 1 rect. face on base
                + 2*(Q+1)*(R+1)              // other 2 rect. faces
                + 2*(R+1) + P*(1 + 2*R - P); // 2 tri. faces
        }

        int StdPrismExp::v_GetEdgeNcoeffs(const int i) const
        {
            ASSERTL2(i >= 0 && i <= 8, "edge id is out of range");
            
            if (i == 0 || i == 2)
            {
                return GetBasisNumModes(0);
            }
            else if (i == 1 || i == 3 || i == 8)
            {
                return GetBasisNumModes(1);
            }
            else
            {
                return GetBasisNumModes(2);
            }
        }
        
        int StdPrismExp::v_GetFaceNcoeffs(const int i) const
        {
            ASSERTL2(i >= 0 && i <= 4, "face id is out of range");
            if (i == 0)
            {
                return GetBasisNumModes(0)*GetBasisNumModes(1);
            }
            else if (i == 1 || i == 3)
            {
                int P = GetBasisNumModes(0)-1, Q = GetBasisNumModes(2)-1;
                return Q+1 + (P*(1 + 2*Q - P))/2;
            }
            else
            {
                return GetBasisNumModes(1)*GetBasisNumModes(2);
            }
        }

        int StdPrismExp::v_GetFaceIntNcoeffs(const int i) const
        {
            ASSERTL2(i >= 0 && i <= 4, "face id is out of range");
            
            int Pi = GetBasisNumModes(0) - 2;
            int Qi = GetBasisNumModes(1) - 2;
            int Ri = GetBasisNumModes(2) - 2;
            
            if (i == 0)
            {
                return Pi * Qi;
            }
            else if (i == 1 || i == 3)
            {
                return Pi * (2*Ri - Pi - 1) / 2;
            }
            else
            {
                return Qi * Ri;
            }
        }
        
        int StdPrismExp::v_GetFaceNumPoints(const int i) const
        {
            ASSERTL2(i >= 0 && i <= 4, "face id is out of range");
            
            if (i == 0)
            {
                return m_base[0]->GetNumPoints()*
                       m_base[1]->GetNumPoints();
            }
            else if (i == 1 || i == 3)
            {
                return m_base[0]->GetNumPoints()*
                       m_base[2]->GetNumPoints();
            }
            else
            {
                return m_base[1]->GetNumPoints()*
                       m_base[2]->GetNumPoints();
            }
        }

        LibUtilities::PointsKey StdPrismExp::v_GetFacePointsKey(
            const int i, const int j) const
        {
            ASSERTL2(i >= 0 && i <= 4, "face id is out of range");
            ASSERTL2(j == 0 || j == 1, "face direction is out of range");
            
            if (i == 0)
            {
                return m_base[j]->GetPointsKey();
            }
            else if (i == 1 || i == 3)
            {
                return m_base[2*j]->GetPointsKey();
            }
            else
            {
                return m_base[j+1]->GetPointsKey();
            }
        }

        int StdPrismExp::v_CalcNumberOfCoefficients(const std::vector<unsigned int> &nummodes, 
                                                    int &modes_offset)
        {
            int nmodes = StdPrismData::getNumberOfCoefficients(
                nummodes[modes_offset],
                nummodes[modes_offset+1],
                nummodes[modes_offset+2]);

            modes_offset += 3;
            return nmodes;
        }

        LibUtilities::BasisType StdPrismExp::v_GetEdgeBasisType(const int i) const
        {
            ASSERTL2(i >= 0 && i <= 8, "edge id is out of range");
            if (i == 0 || i == 2)
            {
                return GetBasisType(0);
            }
            else if (i == 1 || i == 3 || i == 8)
            {
                return GetBasisType(1);
            }
            else
            {
                return GetBasisType(2);
            }
        }

        void StdPrismExp::v_WriteToFile(std::ofstream &outfile, 
                                        OutputFormat   format, 
                                        const bool     dumpVar, 
                                        std::string    var)
        {
            if (format == eTecplot)
            {
                int  Qx = m_base[0]->GetNumPoints();
                int  Qy = m_base[1]->GetNumPoints();
                int  Qz = m_base[2]->GetNumPoints();
                
                Array<OneD, const NekDouble> eta_x, eta_y, eta_z;
                eta_x = m_base[0]->GetZ();
                eta_y = m_base[1]->GetZ();
                eta_z = m_base[2]->GetZ();
                
                if(dumpVar)
                {
                    outfile << "Variables = z1,  z2,  z3"; 
                    outfile << ", "<< var << std::endl << std::endl;
                }      
                outfile << "Zone, I=" << Qx <<", J=" << Qy <<", K=" << Qz <<", F=Point" << std::endl;
                
                for (int k = 0; k < Qz; ++k) 
                {
                    for (int j = 0; j < Qy; ++j)
                    {
                        for (int i = 0; i < Qx; ++i)
                        {
                            outfile << 0.5*(1.0+eta_x[i])*(1.0-eta_z[k])-1.0 << " "
                                    << eta_y[j]                              << " "
                                    << eta_z[k]                              << " "
                                    << m_phys[i + Qx*(j + Qy*k)]             << std::endl;
                        }
                    }
                }
            }
            else if (format == eGnuplot)
            {
                Array<OneD, const NekDouble> eta_x, eta_y, eta_z;
                int  Qx = m_base[0]->GetNumPoints();
                int  Qy = m_base[1]->GetNumPoints();
                int  Qz = m_base[2]->GetNumPoints();
                eta_x   = m_base[0]->GetZ();
                eta_y   = m_base[1]->GetZ();
                eta_z   = m_base[2]->GetZ();
                
                for (int k = 0; k < Qz; ++k) 
                {
                    for (int j = 0; j < Qy; ++j)
                    {
                        for (int i = 0; i < Qx; ++i)
                        {
                            outfile << 0.5*(1.0+eta_x[i])*(1.0-eta_z[k])-1.0 << " "
                                    << eta_y[j]                              << " "
                                    << eta_z[k]                              << " "
                                    << m_phys[i + Qx*(j + Qy*k)]             << std::endl;
                        }
                    }
                }
            }
            else
            {
                ASSERTL0(false, "Output routine not implemented for requested type of output");
            }            
        }

        bool StdPrismExp::v_IsBoundaryInteriorExpansion()
        {
            return (m_base[0]->GetBasisType() == LibUtilities::eModified_A) &&
                   (m_base[1]->GetBasisType() == LibUtilities::eModified_A) &&
                   (m_base[2]->GetBasisType() == LibUtilities::eModified_B);
        }
        
        //---------------------------------------
        // Mappings
        //---------------------------------------


        void StdPrismExp::v_GetFaceToElementMap(
            const int                  fid,
            const Orientation      faceOrient,
            Array<OneD, unsigned int> &maparray,
            Array<OneD,          int> &signarray,
            int                        nummodesA,
            int                        nummodesB)
        {
            const LibUtilities::BasisType bType0 = GetEdgeBasisType(0);
            const LibUtilities::BasisType bType1 = GetEdgeBasisType(1);
            const LibUtilities::BasisType bType2 = GetEdgeBasisType(4);
            
            ASSERTL1(bType0 == bType1,
                     "Method only implemented if BasisType is identical"
                     "in x and y directions");
            ASSERTL1(bType0 == LibUtilities::eModified_A && 
                     bType2 == LibUtilities::eModified_B,
                     "Method only implemented for Modified_A BasisType"
                     "(x and y direction) and Modified_B BasisType (z direction)");

            int i, j, p, q, r, nFaceCoeffs, idx = 0;

            if (nummodesA == -1)
            {
                switch (fid)
                {
                    case 0:
                        nummodesA = m_base[0]->GetNumModes();
                        nummodesB = m_base[1]->GetNumModes();
                        break;
                    case 1:
                    case 3:
                        nummodesA = m_base[0]->GetNumModes();
                        nummodesB = m_base[2]->GetNumModes();
                        break;
                    case 2:
                    case 4:
                        nummodesA = m_base[1]->GetNumModes();
                        nummodesB = m_base[2]->GetNumModes();
                        break;
                }
                nFaceCoeffs = GetFaceNcoeffs(fid);
            }
            else if (fid == 1 || fid == 3)
            {
                nFaceCoeffs = nummodesB + (nummodesA-1)*(1+2*(nummodesB-1)-(nummodesA-1))/2;
            }
            else
            {
                nFaceCoeffs = nummodesA*nummodesB;
            }

            // Allocate the map array and sign array; set sign array to ones (+)
            if (maparray.num_elements() != nFaceCoeffs)
            {
                maparray = Array<OneD, unsigned int>(nFaceCoeffs);
            }
            
            if (signarray.num_elements() != nFaceCoeffs)
            {
                signarray = Array<OneD, int>(nFaceCoeffs,1);
            }
            else
            {
                fill(signarray.get(), signarray.get() + nFaceCoeffs, 1);
            }
            
            // Set up an array indexing for quads, since the ordering may need
            // to be transposed.
            Array<OneD, int> arrayindx(nFaceCoeffs,-1);
            
            if (fid != 1 && fid != 3) 
            {
                for (i = 0; i < nummodesB; i++)
                {
                    for (j = 0; j < nummodesA; j++)
                    {
                        if (faceOrient < 9)
                        {
                            arrayindx[i*nummodesA+j] = i*nummodesA+j;
                        }
                        else
                        {
                            arrayindx[i*nummodesA+j] = j*nummodesB+i;
                        }
                    }
                }
            }
            
            // Set up ordering inside each 2D face. Also for triangular faces,
            // populate signarray.
            switch (fid) 
            {
                case 0: // Bottom quad
                    for (q = 0; q < nummodesB; ++q)
                    {
                        for (p = 0; p < nummodesA; ++p)
                        {
                            maparray[arrayindx[q*nummodesA+p]] = GetMode(p,q,0);
                        }
                    }
                    break;
                    
                case 1: // Left triangle
                    for (p = 0; p < nummodesA; ++p)
                    {
                        for (r = 0; r < nummodesB-p; ++r)
                        {
                            if ((int)faceOrient == 7 && p > 1)
                            {
                                signarray[idx] = p % 2 ? -1 : 1;
                            }
                            maparray[idx++] = GetMode(p,0,r);
                        }
                    }
                    break;

                case 2: // Slanted quad
                    for (q = 0; q < nummodesA; ++q)
                    {
                        maparray[arrayindx[q]] = GetMode(1,q,0);
                    }
                    for (q = 0; q < nummodesA; ++q)
                    {
                        maparray[arrayindx[nummodesA+q]] = GetMode(0,q,1);
                    }
                    for (r = 1; r < nummodesB-1; ++r)
                    {
                        for (q = 0; q < nummodesA; ++q) 
                        {
                            maparray[arrayindx[(r+1)*nummodesA+q]] = GetMode(1,q,r);
                        }
                    }
                    break;

                case 3: // Right triangle
                    for (p = 0; p < nummodesA; ++p)
                    {
                        for (r = 0; r < nummodesB-p; ++r)
                        {
                            if ((int)faceOrient == 7 && p > 1)
                            {
                                signarray[idx] = p % 2 ? -1 : 1;
                            }
                            maparray[idx++] = GetMode(p, 1, r);
                        }
                    }
                    break;

                case 4: // Rear quad
                    for (r = 0; r < nummodesB; ++r)
                    {
                        for (q = 0; q < nummodesA; ++q)
                        {
                            maparray[arrayindx[r*nummodesA+q]] = GetMode(0, q, r);
                        }
                    }
                    break;
                    
                default:
                    ASSERTL0(false, "Face to element map unavailable.");
            }
            
            if (fid == 1 || fid == 3) 
            {
                // Triangles only have one possible orientation (base
                // direction reversed); swap edge modes.
                if ((int)faceOrient == 7)
                {
                    swap(maparray[0], maparray[nummodesA]);
                    for (i = 1; i < nummodesA-1; ++i)
                    {
                        swap(maparray[i+1], maparray[nummodesA+i]);
                    }
                }
            }
            else
            {
                // The code below is exactly the same as that taken from
                // StdHexExp and reverses the 'b' and 'a' directions as
                // appropriate (1st and 2nd if statements respectively) in
                // quadrilateral faces.
                if (faceOrient == 6 || faceOrient == 8 ||
                    faceOrient == 11 || faceOrient == 12)
                {
                    if (faceOrient < 9)
                    {
                        for (i = 3; i < nummodesB; i += 2)
                        {
                            for (j = 0; j < nummodesA; j++)
                            {
                                signarray[arrayindx[i*nummodesA+j]] *= -1;
                            }
                        }
                        
                        for (i = 0; i < nummodesA; i++)
                        {
                            swap(maparray [i], maparray [i+nummodesA]);
                            swap(signarray[i], signarray[i+nummodesA]);
                        }
                    }
                    else
                    {
                        for (i = 0; i < nummodesB; i++)
                        {
                            for (j = 3; j < nummodesA; j += 2)
                            {
                                signarray[arrayindx[i*nummodesA+j]] *= -1;
                            }
                        }

                        for (i = 0; i < nummodesB; i++)
                        {
                            swap (maparray [i], maparray [i+nummodesB]);
                            swap (signarray[i], signarray[i+nummodesB]);
                        }
                    }
                }

                if (faceOrient == 7 || faceOrient == 8 ||
                    faceOrient == 10 || faceOrient == 12)
                {
                    if (faceOrient < 9)
                    {
                        for (i = 0; i < nummodesB; i++)
                        {
                            for (j = 3; j < nummodesA; j += 2)
                            {
                                signarray[arrayindx[i*nummodesA+j]] *= -1;
                            }
                        }

                        for(i = 0; i < nummodesB; i++)
                        {
                            swap(maparray [i*nummodesA], maparray [i*nummodesA+1]);
                            swap(signarray[i*nummodesA], signarray[i*nummodesA+1]);
                        }
                    }
                    else
                    {
                        for (i = 3; i < nummodesB; i += 2)
                        {
                            for (j = 0; j < nummodesA; j++)
                            {
                                signarray[arrayindx[i*nummodesA+j]] *= -1;
                            }
                        }

                        for (i = 0; i < nummodesA; i++)
                        {
                            swap(maparray [i*nummodesB], maparray [i*nummodesB+1]);
                            swap(signarray[i*nummodesB], signarray[i*nummodesB+1]);
                        }
                    }
                }
            }
        }
        
        int StdPrismExp::v_GetVertexMap(const int vId)
        {
            ASSERTL0(GetEdgeBasisType(vId) == LibUtilities::eModified_A ||
                     GetEdgeBasisType(vId) == LibUtilities::eModified_A ||
                     GetEdgeBasisType(vId) == LibUtilities::eModified_B,
                     "Mapping not defined for this type of basis");
            
            int l = 0;
            
            switch (vId)
            {
                case 0:
                    l = GetMode(0,0,0);
                    break;
                case 1:
                    l = GetMode(1,0,0);
                    break;
                case 2:
                    l = GetMode(1,1,0);
                    break;
                case 3:
                    l = GetMode(0,1,0);
                    break;
                case 4:
                    l = GetMode(0,0,1);
                    break;
                case 5:
                    l = GetMode(0,1,1);
                    break;
                default:
                    ASSERTL0(false, "local vertex id must be between 0 and 5");
            }
            
            return l;
        }

        void StdPrismExp::v_GetEdgeInteriorMap(
            const int                  eid, 
            const Orientation      edgeOrient,
            Array<OneD, unsigned int> &maparray,
            Array<OneD, int>          &signarray)
        {
            int       i;
            bool      signChange;
            const int P              = m_base[0]->GetNumModes() - 1;
            const int Q              = m_base[1]->GetNumModes() - 1;
            const int R              = m_base[2]->GetNumModes() - 1;
            const int nEdgeIntCoeffs = v_GetEdgeNcoeffs(eid)    - 2;
            
            if (maparray.num_elements() != nEdgeIntCoeffs)
            {
                maparray = Array<OneD, unsigned int>(nEdgeIntCoeffs);
            }
            
            if(signarray.num_elements() != nEdgeIntCoeffs)
            {
                signarray = Array<OneD, int>(nEdgeIntCoeffs,1);
            }
            else
            {
                fill(signarray.get(), signarray.get()+nEdgeIntCoeffs, 1);
            }
            
            // If edge is oriented backwards, change sign of modes which have
            // degree 2n+1, n >= 1.
            signChange = edgeOrient == eBackwards;
            
            switch (eid)
            {
                case 0:
                    for (i = 2; i <= P; ++i)
                    {
                        maparray[i-2] = GetMode(i,0,0);
                    }
                    break;

                case 1:
                    for (i = 2; i <= Q; ++i)
                    {
                        maparray[i-2] = GetMode(1,i,0);
                    }
                    break;
                    
                case 2:
                    // Base quad; reverse direction.
                    //signChange = !signChange;
                    for (i = 2; i <= P; ++i)
                    {
                        maparray[i-2] = GetMode(i,1,0);
                    }
                    break;

                case 3:
                    // Base quad; reverse direction.
                    //signChange = !signChange;
                    for (i = 2; i <= Q; ++i)
                    {
                        maparray[i-2] = GetMode(0,i,0);
                    }
                    break;

                case 4:
                    for (i = 2; i <= R; ++i)
                    {
                        maparray[i-2] = GetMode(0,0,i);
                    }
                    break;

                case 5:
                    for (i = 1; i <= R-1; ++i)
                    {
                        maparray[i-1] = GetMode(1,0,i);
                    }
                    break;

                case 6:
                    for (i = 1; i <= R-1; ++i)
                    {
                        maparray[i-1] = GetMode(1,1,i);
                    }
                    break;

                case 7:
                    for (i = 2; i <= R; ++i)
                    {
                        maparray[i-2] = GetMode(0,1,i);
                    }
                    break;
                    
                case 8:
                    for (i = 2; i <= Q; ++i)
                    {
                        maparray[i-2] = GetMode(0,i,1);
                    }
                    break;
                    
                default:
                    ASSERTL0(false, "Edge not defined.");
                    break;
            }

            if (signChange)
            {
                for (i = 1; i < nEdgeIntCoeffs; i += 2)
                {
                    signarray[i] = -1;
                }
            }
        }

        void StdPrismExp::v_GetFaceInteriorMap(
            const int                  fid,
            const Orientation      faceOrient,
            Array<OneD, unsigned int> &maparray,
            Array<OneD, int>          &signarray)
        {
            const int P              = m_base[0]->GetNumModes() - 1;
            const int Q              = m_base[1]->GetNumModes() - 1;
            const int R              = m_base[2]->GetNumModes() - 1;
            const int nFaceIntCoeffs = v_GetFaceIntNcoeffs(fid);
            int       p, q, r, idx   = 0;
            int nummodesA, nummodesB, i, j;

            if (maparray.num_elements() != nFaceIntCoeffs)
            {
                maparray = Array<OneD, unsigned int>(nFaceIntCoeffs);
            }
            
            if (signarray.num_elements() != nFaceIntCoeffs)
            {
                signarray = Array<OneD, int>(nFaceIntCoeffs, 1);
            }
            else
            {
                fill(signarray.get(), signarray.get()+nFaceIntCoeffs, 1);
            }
            
            // Set up an array indexing for quad faces, since the ordering may
            // need to be transposed depending on orientation.
            Array<OneD, int> arrayindx(nFaceIntCoeffs);
            if (fid != 1 && fid != 3) 
            {
                if (fid == 0) // Base quad 
                {
                    nummodesA = P-1;
                    nummodesB = Q-1;
                }
                else if (fid == 2 || fid == 4) // front and back quad
                {
                    nummodesA = Q-1;
                    nummodesB = R-1;
                }
                
                for (i = 0; i < nummodesB; i++)
                {
                    for (j = 0; j < nummodesA; j++)
                    {
                        if (faceOrient < 9)
                        {
                            arrayindx[i*nummodesA+j] = i*nummodesA+j;
                        }
                        else
                        {
                            arrayindx[i*nummodesA+j] = j*nummodesB+i;
                        }
                    }
                }
            }
            
            switch (fid)
            {
                case 0: // Bottom quad
                    for (q = 2; q <= Q; ++q)
                    {
                        for (p = 2; p <= P; ++p)
                        {
                            maparray[arrayindx[(q-2)*nummodesA+(p-2)]] = GetMode(p,q,0);
                        }
                    }
                    break;

                case 1: // Left triangle
                    for (p = 2; p <= P; ++p)
                    {
                        for (r = 1; r <= R-p; ++r)
                        {
                            if ((int)faceOrient == 7)
                            {
                                signarray[idx] = p % 2 ? -1 : 1;
                            }
                            maparray[idx++] = GetMode(p,0,r);
                        }
                    }
                    break;
                    
                case 2: // Slanted quad
                    for (r = 1; r <= R-1; ++r)
                    {
                        for (q = 2; q <= Q; ++q)
                        {
                            maparray[arrayindx[(r-1)*nummodesA+(q-2)]] = GetMode(1, q, r);
                        }
                    }
                    break;
                    
                case 3: // Right triangle
                    for (p = 2; p <= P; ++p)
                    {
                        for (r = 1; r <= R-p; ++r)
                        {
                            if ((int)faceOrient == 7)
                            {
                                signarray[idx] = p % 2 ? -1 : 1;
                            }
                            maparray[idx++] = GetMode(p, 1, r);
                        }
                    }
                    break;

                case 4: // Back quad
                    for (r = 2; r <= R; ++r)
                    {
                        for (q = 2; q <= Q; ++q)
                        {
                            maparray[arrayindx[(r-2)*nummodesA+(q-2)]] = GetMode(0, q, r);
                        }
                    }
                    break;
                    
                default:
                    ASSERTL0(false, "Face interior map not available.");
            }

            // Triangular faces are processed in the above switch loop; for
            // remaining quad faces, set up orientation if necessary.
            if (fid == 1 || fid == 3)
                return;

            if (faceOrient == 6 || faceOrient == 8 ||
                faceOrient == 11 || faceOrient == 12)
            {
                if (faceOrient < 9)
                {
                    for (i = 1; i < nummodesB; i += 2)
                    {
                        for (j = 0; j < nummodesA; j++)
                        {
                            signarray[arrayindx[i*nummodesA+j]] *= -1;
                        }
                    }
                }
                else
                {
                    for (i = 0; i < nummodesB; i++)
                    {
                        for (j = 1; j < nummodesA; j += 2)
                        {
                            signarray[arrayindx[i*nummodesA+j]] *= -1;
                        }
                    }
                }
            }

            if (faceOrient == 7 || faceOrient == 8 ||
                faceOrient == 10 || faceOrient == 12)
            {
                if (faceOrient < 9)
                {
                    for (i = 0; i < nummodesB; i++)
                    {
                        for (j = 1; j < nummodesA; j += 2)
                        {
                            signarray[arrayindx[i*nummodesA+j]] *= -1;
                        }
                    }
                }
                else
                {
                    for (i = 1; i < nummodesB; i += 2)
                    {
                        for (j = 0; j < nummodesA; j++)
                        {
                            signarray[arrayindx[i*nummodesA+j]] *= -1;
                        }
                    }
                }
            }
        }
        
        void StdPrismExp::v_GetInteriorMap(Array<OneD, unsigned int>& outarray)
        {
            ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                     GetBasisType(0) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                     GetBasisType(1) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(2) == LibUtilities::eModified_B ||
                     GetBasisType(2) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");

            int P   = m_base[0]->GetNumModes() - 1, p;
            int Q   = m_base[1]->GetNumModes() - 1, q;
            int R   = m_base[2]->GetNumModes() - 1, r;

            int nIntCoeffs = m_ncoeffs - NumBndryCoeffs();

            if(outarray.num_elements()!=nIntCoeffs)
            {
                outarray = Array<OneD, unsigned int>(nIntCoeffs);
            }

            int idx = 0;
            
            // Loop over all interior modes.
            for (p = 2; p <= P; ++p)
            {
            	for (q = 2; q <= Q; ++q)
            	{
                    for (r = 1; r <= R-p; ++r)
                    {
                        outarray[idx++] = GetMode(p,q,r);
                    }
                }
            }
        }

        void StdPrismExp::v_GetBoundaryMap(Array<OneD, unsigned int> &maparray)
        {
            ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                     GetBasisType(0) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                     GetBasisType(1) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(2) == LibUtilities::eModified_B ||
                     GetBasisType(2) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");

            int P   = m_base[0]->GetNumModes() - 1, p;
            int Q   = m_base[1]->GetNumModes() - 1, q;
            int R   = m_base[2]->GetNumModes() - 1, r;
            int idx = 0;

            // Loop over all boundary modes (in ascending order).
            for (p = 0; p <= P; ++p) 
            {
                // First two q-r planes are entirely boundary modes.
                if (p <= 1)
                {
                    for (q = 0; q <= Q; ++q)
                    {
                        for (r = 0; r <= R-p; ++r)
                        {
                            maparray[idx++] = GetMode(p,q,r);
                        }
                    }
                }
                else
                {
                    // Remaining q-r planes contain boundary modes on the two
                    // left-hand sides and bottom edge.
                    for (q = 0; q <= Q; ++q)
                    {
                        if (q <= 1)
                        {
                            for (r = 0; r <= R-p; ++r)
                            {
                                maparray[idx++] = GetMode(p,q,r);
                            }
                        }
                        else 
                        {
                            maparray[idx++] = GetMode(p,q,0);
                        }
                    }
                }
            }
        }

        
        //---------------------------------------
        // Wrapper functions
        //---------------------------------------
        
        DNekMatSharedPtr StdPrismExp::v_GenMatrix(const StdMatrixKey &mkey)
        {
            return StdExpansion::CreateGeneralMatrix(mkey);
        }
        
        DNekMatSharedPtr StdPrismExp::v_CreateStdMatrix(const StdMatrixKey &mkey)
        {
            return StdExpansion::CreateGeneralMatrix(mkey);
        }
        
        
        //---------------------------------------
        // Private helper functions
        //---------------------------------------
        
        /**
         * @brief Compute the local mode number in the expansion for a
         * particular tensorial combination.
         *
         * Modes are numbered with the r index travelling fastest, followed by
         * q and then p, and each q-r plane is of size (R+1-p). For example,
         * with P=1, Q=2, R=3, the indexing inside each q-r plane (with r
         * increasing upwards and q to the right) is:
         *
         * p = 0:       p = 1:
         * -----------------------
         * 3   7  11
         * 2   6  10    14  17  20
         * 1   5   9    13  16  19
         * 0   4   8    12  15  18
         * 
         * Note that in this element, we must have that \f$ P <= R \f$.
         */
        int StdPrismExp::GetMode(int p, int q, int r)
        {
            int Q = m_base[1]->GetNumModes() - 1;
            int R = m_base[2]->GetNumModes() - 1;
            
            return r +                         // Skip along stacks  (r-direction)
                q*(R+1-p) +                    // Skip along columns (q-direction)
                (Q+1)*(p*R + 1-(p-2)*(p-1)/2); // Skip along rows    (p-direction)
        }

        void StdPrismExp::MultiplyByQuadratureMetric(
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            int i, j;
            int nquad0 = m_base[0]->GetNumPoints();
            int nquad1 = m_base[1]->GetNumPoints();
            int nquad2 = m_base[2]->GetNumPoints();

            const Array<OneD, const NekDouble>& w0 = m_base[0]->GetW();
            const Array<OneD, const NekDouble>& w1 = m_base[1]->GetW();
            const Array<OneD, const NekDouble>& w2 = m_base[2]->GetW();

            const Array<OneD, const NekDouble>& z2 = m_base[2]->GetZ();

            // Multiply by integration constants in x-direction
            for(i = 0; i < nquad1*nquad2; ++i)
            {
                Vmath::Vmul(nquad0, inarray.get()+i*nquad0, 1,
                            w0.get(), 1, outarray.get()+i*nquad0,1);
            }

            // Multiply by integration constants in y-direction
            for(j = 0; j < nquad2; ++j)
            {
                for(i = 0; i < nquad1; ++i)
                {
                    Blas::Dscal(nquad0,w1[i], &outarray[0]+i*nquad0 +
                                j*nquad0*nquad1,1);
                }
            }
            
            // Multiply by integration constants in z-direction; need to
            // incorporate factor (1-eta_3)/2 into weights, but only if using
            // GLL quadrature points.
            switch(m_base[2]->GetPointsType())
            {
                // Legendre inner product.
                case LibUtilities::eGaussLobattoLegendre:
                    for(i = 0; i < nquad2; ++i)
                    {
                        Blas::Dscal(nquad0*nquad1,0.25*(1-z2[i])*w2[i],
                                    &outarray[0]+i*nquad0*nquad1,1);
                    }
                    break;
                
                // (1,0) Jacobi inner product.
                case LibUtilities::eGaussRadauMAlpha1Beta0:
                    for(i = 0; i < nquad2; ++i)
                    {
                        Blas::Dscal(nquad0*nquad1, 0.5*w2[i],
                                    &outarray[0]+i*nquad0*nquad1, 1);
                    }
                    break;
                    
                default:
                    ASSERTL0(false, "Quadrature point type not supported for this element.");
                    break;
            }
        }
    }//end namespace
}//end namespace
