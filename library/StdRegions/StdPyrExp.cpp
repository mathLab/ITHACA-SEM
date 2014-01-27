///////////////////////////////////////////////////////////////////////////////
//
// File StdPyrExp.cpp
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
// Description: pyramadic routines built upon StdExpansion3D
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdPyrExp.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <iomanip>
namespace Nektar
{
    namespace StdRegions
    {
        
        StdPyrExp::StdPyrExp() // Deafult construct of standard expansion directly called. 
        {
        }
        
        StdPyrExp::StdPyrExp(const LibUtilities::BasisKey &Ba,
                             const LibUtilities::BasisKey &Bb,
                             const LibUtilities::BasisKey &Bc) 
            : StdExpansion  (LibUtilities::StdPyrData::getNumberOfCoefficients(
                                 Ba.GetNumModes(),
                                 Bb.GetNumModes(),
                                 Bc.GetNumModes()),
                             3, Ba, Bb, Bc),
              StdExpansion3D(LibUtilities::StdPyrData::getNumberOfCoefficients(
                                 Ba.GetNumModes(),
                                 Bb.GetNumModes(),
                                 Bc.GetNumModes()),
                             Ba, Bb, Bc),
              m_map (m_ncoeffs),
              m_rmap(m_ncoeffs)
        {
            if (Ba.GetNumModes() > Bc.GetNumModes())
            {
                ASSERTL0(false, "order in 'a' direction is higher "
                         "than order in 'c' direction");
            }
            if (Bb.GetNumModes() > Bc.GetNumModes())
            {
                ASSERTL0(false, "order in 'b' direction is higher "
                         "than order in 'c' direction");
            }

            // Set up mode mapping which takes 0\leq i\leq N_coeffs -> (p,q,r)
            // of the 3D tensor product
            const int P = Ba.GetNumModes() - 1;
            const int Q = Bb.GetNumModes() - 1;
            const int R = Bc.GetNumModes() - 1;
            int cnt = 0;

            // Vertices
            m_map [cnt  ] = triple(0, 0, 0);
            m_rmap[cnt++] = 0;
            m_map [cnt  ] = triple(1, 0, 0);
            m_rmap[cnt++] = 0;
            m_map [cnt  ] = triple(1, 1, 0);
            m_rmap[cnt++] = 0;
            m_map [cnt  ] = triple(0, 1, 0);
            m_rmap[cnt++] = 0;
            m_map [cnt  ] = triple(0, 0, 1);
            m_rmap[cnt++] = 1;

            // Edge 0
            for (int i = 2; i <= P; ++i)
            {
                m_map [cnt  ] = triple    (i, 0, 0);
                m_rmap[cnt++] = GetTetMode(i, 0, 0);
            }

            // Edge 1
            for (int i = 2; i <= Q; ++i)
            {
                m_map [cnt  ] = triple    (1, i, 0);
                m_rmap[cnt++] = GetTetMode(0, i, 0);
            }

            // Edge 2
            for (int i = 2; i <= P; ++i)
            {
                m_map [cnt  ] = triple    (i, 1, 0);
                m_rmap[cnt++] = GetTetMode(i, 0, 0);
            }

            // Edge 3
            for (int i = 2; i <= Q; ++i)
            {
                m_map [cnt  ] = triple    (0, i, 0);
                m_rmap[cnt++] = GetTetMode(0, i, 0);
            }

            // Edge 4
            for (int i = 2; i <= R; ++i)
            {
                m_map [cnt  ] = triple(0, 0, i);
                m_rmap[cnt++] = i;
            }

            // Edge 5
            for (int i = 2; i <= R; ++i)
            {
                m_map [cnt  ] = triple(1, 0, i);
                m_rmap[cnt++] = i;
            }

            // Edge 6
            for (int i = 2; i <= R; ++i)
            {
                m_map [cnt  ] = triple(1, 1, i);
                m_rmap[cnt++] = i;
            }

            // Edge 7
            for (int i = 2; i <= R; ++i)
            {
                m_map [cnt  ] = triple(0, 1, i);
                m_rmap[cnt++] = i;
            }

            // Face 0 - TODO check this
            for (int j = 2; j <= Q; ++j)
            {
                for (int i = 2; i <= P; ++i)
                {
                    m_map [cnt  ] = triple    (i, j, 0);
                    m_rmap[cnt++] = GetTetMode(2, 0, 0);
                }
            }

            // Face 1
            for (int i = 2; i <= P; ++i)
            {
                for (int j = 1; j <= R-i; ++j)
                {
                    m_map [cnt  ] = triple    (i, 0, j);
                    m_rmap[cnt++] = GetTetMode(i, 0, j);
                }
            }

            // Face 2
            for (int i = 2; i <= Q; ++i)
            {
                for (int j = 1; j <= R-i; ++j)
                {
                    m_map [cnt  ] = triple    (1, i, j);
                    m_rmap[cnt++] = GetTetMode(0, i, j);
                }
            }

            // Face 3
            for (int i = 2; i <= P; ++i)
            {
                for (int j = 1; j <= R-i; ++j)
                {
                    m_map [cnt  ] = triple    (i, 1, j);
                    m_rmap[cnt++] = GetTetMode(i, 0, j);
                }
            }

            // Face 4
            for (int i = 2; i <= Q; ++i)
            {
                for (int j = 1; j <= R-i; ++j)
                {
                    m_map [cnt  ] = triple    (0, i, j);
                    m_rmap[cnt++] = GetTetMode(0, i, j);
                }
            }

            // Interior (tetrahedral modes)
            for (int i = 2; i <= P+1; ++i)
            {
            	for (int j = 1; j <= Q-i+1; ++j)
            	{
                    for (int k = 1; k <= R-i-j+1; ++k)
                    {
                        // need to go to j+1-th mode in the 'b' direction to
                        // select correct modified_a mode
                        m_map [cnt  ] = triple    (i,   j+1, k);
                        m_rmap[cnt++] = GetTetMode(i-1, j,   k);
                    }
            	}
            }
        }

        StdPyrExp::StdPyrExp(const StdPyrExp &T)
            : StdExpansion  (T),
              StdExpansion3D(T)
        {
        }


        // Destructor
        StdPyrExp::~StdPyrExp()
        {   
        } 


        //---------------------------------------
        // Differentiation/integration Methods
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
        void StdPyrExp::v_PhysDeriv(
            const Array<OneD, const NekDouble> &u_physical,
                  Array<OneD,       NekDouble> &out_dxi1,
                  Array<OneD,       NekDouble> &out_dxi2,
                  Array<OneD,       NekDouble> &out_dxi3)
        {
            // PhysDerivative implementation based on Spen's book page 152.
            int    Qx = m_base[0]->GetNumPoints();
            int    Qy = m_base[1]->GetNumPoints();
            int    Qz = m_base[2]->GetNumPoints();

            Array<OneD, NekDouble> dEta_bar1(Qx*Qy*Qz,0.0);
            Array<OneD, NekDouble> dXi2     (Qx*Qy*Qz,0.0);
            Array<OneD, NekDouble> dEta3    (Qx*Qy*Qz,0.0);
            PhysTensorDeriv(u_physical, dEta_bar1, dXi2, dEta3);

            Array<OneD, const NekDouble> eta_x, eta_y, eta_z;
            eta_x = m_base[0]->GetZ();
            eta_y = m_base[1]->GetZ();
            eta_z = m_base[2]->GetZ();

            int i, j, k, n;

            for (k = 0, n = 0; k < Qz; ++k)
            {
                for (j = 0; j < Qy; ++j)
                {
                    for (i = 0; i < Qx; ++i, ++n)
                    {
                        if (out_dxi1.num_elements() > 0)
                            out_dxi1[n] = 2.0/(1.0 - eta_z[k]) * dEta_bar1[n];
                        if (out_dxi2.num_elements() > 0)
                            out_dxi2[n] = 2.0/(1.0 - eta_z[k]) * dXi2[n];
                        if (out_dxi3.num_elements() > 0)
                            out_dxi3[n] = (1.0+eta_x[i])/(1.0-eta_z[k])*dEta_bar1[n] +
                                (1.0+eta_y[j])/(1.0-eta_z[k])*dXi2[n] + dEta3[n];
                    } 
                }
            }
        }

        void StdPyrExp::v_PhysDeriv(const int dir,
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

        void StdPyrExp::v_StdPhysDeriv(const Array<OneD, const NekDouble> &inarray, 
                                             Array<OneD,       NekDouble> &out_d0,
                                             Array<OneD,       NekDouble> &out_d1,
                                             Array<OneD,       NekDouble> &out_d2)
        {
            StdPyrExp::v_PhysDeriv(inarray, out_d0, out_d1, out_d2);
        }
        
        //---------------------------------------
        // Transforms
        //---------------------------------------
        
	/** 
         * \brief Backward transformation is evaluated at the quadrature
         * points.
         *
         * \f$ u^{\delta} (\xi_{1i}, \xi_{2j}, \xi_{3k}) = \sum_{m(pqr)} \hat
         * u_{pqr} \phi_{pqr} (\xi_{1i}, \xi_{2j}, \xi_{3k})\f$
         * 
         * Backward transformation is three dimensional tensorial expansion
         *
         * \f$ u (\xi_{1i}, \xi_{2j}, \xi_{3k}) = \sum_{p=0}^{Q_x} \psi_p^a
	 *  (\xi_{1i}) \lbrace { \sum_{q=0}^{Q_y} \psi_{q}^a (\xi_{2j})
	 *  \lbrace { \sum_{r=0}^{Q_z} \hat u_{pqr} \psi_{pqr}^c (\xi_{3k})
	 *  \rbrace} \rbrace}. \f$
	 *
         * And sumfactorizing step of the form is as:\ \ \f$ f_{pqr}
         * (\xi_{3k}) = \sum_{r=0}^{Q_z} \hat u_{pqr} \psi_{pqr}^c
         * (\xi_{3k}),\\ g_{p} (\xi_{2j}, \xi_{3k}) = \sum_{r=0}^{Q_y}
         * \psi_{p}^a (\xi_{2j}) f_{pqr} (\xi_{3k}),\\ u(\xi_{1i}, \xi_{2j},
         * \xi_{3k}) = \sum_{p=0}^{Q_x} \psi_{p}^a (\xi_{1i}) g_{p}
         * (\xi_{2j}, \xi_{3k}).  \f$
         **/
        void StdPyrExp::v_BwdTrans(const Array<OneD, const NekDouble> &inarray,
                                         Array<OneD,       NekDouble> &outarray)
        {
            ASSERTL1(m_base[1]->GetBasisType() != LibUtilities::eOrtho_B   ||
                     m_base[1]->GetBasisType() != LibUtilities::eModified_B,
                     "Basis[1] is not a general tensor type");

            ASSERTL1(m_base[2]->GetBasisType() != LibUtilities::eOrtho_C   ||
                     m_base[2]->GetBasisType() != LibUtilities::eModified_C,
                     "Basis[2] is not a general tensor type");

            const int Qx = m_base[0]->GetNumPoints();
            const int Qy = m_base[1]->GetNumPoints();
            const int Qz = m_base[2]->GetNumPoints();
            int s = 0;

            const Array<OneD, const NekDouble> &bx = m_base[0]->GetBdata();
            const Array<OneD, const NekDouble> &by = m_base[1]->GetBdata();
            const Array<OneD, const NekDouble> &bz = m_base[2]->GetBdata();

            int Q = m_base[2]->GetNumModes()-1;

            for (int k = 0; k < Qz; ++k)
            {
                for (int j = 0; j < Qy; ++j)
                {
                    for (int i = 0; i < Qx; ++i, ++s)
                    {
                        NekDouble sum = 0.0;

                        for (int cnt = 0; cnt < m_ncoeffs; ++cnt)
                        {
                            triple &idx = m_map[cnt];
                            const int p = boost::get<0>(idx);
                            const int q = boost::get<1>(idx);
                            const int r = boost::get<2>(idx);
                            NekDouble tmp;

                            tmp = inarray[cnt]*
                                bx[i + Qx*p]*
                                by[j + Qy*q]*
                                bz[k + Qz*m_rmap[cnt]];

                            if (r == 0 && p >= 2 && q >= 2)
                            {
                                int blah = (p-2+q-2) % (Q-1) + 2;

                                if (blah - 3 >= 0)
                                {
                                    tmp *= bz[k + Qz*GetTetMode(blah-2,0,0)];
                                }
                            }
                            sum += tmp;
                        }

                        // Add in contributions from singular vertices;
                        // i.e. (p,q,r) = (1,1,1),(0,1,1),(1,0,1)
                        int m = 4;
                        sum += inarray[m]*bz[k+Qz]*(bx[i+Qx]*by[j+Qy]+
                                                    bx[i   ]*by[j+Qy]+
                                                    bx[i+Qx]*by[j   ]);
                        outarray[s] = sum;
                    }
                }
            }
        }

        void StdPyrExp::v_BwdTrans_SumFacKernel(
                    const Array<OneD, const NekDouble>& base0,
                    const Array<OneD, const NekDouble>& base1,
                    const Array<OneD, const NekDouble>& base2,
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> &outarray,
                          Array<OneD, NekDouble> &wsp,
                    bool doCheckCollDir0,
                    bool doCheckCollDir1,
                    bool doCheckCollDir2)
        {
            ASSERTL0(false, "BwdTrans_SumFacKernel not yet implemented.");
        }

        void StdPyrExp::v_IProductWRTBase_SumFacKernel(
                const Array<OneD, const NekDouble>& base0,
                const Array<OneD, const NekDouble>& base1,
                const Array<OneD, const NekDouble>& base2,
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD, NekDouble> &outarray,
                      Array<OneD, NekDouble> &wsp,
                bool doCheckCollDir0,
                bool doCheckCollDir1,
                bool doCheckCollDir2)
        {
            ASSERTL0(false, "IProductWRTBase_SumFacKernel not yet implemented.");
        }

	/** \brief Forward transform from physical quadrature space
            stored in \a inarray and evaluate the expansion coefficients and
            store in \a outarray
            
            Inputs:\n
            
            - \a inarray: array of physical quadrature points to be transformed
            
            Outputs:\n
            
            - \a outarray: updated array of expansion coefficients. 
            
        */    
        void StdPyrExp::v_FwdTrans(const Array<OneD, const NekDouble> &inarray,
                                         Array<OneD,       NekDouble> &outarray)
        {
            v_IProductWRTBase(inarray,outarray);

            // get Mass matrix inverse
            StdMatrixKey      masskey(eInvMass,DetShapeType(),*this);
            DNekMatSharedPtr  matsys = GetStdMatrix(masskey);

            // copy inarray in case inarray == outarray
            DNekVec in (m_ncoeffs, outarray);
            DNekVec out(m_ncoeffs, outarray, eWrapper);

            out = (*matsys)*in;
        }
        
        
        //---------------------------------------
        // Inner product functions
        //---------------------------------------

        /** \brief  Inner product of \a inarray over region with respect to the 
            expansion basis m_base[0]->GetBdata(),m_base[1]->GetBdata(), m_base[2]->GetBdata() and return in \a outarray 
            
            Wrapper call to StdPyrExp::IProductWRTBase
            
            Input:\n
            
            - \a inarray: array of function evaluated at the physical collocation points
            
            Output:\n
            
            - \a outarray: array of inner product with respect to each basis over region
            
        */
        void StdPyrExp::v_IProductWRTBase(
            const Array<OneD, const NekDouble> &inarray, 
                  Array<OneD,       NekDouble> &outarray)
        {
            int Qx = m_base[0]->GetNumPoints();
            int Qy = m_base[1]->GetNumPoints();
            int Qz = m_base[2]->GetNumPoints();

            Array<OneD, NekDouble> tmp(Qx*Qy*Qz);

            MultiplyByStdQuadratureMetric(inarray, tmp);

            const Array<OneD, const NekDouble> &bx = m_base[0]->GetBdata();
            const Array<OneD, const NekDouble> &by = m_base[1]->GetBdata();
            const Array<OneD, const NekDouble> &bz = m_base[2]->GetBdata();
            
            int Q = m_base[1]->GetNumModes()-1;

            // Initial pyramid implementation. We need to iterate over vertices,
            // edge int, face int and then interior.
            for (int cnt = 0; cnt < m_ncoeffs; ++cnt)
            {
                // Get triple (p,q,r) which corresponds to this coefficient.
                triple &idx = m_map[cnt];
                const int  p  = boost::get<0>(idx);
                const int  q  = boost::get<1>(idx);
                const int  r  = boost::get<2>(idx);
                NekDouble tmp = 0.0;
                int s = 0;
                
                for (int k = 0; k < Qz; ++k)
                {
                    for (int j = 0; j < Qy; ++j)
                    {
                        for (int i = 0; i < Qx; ++i, ++s)
                        {
                            NekDouble tmp2;
                            tmp2 = inarray[s]*
                                bx[i + Qx*p]*
                                by[j + Qy*q]*
                                bz[k + Qz*m_rmap[cnt]];

                            if (r == 0 && p >= 2 && q >= 2)
                            {
                                int blah = (p-2+q-2) % (Q-1) + 2;

                                if (blah - 3 >= 0)
                                {
                                    tmp2 *= bz[k + Qz*GetTetMode(blah-2,0,0)];
                                }
                            }

                            tmp += tmp2;

                            // Add missing contributions from top vertex mode.
                            if (p == 0 && q == 0 && r == 1)
                            {
                                tmp2 += inarray[s] * bz[k+Qz]*(
                                    bx[i+Qx]*by[j+Qy] + 
                                    bx[i+Qx]*by[j   ] + 
                                    bx[i   ]*by[j+Qy]);
                            }
                        }
                    }
                }

                outarray[cnt] = tmp;
            }
        }

        //---------------------------------------
        // Evaluation functions
        //---------------------------------------
        
        NekDouble StdPyrExp::v_PhysEvaluate(
            const Array<OneD, const NekDouble>& xi,
            const Array<OneD, const NekDouble>& physvals)
        {
            Array<OneD, NekDouble> eta = Array<OneD, NekDouble>(3);

            if (fabs(xi[2]-1.0) < NekConstants::kNekZeroTol)
            {
                // Very top point of the pyramid
                eta[0] = -1.0;
                eta[1] = -1.0;
                eta[2] = xi[2];
            }
            else  
            {
                // Below the line-singularity -- Common case
                eta[2] = xi[2]; // eta_z = xi_z
                eta[1] = 2.0*(1.0 + xi[1])/(1.0 - xi[2]) - 1.0; 
                eta[0] = 2.0*(1.0 + xi[0])/(1.0 - xi[2]) - 1.0;
            } 
            
            return StdExpansion3D::v_PhysEvaluate(eta, physvals);
        }

        void StdPyrExp::v_GetCoords(Array<OneD, NekDouble> &xi_x, 
                                    Array<OneD, NekDouble> &xi_y,
                                    Array<OneD, NekDouble> &xi_z)
        {
            Array<OneD, const NekDouble> etaBar_x = m_base[0]->GetZ();
            Array<OneD, const NekDouble> eta_y    = m_base[1]->GetZ();
            Array<OneD, const NekDouble> eta_z    = m_base[2]->GetZ();
            int Qx = GetNumPoints(0);
            int Qy = GetNumPoints(1);
            int Qz = GetNumPoints(2);

            // Convert collapsed coordinates into cartesian coordinates: eta --> xi
            for (int k = 0; k < Qz; ++k )
            {
                for (int j = 0; j < Qy; ++j) 
                {
                    for (int i = 0; i < Qx; ++i) 
                    {
                        int s = i + Qx*(j + Qy*k);

                        xi_z[s] = eta_z[k];
                        xi_y[s] = (1.0 + eta_y[j]) * (1.0 - eta_z[k]) / 2.0  -  1.0;
                        xi_x[s] = (1.0 + etaBar_x[i]) * (1.0 - eta_z[k]) / 2.0  -  1.0;
                    }
                }
            }
        }

        void StdPyrExp::v_FillMode(const int mode, Array<OneD, NekDouble> &outarray)
        {
            Array<OneD, const NekDouble> bx = m_base[0]->GetBdata();
            Array<OneD, const NekDouble> by = m_base[1]->GetBdata();
            Array<OneD, const NekDouble> bz = m_base[2]->GetBdata();

            int Qx = m_base[0]->GetNumPoints();
            int Qy = m_base[1]->GetNumPoints();
            int Qz = m_base[2]->GetNumPoints();

            int p = boost::get<0>(m_map[mode]);
            int q = boost::get<1>(m_map[mode]);
            int r = m_rmap[mode];

            // Compute tensor product of inarray with the 3 basis functions
            for (int k = 0; k < Qz; ++k)
            {
                for (int j = 0; j < Qy; ++j)
                {
                    for (int i = 0; i < Qx; ++i)
                    {
                        int s = i + Qx*(j + Qy*k);
                        outarray[s] =
                            bx[i + Qx*p] * 
                            by[j + Qy*q] * 
                            bz[k + Qz*r];
                    }
                }
            }
        }
        
        
        //---------------------------------------
        // Helper functions
        //---------------------------------------

        int StdPyrExp::v_GetNverts() const
        {
            return 5;
        }

        int StdPyrExp::v_GetNedges() const
        {
            return 8;
        }

        int StdPyrExp::v_GetNfaces() const
        {
            return 5;
        }
        
        LibUtilities::ShapeType StdPyrExp::v_DetShapeType() const
        {
            return LibUtilities::ePyramid;
        }

        int StdPyrExp::v_NumBndryCoeffs() const
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
            
            return (P+1)*(Q+1)              // 1 rect. face in p-q plane
                + 2*(R+1) + P*(1+2*R-P)     // 2 tri. faces in p-r plane
                + 2*(R+1) + Q*(1+2*R-Q)     // 2 tri. faces in q-r plane
                - 2*(P+1)-2*(Q+1)-4*(R+1)   // subtract double counted edges
                + 5;                        // add vertices
        }

        int StdPyrExp::v_GetEdgeNcoeffs(const int i) const
        {
            ASSERTL2(i >= 0 && i <= 7, "edge id is out of range");
            
            if (i == 0 || i == 2)
            {
                return GetBasisNumModes(0);
            }
            else if (i == 1 || i == 3)
            {
                return GetBasisNumModes(1);
            }
            else
            {
                return GetBasisNumModes(2);
            }
        }

        int StdPyrExp::v_GetFaceNcoeffs(const int i) const
        {
            ASSERTL2(i >= 0 && i <= 4, "face id is out of range");
            
            if (i == 0)
            {
                return GetBasisNumModes(0)*GetBasisNumModes(1);
            }
            else if (i == 1 || i == 3)
            {
                return GetBasisNumModes(0)*GetBasisNumModes(2);
            }
            else
            {
                return GetBasisNumModes(1)*GetBasisNumModes(2);
            }
        }
        
        int StdPyrExp::v_GetFaceIntNcoeffs(const int i) const
        {
            ASSERTL2(i >= 0 && i <= 4, "face id is out of range");

            int P = m_base[0]->GetNumModes()-1;
            int Q = m_base[1]->GetNumModes()-1;
            int R = m_base[2]->GetNumModes()-1;

            if (i == 0)
            {
                return (P-1)*(Q-1);
            }
            else if (i == 1 || i == 3)
            {
                return (P-1) * (2*(R-1) - (P-1) - 1) / 2;
            }
            else
            {
                return (Q-1) * (2*(R-1) - (Q-1) - 1) / 2;
            }
        }

        int StdPyrExp::v_CalcNumberOfCoefficients(
            const std::vector<unsigned int> &nummodes, 
            int &modes_offset)
        {
            int nmodes = LibUtilities::StdPyrData::getNumberOfCoefficients(
                nummodes[modes_offset],
                nummodes[modes_offset+1],
                nummodes[modes_offset+2]);
            
            modes_offset += 3;
            return nmodes;
        }
        
        LibUtilities::BasisType StdPyrExp::v_GetEdgeBasisType(const int i) const
        {
            ASSERTL2(i >= 0 && i <= 7, "edge id is out of range");
            if (i == 0 || i == 2)
            {
                return GetBasisType(0);
            }
            else if (i == 1 || i == 3)
            {
                return GetBasisType(1);
            }
            else
            {
                return GetBasisType(2);
            }
        }


        //---------------------------------------
        // Mappings
        //---------------------------------------
        
        void StdPyrExp::v_GetFaceToElementMap(
            const int                  fid, 
            const Orientation          faceOrient,
            Array<OneD, unsigned int> &maparray,
            Array<OneD,          int> &signarray,
            int                        nummodesA, 
            int                        nummodesB)
        {
#if 0
            int i,j,P,Q;
            const int nummodes0 = m_base[0]->GetNumModes();
            const int nummodes1 = m_base[1]->GetNumModes();
            const int nummodes2 = m_base[2]->GetNumModes();
            //int nummodesA, nummodesB, P, Q;

            ASSERTL1( GetEdgeBasisType(0) == GetEdgeBasisType(1),
                      "Method only implemented if BasisType is indentical in "
                      "x and y directions");
            ASSERTL1( GetEdgeBasisType(0) == LibUtilities::eModified_A &&
                      GetEdgeBasisType(1) == LibUtilities::eModified_A &&
                      GetEdgeBasisType(4) == LibUtilities::eModified_C,
                      "Method only implemented for Modified_A BasisType (x "
                      "and y direction) and Modified_C BasisType (z "
                      "direction)");

            bool isQuad = true;

            int nFaceCoeffs = 0;
            if( fid == 0 ) // Base quad 
            {
                nummodesA = nummodes0;
                nummodesB = nummodes1;
                P = nummodesA-1;
                Q = nummodesB-1;
                nFaceCoeffs = nummodesA*nummodesB;
            }
            else if((fid == 2) || (fid == 4))
            {
                nummodesA = nummodes1;
                nummodesB = nummodes2;
                P = nummodesA-1;
                Q = nummodesB-1;
                nFaceCoeffs = R+1 + (P*(1 + 2*R - P))/2;
                isQuad = false;
            }
            else  // left and right triangles
            {
                nummodesA = nummodes0;
                nummodesB = nummodes2;
                P = nummodesA-1;
                Q = nummodesB-1;
                nFaceCoeffs = Q+1 + (P*(1 + 2*Q - P))/2;
                isQuad = false;
            }

            // Allocate the map array and sign array; set sign array to ones (+)
            if(maparray.num_elements() != nFaceCoeffs)
            {
                maparray = Array<OneD, unsigned int>(nFaceCoeffs,1);
            }
            
            if(signarray.num_elements() != nFaceCoeffs)
            {
                signarray = Array<OneD, int>(nFaceCoeffs,1);
            }
            else
            {
                fill( signarray.get() , signarray.get()+nFaceCoeffs, 1 );
            }

            Array<OneD, int> arrayindex(nFaceCoeffs,-1);

            for(int a = 0; a < nummodesA; ++a)
            {
                for(int b = 0; isQuad ? (b <  nummodesB) : (b < nummodesB - a); ++b)
                {
                    if( faceOrient < 9 ) // Not transposed
                    {
                        arrayindex[b + nummodesB*a] = b + nummodesB*a;
                    }
                    else // Transposed
                    {
                        arrayindex[b + nummodesB*a] = a + nummodesA*b;
                    }
                }
            }


            int baseCoefficient = 0;
            
            switch(fid)
            {
                // Base quad
                case 0: 
                    for(int a = 0; a < nummodesA; ++a) {
                        for(int b = 0; b < nummodesB; ++b) {
                            ASSERTL0(arrayindex[b + nummodesB*a] != -1, "arrayindex is not set up properly.");
                            maparray[ arrayindex[b + nummodesB*a] ] = b + nummodesB*a;
                        }
                }
            break;
            
            // Rear triangle
            case 3:
                baseCoefficient = (nummodes1 - 1) * nummodes2;
                for(int a = 0; a < nummodesA; ++a) {
                    for(int b = 0; b <  nummodesB - a; ++b) {
                        ASSERTL0(arrayindex[b + nummodesB*a] != -1, "arrayindex is not set up properly.");
                        maparray[ arrayindex[b + nummodesB*a] ] = baseCoefficient + b;
                    }
                    baseCoefficient += nummodes1*(nummodesB-1 - a)  +  1;
                }
            break;

            // Front triangle
            case 1: 
                for(int a = 0; a < nummodesA; ++a) {
                    for(int b = 0; b <  nummodesB - a; ++b) {
                        ASSERTL0(arrayindex[b + nummodesB*a] != -1, "arrayindex is not set up properly.");
                        maparray[ arrayindex[b + nummodesB*a] ] = baseCoefficient + b;
                    }
                    baseCoefficient += nummodes1 * (nummodes2 - a);
                }
            break;


            // Vertical triangle
            case 4: 
                for(int a = 0, n = 0; a < nummodesA; ++a) {
                    for(int b = 0; b < nummodesB - a; ++b, ++n) {
                        ASSERTL0(arrayindex[b + nummodesB*a] != -1, "arrayindex is not set up properly.");
                        maparray[ arrayindex[b + nummodesB*a] ] = n;   
                    }
                }
            break;
            

            // Slanted triangle
            case 2:
                for(int b = nummodesB-1; b >= 0; --b) {
                    for(int a = 0; a < nummodesA - b; ++a) {
                        ASSERTL0(arrayindex[b + nummodesB*a] != -1, "arrayindex is not set up properly.");
                        maparray[ arrayindex[b + nummodesB*a] ] = baseCoefficient + (a+1)*(b+1) - 1;
                    }
                    baseCoefficient += nummodesA*(b+1);
                }
            break;
            
            }

            if( (faceOrient==6) || (faceOrient==8) ||
                (faceOrient==11) || (faceOrient==12) )
            {    

                if(faceOrient<9)
                {
                    for(i = 3; i < nummodesB; i+=2)
                    {
                        for(j = 0; j < nummodesA; j++)
                        {
                            if( arrayindex[i*nummodesA+j] >= 0 )
                                signarray[ arrayindex[i*nummodesA+j] ] *= -1;
                        }
                    }
                        
                    for(i = 0; i < nummodesA; i++)
                    {
                        swap( maparray[i] , maparray[i+nummodesA] );
                        swap( signarray[i] , signarray[i+nummodesA] );
                    }
                }
                else
                {  
                    for(i = 0; i < nummodesB; i++)
                    {
                        for(j = 3; j < nummodesA; j+=2)
                        {
                            if( arrayindex[i*nummodesA+j] >= 0 )
                                signarray[ arrayindex[i*nummodesA+j] ] *= -1;
                        }
                    } 
                        
                    for(i = 0; i < nummodesB; i++)
                    {
                        swap( maparray[i] , maparray[i+nummodesB] );
                        swap( signarray[i] , signarray[i+nummodesB] );
                    }
                }
            }
                
            if( (faceOrient==7) || (faceOrient==8) ||
                (faceOrient==10) || (faceOrient==12) )
            {  
                if(faceOrient<9)
                {                                   
                    for(i = 0; i < nummodesB; i++)
                    {
                        for(j = 3; j < nummodesA; j+=2)
                        {
                            if( arrayindex[i*nummodesA+j] >= 0 )
                                signarray[ arrayindex[i*nummodesA+j] ] *= -1;
                        }
                    }                 
                        
                    for(i = 0; i < nummodesB; i++)
                    {
                        swap( maparray[i*nummodesA] , maparray[i*nummodesA+1] );
                        swap( signarray[i*nummodesA] , signarray[i*nummodesA+1] );
                    }
                }
                else
                { 
                    for(i = 3; i < nummodesB; i+=2)
                    {
                        for(j = 0; j < nummodesA; j++)
                        {
                            if( arrayindex[i*nummodesA+j] >= 0 )
                                signarray[ arrayindex[i*nummodesA+j] ] *= -1;
                        }
                    }                
                        
                    for(i = 0; i < nummodesA; i++)
                    {
                        swap( maparray[i*nummodesB] , maparray[i*nummodesB+1] );
                        swap( signarray[i*nummodesB] , signarray[i*nummodesB+1] );
                    }
                }
            }      
#endif
        }

        int StdPyrExp::v_GetVertexMap(int vId, bool useCoeffPacking)
        {
            ASSERTL0(GetEdgeBasisType(vId) == LibUtilities::eModified_A ||
                     GetEdgeBasisType(vId) == LibUtilities::eModified_A ||
                     GetEdgeBasisType(vId) == LibUtilities::eModified_B,
                     "Mapping not defined for this type of basis");
            return vId;
        }


        //---------------------------------------
        // Wrapper functions
        //---------------------------------------
        
        DNekMatSharedPtr StdPyrExp::v_GenMatrix(const StdMatrixKey &mkey)
        {
            return CreateGeneralMatrix(mkey);
        }
        
        DNekMatSharedPtr StdPyrExp::v_CreateStdMatrix(const StdMatrixKey &mkey)
        {
            return v_GenMatrix(mkey);
        }

        /**
         * 3
         * 2  6          12
         * 1  5  8       11 14      17
         * 0  4  7  9    10 13 15   16  18   19
         */

        int StdPyrExp::GetTetMode(const int I, const int J, const int K)
        {
            const int R = m_base[2]->GetNumModes();
            int i, j, cnt = 0;
            for (i = 0; i < I; ++i)
            {
                cnt += (R-i)*(R-i+1)/2;
            }

            i = R-I;
            for (j = 0; j < J; ++j)
            {
                cnt += i;
                i--;
            }

            return cnt + K;
        }

        void StdPyrExp::v_MultiplyByStdQuadratureMetric(
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            int i, j;
            
            int  nquad0 = m_base[0]->GetNumPoints();
            int  nquad1 = m_base[1]->GetNumPoints();
            int  nquad2 = m_base[2]->GetNumPoints();
            
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
            // incorporate factor [(1-eta_3)/2]^2 into weights, but only if
            // using GLL quadrature points.
            switch(m_base[2]->GetPointsType())
            {
                // Legendre inner product.
                case LibUtilities::eGaussLobattoLegendre:
                    for(i = 0; i < nquad2; ++i)
                    {
                        Blas::Dscal(nquad0*nquad1,0.125*(1-z2[i])*(1-z2[i])*w2[i],
                                    &outarray[0]+i*nquad0*nquad1,1);
                    }
                    break;
                
                // (2,0) Jacobi inner product.
                case LibUtilities::eGaussRadauMAlpha2Beta0:
                    for(i = 0; i < nquad2; ++i)
                    {
                        Blas::Dscal(nquad0*nquad1, 0.25*w2[i],
                                    &outarray[0]+i*nquad0*nquad1, 1);
                    }
                    break;
                
                default:
                    ASSERTL0(false, "Quadrature point type not supported for this element.");
                    break;
            }
        }
    }
}
