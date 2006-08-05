///////////////////////////////////////////////////////////////////////////////
//
// File StdBasis.cpp
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
// Description: Basis definition 
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <LibUtilities/Polylib/Polylib.h>
#include <StdRegions/StdBasis.h>


namespace Nektar
{
    namespace StdRegions 
    {

        /** \brief Generate a basis array stored in \a _bdata[m][i] of type \a
        init_ptype of order \a init_order at \a init_zorder points
        distributed according to \a init_z[i].

        \b init_pointstype can take the following arguments: 

        \li \b Ortho_A, \b Legendre: Orthogonal cardinal function A 
        (Legendre polynomials) \f$\widetilde{\psi}^a_p(z) = L_p(z)\f$\n\n
        %
        \f$ \mbox{\bf \_bdata}[p][i] = \widetilde{\psi}^a_p(z_i) 
        = P^{0,0}_p(z_i) \f$ 
        normalised so that 
        \f$ (\widetilde{\psi}^a_{p},\widetilde{\psi}^a_{p})=1\f$

        \li \b Ortho_B: Orthogonal cardinal function B 
        \f$\widetilde{\psi}^b_{pq}(z)\f$\n\n
        %	
        \f$ \mbox{\bf \_bdata}[m(p,q)][i]=\widetilde{\psi}^b_{pq}(z_i) = 
        [(1-z_i)]^p P^{2p+1,0}_q(z_i) \f$
        %     
        normalised so that 
        \f$(\widetilde{\psi}^b_{pq},\widetilde{\psi}^b_{pq})=1\f$\n\n
        %
        where \f$ m(p,q)=\frac{p}{2}(2\ \mbox{init\_order}-p+1)+q\f$ 
        and  \f$0 \leq p, p+q < \mbox{init\_order}, \f$

        \li \b Ortho_C: Orthogonal cadinal function C:
        \f$ \widetilde{\psi}^b_{pqr}(z)\f$ \n
        %	
        \f$ \widetilde{\psi}^c_{pqr}(z_i)
        =  [(1-z_i)]^{p+q} P^{2p+2q+1,0}_q(z_i) \f$
        %     
        normalised so that 
        \f$(\widetilde{\psi}^c_{pqr},\widetilde{\psi}^c_{pqr})=1\f$\n
        %
        Due symmetry we note that \n
        \f$ \widetilde{\psi}^c_{pqr}(z_i) = 
        \widetilde{\psi}^{c\star}_{(p+q)r}(z_i) = 
        \widetilde{\psi}^{c\star}_{sr}(z_i) \f$ \n
        %
        and so we only need calculate \n\n
        %
        \f$ \mbox{\bf \_bdata}[m(s,r)][i]=
        \widetilde{\psi}^{c\star}_{sr}(z_i)\
        = [(1-z)]^{s} P^{2s+2,0}_q(z)\f$
        %
        (Which is very similar to the cardinal function B)\n\n
        %
        where \f$ s = p+q, m(s,r)=\frac{s}{2}
        (2\ \mbox{init\_order}-s+1)+r\f$		    
        and  \f$0 \leq s, s+r < \mbox{init\_order}, \f$


        \li \b Modified_A:  Modified cardial functions A, \f$ \psi^a_p(z)\f$\n\n
        %
        \f$ \mbox{\bf \_bdata}[p][i] = \left \{ \begin{array}{ll}
        (1-z_i)/2    & p=0          \\          
        (1+z_i)/2    & p=1          \\         
        (1-z_i)/2 (1+z_i)/2 & p=2   \\
        (1-z_i)/2 (1+z_i)/2 P_1^{1,1}(z_i)  & p=3\\
        \vdots &   \vdots \\
        (1-z_i)/2 (1+z_i)/2 P_{j-2}^{1,1}(z_i) & p=j \end{array} \right . 
        \f$\n\n
        %
        Note this has been changed from Nektar where (1+z)/2 was
        stored first in list rather than (1-z)/2. The new format is
        more consistent with a left to right convention of
        interpreting the basis index \a p.\n
        %
        The vertex modes are ordered in a hierarchical fashion rather
        than the nodal format as used in the book.


        \li \b Modified_B, \b Modified_C Modified cardinal functions B and C
        \f$ \psi^b_{pq}(z)\f$, \f$ \psi^c_{pqr}(z)\f$\n
        %
        \b Note: Modified_B and  Modified_C can use the same definition
        since \f$ \psi^c_{pqr}(z) =  \psi^b_{(p+q)r}(z) \f$	 
        %
        \f$ \mbox{\bf \_bdata[m(p,q)][i]} = 
        \left \{ \begin{array}{ll}
        (1-z_i)/2            & m = 0,\ \  (p=0,q=0) \\
        (1+z_i)/2            & m = 1,\ \  (p=0,q=1) \\
        (1-z_i)/2 (1+z_i)/2  & m = 2,\ \  (p=0,q=2) \\
        \vdots       &       \vdots         \\
        (1-z_i)/2 (1+z_i)/2 P_{j-2}^{1,1}(z_i) & m = j,\ \ (p=0,q=j) \\  
        & \\
        (1-z_i)/2            & m = \mbox{order},  \ \  (p=1,q=0) \\
        (1-z_i)/2 (1+z_i)/2  & m = \mbox{order}+1,\ \  (p=1,q=2) \\
        \vdots       &       \vdots         \\
        (1-z_i)/2 (1+z_i)/2 P_{j-2}^{1,1}(z_i) & m = \mbox{order}+j,\ \ (p=1,q=j) \\  	   & \\
        \left [(1-z_i)/2\right]^2           & m=\mbox{2 order-1}\ \ 
        (p=2,q=0) \\
        \left [(1-z_i)/2\right]^2 (1+z_i)/2 & m=\mbox{2 order}\ \ 
        (p=2,q=1) \\
        \vdots	           &         \vdots                         \\
        \left [(1-z_i)/2\right]^2 (1+z_i)/2 P_{j-1}^{3,1}(z_i) &  
        m=\mbox{2 order}+j-1\ \  (p=2,q=j)       \\
        &        \\
        \left [(1-z_i)/2\right]^{k+1}          & m(p,q)  \ \ 
        (p=3,q=0) \\
        \left [(1-z_i)/2\right]^{k+1}(1+z_i)/2 & m(p,q).\ \ 
        (p=3.q=1) \\
\vdots	              &         \vdots                      \\
        \left [(1-z_i)/2\right]^{k+1} (1+z_i)/2 P_{j-1}^{2k+1,1}(z_i) & 
        m(p,q)\ \  (p=k.q=j) 
        \end{array} \right . \f$ \n \n
        %
        Note: this has been changed from Nektar where (1+z)/2 was
        stored first in list rather than (1-z)/2. The new format is
        more consistent with a left to right convention. The modes
        used in edge 0 are now also interlaced with the interior modes
        which is convenient for matrix operations in the  
        sumfactorizations\n\n
        %
        The vertex modes are ordered in a hierarchical fashion
        rather than the nodal format as used in the book.\n\n
        %
        Finally the second row which is almost identical to the first 
        row but does not include the mode (1+z)/2. This is to make the block 
        consistent with the size of the orthonormal basis and means that
        similar sumfactorisation routines can be called for both basis althoug
        the top vertex mode will remain a special case in the Modified_B
        bases.

        \li \b Fourier: Fourier bases stored in real and then imaginary parts, 
        in the interval [-1,1] i.e.\n \n
        %
        \f$ \mbox{\bf \_bdata}[p][i] = \left \{ \begin{array}{ll}
        1 & p=0 \\ 
        \sin((m+1) \pi z) & p=1 \\ 
        \cos(\pi z_i) & p=2 \\
        \sin(\pi z_i) & p=3 \\
        \cos(2 \pi z_i) & p=4 \\
        \sin(2 \pi z_i) & p=5 \\
        \vdots & \vdots \\
        \cos(m \pi z) & p = 2m \\
        \sin(m \pi z) & p = 2m+1  \end{array} \right . \f$\n \n
        %	 
        Note: Should alway be called with a factor of 2 modes
        \li \b Chebychev: Chebychev polynomials \f$ T_p(z) \f$\n\n
        %
        \f$\mbox{\bf \_bdata}[p][i] = T_p(z_i) = 
        2^{2p} (p!)^2/(2n!) P^{-1/2,-1/2}_p(z_i) \f$\n

        */

        void Basis::GenBasis(){

            ASSERTL0(m_basisorder>0, "Cannot call Basis initialisation with zero or negative order");

            int i,p,q;
            double scal,*mode;
            const double *z,*w,*D;

            PolyManagerSingleton::Instance().GetZW(m_pointstype, m_pointsorder, 
                z, w, m_alpha, m_beta);

            if(m_basistype != eFourier)
            {
                PolyManagerSingleton::Instance().GetD (m_pointstype,m_pointsorder,
                    D,m_alpha,m_beta);
            }

            switch(m_basistype)
            {
            case eOrtho_A: case eLegendre:
                mode = m_bdata;

                for (p=0; p<m_basisorder; ++p, mode += m_pointsorder)
                {
                    Polylib::jacobfd(m_pointsorder, z, mode, NULL, p, 0.0, 0.0);

                    // normalise 
                    scal = sqrt(0.5*(2.0*p+1.0));
                    for(i = 0; i < m_pointsorder; ++i)
                    {
                        mode[i] *= scal;
                    }
                }

                // define derivative basis
                Blas::Dgemm('t','n',m_pointsorder,m_basisorder,m_pointsorder,1.0,D,
                    m_pointsorder,m_bdata,m_pointsorder,0.0,
                    m_dbdata,m_pointsorder);
                break;


            case eOrtho_B:
                {
                    double *one_m_z_pow;

                    // bdata should be of size order*(order+1)/2*zorder;

                    mode = m_bdata;

                    for(i = 0; i < m_pointsorder; ++i)
                        mode[i] = 1.0;

                    mode += m_pointsorder;

                    for (q = 1; q < m_basisorder; ++q, mode += m_pointsorder)
                    {
                        Polylib::jacobfd(m_pointsorder, z, mode, NULL, q, 1.0, 0.0);
                    }

                    one_m_z_pow = m_bdata;

                    for(p = 1; p < m_basisorder; ++p)
                    {

                        for(i = 0; i < m_pointsorder; ++i)
                        {
                            mode[i] = 0.5*(1-z[i])*one_m_z_pow[i];
                        }

                        one_m_z_pow = mode;
                        mode += m_pointsorder;

                        for(q = 1; q < m_basisorder-p; ++q, mode+=m_pointsorder){
                            Polylib::jacobfd(m_pointsorder, z, mode, NULL, q, 2.0*p+1.0, 0.0);

                            for(i = 0; i < m_pointsorder; ++i)
                            {
                                mode[i] *= one_m_z_pow[i];
                            }
                        }
                    }

                    // normalise (recalling factor of 2 for weight (1-b)/2) 
                    for(p = 0, mode=m_bdata; p < m_basisorder; ++p)
                    {
                        for(q = 0; q < m_basisorder-p; ++q,mode+=m_pointsorder)
                        {
                            scal = sqrt((p+q+1.0));
                            for(i = 0; i < m_pointsorder; ++i)
                            {
                                mode[i] *= scal;
                            }
                        }
                    }

                    // define derivative basis 
                    Blas::Dgemm('t','n',m_pointsorder,m_basisorder*(m_basisorder+1)/2,m_pointsorder,1.0,D,m_pointsorder,
                        m_bdata,m_pointsorder,0.0,m_dbdata,m_pointsorder);
                }
                break; 
            case eOrtho_C:
                {

                    double *one_m_z_pow,scal;

                    // bdata should be of size _order*(order+1)*(order+2)/6*zorder;

                    mode = m_bdata;
                    for(i = 0; i < m_pointsorder; ++i)
                    {
                        mode[0] = 1.0;
                    }

                    mode += m_pointsorder;
                    for (q = 1; q < m_basisorder; ++q, mode += m_pointsorder)
                    {
                        Polylib::jacobfd(m_pointsorder, z, mode, NULL, q, 2.0, 0.0);
                    }

                    one_m_z_pow = m_bdata;
                    for(p = 1; p < m_basisorder; ++p)
                    {
                        for(i = 0; i < m_pointsorder; ++i)
                        {
                            mode[i] = 0.5*(1-z[i])*one_m_z_pow[i];
                        }

                        one_m_z_pow = mode;
                        mode += m_pointsorder;

                        for(q = 1; q < m_basisorder-p; ++q,mode+=m_pointsorder)
                        {
                            Polylib::jacobfd(m_pointsorder, z, mode, NULL, q, 2.0*p+2.0, 0.0);

                            for(i = 0; i < m_pointsorder; ++i)
                            {
                                mode[i] *= one_m_z_pow[i];
                            }
                        }
                    }

                    ASSERTL2(false, "Normalisation might need fixing");

                    for(p = 0, mode=m_bdata; p < m_basisorder; ++p)
                    {
                        for(q = 0; q < m_basisorder-p; ++q,mode+=m_pointsorder)
                        {
                            scal = sqrt((p+q+1.5));
                            for(i = 0; i < m_pointsorder; ++i)
                            {
                                mode[i] *= scal;
                            }
                        }
                    }

                    // define derivative basis 
                    Blas::Dgemm('t','n',m_pointsorder,m_basisorder*(m_basisorder+1)*
                        (m_basisorder+2)/6,m_pointsorder,1.0,D,m_pointsorder,
                        m_bdata,m_pointsorder,0.0,m_dbdata,m_pointsorder);
                }       
                break;

            case eModified_A:

                for(i = 0; i < m_pointsorder; ++i)
                {
                    m_bdata[i] = 0.5*(1-z[i]);
                    m_bdata[m_pointsorder + i] = 0.5*(1+z[i]);
                }

                mode = m_bdata + 2*m_pointsorder;

                for(p = 2; p < m_basisorder; ++p, mode += m_pointsorder)
                {
                    Polylib::jacobfd(m_pointsorder, z, mode, NULL, p-2,1.0,1.0);

                    for(i = 0; i < m_pointsorder; ++i)
                    {
                        mode[i] *= m_bdata[i]*m_bdata[m_pointsorder+i];
                    }
                }

                // define derivative basis 
                Blas::Dgemm('t','n',m_pointsorder,m_basisorder,m_pointsorder,1.0,D,
                    m_pointsorder,m_bdata,m_pointsorder,0.0,m_dbdata,
                    m_pointsorder);
                break;
            case eModified_B: case eModified_C:
                {

                    double *one_m_z_pow,*one_p_z;

                    // bdata should be of size order*(order+1)/2*zorder

                    // first fow 
                    for(i = 0; i < m_pointsorder; ++i)
                    {
                        m_bdata[0*m_pointsorder + i] = 0.5*(1-z[i]);
                        m_bdata[1*m_pointsorder + i] = 0.5*(1+z[i]);
                    }

                    mode = m_bdata + 2*m_pointsorder;

                    for(q = 2; q < m_basisorder; ++q, mode+=m_pointsorder)
                    {
                        Polylib::jacobfd(m_pointsorder, z, mode, NULL, q-2,1.0,1.0);

                        for(i = 0; i < m_pointsorder; ++i)
                        {
                            mode[i] *= m_bdata[i]*m_bdata[m_pointsorder+i];
                        }
                    }

                    // second row
                    for(i = 0; i < m_pointsorder; ++i)
                    {
                        mode[i] = 0.5*(1-z[i]);
                    }

                    mode += m_pointsorder;

                    for(q = 2; q < m_basisorder; ++q, mode+=m_pointsorder)
                    {
                        Polylib::jacobfd(m_pointsorder, z, mode, NULL, q-2,1.0,1.0);

                        for(i = 0; i < m_pointsorder; ++i)
                        {
                            mode[i] *= m_bdata[i]*m_bdata[m_pointsorder+i];
                        }
                    }

                    // third and higher rows 
                    one_m_z_pow = m_bdata;
                    one_p_z     = m_bdata+m_pointsorder;

                    for(p = 2; p < m_basisorder; ++p)
                    {
                        for(i = 0; i < m_pointsorder; ++i)
                        {
                            mode[i] = m_bdata[i]*one_m_z_pow[i];
                        }

                        one_m_z_pow  = mode;
                        mode        += m_pointsorder;

                        for(q = 1; q < m_basisorder-p; ++q, mode+=m_pointsorder)
                        {
                            Polylib::jacobfd(m_pointsorder,z,mode,NULL, q-1,2*p+1,1.0);

                            for(i = 0; i <  m_pointsorder; ++i)
                            {
                                mode[i] *= one_m_z_pow[i]*one_p_z[i];
                            }
                        }
                    }

                    Blas::Dgemm('t','n',m_pointsorder,m_basisorder*(m_basisorder+1)/2,
                        m_pointsorder,1.0,D,m_pointsorder,
                        m_bdata,m_pointsorder,0.0,m_dbdata,m_pointsorder);
                }
                break;

            case eGLL_Lagrange: 
                {
                    const double *zp,*wp;

                    mode = m_bdata;

                    // get zeros and weights  of GLL points
                    PolyManagerSingleton::Instance().GetZW((PointsType)eLobatto,m_basisorder,
                        zp, wp, 0.0 ,0.0);

                    for (p=0; p<m_basisorder; ++p,mode += m_pointsorder)
                    {
                        for(q = 0; q < m_pointsorder; ++q)
                        {
                            mode[q] = Polylib::hglj(p,z[q],zp,m_basisorder,0.0,0.0);
                        }
                    }

                    // define derivative basis 
                    Blas::Dgemm('t','n',m_pointsorder,m_basisorder,m_pointsorder,1.0,D,
                        m_pointsorder, m_bdata,m_pointsorder,0.0,
                        m_dbdata,m_pointsorder);

                }//end scope
                break;
            case eFourier:

                ASSERTL0(m_basisorder%2==0, "Fourier modes should be a factor of 2");

                for(i = 0; i < m_pointsorder; ++i)
                {
                    m_bdata[i] = 1.0;
                    m_bdata[m_pointsorder+i] = 0.0; 

                    m_dbdata[i] = m_dbdata[m_pointsorder+i] = 0.0; 
                }

                for (p=1; p < m_basisorder/2; ++p)
                {
                    for(i = 0; i < m_pointsorder; ++i)
                    {
                        m_bdata[ 2*p   *m_pointsorder+i] = cos(p*M_PI*z[i]);
                        m_bdata[(2*p+1)*m_pointsorder+i] = sin(p*M_PI*z[i]);

                        m_dbdata[ 2*p   *m_pointsorder+i] = -p*M_PI*sin(p*M_PI*z[i]);
                        m_dbdata[(2*p+1)*m_pointsorder+i] =  p*M_PI*cos(p*M_PI*z[i]);
                    }
                }

                break;

            case eChebyshev:

                mode = m_bdata;

                for (p=0,scal = 1; p<m_basisorder; ++p,mode += m_pointsorder)
                {
                    Polylib::jacobfd(m_pointsorder, z, mode, NULL, p, -0.5, -0.5);

                    for(i = 0; i < m_pointsorder; ++i)
                    {
                        mode[i] *= scal;
                    }

                    scal *= 4*(p+1)*(p+1)/(2*p+2)/(2*p+1);	
                }

                // define derivative basis 
                Blas::Dgemm('t','n',m_pointsorder,m_basisorder,m_pointsorder,1.0,D,
                    m_pointsorder, m_bdata,m_pointsorder,0.0,m_dbdata,
                    m_pointsorder);
                break;

            default:
                ASSERTL0(false, "Basis Type not known or "
                    "not implemented at this time.");
            }
        }

        /// \brief Determine if polynomial basis can be eactly integrated
        /// with itself
        int BasisKey::ExactIprodInt(void) const 
        {

            switch(m_pointstype)
            {
            case eGauss:
                if(m_pointsorder >= m_basisorder)
                {
                    return true;
                }
                break;
            case eRadauM: case eRadauP:
                if(m_pointsorder >= m_basisorder)
                {
                    return true;
                }
                break;
            case eLobatto:
                if(m_pointsorder > m_basisorder)
                {
                    return true;
                }
                break;
            }

            return false;
        }


        /// \brief Determine if basis has collocation proeprty,
        /// i.e. GLL_Lagrange with Lobatto integration of appropriate order.

        int BasisKey::Collocation() const 
        {

            if(m_basistype == eGLL_Lagrange)
            {
                if((m_pointstype == eLobatto)&&(m_alpha == 0.0)&&(m_beta == 0.0))
                {
                    if(m_basisorder == m_pointsorder)
                    {
                        return true;
                    }
                }
            }

            return false;
        }

        void BasisKey::GetInterpVec(const double zi, double *I) const
        {
            const double *z,*w;

            PolyManagerSingleton::Instance().GetZW(m_pointstype, m_pointsorder, z, w, 
                m_alpha, m_beta);
            PolyManagerSingleton::Instance().GetInterpVec(zi, m_pointstype, z,
                m_pointsorder, m_alpha, m_beta, I);
        }

        // BasisKey compared to BasisKey
        bool operator  == (const BasisKey& x, const BasisKey& y)
        {

            ASSERTL0(x.m_basistype == y.m_basistype, "Basis type not the same");

            if((x.m_pointsorder == y.m_pointsorder)&&(x.m_pointstype == y.m_pointstype) &&
                (x.m_alpha == y.m_alpha)&&(x.m_beta  == y.m_beta) && 
                (x.m_basisorder == y.m_basisorder))
            {
                return true;
            }
            else
            {
                return false;
            }
        }


        // BasisKey* compared to BasisKey
        bool operator  == (const BasisKey* x, const BasisKey& y)
        {
            ASSERTL0((*x).m_basistype == y.m_basistype, "Basis type not the same");

            if(((*x).m_pointsorder == y.m_pointsorder)&&((*x).m_pointstype == y.m_pointstype) &&
                ((*x).m_alpha == y.m_alpha)&&((*x).m_beta  == y.m_beta) && 
                ((*x).m_basisorder == y.m_basisorder))
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        // \brief BasisKey compared to BasisKey*
        bool operator  == (const BasisKey& x, const BasisKey *y)
        {
            ASSERTL0(x.m_basistype == (*y).m_basistype,"Basis type not the same");

            if((x.m_pointsorder == (*y).m_pointsorder)&&(x.m_pointstype == (*y).m_pointstype) &&
                (x.m_alpha == (*y).m_alpha)&&(x.m_beta  == (*y).m_beta)  && 
                (x.m_basisorder == (*y).m_basisorder))
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        // \brief BasisKey compared to BasisKey
        bool operator  != (const BasisKey& x, const BasisKey& y)
        {
            ASSERTL0(x.m_basistype == y.m_basistype, "Basis type not the same");

            if((x.m_pointsorder == y.m_pointsorder)&&(x.m_pointstype == y.m_pointstype) &&
                (x.m_alpha == y.m_alpha)&&(x.m_beta  ==  y.m_beta) && 
                (x.m_basisorder == y.m_basisorder))
            {
                return false;
            }
            else
            {
                return true;
            }
        }


        // BasisKey* compared to BasisKey
        bool operator  != (const BasisKey* x, const BasisKey& y)
        {

            ASSERTL0((*x).m_basistype == y.m_basistype, "Basis type not the same");

            if(((*x).m_pointsorder == y.m_pointsorder)&&((*x).m_pointstype == y.m_pointstype) &&
                ((*x).m_alpha == y.m_alpha)&&((*x).m_beta  == y.m_beta)  && 
                ((*x).m_basisorder == y.m_basisorder))
            {
                return false;
            }
            else
            {
                return true;
            }
        }


        // BasisKey compared to BasisKey*
        bool operator  != (const BasisKey& x, const BasisKey* y)
        {
            ASSERTL0(x.m_basistype == y->m_basistype, "Basis type not the same");

            if((x.m_pointsorder == y->m_pointsorder)&&(x.m_pointstype == y->m_pointstype) &&
                (x.m_alpha == y->m_alpha)&&(x.m_beta  == y->m_beta)  && 
                (x.m_basisorder == y->m_basisorder))
            {
                return false;
            }
            else
            {
                return true;
            }
        }

    } // end of namespace stdregion
} // end of namespace stdregion

/** 
* $Log: StdBasis.cpp,v $
* Revision 1.4  2006/06/13 18:05:02  sherwin
* Modifications to make MultiRegions demo ProjectLoc2D execute properly.
*
* Revision 1.3  2006/06/06 15:25:21  jfrazier
* Removed unreferenced variables and replaced ASSERTL0(false, ....) with
* NEKERROR.
*
* Revision 1.2  2006/06/01 13:43:19  kirby
* *** empty log message ***
*
* Revision 1.1  2006/05/04 18:58:30  kirby
* *** empty log message ***
*
* Revision 1.30  2006/04/25 20:23:33  jfrazier
* Various fixes to correct bugs, calls to ASSERT, etc.
*
* Revision 1.29  2006/03/05 22:11:02  sherwin
*
* Sorted out Project1D, Project2D and Project_Diff2D as well as some test scripts
*
* Revision 1.28  2006/03/03 23:04:54  sherwin
*
* Corrected Mistake in StdBasis.cpp to do with eModified_B
*
* Revision 1.27  2006/03/01 17:07:33  sherwin
*
* Added new location of polylib in header definitions
*
* Revision 1.26  2006/02/27 23:47:23  sherwin
*
* Standard coding update upto compilation of StdHexExp.cpp
*
* Revision 1.25  2006/02/26 23:37:29  sherwin
*
* Updates and compiling checks upto StdExpansions1D
*
* Revision 1.24  2006/02/19 13:26:13  sherwin
*
* Coding standard revisions so that libraries compile
*
**/ 




