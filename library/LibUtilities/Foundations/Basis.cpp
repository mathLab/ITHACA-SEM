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


#include <LibUtilities/Foundations/Basis.h>
#include <LibUtilities/Polylib/Polylib.h>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>

namespace Nektar
{
    namespace LibUtilities 
    {
        bool operator<(const BasisKey &lhs, const BasisKey &rhs)
        {
            //if (lhs.m_pointsKey < rhs.m_pointsKey) return true;
            //if (lhs.m_pointsKey > rhs.m_pointsKey) return false;

            //if (lhs.m_nummodes < rhs.m_nummodes) return true;
            //if (lhs.m_nummodes > rhs.m_nummodes) return false;

            //if (lhs.m_basistype < rhs.m_basistype) return true;
            //if (lhs.m_basistype > rhs.m_basistype) return false;
        }

        bool BasisKey::opLess::operator()(const BasisKey &lhs, const BasisKey &rhs)
        {
            /*if (lhs.m_pointsKey < rhs.m_pointsKey) return true;
            if (lhs.m_pointsKey > rhs.m_pointsKey) return false;

            if (lhs.m_pointstype < rhs.m_pointstype) return true;
            if (lhs.m_pointstype > rhs.m_pointstype) return false;

            return (lhs.m_pointsid < rhs.m_pointsid);*/
        }

        std::ostream& operator<<(std::ostream& os, const BasisKey& rhs)
        {
            return os;
        }


        int Basis::BasisMem()
        {
	        switch(m_basistype)
            {
            case eOrtho_B: case eModified_B: case eModified_C:
                return  m_basisorder*(m_basisorder+1)/2*m_pointsorder;
                break;
	        case eOrtho_C:
                return  m_basisorder*(m_basisorder+1)*(m_basisorder+2)/6*m_pointsorder;
                break;
        	default:
                return m_basisorder*m_pointsorder;
        	}       
        }  
      
      void Initialize()
      {
	// Allocate Memory
	int size = BasisMem();
	m_bdata  = new double [size];
	m_dbdata = new double [size];
	
	GenBasis();
      };

      // Method used to generate appropriate basis
      void GenBasis();



        void Basis::GenBasis(){

            ASSERTL0(m_basisorder>0, "Cannot call Basis initialisation with zero or negative order");

            int i,p,q;
            double scal,*mode;
            const double *z,*w,*D;
            int localPManager = 0;

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
                    double *one_m_z_pow,fac;

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

                    double *one_m_z_pow,scal,fac;

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

        /** \brief Determine if polynomial basis can be eactly integrated
	 *  with itself
	 */
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


        /** \brief Determine if basis has collocation property,
	 *  i.e. GLL_Lagrange with Lobatto integration of appropriate order.
	 */
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
* $Log: Basis.cpp,v $
* Revision 1.2  2007/01/17 11:35:52  pvos
* updating doxygen documentation
*
* Revision 1.1  2006/06/06 18:45:14  kirby
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




